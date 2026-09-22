use crate::allele::CanonicalAllele;
use crate::coords::TranscriptPos;
use crate::data::{DataProvider, IdentifierKind, IdentifierType, TranscriptData, TranscriptSearch};
use crate::error::HgvsError;
use crate::normalize::{self, PlacedEdit};
use crate::reference::ReferenceStore;
use crate::structs::Anchor;
use crate::structs::Variant;
use crate::structs::{
    BaseOffsetInterval, BaseOffsetPosition, CVariant, GVariant, GenomicPos, LinearVariant,
    NVariant, PVariant, RVariant, SimpleInterval, SimplePosition, TranscriptVariant,
};
use crate::transcript_mapper::TranscriptMapper;
use crate::vrs::{
    vrs_type, VrsAllele, VrsBound, VrsCisPhasedBlock, VrsCopyChange, VrsCopyNumberChange,
    VrsCopyNumberCount, VrsMolecule, VrsSequenceLocation, VrsState, VrsVariation,
};

fn make_base_offset_position(
    base: crate::coords::HgvsTranscriptPos,
    combined_offset: i32,
    anchor: crate::coords::Anchor,
) -> BaseOffsetPosition {
    BaseOffsetPosition {
        base,
        offset: if combined_offset != 0 {
            Some(crate::structs::IntronicOffset(combined_offset))
        } else {
            None
        },
        anchor,
        uncertain: false,
    }
}

fn make_simple_position(base: crate::coords::HgvsGenomicPos) -> crate::structs::SimplePosition {
    crate::structs::SimplePosition {
        base,
        end: None,
        uncertain: false,
    }
}

/// `edit` rewritten against `actual`, the bases the target sequence holds
/// over the projected range. A transcript record and its genome can differ
/// at a base (RefSeq transcripts are curated against submitted mRNAs), and a
/// projection that copies the stated bases across then names bases the
/// target does not have: MUC2 c.12468C>A over a genomic G read g.…C>A where
/// the genome says G>A. So the stated reference becomes the target's, and a
/// change whose alternate is what the target already holds collapses to
/// `=`, as biocommons hgvs's replace_reference does. Edits that state no
/// bases, or state a length, are left as given.
fn replace_reference(edit: crate::edits::NaEdit, actual: &str) -> crate::edits::NaEdit {
    use crate::edits::{is_length, NaEdit};
    let stated = |r: &Option<String>| r.as_deref().is_some_and(|r| !is_length(r));
    match edit {
        NaEdit::RefAlt {
            ref_,
            alt,
            uncertain,
        } if stated(&ref_) => {
            if alt.as_deref() == Some(actual) {
                NaEdit::RefAlt {
                    ref_: None,
                    alt: None,
                    uncertain,
                }
            } else {
                NaEdit::RefAlt {
                    ref_: Some(actual.to_string()),
                    alt,
                    uncertain,
                }
            }
        }
        NaEdit::Del { ref_, uncertain } if stated(&ref_) => NaEdit::Del {
            ref_: Some(actual.to_string()),
            uncertain,
        },
        NaEdit::Dup { ref_, uncertain } if stated(&ref_) => NaEdit::Dup {
            ref_: Some(actual.to_string()),
            uncertain,
        },
        NaEdit::Inv { ref_, uncertain } if stated(&ref_) => NaEdit::Inv {
            ref_: Some(actual.to_string()),
            uncertain,
        },
        other => other,
    }
}

/// Whether `replace_reference` would change anything: the edit states bases.
fn states_bases(edit: &crate::edits::NaEdit) -> bool {
    use crate::edits::{is_length, NaEdit};
    match edit {
        NaEdit::RefAlt { ref_: Some(r), .. }
        | NaEdit::Del { ref_: Some(r), .. }
        | NaEdit::Dup { ref_: Some(r), .. }
        | NaEdit::Inv { ref_: Some(r), .. } => !is_length(r),
        _ => false,
    }
}

fn apply_strand_complement(
    edit: crate::edits::NaEdit,
    strand: crate::data::Strand,
) -> crate::edits::NaEdit {
    if strand == crate::data::Strand::Minus {
        edit.map_sequence(crate::utils::reverse_complement)
    } else {
        edit
    }
}

/// True unless the edit states literal reference bases that differ from `actual`.
fn stated_ref_matches(edit: &crate::edits::NaEdit, actual: &str) -> bool {
    match edit {
        crate::edits::NaEdit::RefAlt { .. } => edit.stated_ref().is_none_or(|r| r == actual),
        _ => true,
    }
}

/// The 0-based half-open residue range a p. interval names.
pub(crate) fn aa_interval_range(
    pos: &crate::structs::AaInterval,
) -> Result<(usize, usize), HgvsError> {
    let start = pos.start.base.to_index().0;
    let last = pos.end.as_ref().map_or(start, |e| e.base.to_index().0);
    if start < 0 || last < start {
        return Err(HgvsError::ValidationError(format!(
            "Protein interval {}..{} is not a valid range",
            start + 1,
            last + 1
        )));
    }
    Ok((start as usize, last as usize + 1))
}

/// The VRS bounds of a g. interval whose breakpoints are uncertain, HGVS
/// `(a_b)_(c_d)` with `?` for an unknown side; `None` when every position is
/// exact. A deletion `(a_b)_(c_d)` starts at an interbase coordinate in
/// `[a-1, b-1]` and ends at one in `[c, d]`; a single `(a_b)` removes one
/// base somewhere in `a..=b`.
fn uncertain_bounds(pos: &SimpleInterval) -> Option<(VrsBound, VrsBound)> {
    let ranged = |p: &SimplePosition| p.end.is_some() || p.base.is_unknown();
    let last = pos.end.as_ref().unwrap_or(&pos.start);
    if !ranged(&pos.start) && !ranged(last) {
        return None;
    }
    let known = |b: crate::coords::HgvsGenomicPos| (b.0 > 0).then_some(b.0 as usize);
    let bound = |p: &SimplePosition, interbase: fn(usize) -> usize| match (p.end, known(p.base)) {
        (Some(e), lo) => VrsBound::Range(lo.map(interbase), known(e).map(interbase)),
        (None, Some(b)) => VrsBound::Exact(interbase(b)),
        (None, None) => VrsBound::Range(None, None),
    };
    Some((bound(&pos.start, |b| b - 1), bound(last, |b| b)))
}

/// The VRS bounds of a g. interval: Range bounds when a breakpoint is
/// uncertain, else the exact interbase `[start-1, end)`.
fn interval_bounds(pos: &SimpleInterval) -> Result<(VrsBound, VrsBound), HgvsError> {
    Ok(match uncertain_bounds(pos) {
        Some(bounds) => bounds,
        None => {
            let (s, e) = simple_interval_range(pos)?;
            (VrsBound::Exact(s), VrsBound::Exact(e))
        }
    })
}

/// The g. interval of a VRS location over a non-empty range: `(a_b)_(c_d)`
/// for Range bounds, `<start+1>_<end>` for exact ones.
fn location_interval(start: VrsBound, end: VrsBound) -> Result<SimpleInterval, HgvsError> {
    Ok(match (start, end) {
        (VrsBound::Exact(s), VrsBound::Exact(e)) => {
            if e <= s {
                return Err(HgvsError::ValidationError(format!(
                    "Location end {e} is not after start {s}"
                )));
            }
            interval(s, e)
        }
        (start, end) => bounded_interval(start, end),
    })
}

/// The g. interval over the 0-based half-open range `[s, e)`.
fn interval(s: usize, e: usize) -> SimpleInterval {
    let position = |i: usize| SimplePosition {
        base: GenomicPos(i as i32).to_hgvs(),
        end: None,
        uncertain: false,
    };
    SimpleInterval {
        start: position(s),
        end: (e > s + 1).then(|| position(e - 1)),
        uncertain: false,
    }
}

/// The g. interval `(a_b)_(c_d)` from VRS bounds: the inverse of
/// `uncertain_bounds`. Exact bounds come out as parenthesis-free positions.
fn bounded_interval(start: VrsBound, end: VrsBound) -> SimpleInterval {
    use crate::coords::HgvsGenomicPos;
    let position = |b: VrsBound, to_hgvs: fn(usize) -> i32| {
        let known =
            |v: Option<usize>| v.map_or(HgvsGenomicPos::UNKNOWN, |n| HgvsGenomicPos(to_hgvs(n)));
        match b {
            VrsBound::Exact(n) => SimplePosition {
                base: HgvsGenomicPos(to_hgvs(n)),
                end: None,
                uncertain: false,
            },
            VrsBound::Range(lo, hi) => SimplePosition {
                base: known(lo),
                end: Some(known(hi)),
                uncertain: true,
            },
        }
    };
    // A single base somewhere in `a..=b` came out as `[a-1, b-1]`, `[a, b]`.
    let single = matches!((start, end), (VrsBound::Range(lo, hi), VrsBound::Range(lo2, hi2))
        if lo.map(|n| n + 1) == lo2 && hi.map(|n| n + 1) == hi2);
    SimpleInterval {
        start: position(start, |n| n as i32 + 1),
        end: (!single).then(|| position(end, |n| n as i32)),
        uncertain: false,
    }
}

/// The `g.` variant with edit `edit` at `pos` on `ac`.
fn g_variant(ac: &str, pos: SimpleInterval, edit: crate::edits::NaEdit) -> GVariant {
    GVariant::from_parts(
        ac.to_string(),
        None,
        crate::structs::PosEdit {
            pos: Some(pos),
            edit,
            uncertain: false,
            predicted: false,
        },
    )
}

/// The nucleotide edit of `var`; `None` for a protein variant.
fn na_edit(var: &crate::SequenceVariant) -> Option<&crate::edits::NaEdit> {
    use crate::SequenceVariant as SV;
    Some(match var {
        SV::Genomic(v) => &v.posedit.edit,
        SV::Mitochondrial(v) => &v.posedit.edit,
        SV::Coding(v) => &v.posedit.edit,
        SV::NonCoding(v) => &v.posedit.edit,
        SV::Rna(v) => &v.posedit.edit,
        SV::Protein(_) | SV::CisPhased(_) => return None,
    })
}

/// Whether `var` is a copy-number edit, `copyN`, which VRS renders as a
/// `CopyNumberCount` rather than an `Allele`.
fn is_copy_number(var: &crate::SequenceVariant) -> bool {
    matches!(na_edit(var), Some(crate::edits::NaEdit::NACopy { .. }))
}

/// Whether `var` is a `g.` or `m.` duplication with uncertain breakpoints,
/// `g.(a_b)_(c_d)dup`. Its bases are not known, so it has no Allele; VRS
/// renders it as a `CopyNumberChange`.
fn is_imprecise_duplication(var: &crate::SequenceVariant) -> bool {
    use crate::SequenceVariant as SV;
    let posedit = match var {
        SV::Genomic(v) => &v.posedit,
        SV::Mitochondrial(v) => &v.posedit,
        _ => return false,
    };
    matches!(posedit.edit, crate::edits::NaEdit::Dup { .. })
        && posedit.pos.as_ref().and_then(uncertain_bounds).is_some()
}

/// The number of unspecified bases an edit inserts, `insN[20]` or
/// `delinsN[(20_30)]`, as a VRS bound; `None` for an edit that spells its
/// bases or has none.
fn inserted_length(edit: &crate::edits::NaEdit) -> Option<VrsBound> {
    use crate::edits::NaEdit;
    match edit {
        NaEdit::InsLength { min, max, .. } | NaEdit::DelInsLength { min, max, .. } => {
            Some(if min == max {
                VrsBound::Exact(*min)
            } else {
                VrsBound::Range(Some(*min), Some(*max))
            })
        }
        _ => None,
    }
}

/// The p. variant for a normalised edit over residues `[s, e)` of a protein.
fn protein_variant(
    ac: &str,
    reference: &crate::reference::Reference<'_, '_>,
    s: usize,
    e: usize,
    edit: crate::edits::NaEdit,
) -> Result<PVariant, HgvsError> {
    use crate::edits::{AaEdit, NaEdit};
    use crate::structs::{AAPosition, AaInterval, PosEdit, ProteinPos};
    use crate::utils::seq1_to_aa3;
    let residue = |i: usize| -> Result<AAPosition, HgvsError> {
        Ok(AAPosition {
            base: ProteinPos(i as i32).to_hgvs(),
            aa: seq1_to_aa3(&reference.slice(i, i + 1)?),
            uncertain: false,
        })
    };
    let aa_edit = match edit {
        NaEdit::RefAlt {
            ref_: None,
            alt: None,
            ..
        } => AaEdit::Identity { uncertain: false },
        NaEdit::RefAlt {
            ref_: Some(r),
            alt: Some(a),
            ..
        } if r.len() == 1 && a.len() == 1 => AaEdit::Subst {
            ref_: seq1_to_aa3(&r),
            alt: seq1_to_aa3(&a),
            uncertain: false,
        },
        NaEdit::RefAlt { alt: Some(a), .. } => AaEdit::DelIns {
            ref_: String::new(),
            alt: seq1_to_aa3(&a),
            uncertain: false,
        },
        NaEdit::Del { .. } => AaEdit::Del {
            ref_: String::new(),
            uncertain: false,
        },
        NaEdit::Ins { alt: Some(a), .. } => AaEdit::Ins {
            alt: seq1_to_aa3(&a),
            uncertain: false,
        },
        NaEdit::Dup { .. } => AaEdit::Dup {
            ref_: None,
            uncertain: false,
        },
        other => {
            return Err(HgvsError::UnsupportedOperation(format!(
                "{other:?} has no protein form"
            )))
        }
    };
    Ok(PVariant {
        ac: ac.to_string(),
        gene: None,
        posedit: PosEdit {
            pos: Some(AaInterval {
                start: residue(s)?,
                end: if e > s + 1 {
                    Some(residue(e - 1)?)
                } else {
                    None
                },
                uncertain: false,
            }),
            edit: aa_edit,
            uncertain: false,
            predicted: false,
        },
    })
}

/// What a coding variant does to its protein.
enum CodingOutcome {
    /// A statement rather than a change: `p.?`, `p.Met1?`, `p.0?`.
    Statement(PVariant),
    Change(crate::protein::CodingChange),
}

/// `p.?`, `p.Met1?` or `p.0?` for an edit whose codons cannot be read: an
/// intronic position, or a start in the 5'UTR (deleting the whole CDS predicts
/// no protein, reaching into it disrupts the start codon, staying upstream
/// says nothing). `None` when the edit is a change to read.
fn statement_about_transcript(
    var_c: &CVariant,
    transcript: &TranscriptData,
    protein_ac: &str,
) -> Option<PVariant> {
    use crate::structs::{AAPosition, AaInterval, PosEdit, ProteinPos};
    let pos = var_c.posedit.pos.as_ref()?;
    let statement = |pos: Option<AaInterval>, value: &str| PVariant {
        ac: protein_ac.to_string(),
        gene: var_c.gene.clone(),
        posedit: PosEdit {
            pos,
            edit: crate::edits::AaEdit::Special {
                value: value.to_string(),
                uncertain: false,
            },
            uncertain: false,
            predicted: false,
        },
    };
    if has_intronic_offset(pos) {
        return Some(statement(None, "?"));
    }
    if pos.start.anchor != Anchor::CdsStart || pos.start.base.0 >= 0 {
        return None;
    }
    let cds_len = transcript
        .cds_start_index
        .zip(transcript.cds_end_index)
        .map(|(s, e)| e.0 - s.0 + 1);
    let end = pos.end.as_ref();
    let reaches_cds = end.is_some_and(|e| e.anchor == Anchor::CdsEnd || e.base.0 > 0);
    let covers_cds = end.is_some_and(|e| {
        e.anchor == Anchor::CdsEnd
            || (e.anchor == Anchor::CdsStart && cds_len.is_some_and(|n| e.base.0 >= n))
    });
    let deletes = matches!(var_c.posedit.edit, crate::edits::NaEdit::Del { .. });
    Some(if covers_cds && deletes {
        statement(None, "0?")
    } else if reaches_cds {
        let met1 = AaInterval {
            start: AAPosition {
                base: ProteinPos(0).to_hgvs(),
                aa: "Met".to_string(),
                uncertain: false,
            },
            end: None,
            uncertain: false,
        };
        statement(Some(met1), "?")
    } else {
        statement(None, "?")
    })
}

/// The CDS as 0-based transcript indices (start, last base of the stop),
/// checked against the transcript's length.
fn cds_bounds(transcript: &TranscriptData, len: usize) -> Result<(usize, usize), HgvsError> {
    let cds_start_tx = transcript
        .cds_start_index
        .ok_or_else(|| HgvsError::ValidationError("Missing CDS start".into()))?;
    let cds_end_tx = transcript
        .cds_end_index
        .ok_or_else(|| HgvsError::ValidationError("Missing CDS end".into()))?;
    let cds_start = checked_usize(cds_start_tx.0, "CDS start")?;
    let cds_end = checked_usize(cds_end_tx.0, "CDS end")?;
    if len < cds_end {
        return Err(HgvsError::ValidationError(format!(
            "Transcript sequence too short (len={}, expected at least {})",
            len, cds_end
        )));
    }
    if cds_start > len {
        return Err(HgvsError::ValidationError(format!(
            "CDS start {} out of sequence bounds {}",
            cds_start, len
        )));
    }
    Ok((cds_start, cds_end))
}

/// The change `var_c` names on the transcript, over 0-based indices: the range
/// checked against the sequence, and the stated bases against what is there.
fn resolve_in_transcript(
    var_c: &CVariant,
    am: &TranscriptMapper,
    ref_seq: &str,
) -> Result<crate::edits::ResolvedEdit, HgvsError> {
    let pos = var_c
        .posedit
        .pos
        .as_ref()
        .ok_or_else(|| HgvsError::ValidationError("Missing position".into()))?;
    let (n_start, n_end) = am.interval_to_n(pos)?;
    if n_start.0 < 0 {
        return Err(HgvsError::ValidationError(format!(
            "Position {} before transcript start",
            n_start.0
        )));
    }
    let (start, end) = (n_start.0 as usize, n_end.0 as usize);
    if end > ref_seq.len() {
        let first_bad = if start >= ref_seq.len() { start } else { end };
        return Err(HgvsError::ValidationError(format!(
            "Coordinate out of bounds: index {} is beyond transcript length {}",
            first_bad,
            ref_seq.len()
        )));
    }
    let window = |s: usize, e: usize| -> &str {
        let s = s.min(ref_seq.len());
        &ref_seq[s..e.min(ref_seq.len()).max(s)]
    };
    if let Some(stated) = var_c.posedit.edit.stated_ref() {
        let actual = window(start, end);
        if actual != stated {
            return Err(HgvsError::TranscriptMismatch {
                expected: stated.to_string(),
                found: actual.to_string(),
                start,
                end,
            });
        }
    }
    var_c
        .posedit
        .edit
        .resolve_with(start, end, |s, e| Ok(window(s, e).to_string()))
}

/// The residue a `p.` substitution puts at which 1-based position.
fn substituted_residue(var_p: &PVariant) -> Result<(usize, char), HgvsError> {
    let pos = var_p
        .posedit
        .pos
        .as_ref()
        .ok_or_else(|| HgvsError::ValidationError("Missing protein position".into()))?;
    let alt = match &var_p.posedit.edit {
        crate::edits::AaEdit::Subst { alt, .. } => alt,
        _ => {
            return Err(HgvsError::UnsupportedOperation(
                "p_to_c only supports single amino acid substitutions".into(),
            ))
        }
    };
    let raw_pos = pos.start.base.0;
    if raw_pos <= 0 {
        return Err(HgvsError::ValidationError(format!(
            "Protein position {raw_pos} is not valid (must be >= 1)"
        )));
    }
    let aa_pos = checked_usize(raw_pos, "protein position")?;
    crate::utils::aa3_to_aa1(&pos.start.aa)
        .chars()
        .next()
        .ok_or_else(|| HgvsError::ValidationError("Invalid reference AA".into()))?;
    let alt_aa = crate::utils::aa3_to_aa1(alt)
        .chars()
        .next()
        .ok_or_else(|| HgvsError::ValidationError("Invalid alternate AA".into()))?;
    Ok((aa_pos, alt_aa))
}

/// The codon for `alt_aa` fewest bases away from `ref_codon`, and whether it
/// is the only one that close (ties go to the alphabetically first).
fn closest_codon(ref_codon: &str, alt_aa: char) -> Result<(&'static str, bool), HgvsError> {
    let candidates = crate::utils::codons_for_aa(alt_aa);
    if candidates.is_empty() {
        return Err(HgvsError::ValidationError(format!(
            "No codons found for amino acid '{}'",
            alt_aa
        )));
    }
    let distance = |codon: &str| {
        codon
            .bytes()
            .zip(ref_codon.bytes())
            .filter(|(a, b)| a != b)
            .count()
    };
    let mut scored: Vec<(usize, &'static str)> =
        candidates.iter().map(|c| (distance(c), *c)).collect();
    scored.sort();
    let best = scored[0].0;
    let ties = scored.iter().filter(|(d, _)| *d == best).count();
    Ok((scored[0].1, ties == 1))
}

/// The c. edit that turns `ref_codon` (residue `aa_pos`, 1-based) into
/// `alt_codon`: a substitution for one changed base, a delins for more.
fn codon_change(
    aa_pos: usize,
    ref_codon: &str,
    alt_codon: &str,
) -> crate::structs::PosEdit<BaseOffsetInterval, crate::edits::NaEdit> {
    use crate::coords::HgvsTranscriptPos;
    let changes: Vec<(usize, u8, u8)> = ref_codon
        .bytes()
        .zip(alt_codon.bytes())
        .enumerate()
        .filter(|(_, (r, a))| r != a)
        .map(|(i, (r, a))| ((aa_pos - 1) * 3 + i + 1, r, a))
        .collect();
    let edit = crate::edits::NaEdit::RefAlt {
        ref_: Some(changes.iter().map(|(_, r, _)| *r as char).collect()),
        alt: Some(changes.iter().map(|(_, _, a)| *a as char).collect()),
        uncertain: false,
    };
    let start = changes.first().map_or(1, |(p, _, _)| *p);
    let end = changes.last().map_or(start, |(p, _, _)| *p);
    let position = |c: usize| BaseOffsetPosition {
        base: HgvsTranscriptPos(c as i32),
        offset: None,
        anchor: Anchor::CdsStart,
        uncertain: false,
    };
    crate::structs::PosEdit {
        pos: Some(BaseOffsetInterval {
            start: position(start),
            end: (end != start).then(|| position(end)),
            uncertain: false,
        }),
        edit,
        uncertain: false,
        predicted: false,
    }
}

/// `ref_` and `alt` with what they share at both ends removed, prefix first,
/// and the index the remainder starts at. Prefix first so that a fully
/// justified insertion or deletion lands at the 3' end of its run, where HGVS
/// writes it.
fn trim_ends(start: usize, ref_: &str, alt: &str) -> (usize, String, String) {
    let prefix = ref_
        .bytes()
        .zip(alt.bytes())
        .take_while(|(x, y)| x == y)
        .count();
    let (r, a) = (&ref_[prefix..], &alt[prefix..]);
    let suffix = r
        .bytes()
        .rev()
        .zip(a.bytes().rev())
        .take_while(|(x, y)| x == y)
        .count();
    (
        start + prefix,
        r[..r.len() - suffix].to_string(),
        a[..a.len() - suffix].to_string(),
    )
}

/// The HGVS edit, and the range it is written over before normalisation, for
/// "the bases over `[start, end)` (which are `ref_`) become `alt`": identity,
/// substitution, inversion, deletion, insertion or delins.
fn hgvs_edit_for(
    reference: &crate::reference::Reference<'_, '_>,
    kind: IdentifierType,
    start: usize,
    end: usize,
    ref_: &str,
    alt: &str,
) -> Result<(crate::edits::NaEdit, usize, usize), HgvsError> {
    use crate::edits::NaEdit;
    let (at, r, a) = trim_ends(start, ref_, alt);
    let uncertain = false;
    Ok(if r.is_empty() && a.is_empty() {
        let edit = NaEdit::RefAlt {
            ref_: None,
            alt: None,
            uncertain,
        };
        if end > start {
            (edit, start, end)
        } else {
            // Nothing over an empty range: say so of the base before it.
            let s = start.saturating_sub(1);
            (edit, s, s + 1)
        }
    } else if r.is_empty() && at == 0 {
        // HGVS has no insertion before the first base; it is written as a
        // delins of that base.
        let first = reference.slice(0, 1)?;
        let edit = NaEdit::RefAlt {
            ref_: Some(String::new()),
            alt: Some(format!("{a}{first}")),
            uncertain,
        };
        (edit, 0, 1)
    } else if r.is_empty() {
        let edit = NaEdit::Ins {
            alt: Some(a),
            uncertain,
        };
        (edit, at - 1, at + 1)
    } else if a.is_empty() {
        let edit = NaEdit::Del {
            ref_: None,
            uncertain,
        };
        (edit, at, at + r.len())
    } else if r.len() == 1 && a.len() == 1 {
        let edit = NaEdit::RefAlt {
            ref_: Some(r),
            alt: Some(a),
            uncertain,
        };
        (edit, at, at + 1)
    } else if kind != IdentifierType::ProteinAccession
        && r.len() > 1
        && a == crate::utils::reverse_complement(&r)
    {
        let edit = NaEdit::Inv {
            ref_: None,
            uncertain,
        };
        (edit, at, at + r.len())
    } else {
        // A delins, written without the deleted bases.
        let edit = NaEdit::RefAlt {
            ref_: Some(String::new()),
            alt: Some(a),
            uncertain,
        };
        (edit, at, at + r.len())
    })
}

/// The alphabet a transcript edit is written in.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Letters {
    /// Uppercase, T.
    Dna,
    /// Lowercase, u.
    Rna,
}

/// `edit` with its bases written in `letters`; stated lengths are untouched.
fn relettered(edit: &crate::edits::NaEdit, letters: Letters) -> crate::edits::NaEdit {
    use crate::edits::NaEdit;
    let word = |s: &String| -> String {
        if crate::edits::is_length(s) {
            s.clone()
        } else {
            match letters {
                Letters::Dna => s.to_uppercase().replace('U', "T"),
                Letters::Rna => s.to_lowercase().replace('t', "u"),
            }
        }
    };
    let opt = |o: &Option<String>| o.as_ref().map(word);
    match edit {
        NaEdit::RefAlt {
            ref_,
            alt,
            uncertain,
        } => NaEdit::RefAlt {
            ref_: opt(ref_),
            alt: opt(alt),
            uncertain: *uncertain,
        },
        NaEdit::Del { ref_, uncertain } => NaEdit::Del {
            ref_: opt(ref_),
            uncertain: *uncertain,
        },
        NaEdit::Ins { alt, uncertain } => NaEdit::Ins {
            alt: opt(alt),
            uncertain: *uncertain,
        },
        NaEdit::Dup { ref_, uncertain } => NaEdit::Dup {
            ref_: opt(ref_),
            uncertain: *uncertain,
        },
        NaEdit::Inv { ref_, uncertain } => NaEdit::Inv {
            ref_: opt(ref_),
            uncertain: *uncertain,
        },
        NaEdit::DelInsLength {
            ref_,
            min,
            max,
            uncertain,
        } => NaEdit::DelInsLength {
            ref_: opt(ref_),
            min: *min,
            max: *max,
            uncertain: *uncertain,
        },
        NaEdit::Repeat {
            ref_,
            min,
            max,
            uncertain,
        } => NaEdit::Repeat {
            ref_: opt(ref_),
            min: *min,
            max: *max,
            uncertain: *uncertain,
        },
        other => other.clone(),
    }
}

/// `posedit` re-anchored through `anchor` and re-lettered, for moving between
/// r. and c./n. Statements about the transcript (r.0, r.spl) have no other
/// spelling.
fn relettered_posedit(
    var: &dyn std::fmt::Display,
    posedit: &crate::structs::PosEdit<BaseOffsetInterval, crate::edits::NaEdit>,
    letters: Letters,
    anchor: impl Fn(Anchor) -> Anchor,
) -> Result<crate::structs::PosEdit<BaseOffsetInterval, crate::edits::NaEdit>, HgvsError> {
    if matches!(posedit.edit, crate::edits::NaEdit::Special { .. }) {
        return Err(HgvsError::UnsupportedOperation(format!(
            "{var} describes the transcript as a whole and has no other spelling"
        )));
    }
    let mut out = posedit.clone();
    if let Some(pos) = &mut out.pos {
        pos.start.anchor = anchor(pos.start.anchor);
        if let Some(end) = &mut pos.end {
            end.anchor = anchor(end.anchor);
        }
    }
    out.edit = relettered(&posedit.edit, letters);
    // Only r. has a predicted spelling, r.(123a>g); c. and n. cannot carry
    // the flag, so it is dropped on the way to DNA letters.
    if matches!(letters, Letters::Dna) {
        out.predicted = false;
    }
    Ok(out)
}

/// The 0-based half-open index range a g. interval names.
fn simple_interval_range(pos: &SimpleInterval) -> Result<(usize, usize), HgvsError> {
    let start_i = pos.start.base.to_index().0;
    if start_i < 0 {
        return Err(HgvsError::ValidationError(format!(
            "Genomic start position {} is negative or uncertain",
            start_i
        )));
    }
    let start = start_i as usize;
    let last = match &pos.end {
        Some(e) => {
            let end_i = e.base.to_index().0;
            if end_i < 0 {
                return Err(HgvsError::ValidationError(format!(
                    "Genomic end position {} is negative or uncertain",
                    end_i
                )));
            }
            end_i as usize
        }
        None => start,
    };
    let end = last
        .checked_add(1)
        .ok_or_else(|| HgvsError::ValidationError("Genomic end position overflow".into()))?;
    Ok((start, end))
}

/// Whether `after` names different bases than `before`, so HGVS positions must
/// be rewritten.
fn placement_changed(before: &PlacedEdit, after: &PlacedEdit) -> bool {
    (before.start, before.end) != (after.start, after.end)
        || before.is_insertion() != after.is_insertion()
}

/// The HGVS start index and optional end index to write for `after`.
///
/// An insertion that became a duplication names the duplicated bases, with an
/// end only when there is more than one. Anything else keeps an end position
/// exactly when the input had one.
fn hgvs_positions(
    before: &PlacedEdit,
    after: &PlacedEdit,
    had_end: bool,
) -> (usize, Option<usize>) {
    let (s, e) = after.hgvs_range();
    let last = e - 1;
    let converted = before.is_insertion() && !after.is_insertion();
    let end = if converted {
        (last > s).then_some(last)
    } else if had_end {
        Some(last)
    } else {
        None
    };
    (s, end)
}

fn has_intronic_offset(pos: &BaseOffsetInterval) -> bool {
    pos.start.offset.is_some_and(|o| o.0 != 0)
        || pos
            .end
            .as_ref()
            .is_some_and(|e| e.offset.is_some_and(|o| o.0 != 0))
}

fn checked_usize(val: i32, context: &str) -> Result<usize, HgvsError> {
    if val < 0 {
        Err(HgvsError::ValidationError(format!(
            "Invalid negative value for {}: {}",
            context, val
        )))
    } else {
        Ok(val as usize)
    }
}

/// High-level mapper for transforming variants between coordinate systems.
pub struct VariantMapper<'a> {
    /// Cached, random-access view of every sequence the provider serves; the
    /// provider itself is behind it.
    pub refs: ReferenceStore<'a>,
}

impl<'a> VariantMapper<'a> {
    /// Creates a new `VariantMapper` with the given data provider.
    pub fn new(hdp: &'a dyn DataProvider) -> Self {
        Self::from_store(ReferenceStore::new(hdp))
    }

    /// A mapper whose refget accessions come from `refget` (and which can
    /// look a refget accession up); without one they are computed from the
    /// whole sequence and `from_vrs` needs the accession passed.
    pub fn with_refget(hdp: &'a dyn DataProvider, refget: &'a dyn crate::refget::Refget) -> Self {
        Self::from_store(ReferenceStore::with_refget(hdp, refget))
    }

    /// A mapper over an existing store, for callers that keep a cache alive
    /// across mappers.
    pub fn from_store(refs: ReferenceStore<'a>) -> Self {
        VariantMapper { refs }
    }

    /// The provider behind this mapper, for transcript models and symbols.
    pub fn provider(&self) -> &'a dyn DataProvider {
        self.refs.provider()
    }

    /// Transforms a genomic variant (`g.`) to a coding cDNA variant (`c.`).
    pub fn g_to_c(&self, var_g: &GVariant, transcript_ac: &str) -> Result<CVariant, HgvsError> {
        let transcript = self
            .provider()
            .get_transcript(transcript_ac, Some(&var_g.ac))?;
        let am = TranscriptMapper::new(transcript)?;

        let pos = var_g
            .posedit
            .pos
            .as_ref()
            .ok_or_else(|| HgvsError::ValidationError("Missing genomic position".into()))?;
        let run;
        let pos = match self.repeat_run_g(var_g, pos)? {
            Some(iv) => {
                run = iv;
                &run
            }
            None => pos,
        };
        let (mut n_lo, mut off_lo) = am.g_to_n(pos.start.base.to_index())?;
        let (mut n_hi, mut off_hi) = match &pos.end {
            Some(end) => am.g_to_n(end.base.to_index())?,
            None => (n_lo, off_lo),
        };
        if n_lo.0 > n_hi.0 {
            std::mem::swap(&mut n_lo, &mut n_hi);
            std::mem::swap(&mut off_lo, &mut off_hi);
        }
        let position = |n: TranscriptPos, off: crate::structs::IntronicOffset| {
            let (c, c_off, anchor) = am.n_to_c(n)?;
            Ok::<_, HgvsError>(make_base_offset_position(
                c.to_hgvs(),
                c_off.0 + off.0,
                anchor,
            ))
        };
        let start = position(n_lo, off_lo)?;
        let end = match &pos.end {
            Some(_) => Some(position(n_hi, off_hi)?),
            None => None,
        };

        // The edit in transcript orientation, then against the transcript's
        // bases: an intronic position has none, so it is left as given.
        let mut edit = apply_strand_complement(var_g.posedit.edit.clone(), am.transcript.strand);
        let exonic = off_lo.0 == 0 && off_hi.0 == 0 && n_lo.0 >= 0;
        if states_bases(&edit) && exonic {
            let actual = self.target_bases(
                transcript_ac,
                IdentifierType::TranscriptAccession,
                n_lo.0 as usize,
                n_hi.0 as usize + 1,
            )?;
            edit = replace_reference(edit, &actual);
        }

        Ok(CVariant {
            ac: transcript_ac.to_string(),
            gene: var_g.gene.clone(),
            posedit: crate::structs::PosEdit {
                pos: Some(crate::structs::BaseOffsetInterval {
                    start,
                    end,
                    uncertain: false,
                }),
                edit,
                uncertain: var_g.posedit.uncertain,
                predicted: var_g.posedit.predicted,
            },
        })
    }

    /// A repeat names its first unit by its first base, and its run extends 3'
    /// on the sequence it is written on. Projected to the other strand, that
    /// base is the unit's last and the run lies below it, so a run detected from
    /// the projected start would be wrong. Before projecting, widen the interval
    /// to the whole run here; the projected range then covers every copy and
    /// reads as a run from its lower end on either strand.
    fn repeat_run_tx<V: TranscriptVariant>(
        &self,
        am: &TranscriptMapper,
        var: &V,
        pos: &BaseOffsetInterval,
    ) -> Result<Option<BaseOffsetInterval>, HgvsError> {
        let edit = &var.posedit().edit;
        if !matches!(edit, crate::edits::NaEdit::Repeat { .. }) || has_intronic_offset(pos) {
            return Ok(None);
        }
        let (n_start, n_end) = am.interval_to_n(pos)?;
        if n_start.0 < 0 {
            return Ok(None);
        }
        let reference = self
            .refs
            .reference(var.ac(), IdentifierType::TranscriptAccession);
        let resolved = edit.resolve(&reference, n_start.0 as usize, n_end.0 as usize)?;
        if resolved.end <= resolved.start {
            return Ok(None);
        }
        Ok(Some(BaseOffsetInterval {
            start: pos.start,
            end: Some(V::position_from_index(am, (resolved.end - 1) as i32)?),
            uncertain: pos.uncertain,
        }))
    }

    /// The genomic counterpart of [`repeat_run_tx`](Self::repeat_run_tx).
    fn repeat_run_g(
        &self,
        var_g: &GVariant,
        pos: &SimpleInterval,
    ) -> Result<Option<SimpleInterval>, HgvsError> {
        let edit = &var_g.posedit.edit;
        if !matches!(edit, crate::edits::NaEdit::Repeat { .. }) {
            return Ok(None);
        }
        let (start, end) = simple_interval_range(pos)?;
        let reference = self
            .refs
            .reference(&var_g.ac, IdentifierType::GenomicAccession);
        let resolved = edit.resolve(&reference, start, end)?;
        if resolved.end <= resolved.start {
            return Ok(None);
        }
        Ok(Some(SimpleInterval {
            start: pos.start,
            end: Some(SimplePosition {
                base: GenomicPos((resolved.end - 1) as i32).to_hgvs(),
                end: None,
                uncertain: false,
            }),
            uncertain: pos.uncertain,
        }))
    }

    /// Transforms a coding cDNA variant (`c.`) to a genomic variant (`g.`).
    pub fn c_to_g(
        &self,
        var_c: &CVariant,
        reference_ac: Option<&str>,
    ) -> Result<GVariant, HgvsError> {
        self.tx_to_g(var_c, reference_ac)
    }

    /// Transforms a non-coding cDNA variant (`n.`) to a genomic variant (`g.`).
    pub fn n_to_g(
        &self,
        var_n: &NVariant,
        reference_ac: Option<&str>,
    ) -> Result<GVariant, HgvsError> {
        self.tx_to_g(var_n, reference_ac)
    }

    /// Transforms any transcript-space variant (`c.` or `n.`) to a genomic variant (`g.`).
    pub fn tx_to_g<V: TranscriptVariant>(
        &self,
        var_c: &V,
        reference_ac: Option<&str>,
    ) -> Result<GVariant, HgvsError> {
        let transcript = self.provider().get_transcript(var_c.ac(), reference_ac)?;
        let am = TranscriptMapper::new(transcript)?;
        let target_ac = reference_ac
            .unwrap_or(am.transcript.reference_accession.as_str())
            .to_string();

        let pos = var_c
            .posedit()
            .pos
            .as_ref()
            .ok_or_else(|| HgvsError::ValidationError("Missing cDNA position".into()))?;
        let run;
        let pos = match self.repeat_run_tx(&am, var_c, pos)? {
            Some(iv) => {
                run = iv;
                &run
            }
            None => pos,
        };
        let project = |p: &BaseOffsetPosition| -> Result<GenomicPos, HgvsError> {
            let n = am.c_to_n(p.base.to_index(), p.anchor)?;
            am.n_to_g(n, p.offset.unwrap_or(crate::structs::IntronicOffset(0)))
        };
        let mut lo = project(&pos.start)?;
        let mut hi = match &pos.end {
            Some(end) => project(end)?,
            None => lo,
        };
        if lo.0 > hi.0 {
            std::mem::swap(&mut lo, &mut hi);
        }

        // The edit in genome orientation, then against the genome's bases:
        // the transcript record may not agree with the genome here.
        let mut edit = apply_strand_complement(var_c.posedit().edit.clone(), am.transcript.strand);
        if states_bases(&edit) && lo.0 >= 0 {
            let actual = self.target_bases(
                &target_ac,
                IdentifierType::GenomicAccession,
                lo.0 as usize,
                hi.0 as usize + 1,
            )?;
            edit = replace_reference(edit, &actual);
        }

        Ok(GVariant {
            ac: target_ac,
            gene: var_c.gene().map(str::to_string),
            posedit: crate::structs::PosEdit {
                pos: Some(crate::structs::SimpleInterval {
                    start: make_simple_position(lo.to_hgvs()),
                    end: pos.end.as_ref().map(|_| make_simple_position(hi.to_hgvs())),
                    uncertain: false,
                }),
                edit,
                uncertain: var_c.posedit().uncertain,
                predicted: var_c.posedit().predicted,
            },
        })
    }

    /// The bases of `ac` over `[start, end)`, which must all exist: a projected
    /// range past the end of the target is a data error, not a shorter edit.
    fn target_bases(
        &self,
        ac: &str,
        kind: IdentifierType,
        start: usize,
        end: usize,
    ) -> Result<String, HgvsError> {
        let bases = self.refs.reference(ac, kind).slice(start, end)?;
        if bases.len() != end - start {
            return Err(HgvsError::ValidationError(format!(
                "{ac} has {} bases over [{start}, {end}); the projected range runs past its end",
                bases.len()
            )));
        }
        Ok(bases)
    }

    /// Discovers all possible cDNA consequences for a genomic variant.
    pub fn g_to_c_all(
        &self,
        var_g: &GVariant,
        searcher: &dyn TranscriptSearch,
    ) -> Result<Vec<CVariant>, HgvsError> {
        let pos = var_g
            .posedit
            .pos
            .as_ref()
            .ok_or_else(|| HgvsError::ValidationError("Missing position".into()))?;
        let start_0 = pos.start.base.to_index().0;
        let end_0 = pos
            .end
            .as_ref()
            .map_or(start_0 + 1, |e| e.base.to_index().0 + 1);

        let transcripts = searcher.get_transcripts_for_region(&var_g.ac, start_0, end_0)?;
        if transcripts.is_empty() {
            return Err(HgvsError::ValidationError(format!(
                "No transcripts found for region {}:{}-{}",
                var_g.ac, start_0, end_0
            )));
        }

        let mut results = Vec::new();
        let mut errors = Vec::new();

        for tx_ac in transcripts {
            match self.g_to_c(var_g, &tx_ac) {
                Ok(vc) => results.push(vc),
                Err(e) => errors.push(format!("{}: {}", tx_ac, e)),
            }
        }

        if results.is_empty() && !errors.is_empty() {
            return Err(HgvsError::ValidationError(format!(
                "All mapping attempts failed: {}",
                errors.join("; ")
            )));
        }

        Ok(results)
    }

    /// `protein_ac` if given, else the protein the provider maps `transcript_ac` to.
    fn protein_accession(
        &self,
        transcript_ac: &str,
        protein_ac: Option<&str>,
    ) -> Result<String, HgvsError> {
        if let Some(ac) = protein_ac {
            return Ok(ac.to_string());
        }
        Ok(self
            .provider()
            .get_symbol_accessions(
                transcript_ac,
                IdentifierKind::Transcript,
                IdentifierKind::Protein,
            )?
            .first()
            .ok_or_else(|| {
                HgvsError::ValidationError(format!(
                    "No protein accession found for {}",
                    transcript_ac
                ))
            })?
            .1
            .clone())
    }

    /// The protein consequence of an r. variant: a statement about the
    /// transcript (`r.0`, `r.spl`, `r.?`, `r.=`) becomes the matching statement
    /// about the protein; anything else is predicted from its c. spelling.
    pub fn r_to_p(&self, r: &RVariant, protein_ac: Option<&str>) -> Result<PVariant, HgvsError> {
        if let crate::edits::NaEdit::Special { value, .. } = &r.posedit.edit {
            let (p, predicted) = match value.as_str() {
                "0" => ("0", false),
                "0?" => ("0?", false),
                "=" => ("=", true),
                _ => ("?", false), // r.?, r.spl, r.spl?
            };
            return Ok(PVariant {
                ac: self.protein_accession(&r.ac, protein_ac)?,
                gene: r.gene.clone(),
                posedit: crate::structs::PosEdit {
                    pos: None,
                    edit: crate::edits::AaEdit::Special {
                        value: p.to_string(),
                        uncertain: false,
                    },
                    uncertain: false,
                    predicted,
                },
            });
        }
        self.c_to_p(&self.r_to_c(r)?, protein_ac)
    }

    /// Transforms a coding cDNA variant (`c.`) to a protein variant (`p.`).
    pub fn c_to_p(
        &self,
        var_c: &CVariant,
        protein_ac: Option<&str>,
    ) -> Result<PVariant, HgvsError> {
        match self.coding_outcome(var_c, protein_ac)? {
            CodingOutcome::Statement(p) => Ok(p),
            CodingOutcome::Change(change) => {
                let mut var_p = crate::protein::describe(&change)?;
                var_p.posedit.predicted = true;
                Ok(var_p)
            }
        }
    }

    /// The protein allele of a coding variant: on the protein the transcript
    /// is paired with, the residues from the first change to the end become
    /// the residues the edited transcript encodes, up to its new stop. Unlike
    /// the allele of a p. variant this covers frameshifts, extensions and
    /// stop losses, which are computed rather than described. The translated
    /// CDS must be the protein the provider serves; a difference is an
    /// annotation error and is reported as one.
    pub fn protein_allele(
        &self,
        var_c: &CVariant,
        protein_ac: Option<&str>,
    ) -> Result<CanonicalAllele, HgvsError> {
        let strip_stop = |s: String| s.trim_end_matches('*').to_string();
        // For a silent change the allele is the reference over the codons the
        // edit touched, as the p. description names them.
        let mut touched: Option<(usize, usize)> = None;
        let (protein_ac, reference, alternate) = match self.coding_outcome(var_c, protein_ac)? {
            CodingOutcome::Statement(p) => match &p.posedit.edit {
                // A deleted CDS: no protein at all.
                crate::edits::AaEdit::Special { value, .. } if value == "0?" => {
                    let whole = self
                        .refs
                        .reference(&p.ac, IdentifierType::ProteinAccession)
                        .whole()?;
                    (p.ac.clone(), strip_stop(whole), String::new())
                }
                _ => {
                    return Err(HgvsError::UnsupportedOperation(format!(
                        "{p} names no protein sequence, so it has no allele"
                    )))
                }
            },
            CodingOutcome::Change(change) => {
                let (r, a) = crate::protein::proteins(&change)?;
                if r == a {
                    let first = change.edit.start / 3;
                    let last = change.edit.end.max(change.edit.start + 1) - 1;
                    touched = Some((
                        first.min(r.len()),
                        (last / 3 + 1).clamp(first.min(r.len()), r.len()),
                    ));
                }
                (change.protein_ac, r, a)
            }
        };
        let np = self
            .refs
            .reference(&protein_ac, IdentifierType::ProteinAccession);
        let actual = strip_stop(np.whole()?);
        if actual != reference {
            let at = actual
                .bytes()
                .zip(reference.bytes())
                .position(|(x, y)| x != y)
                .unwrap_or(actual.len().min(reference.len()));
            return Err(HgvsError::ValidationError(format!(
                "The CDS of {} does not translate to {protein_ac}: they differ from residue {} ({} vs {}); the annotation pairs them wrongly",
                var_c.ac,
                at + 1,
                &actual[at.min(actual.len())..(at + 10).min(actual.len())],
                &reference[at.min(reference.len())..(at + 10).min(reference.len())],
            )));
        }
        let edit = match touched {
            Some((s, e)) => crate::edits::ResolvedEdit {
                start: s,
                end: e,
                ref_: reference[s..e].to_string(),
                alt: reference[s..e].to_string(),
            },
            None => crate::edits::ResolvedEdit {
                start: 0,
                end: reference.len(),
                ref_: reference,
                alt: alternate,
            },
        };
        CanonicalAllele::canonicalize(&np, &protein_ac, &edit)
    }

    /// The protein a coding variant leaves, in 1-letter code up to its stop:
    /// empty when the CDS is deleted, `None` when nothing can be said (an
    /// intronic position, a disrupted start codon).
    pub fn predicted_protein(
        &self,
        var_c: &CVariant,
        protein_ac: Option<&str>,
    ) -> Result<Option<String>, HgvsError> {
        Ok(match self.coding_outcome(var_c, protein_ac)? {
            CodingOutcome::Statement(p) => match &p.posedit.edit {
                crate::edits::AaEdit::Special { value, .. } if value.starts_with('0') => {
                    Some(String::new())
                }
                _ => None,
            },
            CodingOutcome::Change(change) => Some(crate::protein::proteins(&change)?.1),
        })
    }

    /// The GA4GH VRS 2.0 Allele of a coding variant's protein consequence, on
    /// the protein sequence, with the predicted p. description as its
    /// expression when there is one.
    pub fn protein_vrs(
        &self,
        var_c: &CVariant,
        protein_ac: Option<&str>,
    ) -> Result<VrsAllele, HgvsError> {
        let allele = self.protein_allele(var_c, protein_ac)?;
        let refget = self
            .refs
            .reference(&allele.accession, IdentifierType::ProteinAccession)
            .refget_accession()?;
        let hgvs = self.c_to_p(var_c, protein_ac).ok().map(|p| p.to_string());
        Ok(VrsAllele::new(
            &allele,
            &refget,
            VrsMolecule::Protein,
            hgvs.as_deref().map(|h| ("hgvs.p", h)),
        ))
    }

    /// What a coding variant does to the protein: a statement (`p.?`,
    /// `p.Met1?`, `p.0?`) when the edit lies outside or across the CDS
    /// boundary, else the resolved change to describe or read.
    fn coding_outcome(
        &self,
        var_c: &CVariant,
        protein_ac: Option<&str>,
    ) -> Result<CodingOutcome, HgvsError> {
        let protein_ac = self.protein_accession(&var_c.ac, protein_ac)?;
        let transcript = self.provider().get_transcript(&var_c.ac, None)?;
        if let Some(statement) = statement_about_transcript(var_c, &transcript, &protein_ac) {
            return Ok(CodingOutcome::Statement(statement));
        }
        let ref_seq = self
            .refs
            .reference(&var_c.ac, IdentifierType::TranscriptAccession)
            .whole()?;
        let (cds_start, cds_end) = cds_bounds(&transcript, ref_seq.len())?;
        let am = TranscriptMapper::new(transcript)?;
        let resolved = resolve_in_transcript(var_c, &am, &ref_seq)?;
        let rel = |i: usize| {
            i.checked_sub(cds_start).ok_or_else(|| {
                HgvsError::ValidationError(format!("Position {} before the CDS start", i))
            })
        };
        Ok(CodingOutcome::Change(crate::protein::CodingChange {
            coding: ref_seq[cds_start..].to_string(),
            cds_len: cds_end + 1 - cds_start,
            edit: crate::edits::ResolvedEdit {
                start: rel(resolved.start)?,
                end: rel(resolved.end)?,
                ref_: resolved.ref_,
                alt: resolved.alt,
            },
            protein_ac,
        }))
    }

    /// Back-converts a single amino acid substitution (p.Xxx###Yyy) to a coding variant (c.).
    ///
    /// Only handles `AaEdit::Subst` currently. For ambiguous back-conversions (multiple
    /// codons with the same minimum nucleotide distance), the first alphabetically is chosen
    /// and `is_unique` is set to false.
    ///
    /// Returns `(CVariant, is_unique)`.
    pub fn p_to_c(
        &self,
        var_p: &PVariant,
        transcript_ac: Option<&str>,
    ) -> Result<(CVariant, bool), HgvsError> {
        let (aa_pos, alt_aa) = substituted_residue(var_p)?;
        let tx_ac = self.transcript_for_protein(&var_p.ac, transcript_ac)?;
        let ref_codon = self.codon_at(&tx_ac, aa_pos)?;
        let (alt_codon, is_unique) = closest_codon(&ref_codon, alt_aa)?;
        let c_variant = CVariant {
            ac: tx_ac,
            gene: var_p.gene.clone(),
            posedit: codon_change(aa_pos, &ref_codon, alt_codon),
        };
        Ok((c_variant, is_unique))
    }

    /// `transcript_ac` if given, else the transcript the provider maps `protein_ac` to.
    fn transcript_for_protein(
        &self,
        protein_ac: &str,
        transcript_ac: Option<&str>,
    ) -> Result<String, HgvsError> {
        if let Some(ac) = transcript_ac {
            return Ok(ac.to_string());
        }
        Ok(self
            .provider()
            .get_symbol_accessions(
                protein_ac,
                IdentifierKind::Protein,
                IdentifierKind::Transcript,
            )?
            .first()
            .ok_or_else(|| {
                HgvsError::ValidationError(format!(
                    "No transcript accession found for {}",
                    protein_ac
                ))
            })?
            .1
            .clone())
    }

    /// The codon of residue `aa_pos` (1-based) on `tx_ac`, uppercase.
    fn codon_at(&self, tx_ac: &str, aa_pos: usize) -> Result<String, HgvsError> {
        let transcript = self.provider().get_transcript(tx_ac, None)?;
        let cds_start = transcript
            .cds_start_index
            .ok_or_else(|| HgvsError::ValidationError("No CDS start for transcript".into()))?
            .0 as usize;
        let tx_seq = self
            .refs
            .reference(tx_ac, IdentifierType::TranscriptAccession)
            .whole()?;
        let codon_start = cds_start + (aa_pos - 1) * 3;
        let codon_end = codon_start + 3;
        if codon_end > tx_seq.len() {
            return Err(HgvsError::ValidationError(format!(
                "Codon position {}-{} out of range for sequence length {}",
                codon_start,
                codon_end,
                tx_seq.len()
            )));
        }
        Ok(tx_seq[codon_start..codon_end].to_uppercase())
    }

    /// Normalizes a variant to its 3' most position.
    /// Checks a variant's stated reference bases against the reference sequence.
    ///
    /// Returns `Ok(true)` when the stated reference matches, or when the edit
    /// states no reference (or only a length) so there is nothing to check.
    /// Intronic c. positions are accepted unchecked. Only g. and c. variants
    /// are supported.
    pub fn validate(&self, var: &crate::SequenceVariant) -> Result<bool, HgvsError> {
        match var {
            crate::SequenceVariant::Genomic(v) => self.validate_linear(v),
            crate::SequenceVariant::Mitochondrial(v) => self.validate_linear(v),
            crate::SequenceVariant::Coding(v) => self.validate_transcript(v),
            crate::SequenceVariant::NonCoding(v) => self.validate_transcript(v),
            crate::SequenceVariant::Protein(v) => self.validate_protein(v),
            crate::SequenceVariant::Rna(r) => self.validate(&self.r_as_transcript(r)?),
            // A cis allele holds when every member does.
            crate::SequenceVariant::CisPhased(cis) => {
                for m in &cis.members {
                    if !self.validate(m)? {
                        return Ok(false);
                    }
                }
                Ok(true)
            }
        }
    }

    /// Whether the residues a p. variant names at its positions, and any it
    /// states as reference, are what the protein sequence holds.
    fn validate_protein(&self, v: &PVariant) -> Result<bool, HgvsError> {
        let pos = v
            .posedit
            .pos
            .as_ref()
            .ok_or_else(|| HgvsError::ValidationError("Missing position".into()))?;
        let (start, end) = aa_interval_range(pos)?;
        let actual = self
            .refs
            .reference(&v.ac, IdentifierType::ProteinAccession)
            .slice(start, end)?;
        if actual.len() != end - start {
            return Ok(false); // the range runs past the end of the protein
        }
        let named = |aa: &str, at: usize| -> Result<bool, HgvsError> {
            Ok(aa.is_empty() || crate::utils::residues_1(aa)? == actual[at..=at])
        };
        if !named(&pos.start.aa, 0)? {
            return Ok(false);
        }
        if let Some(e) = &pos.end {
            if !named(&e.aa, actual.len() - 1)? {
                return Ok(false);
            }
        }
        use crate::edits::AaEdit;
        let stated = match &v.posedit.edit {
            AaEdit::Subst { ref_, .. } | AaEdit::DelIns { ref_, .. } | AaEdit::Del { ref_, .. } => {
                Some(ref_.as_str())
            }
            AaEdit::Dup { ref_, .. } | AaEdit::RefAlt { ref_, .. } => ref_.as_deref(),
            AaEdit::Repeat { ref_, .. } => ref_
                .as_deref()
                .filter(|u| !u.chars().all(|c| c.is_ascii_digit())),
            _ => None,
        }
        .filter(|r| !r.is_empty());
        match stated {
            Some(r) => Ok(crate::utils::residues_1(r)? == actual),
            None => Ok(true),
        }
    }

    fn validate_linear<L: LinearVariant>(&self, v: &L) -> Result<bool, HgvsError> {
        let pos = v
            .posedit()
            .pos
            .as_ref()
            .ok_or_else(|| HgvsError::ValidationError("Missing position".into()))?;
        let (start, end) = simple_interval_range(pos)?;
        let ref_seq = self
            .refs
            .reference(v.ac(), IdentifierType::GenomicAccession)
            .slice(start, end)?;
        Ok(stated_ref_matches(&v.posedit().edit, &ref_seq))
    }

    fn validate_transcript<V: TranscriptVariant>(&self, v: &V) -> Result<bool, HgvsError> {
        let transcript = self.provider().get_transcript(v.ac(), None)?;
        let pos = v
            .posedit()
            .pos
            .as_ref()
            .ok_or_else(|| HgvsError::ValidationError("Missing position".into()))?;
        if has_intronic_offset(pos) {
            return Ok(true);
        }
        let am = TranscriptMapper::new(transcript)?;
        let (n_start, n_end) = am.interval_to_n(pos)?;
        let start_idx = checked_usize(n_start.0, "transcript start index")?;
        let end_idx = checked_usize(n_end.0, "transcript end index")?;

        let ref_seq = self
            .refs
            .reference(v.ac(), IdentifierType::TranscriptAccession)
            .whole()?;
        if start_idx >= ref_seq.len() || end_idx > ref_seq.len() {
            return Err(HgvsError::ValidationError(
                "Transcript sequence too short".into(),
            ));
        }
        Ok(stated_ref_matches(
            &v.posedit().edit,
            &ref_seq[start_idx..end_idx],
        ))
    }

    pub fn normalize_variant(
        &self,
        var: crate::SequenceVariant,
    ) -> Result<crate::SequenceVariant, HgvsError> {
        use crate::SequenceVariant as SV;
        Ok(match var {
            SV::Genomic(v) => SV::Genomic(self.normalize_linear(v)?),
            SV::Mitochondrial(v) => SV::Mitochondrial(self.normalize_linear(v)?),
            SV::Coding(v) => SV::Coding(self.normalize_transcript(v)?),
            SV::NonCoding(v) => SV::NonCoding(self.normalize_transcript(v)?),
            // r. normalises as the c. or n. variant it is spelled from, and
            // comes back in RNA letters. A statement about the transcript
            // (r.0, r.spl) has nothing to normalise.
            SV::Rna(r) if matches!(r.posedit.edit, crate::edits::NaEdit::Special { .. }) => {
                SV::Rna(r)
            }
            SV::Rna(r) => {
                let normalised = self.normalize_variant(self.r_as_transcript(&r)?)?;
                SV::Rna(self.tx_to_r(&normalised)?)
            }
            // Each member normalises on its own; the allele keeps their order.
            SV::CisPhased(cis) => {
                let members = cis
                    .members
                    .into_iter()
                    .map(|m| self.normalize_variant(m))
                    .collect::<Result<Vec<_>, _>>()?;
                SV::CisPhased(crate::structs::CisPhasedVariant {
                    ac: cis.ac,
                    gene: cis.gene,
                    members,
                })
            }
            other => other,
        })
    }

    fn normalize_linear<L: LinearVariant>(&self, mut v: L) -> Result<L, HgvsError> {
        let ac = v.ac().to_string();
        let Some(pos) = &mut v.posedit_mut().pos else {
            return Ok(v);
        };
        let (start, end) = simple_interval_range(pos)?;
        let before = PlacedEdit::from_hgvs_range(start, end, v.posedit().edit.clone());
        let reference = self.refs.reference(&ac, IdentifierType::GenomicAccession);
        let after = normalize::normalize(&reference, before.clone())?;
        let posedit = v.posedit_mut();
        let pos = posedit.pos.as_mut().expect("checked above");
        if placement_changed(&before, &after) {
            let (s, e) = hgvs_positions(&before, &after, pos.end.is_some());
            pos.start.base = GenomicPos(s as i32).to_hgvs();
            pos.end = e.map(|last| SimplePosition {
                base: GenomicPos(last as i32).to_hgvs(),
                end: None,
                uncertain: false,
            });
        }
        posedit.edit = after.edit;
        Ok(v)
    }

    fn normalize_transcript<V: TranscriptVariant>(&self, mut v: V) -> Result<V, HgvsError> {
        let ac = v.ac().to_string();
        let am = TranscriptMapper::new(self.provider().get_transcript(&ac, None)?)?;
        let Some(pos) = &v.posedit().pos else {
            return Ok(v);
        };
        if has_intronic_offset(pos) {
            // An intronic base has no transcript index; there is nothing to
            // normalise against in transcript space. Leave the variant as written.
            return Ok(v);
        }
        let (start, end) = self.get_c_indices(pos, &am)?;
        let before = PlacedEdit::from_hgvs_range(start, end, v.posedit().edit.clone());
        let reference = self
            .refs
            .reference(&ac, IdentifierType::TranscriptAccession);
        let after = normalize::normalize(&reference, before.clone())?;
        let posedit = v.posedit_mut();
        let pos = posedit.pos.as_mut().expect("checked above");
        if placement_changed(&before, &after) {
            // Re-derive positions from indices so that, for c., the c.0 gap and
            // the CDS anchors come out right when a shift crosses the CDS bounds.
            let (s, e) = hgvs_positions(&before, &after, pos.end.is_some());
            pos.start = V::position_from_index(&am, s as i32)?;
            pos.end = e
                .map(|last| V::position_from_index(&am, last as i32))
                .transpose()?;
        }
        posedit.edit = after.edit;
        Ok(v)
    }

    /// The 0-based half-open transcript index range a c./n. interval names.
    pub fn get_c_indices(
        &self,
        pos: &BaseOffsetInterval,
        am: &TranscriptMapper,
    ) -> Result<(usize, usize), HgvsError> {
        let (n_start, n_end) = am.interval_to_n(pos)?;
        if n_start.0 < 0 {
            return Err(HgvsError::ValidationError(format!(
                "Transcript start position {:?} maps to a negative index; \
                 normalization of 5'UTR-spanning variants is not yet supported",
                pos.start.base
            )));
        }
        Ok((
            checked_usize(n_start.0, "transcript start index")?,
            checked_usize(n_end.0, "transcript end index")?,
        ))
    }

    /// The canonical allele of a nucleotide variant on its genomic reference:
    /// the one unambiguous statement of the change that SPDI, VRS and
    /// equivalence all derive from. See [`CanonicalAllele`].
    pub fn canonical_allele(
        &self,
        var: &crate::SequenceVariant,
    ) -> Result<CanonicalAllele, HgvsError> {
        if let crate::SequenceVariant::CisPhased(cis) = var {
            return Err(HgvsError::UnsupportedOperation(format!(
                "A cis allele has one canonical allele per member, not one of its own; \
                 take canonical_allele of each of the {} members of {var}, or to_vrs_variation",
                cis.members.len()
            )));
        }
        if let crate::SequenceVariant::Protein(vp) = var {
            let pos = vp
                .posedit
                .pos
                .as_ref()
                .ok_or_else(|| HgvsError::ValidationError("Missing position".into()))?;
            let (start, end) = aa_interval_range(pos)?;
            let reference = self
                .refs
                .reference(&vp.ac, IdentifierType::ProteinAccession);
            let mut resolved = vp.posedit.edit.resolve(&reference, start, end)?;
            // A stop among the new residues ends the protein there: everything
            // from it to the end of the reference goes too. p.Tyr165Ter and
            // p.Ala164_Tyr165insTer are then the same allele.
            if let Some(k) = resolved.alt.find('*') {
                let len = reference.whole()?.trim_end_matches('*').len();
                resolved.alt.truncate(k);
                resolved.end = len.max(resolved.start);
                resolved.ref_ = reference.slice(resolved.start, resolved.end)?;
            }
            return CanonicalAllele::canonicalize(&reference, &vp.ac, &resolved);
        }
        let g = self.as_genomic(var).ok_or_else(|| {
            HgvsError::UnsupportedOperation(
                "Canonical alleles exist for genomic, mitochondrial, coding, non-coding, RNA and protein variants only".into(),
            )
        })??;
        let pos = g
            .posedit
            .pos
            .as_ref()
            .ok_or_else(|| HgvsError::ValidationError("Missing position".into()))?;
        let edit = &g.posedit.edit;
        if matches!(
            edit,
            crate::edits::NaEdit::None
                | crate::edits::NaEdit::Con { .. }
                | crate::edits::NaEdit::NACopy { .. }
                | crate::edits::NaEdit::Special { .. }
                | crate::edits::NaEdit::InsLength { .. }
                | crate::edits::NaEdit::DelInsLength { .. }
        ) {
            return Err(HgvsError::UnsupportedOperation(format!(
                "Edit type {:?} has no canonical allele",
                edit
            )));
        }
        let (hgvs_start, hgvs_end) = simple_interval_range(pos)?;
        let reference = self.refs.reference(&g.ac, IdentifierType::GenomicAccession);
        let resolved = edit.resolve(&reference, hgvs_start, hgvs_end)?;
        CanonicalAllele::canonicalize(&reference, &g.ac, &resolved)
    }

    /// The unambiguous SPDI of a variant: its canonical allele, rendered.
    pub fn to_spdi_unambiguous(&self, var: &crate::SequenceVariant) -> Result<String, HgvsError> {
        Ok(self.canonical_allele(var)?.spdi())
    }

    /// The GA4GH VRS 2.0 Allele of a variant, with computed identifiers. The
    /// input HGVS is carried as an expression.
    ///
    /// An insertion of bases known only by number, `insN[20]` or
    /// `delinsN[20]`, has a `LengthExpression` state over the insertion point
    /// or the deleted range, unnormalised; a deletion with uncertain
    /// breakpoints has Range bounds. A copy-number edit or an imprecise
    /// duplication has no Allele: see [`to_vrs_variation`](Self::to_vrs_variation).
    pub fn to_vrs(&self, var: &crate::SequenceVariant) -> Result<VrsAllele, HgvsError> {
        if let crate::SequenceVariant::CisPhased(_) = var {
            return Err(HgvsError::UnsupportedOperation(format!(
                "{var} is an allele in cis, which is a VRS CisPhasedBlock, not an Allele; \
                 use to_vrs_variation"
            )));
        }
        let syntax = format!("hgvs.{}", var.coordinate_type());
        let hgvs = var.to_string();
        // Breakpoints known only to ranges cannot be normalised; VRS carries
        // them as Range bounds, for deletions.
        let linear = match var {
            crate::SequenceVariant::Genomic(v) => Some(v.clone()),
            crate::SequenceVariant::Mitochondrial(v) => Some(v.to_genomic()),
            _ => None,
        };
        if let Some(g) = linear {
            if let Some((start, end)) = g.posedit.pos.as_ref().and_then(uncertain_bounds) {
                if matches!(g.posedit.edit, crate::edits::NaEdit::Dup { .. }) {
                    return Err(HgvsError::UnsupportedOperation(
                        "A duplication with uncertain breakpoints has no Allele; \
                         to_vrs_copy_number_change or to_vrs_variation renders it as a \
                         CopyNumberChange"
                            .into(),
                    ));
                }
                if !matches!(g.posedit.edit, crate::edits::NaEdit::Del { .. }) {
                    return Err(HgvsError::UnsupportedOperation(format!(
                        "Only a deletion can have uncertain breakpoints in an Allele, not {:?}",
                        g.posedit.edit
                    )));
                }
                let refget = self
                    .refs
                    .reference(&g.ac, IdentifierType::GenomicAccession)
                    .refget_accession()?;
                return Ok(VrsAllele::imprecise_deletion(
                    &refget,
                    start,
                    end,
                    VrsMolecule::Genomic,
                    Some((&syntax, &hgvs)),
                ));
            }
        }
        if let Some(length) = na_edit(var).and_then(inserted_length) {
            return self.length_allele(var, length, &syntax, &hgvs);
        }
        let allele = self.canonical_allele(var)?;
        let (kind, molecule) = match var {
            crate::SequenceVariant::Protein(_) => {
                (IdentifierType::ProteinAccession, VrsMolecule::Protein)
            }
            _ => (IdentifierType::GenomicAccession, VrsMolecule::Genomic),
        };
        let refget = self
            .refs
            .reference(&allele.accession, kind)
            .refget_accession()?;
        Ok(VrsAllele::new(
            &allele,
            &refget,
            molecule,
            Some((&syntax, &hgvs)),
        ))
    }

    /// The Allele of an insertion of unspecified bases on the genomic
    /// reference: for `insN[20]` the insertion point, the empty interbase
    /// range in front of the second flanking base (where `NaEdit::Ins` is
    /// anchored too); for `delinsN[20]` the deleted range. The state is a
    /// `LengthExpression`. Unknown bases cannot slide, so nothing is
    /// normalised.
    fn length_allele(
        &self,
        var: &crate::SequenceVariant,
        length: VrsBound,
        syntax: &str,
        hgvs: &str,
    ) -> Result<VrsAllele, HgvsError> {
        let g = self.as_genomic(var).ok_or_else(|| {
            HgvsError::UnsupportedOperation(
                "An insertion of a stated length is a nucleotide variant".into(),
            )
        })??;
        let pos = g
            .posedit
            .pos
            .as_ref()
            .ok_or_else(|| HgvsError::ValidationError("Missing position".into()))?;
        let (start, end) = simple_interval_range(pos)?;
        let (start, end) = if matches!(g.posedit.edit, crate::edits::NaEdit::InsLength { .. }) {
            let anchor = if end > start { end - 1 } else { start };
            (anchor, anchor)
        } else {
            (start, end)
        };
        let refget = self
            .refs
            .reference(&g.ac, IdentifierType::GenomicAccession)
            .refget_accession()?;
        Ok(VrsAllele::length_expression(
            &refget,
            VrsBound::Exact(start),
            VrsBound::Exact(end),
            length,
            VrsMolecule::Genomic,
            Some((syntax, hgvs)),
        ))
    }

    /// The GA4GH VRS 2.0 `CopyNumberCount` of a `g.` or `m.` copy-number
    /// variant, `g.1000_2000copy3`: the location over the range, as interbase
    /// `[start-1, end)` or as Range bounds when the breakpoints are uncertain,
    /// and the count. Nothing is normalised: a count has no placement to
    /// shift. The input HGVS is carried as an expression.
    pub fn to_vrs_copy_number(
        &self,
        var: &crate::SequenceVariant,
    ) -> Result<VrsCopyNumberCount, HgvsError> {
        let g = match var {
            crate::SequenceVariant::Genomic(v) => v.clone(),
            crate::SequenceVariant::Mitochondrial(v) => v.to_genomic(),
            _ => {
                return Err(HgvsError::UnsupportedOperation(
                    "Copy number counts exist for genomic and mitochondrial variants only".into(),
                ))
            }
        };
        let crate::edits::NaEdit::NACopy { copy, .. } = g.posedit.edit else {
            return Err(HgvsError::UnsupportedOperation(format!(
                "Edit type {:?} is not a copy number change",
                g.posedit.edit
            )));
        };
        let copies = usize::try_from(copy)
            .map_err(|_| HgvsError::ValidationError(format!("Copy number {copy} is negative")))?;
        let pos = g
            .posedit
            .pos
            .as_ref()
            .ok_or_else(|| HgvsError::ValidationError("Missing position".into()))?;
        let (start, end) = interval_bounds(pos)?;
        let refget = self
            .refs
            .reference(&g.ac, IdentifierType::GenomicAccession)
            .refget_accession()?;
        Ok(VrsCopyNumberCount::new(
            &refget,
            start,
            end,
            VrsBound::Exact(copies),
            VrsMolecule::Genomic,
            Some((&format!("hgvs.{}", var.coordinate_type()), &var.to_string())),
        ))
    }

    /// The GA4GH VRS 2.0 `CopyNumberChange` of a `g.` or `m.` duplication or
    /// deletion: the location over the range, as interbase `[start-1, end)`
    /// or as Range bounds when the breakpoints are uncertain, and the
    /// direction of the change as an EFO term, "gain" (EFO:0030070, copy
    /// number gain) for a duplication and "loss" (EFO:0030067, copy number
    /// loss) for a deletion. Nothing is normalised: a change of copies has no
    /// placement to shift. The input HGVS is carried as an expression.
    pub fn to_vrs_copy_number_change(
        &self,
        var: &crate::SequenceVariant,
    ) -> Result<VrsCopyNumberChange, HgvsError> {
        let g = match var {
            crate::SequenceVariant::Genomic(v) => v.clone(),
            crate::SequenceVariant::Mitochondrial(v) => v.to_genomic(),
            _ => {
                return Err(HgvsError::UnsupportedOperation(
                    "Copy number changes exist for genomic and mitochondrial variants only".into(),
                ))
            }
        };
        let copy_change = match &g.posedit.edit {
            crate::edits::NaEdit::Dup { .. } => VrsCopyChange::Gain,
            crate::edits::NaEdit::Del { .. } => VrsCopyChange::Loss,
            other => {
                return Err(HgvsError::UnsupportedOperation(format!(
                    "Edit type {other:?} is neither a duplication nor a deletion"
                )))
            }
        };
        let pos = g
            .posedit
            .pos
            .as_ref()
            .ok_or_else(|| HgvsError::ValidationError("Missing position".into()))?;
        let (start, end) = interval_bounds(pos)?;
        let refget = self
            .refs
            .reference(&g.ac, IdentifierType::GenomicAccession)
            .refget_accession()?;
        Ok(VrsCopyNumberChange::new(
            &refget,
            start,
            end,
            copy_change,
            VrsMolecule::Genomic,
            Some((&format!("hgvs.{}", var.coordinate_type()), &var.to_string())),
        ))
    }

    /// The GA4GH VRS 2.0 `CisPhasedBlock` of a cis allele, `c.[145C>T;147C>G]`:
    /// the `Allele` of each member (`to_vrs`, so each is canonicalised on its
    /// own), which must all lie on one sequence, that sequence as
    /// `sequenceReference`, and the input HGVS as an expression. The
    /// identifier does not depend on the order of the members.
    pub fn to_vrs_cis_phased(
        &self,
        var: &crate::SequenceVariant,
    ) -> Result<VrsCisPhasedBlock, HgvsError> {
        let crate::SequenceVariant::CisPhased(cis) = var else {
            return Err(HgvsError::UnsupportedOperation(format!(
                "Only an allele in cis, ac:c.[a;b], renders as a CisPhasedBlock, not {var}"
            )));
        };
        let members = cis
            .members
            .iter()
            .map(|m| self.to_vrs(m))
            .collect::<Result<Vec<_>, _>>()?;
        let Some(first) = members.first() else {
            return Err(HgvsError::ValidationError(
                "A cis allele needs at least one member".into(),
            ));
        };
        let reference = first.location.sequence_reference.clone();
        if let Some(other) = members
            .iter()
            .find(|m| m.location.sequence_reference.refget_accession != reference.refget_accession)
        {
            return Err(HgvsError::ValidationError(format!(
                "The members of {var} lie on different sequences, {} and {}",
                reference.refget_accession, other.location.sequence_reference.refget_accession
            )));
        }
        Ok(VrsCisPhasedBlock::new(
            members,
            Some(reference),
            Some((&format!("hgvs.{}", var.coordinate_type()), &var.to_string())),
        ))
    }

    /// The VRS 2.0 object of a variant: the `CisPhasedBlock` of a cis allele,
    /// the `CopyNumberCount` of a copy-number edit, the `CopyNumberChange` of a
    /// duplication with uncertain breakpoints (which has no Allele), else the
    /// `Allele` of `to_vrs` (an imprecise deletion included).
    pub fn to_vrs_variation(
        &self,
        var: &crate::SequenceVariant,
    ) -> Result<VrsVariation, HgvsError> {
        if let crate::SequenceVariant::CisPhased(_) = var {
            return Ok(VrsVariation::CisPhasedBlock(self.to_vrs_cis_phased(var)?));
        }
        Ok(if is_copy_number(var) {
            VrsVariation::CopyNumberCount(self.to_vrs_copy_number(var)?)
        } else if is_imprecise_duplication(var) {
            VrsVariation::CopyNumberChange(self.to_vrs_copy_number_change(var)?)
        } else {
            VrsVariation::Allele(self.to_vrs(var)?)
        })
    }

    /// The variant an SPDI string names (`accession:position:deletion:insertion`,
    /// interbase position, the deletion as bases or as a length), written in
    /// HGVS on its own sequence and 3'-normalised: `g.` for a nucleotide
    /// accession, `p.` for a protein.
    pub fn from_spdi(&self, spdi: &str) -> Result<crate::SequenceVariant, HgvsError> {
        let parts: Vec<&str> = spdi.trim().split(':').collect();
        let [ac, pos, del, ins] = parts[..] else {
            return Err(HgvsError::ValidationError(format!(
                "SPDI is accession:position:deletion:insertion, not {spdi:?}"
            )));
        };
        let start: usize = pos.parse().map_err(|_| {
            HgvsError::ValidationError(format!("SPDI position {pos:?} is not a number"))
        })?;
        let kind = self.sequence_kind(ac)?;
        let end = if !del.is_empty() && del.chars().all(|c| c.is_ascii_digit()) {
            start + del.parse::<usize>().unwrap_or(0)
        } else {
            let actual = self
                .refs
                .reference(ac, kind)
                .slice(start, start + del.len())?;
            if actual != del {
                return Err(HgvsError::ValidationError(format!(
                    "SPDI deletion {del:?} is not what {ac} holds at {start}: {actual:?}"
                )));
            }
            start + del.len()
        };
        self.allele_to_variant(ac, kind, start, end, ins.to_string())
    }

    /// The variant a GA4GH VRS 2.0 object (as JSON) names, written in HGVS on
    /// its own sequence. An `Allele` comes back 3'-normalised, with Range
    /// bounds accepted for a deletion, `g.(a_b)_(c_d)del`, and a
    /// `LengthExpression` state as `g.<a>_<a+1>insN[n]` or
    /// `g.<start+1>_<end>delinsN[n]` (`N[(min_max)]` for a range of lengths);
    /// a `CopyNumberCount` comes back as `g.<start+1>_<end>copyN`, and a
    /// `CopyNumberChange` as `g.<start+1>_<end>dup` for a gain or `del` for a
    /// loss (Range bounds as `(a_b)_(c_d)`); a `CisPhasedBlock` comes back as the
    /// cis allele of its members, `g.[a;b]` (or `p.`). The sequence is identified by
    /// its refget accession: `accession` names it when given, else the
    /// mapper's `Refget` lookup must; the digest is checked against the
    /// sequence either way.
    pub fn from_vrs(
        &self,
        json: &str,
        accession: Option<&str>,
    ) -> Result<crate::SequenceVariant, HgvsError> {
        if vrs_type(json)? == "CopyNumberCount" {
            return self.copy_number_from_vrs(json, accession);
        }
        if vrs_type(json)? == "CopyNumberChange" {
            return self.copy_number_change_from_vrs(json, accession);
        }
        if vrs_type(json)? == "CisPhasedBlock" {
            return self.cis_phased_from_vrs(json, accession);
        }
        let allele = VrsAllele::from_json(json)?;
        let (ac, kind) = self.located_sequence(&allele.location, accession)?;
        let reference = self.refs.reference(&ac, kind);
        match (allele.location.start, allele.location.end) {
            (VrsBound::Exact(start), VrsBound::Exact(end)) => {
                if end < start {
                    return Err(HgvsError::ValidationError(format!(
                        "Location end {end} is before start {start}"
                    )));
                }
                let alt = match &allele.state {
                    VrsState::Literal { sequence, .. } => sequence.clone(),
                    VrsState::Length { length, .. } => {
                        return self.length_variant(&ac, kind, start, end, *length);
                    }
                    VrsState::ReferenceLength {
                        length,
                        repeat_subunit_length,
                        ..
                    } => {
                        let r = reference.slice(start, end)?;
                        let unit = *repeat_subunit_length;
                        if unit == 0 || unit > r.len() {
                            return Err(HgvsError::ValidationError(format!(
                                "repeatSubunitLength {unit} does not fit a location of {} bases",
                                r.len()
                            )));
                        }
                        r[..unit].chars().cycle().take(*length).collect()
                    }
                };
                self.allele_to_variant(&ac, kind, start, end, alt)
            }
            (start, end) => {
                let empty = matches!(&allele.state, VrsState::Literal { sequence, .. } if sequence.is_empty());
                if kind == IdentifierType::ProteinAccession || !empty {
                    return Err(HgvsError::UnsupportedOperation(
                        "Range bounds are supported for a deletion on a nucleotide sequence only"
                            .into(),
                    ));
                }
                let del = crate::edits::NaEdit::Del {
                    ref_: None,
                    uncertain: false,
                };
                Ok(crate::SequenceVariant::Genomic(g_variant(
                    &ac,
                    bounded_interval(start, end),
                    del,
                )))
            }
        }
    }

    /// The HGVS variant for "`length` unspecified bases replace interbase
    /// `[start, end)` of `ac`": `g.<start>_<start+1>insN[n]` when the range
    /// is empty, `g.<start+1>_<end>delinsN[n]` otherwise, as written (unknown
    /// bases cannot be normalised).
    fn length_variant(
        &self,
        ac: &str,
        kind: IdentifierType,
        start: usize,
        end: usize,
        length: VrsBound,
    ) -> Result<crate::SequenceVariant, HgvsError> {
        use crate::edits::NaEdit;
        if kind == IdentifierType::ProteinAccession {
            return Err(HgvsError::UnsupportedOperation(
                "A LengthExpression is read on a nucleotide sequence only".into(),
            ));
        }
        let (min, max) = match length {
            VrsBound::Exact(n) => (n, n),
            VrsBound::Range(Some(lo), Some(hi)) if lo <= hi => (lo, hi),
            _ => return Err(HgvsError::UnsupportedOperation(
                "HGVS states an insertion's length as a number or a closed range, insN[(20_30)]"
                    .into(),
            )),
        };
        let uncertain = false;
        if start == end {
            if start == 0 {
                return Err(HgvsError::UnsupportedOperation(
                    "HGVS has no insertion before the first base".into(),
                ));
            }
            // The insertion point must have a base on each side.
            if self
                .refs
                .reference(ac, kind)
                .slice(start, start + 1)?
                .is_empty()
            {
                return Err(HgvsError::ValidationError(format!(
                    "{ac} is shorter than position {}",
                    start + 1
                )));
            }
            let edit = NaEdit::InsLength {
                min,
                max,
                uncertain,
            };
            return Ok(crate::SequenceVariant::Genomic(g_variant(
                ac,
                interval(start - 1, start + 1),
                edit,
            )));
        }
        if self.refs.reference(ac, kind).slice(start, end)?.len() != end - start {
            return Err(HgvsError::ValidationError(format!(
                "{ac} is shorter than position {end}"
            )));
        }
        let edit = NaEdit::DelInsLength {
            ref_: None,
            min,
            max,
            uncertain,
        };
        Ok(crate::SequenceVariant::Genomic(g_variant(
            ac,
            interval(start, end),
            edit,
        )))
    }

    /// `ac:g.[a;b]` (or `p.`) from a VRS `CisPhasedBlock`: each member read
    /// back as the `Allele` it is, on the sequence `accession` or the first
    /// member's refget accession names; every member must lie on it.
    fn cis_phased_from_vrs(
        &self,
        json: &str,
        accession: Option<&str>,
    ) -> Result<crate::SequenceVariant, HgvsError> {
        let block = VrsCisPhasedBlock::from_json(json)?;
        if let (Some(reference), Some(first)) = (&block.sequence_reference, block.members.first()) {
            let stated = &first.location.sequence_reference.refget_accession;
            if reference.refget_accession != *stated {
                return Err(HgvsError::ValidationError(format!(
                    "The CisPhasedBlock is on {} but its members on {stated}",
                    reference.refget_accession
                )));
            }
        }
        let mut ac = accession.map(str::to_string);
        let mut members = Vec::with_capacity(block.members.len());
        for member in &block.members {
            let var = self.from_vrs(&member.to_json(), ac.as_deref())?;
            ac.get_or_insert_with(|| var.ac().to_string());
            members.push(var);
        }
        let ac = ac.ok_or_else(|| {
            HgvsError::ValidationError("A CisPhasedBlock has at least one member".into())
        })?;
        Ok(crate::SequenceVariant::CisPhased(
            crate::structs::CisPhasedVariant::new(ac, None, members)?,
        ))
    }

    /// `g.<start+1>_<end>copyN` from a VRS `CopyNumberCount`. HGVS has no
    /// syntax for a range of counts, so `copies` must be exact.
    fn copy_number_from_vrs(
        &self,
        json: &str,
        accession: Option<&str>,
    ) -> Result<crate::SequenceVariant, HgvsError> {
        let count = VrsCopyNumberCount::from_json(json)?;
        let (ac, kind) = self.located_sequence(&count.location, accession)?;
        if kind == IdentifierType::ProteinAccession {
            return Err(HgvsError::UnsupportedOperation(
                "A copy number count is of a nucleotide sequence, not a protein".into(),
            ));
        }
        let VrsBound::Exact(copies) = count.copies else {
            return Err(HgvsError::UnsupportedOperation(
                "HGVS writes an exact copy number, not a range of counts".into(),
            ));
        };
        let copy = i32::try_from(copies).map_err(|_| {
            HgvsError::ValidationError(format!("Copy number {copies} is too large"))
        })?;
        let pos = location_interval(count.location.start, count.location.end)?;
        let edit = crate::edits::NaEdit::NACopy {
            copy,
            uncertain: false,
        };
        Ok(crate::SequenceVariant::Genomic(g_variant(&ac, pos, edit)))
    }

    /// `g.<start+1>_<end>dup` from a VRS `CopyNumberChange` that is a gain
    /// (EFO:0030070, 71, 72), `del` from a loss (EFO:0030067, 68, 69,
    /// EFO:0020073). Regional base ploidy and terms this module does not
    /// know have no HGVS form.
    fn copy_number_change_from_vrs(
        &self,
        json: &str,
        accession: Option<&str>,
    ) -> Result<crate::SequenceVariant, HgvsError> {
        let change = VrsCopyNumberChange::from_json(json)?;
        let (ac, kind) = self.located_sequence(&change.location, accession)?;
        if kind == IdentifierType::ProteinAccession {
            return Err(HgvsError::UnsupportedOperation(
                "A copy number change is of a nucleotide sequence, not a protein".into(),
            ));
        }
        let uncertain = false;
        let edit = if change.copy_change.is_gain() {
            crate::edits::NaEdit::Dup {
                ref_: None,
                uncertain,
            }
        } else if change.copy_change.is_loss() {
            crate::edits::NaEdit::Del {
                ref_: None,
                uncertain,
            }
        } else {
            return Err(HgvsError::UnsupportedOperation(format!(
                "copyChange {:?} is neither a gain nor a loss, so has no HGVS form",
                change.copy_change.label()
            )));
        };
        let pos = location_interval(change.location.start, change.location.end)?;
        Ok(crate::SequenceVariant::Genomic(g_variant(&ac, pos, edit)))
    }

    /// The accession and kind of the sequence a VRS location is on: named by
    /// `accession` when given, else looked up from the refget accession; the
    /// digest is checked against the sequence either way. The kind is read
    /// from the sequence reference when it says, else asked of the provider.
    fn located_sequence(
        &self,
        location: &VrsSequenceLocation,
        accession: Option<&str>,
    ) -> Result<(String, IdentifierType), HgvsError> {
        let sr = &location.sequence_reference;
        let refget = sr.refget_accession.as_str();
        let ac = match accession {
            Some(a) => a.to_string(),
            None => self.refs.accession_for_refget(refget)?.ok_or_else(|| {
                HgvsError::DataProviderError(format!(
                    "No accession is known for {refget}; pass one, or give the mapper a Refget lookup"
                ))
            })?,
        };
        let kind = if sr.residue_alphabet == "aa" || sr.molecule_type == "protein" {
            IdentifierType::ProteinAccession
        } else if sr.residue_alphabet == "na" {
            IdentifierType::GenomicAccession
        } else {
            self.sequence_kind(&ac)?
        };
        let actual = self.refs.reference(&ac, kind).refget_accession()?;
        if actual != refget {
            return Err(HgvsError::ValidationError(format!(
                "{refget} is not the refget accession of {ac}, which is {actual}"
            )));
        }
        Ok((ac, kind))
    }

    /// Whether `ac` names a protein or a nucleotide sequence.
    fn sequence_kind(&self, ac: &str) -> Result<IdentifierType, HgvsError> {
        Ok(match self.provider().get_identifier_type(ac)? {
            IdentifierType::ProteinAccession => IdentifierType::ProteinAccession,
            _ => IdentifierType::GenomicAccession,
        })
    }

    /// The HGVS variant for "the bases over `[start, end)` of `ac` become
    /// `alt`": trimmed to the change, 3'-normalised, written as `g.` or `p.`.
    fn allele_to_variant(
        &self,
        ac: &str,
        kind: IdentifierType,
        start: usize,
        end: usize,
        alt: String,
    ) -> Result<crate::SequenceVariant, HgvsError> {
        let reference = self.refs.reference(ac, kind);
        let ref_ = reference.slice(start, end)?;
        if ref_.len() != end - start {
            return Err(HgvsError::ValidationError(format!(
                "{ac} is shorter than position {end}"
            )));
        }
        let (edit, s, e) = hgvs_edit_for(&reference, kind, start, end, &ref_, &alt)?;
        let placed = normalize::normalize(&reference, PlacedEdit::from_hgvs_range(s, e, edit))?;
        let (s, e) = placed.hgvs_range();
        Ok(match kind {
            IdentifierType::ProteinAccession => {
                crate::SequenceVariant::Protein(protein_variant(ac, &reference, s, e, placed.edit)?)
            }
            _ => crate::SequenceVariant::Genomic(g_variant(ac, interval(s, e), placed.edit)),
        })
    }

    pub fn to_spdi(
        &self,
        var: &crate::SequenceVariant,
        unambiguous: bool,
    ) -> Result<String, HgvsError> {
        // A protein has only the canonical form.
        if unambiguous || matches!(var, crate::SequenceVariant::Protein(_)) {
            self.to_spdi_unambiguous(var)
        } else {
            // 1. Resolve to genomic if possible.
            let g_var_obj = self.as_genomic(var).ok_or_else(|| {
                HgvsError::UnsupportedOperation("SPDI only for genomic/coding/non-coding".into())
            })??;

            // 2. Normalize (3' shift, minimal delins)
            let g_norm_var = self.normalize_variant(crate::SequenceVariant::Genomic(g_var_obj))?;
            let g_norm = match g_norm_var {
                crate::SequenceVariant::Genomic(v) => v,
                _ => unreachable!(),
            };
            g_norm.posedit.to_spdi(&g_norm.ac, &self.refs)
        }
    }

    /// The variant as a `g.` variant on its reference: genomic and mitochondrial
    /// as written, transcript-space variants mapped through their transcript.
    /// `None` for protein and RNA variants.
    pub fn as_genomic(&self, var: &crate::SequenceVariant) -> Option<Result<GVariant, HgvsError>> {
        use crate::SequenceVariant as SV;
        Some(match var {
            SV::Genomic(v) => Ok(v.clone()),
            SV::Mitochondrial(v) => Ok(v.to_genomic()),
            SV::Coding(v) => self.tx_to_g(v, None),
            SV::NonCoding(v) => self.tx_to_g(v, None),
            SV::Rna(r) => self.r_to_g(r, None),
            SV::Protein(_) | SV::CisPhased(_) => return None,
        })
    }

    // --- r.: the transcript in RNA letters ---
    //
    // HGVS numbers r. positions like c. on a coding transcript and like n. on
    // a non-coding one, and writes bases in lowercase with u for T. Every r.
    // operation is therefore a conversion to the c. or n. spelling and back.

    /// Whether `ac` has a CDS, which decides whether r. is numbered like c. or n.
    fn has_cds(&self, ac: &str) -> Result<bool, HgvsError> {
        let t = self.provider().get_transcript(ac, None)?;
        Ok(t.cds_start_index.is_some() && t.cds_end_index.is_some())
    }

    /// An r. variant as the c. or n. variant it is spelled from.
    pub fn r_as_transcript(&self, r: &RVariant) -> Result<crate::SequenceVariant, HgvsError> {
        Ok(if self.has_cds(&r.ac)? {
            crate::SequenceVariant::Coding(self.r_to_c(r)?)
        } else {
            crate::SequenceVariant::NonCoding(self.r_to_n(r)?)
        })
    }

    /// A c. or n. variant in its r. spelling.
    pub fn tx_to_r(&self, var: &crate::SequenceVariant) -> Result<RVariant, HgvsError> {
        match var {
            crate::SequenceVariant::Coding(c) => self.c_to_r(c),
            crate::SequenceVariant::NonCoding(n) => self.n_to_r(n),
            other => Err(HgvsError::UnsupportedOperation(format!(
                "Only c. and n. variants have an r. spelling, not {other}"
            ))),
        }
    }

    /// r. to c.: the same positions, numbered from the CDS start, in DNA letters.
    pub fn r_to_c(&self, r: &RVariant) -> Result<CVariant, HgvsError> {
        if !self.has_cds(&r.ac)? {
            return Err(HgvsError::UnsupportedOperation(format!(
                "{} has no CDS, so its r. positions are n. positions; use r_to_n",
                r.ac
            )));
        }
        let posedit = relettered_posedit(r, &r.posedit, Letters::Dna, |anchor| match anchor {
            Anchor::TranscriptStart => Anchor::CdsStart,
            other => other,
        })?;
        Ok(CVariant::from_parts(r.ac.clone(), r.gene.clone(), posedit))
    }

    /// r. to n. on a non-coding transcript: the same positions in DNA letters.
    pub fn r_to_n(&self, r: &RVariant) -> Result<NVariant, HgvsError> {
        if self.has_cds(&r.ac)? {
            return Err(HgvsError::UnsupportedOperation(format!(
                "{} has a CDS, so its r. positions are c. positions; use r_to_c",
                r.ac
            )));
        }
        let posedit = relettered_posedit(r, &r.posedit, Letters::Dna, |anchor| anchor)?;
        if let Some(pos) = &posedit.pos {
            let cds_anchored = |p: &BaseOffsetPosition| p.anchor == Anchor::CdsEnd || p.base.0 < 1;
            if cds_anchored(&pos.start) || pos.end.as_ref().is_some_and(cds_anchored) {
                return Err(HgvsError::ValidationError(format!(
                    "{r} uses CDS-relative positions on a transcript without a CDS"
                )));
            }
        }
        Ok(NVariant::from_parts(r.ac.clone(), r.gene.clone(), posedit))
    }

    /// c. to r.: the same positions in RNA letters (lowercase, u for T).
    pub fn c_to_r(&self, c: &CVariant) -> Result<RVariant, HgvsError> {
        let posedit = relettered_posedit(c, &c.posedit, Letters::Rna, |anchor| match anchor {
            Anchor::CdsStart => Anchor::TranscriptStart,
            other => other,
        })?;
        Ok(RVariant {
            ac: c.ac.clone(),
            gene: c.gene.clone(),
            posedit,
        })
    }

    /// n. to r.: the same positions in RNA letters.
    pub fn n_to_r(&self, n: &NVariant) -> Result<RVariant, HgvsError> {
        let posedit = relettered_posedit(n, &n.posedit, Letters::Rna, |anchor| anchor)?;
        Ok(RVariant {
            ac: n.ac.clone(),
            gene: n.gene.clone(),
            posedit,
        })
    }

    /// r. to g., for a change within one exon. A change spanning a splice
    /// junction describes the spliced RNA and has no single genomic form.
    pub fn r_to_g(&self, r: &RVariant, reference_ac: Option<&str>) -> Result<GVariant, HgvsError> {
        let tx = self.r_as_transcript(r)?;
        let posedit = match &tx {
            crate::SequenceVariant::Coding(c) => &c.posedit,
            crate::SequenceVariant::NonCoding(n) => &n.posedit,
            _ => unreachable!("r_as_transcript gives c. or n."),
        };
        if let Some(pos) = &posedit.pos {
            let am = TranscriptMapper::new(self.provider().get_transcript(&r.ac, None)?)?;
            // Intronic positions do not resolve to a transcript range; they
            // name the genome directly and need no guard.
            if let Ok((start, end)) = am.interval_to_n(pos) {
                let within_one_exon = am
                    .exons
                    .iter()
                    .any(|x| x.transcript_start.0 <= start.0 && end.0 <= x.transcript_end.0);
                if !within_one_exon {
                    return Err(HgvsError::UnsupportedOperation(format!(
                        "{r} spans a splice junction: the spliced RNA has no single genomic equivalent"
                    )));
                }
            }
        }
        match tx {
            crate::SequenceVariant::Coding(c) => self.tx_to_g(&c, reference_ac),
            crate::SequenceVariant::NonCoding(n) => self.tx_to_g(&n, reference_ac),
            _ => unreachable!("r_as_transcript gives c. or n."),
        }
    }
}
