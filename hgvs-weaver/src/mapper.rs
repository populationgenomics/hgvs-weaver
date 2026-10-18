use crate::allele::CanonicalAllele;
use crate::coords::TranscriptPos;
use crate::data::{DataProvider, IdentifierKind, IdentifierType, TranscriptData, TranscriptSearch};
use crate::error::HgvsError;
use crate::normalize::{self, ExonBounds, PlacedEdit};
use crate::reference::ReferenceStore;
use crate::structs::Anchor;
use crate::structs::Variant;
use crate::structs::{
    BaseOffsetInterval, BaseOffsetPosition, CVariant, GVariant, GenomicPos, LinearVariant,
    NVariant, PVariant, RVariant, SimpleInterval, SimplePosition, TranscriptVariant,
};
use crate::transcript_mapper::{Footprint, TranscriptMapper};
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
        base: base.into(),
        end: None,
        uncertain: false,
    }
}

/// The edits whose meaning is bases over a range: a substitution, delins or
/// identity, a deletion, a duplication, an inversion. An insertion is placed
/// between bases and is not rewritten; a repeat is read against its run, by
/// [`project_repeat`].
fn resolves_over_range(edit: &crate::edits::NaEdit) -> bool {
    use crate::edits::NaEdit;
    matches!(
        edit,
        NaEdit::RefAlt { .. } | NaEdit::Del { .. } | NaEdit::Dup { .. } | NaEdit::Inv { .. }
    )
}

/// The edit to write on the target of a projection, and the target range to
/// write it over, half-open.
///
/// The edit is resolved on the source to the bases it removes and the bases it
/// puts there, and those are what the target must end up holding. Where the
/// target's bases over the projected range are the source's, the edit is
/// carried as given, complemented across strands. Where they differ, the
/// minimal edit that turns the target's bases into the alternate is written:
/// MUC2's `c.12468dup` over a record C and a genomic G is `g.…delinsCC`, as
/// VariantValidator writes it, since the duplicated bases are the record's;
/// `c.12468C>A` is `g.…G>A`; the record's `c.1568=` over a genomic T is
/// `g.…T>C`. This is biocommons hgvs's replace_reference, extended to the
/// edits whose alternate is derived from the source's bases.
///
/// The bases the edit states are checked on the source. With no source, an
/// intronic position, the edit is carried as given, and the caller checks
/// what it states on the target.
fn project_edit(
    edit: &crate::edits::NaEdit,
    source: Option<(&crate::reference::Reference<'_, '_>, (usize, usize))>,
    target: &crate::reference::Reference<'_, '_>,
    target_range: (usize, usize),
    strand: crate::data::Strand,
) -> Result<(crate::edits::NaEdit, usize, usize), HgvsError> {
    let (t_lo, t_hi) = target_range;
    let carried = || (apply_strand_complement(edit.clone(), strand), t_lo, t_hi);
    let Some((source, (s_lo, s_hi))) = source else {
        return Ok(carried());
    };
    edit.check_stated(source.accession(), s_lo, s_hi, |s, e| source.slice(s, e))?;
    if let Some(repeat) = Repeat::of(edit) {
        return project_repeat(
            &repeat,
            (source, (s_lo, s_hi)),
            target,
            target_range,
            strand,
        );
    }
    if !resolves_over_range(edit) {
        return Ok(carried());
    }
    let resolved = edit.resolve(source, s_lo, s_hi)?;
    let orient = |b: &str| match strand {
        crate::data::Strand::Minus => crate::utils::reverse_complement(b),
        crate::data::Strand::Plus => b.to_string(),
    };
    let (source_ref, alt) = (orient(&resolved.ref_), orient(&resolved.alt));
    let target_ref = target.slice(t_lo, t_hi)?;
    if target_ref.len() != t_hi - t_lo {
        return Err(HgvsError::ValidationError(format!(
            "{} has {} bases over [{t_lo}, {t_hi}); the projected range runs past its end",
            target.accession(),
            target_ref.len()
        )));
    }
    if target_ref == source_ref {
        return Ok(carried());
    }
    normalize::minimal_edit(target, t_lo, t_hi, &target_ref, &alt)
}

/// Whether the genome can count a longer run than a transcript repeat's run
/// over indices `start..=last` (`last` `None` for one base): on the minus
/// strand, where the genome reads the run on from its 5' end, and on the plus
/// strand where the run reaches its exon's last base and the intron can
/// continue it, or its exon is aligned with gaps. Otherwise the two runs are
/// one, and need no genome to tell.
pub(crate) fn genome_may_count_more(
    am: &TranscriptMapper,
    start: usize,
    last: Option<usize>,
) -> Result<bool, HgvsError> {
    if am.transcript.strand == crate::data::Strand::Minus {
        return Ok(true);
    }
    let last = last.unwrap_or(start) as i32;
    let exon = am
        .transcript
        .exons
        .iter()
        .find(|x| x.transcript_start.0 <= last && last < x.transcript_end.0)
        .ok_or_else(|| {
            HgvsError::ValidationError(format!(
                "{} has no exon holding transcript index {last}",
                am.transcript.ac
            ))
        })?;
    let gapped = exon.cigar.contains(['I', 'D', 'N', 'S']);
    Ok(gapped || last + 1 == exon.transcript_end.0)
}

/// The fields of an [`NaEdit::Repeat`](crate::edits::NaEdit::Repeat), for
/// projecting it.
struct Repeat {
    unit: Option<String>,
    min: i32,
    max: i32,
    uncertain: bool,
}

impl Repeat {
    fn of(edit: &crate::edits::NaEdit) -> Option<Self> {
        match edit {
            crate::edits::NaEdit::Repeat {
                ref_,
                min,
                max,
                uncertain,
            } => Some(Repeat {
                unit: ref_.clone(),
                min: *min,
                max: *max,
                uncertain: *uncertain,
            }),
            _ => None,
        }
    }

    fn edit(&self) -> crate::edits::NaEdit {
        crate::edits::NaEdit::Repeat {
            ref_: self.unit.clone(),
            min: self.min,
            max: self.max,
            uncertain: self.uncertain,
        }
    }

    /// The repeat of `unit`, stated, with `more` copies added to its count:
    /// the copies a run holds beyond the ones the repeat was written over.
    fn grown(&self, unit: String, more: usize) -> Result<crate::edits::NaEdit, HgvsError> {
        let overflow = || {
            HgvsError::ValidationError(format!(
                "{} more copies of {unit} overflow a repeat count of {} to {}",
                more, self.min, self.max
            ))
        };
        let more = i32::try_from(more).map_err(|_| overflow())?;
        let (min, max) = self
            .min
            .checked_add(more)
            .zip(self.max.checked_add(more))
            .ok_or_else(overflow)?;
        Ok(crate::edits::NaEdit::Repeat {
            ref_: Some(unit),
            min,
            max,
            uncertain: self.uncertain,
        })
    }
}

/// Whether transcript indices `[start, end)` lie in one exon.
fn within_one_exon(am: &TranscriptMapper, start: i32, end: i32) -> bool {
    am.exons
        .iter()
        .any(|x| x.transcript_start.0 <= start && end <= x.transcript_end.0)
}

/// The end of the run of `unit` on `reference` from `start`, half-open.
fn run_up(
    reference: &crate::reference::Reference<'_, '_>,
    unit: &str,
    start: usize,
) -> Result<usize, HgvsError> {
    let mut end = start;
    while reference.slice(end, end + unit.len())? == unit {
        end += unit.len();
    }
    Ok(end)
}

/// The start of the run of `unit` on `reference` that ends at `end`.
fn run_down(
    reference: &crate::reference::Reference<'_, '_>,
    unit: &str,
    end: usize,
) -> Result<usize, HgvsError> {
    let mut start = end;
    while start >= unit.len() && reference.slice(start - unit.len(), start)? == unit {
        start -= unit.len();
    }
    Ok(start)
}

/// The half-open end of the run a repeat counts from `start` on `reference`,
/// and the repeat with its unit stated; `None` for any other edit or an empty
/// run.
fn repeat_over_run(
    edit: &crate::edits::NaEdit,
    reference: &crate::reference::Reference<'_, '_>,
    start: usize,
    end: usize,
) -> Result<Option<(usize, crate::edits::NaEdit)>, HgvsError> {
    let Some(repeat) = Repeat::of(edit) else {
        return Ok(None);
    };
    let unit = crate::edits::repeat_unit(&repeat.unit, start, end, |s, e| reference.slice(s, e))?;
    let stated = repeat.grown(unit, 0)?;
    let run = stated.resolve(reference, start, end)?;
    Ok((run.end > run.start).then_some((run.end, stated)))
}

/// A repeat over the whole run it counts on the source, `[s_lo, s_hi)`,
/// written for the target over `[t_lo, t_hi)`, the same bases.
///
/// The count is the total of the unit in the run read 3' from the repeat's
/// first base, and the target reads its own run from `t_lo`. That run can be
/// longer than the source's: on the minus strand it continues 5' of the
/// transcript's run, and on either strand an intron can continue a run that
/// ends at an exon's edge. The count grows by the copies the target reads past
/// `t_hi`, which the repeat leaves as they are. Where the target's bases over
/// the range are not the source's, the minimal edit is written, as for any
/// other edit; a count known only as a range then has no single alternate,
/// and is refused.
fn project_repeat(
    repeat: &Repeat,
    source: (&crate::reference::Reference<'_, '_>, (usize, usize)),
    target: &crate::reference::Reference<'_, '_>,
    target_range: (usize, usize),
    strand: crate::data::Strand,
) -> Result<(crate::edits::NaEdit, usize, usize), HgvsError> {
    let (source, (s_lo, s_hi)) = source;
    let (t_lo, t_hi) = target_range;
    let resolved = repeat.edit().resolve(source, s_lo, s_hi)?;
    let orient = |b: &str| match strand {
        crate::data::Strand::Minus => crate::utils::reverse_complement(b),
        crate::data::Strand::Plus => b.to_string(),
    };
    let (source_ref, alt) = (orient(&resolved.ref_), orient(&resolved.alt));
    let target_ref = target.slice(t_lo, t_hi)?;
    if target_ref != source_ref {
        if repeat.min != repeat.max {
            return Err(HgvsError::UnsupportedOperation(format!(
                "{} holds {target_ref} over [{t_lo}, {t_hi}), not the {source_ref} the repeat \
                 counts, so a count of {} to {} copies has no single edit there",
                target.accession(),
                repeat.min,
                repeat.max
            )));
        }
        return normalize::minimal_edit(target, t_lo, t_hi, &target_ref, &alt);
    }
    let unit = crate::edits::repeat_unit(&repeat.unit, s_lo, s_hi, |s, e| source.slice(s, e))?;
    let unit = orient(&unit);
    // Short of `t_hi` where the unit is absent from the range altogether.
    let run_end = run_up(target, &unit, t_lo)?;
    if run_end <= t_hi {
        return Ok((apply_strand_complement(repeat.edit(), strand), t_lo, t_hi));
    }
    let more = (run_end - t_hi) / unit.len();
    Ok((repeat.grown(unit, more)?, t_lo, run_end))
}

/// A transcript repeat at an intronic position of a minus-strand transcript,
/// written for the genome, which holds the bases the transcript has none of,
/// over `[g_lo, g_hi)`.
///
/// The transcript reads its run 3' from the repeat's first base, which on the
/// genome is downwards from `g_hi`; the genome reads its own run upwards from
/// the run's lower end, and adds the copies it holds above `g_hi`, as
/// [`project_repeat`] does for an exonic run.
fn intronic_repeat_on_the_minus_strand(
    repeat: &Repeat,
    genome: &crate::reference::Reference<'_, '_>,
    (g_lo, g_hi): (usize, usize),
) -> Result<(crate::edits::NaEdit, usize, usize), HgvsError> {
    let in_transcript_sense = |s, e| {
        genome
            .slice(s, e)
            .map(|b| crate::utils::reverse_complement(&b))
    };
    let unit = crate::edits::repeat_unit(&repeat.unit, g_lo, g_hi, in_transcript_sense)?;
    let unit = crate::utils::reverse_complement(&unit);
    let lo = run_down(genome, &unit, g_hi)?;
    if lo == g_hi {
        let carried = apply_strand_complement(repeat.edit(), crate::data::Strand::Minus);
        return Ok((carried, g_lo, g_hi));
    }
    let run_end = run_up(genome, &unit, lo)?;
    let more = (run_end - g_hi) / unit.len();
    Ok((repeat.grown(unit, more)?, lo, run_end))
}

/// The edit to write on the target where the source range touches a run of
/// bases the target lacks, and the target range to write it over.
///
/// No position inside the run exists on the target, so the run is the unit:
/// the edit is resolved on the source over the range it names, spliced into
/// the whole run, and the target gets the minimal edit that puts the run, as
/// changed, between its flanking bases. `source` carries the range the edit
/// names and the wider one covering it and every run it touches, both
/// half-open.
fn project_edit_over_run(
    edit: &crate::edits::NaEdit,
    source: (
        &crate::reference::Reference<'_, '_>,
        (usize, usize),
        (usize, usize),
    ),
    target: &crate::reference::Reference<'_, '_>,
    target_range: (usize, usize),
    strand: crate::data::Strand,
) -> Result<(crate::edits::NaEdit, usize, usize), HgvsError> {
    let (source, named, wide) = source;
    let resolved = edit.resolve(source, named.0, named.1)?;
    let (w_lo, w_hi) = (wide.0.min(resolved.start), wide.1.max(resolved.end));
    let source_ref = source.slice(w_lo, w_hi)?;
    let alt = format!(
        "{}{}{}",
        source.slice(w_lo, resolved.start)?,
        resolved.alt,
        source.slice(resolved.end, w_hi)?
    );
    let orient = |b: &str| match strand {
        crate::data::Strand::Minus => crate::utils::reverse_complement(b),
        crate::data::Strand::Plus => b.to_string(),
    };
    let (source_ref, alt) = (orient(&source_ref), orient(&alt));
    let (t_lo, t_hi) = target_range;
    let target_ref = target.slice(t_lo, t_hi)?;
    debug_assert_ne!(target_ref, source_ref, "the target lacks the run");
    normalize::minimal_edit(target, t_lo, t_hi, &target_ref, &alt)
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

/// The 0-based half-open residue range a p. interval names.
fn aa_interval_range(pos: &crate::structs::AaInterval) -> Result<(usize, usize), HgvsError> {
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

/// The 0-based half-open residue range a p. variant names on `reference`,
/// once the residues it states are found there: the residue named at each end
/// of its interval (`p.Arg2Leu`'s Arg), and any the edit states as its
/// reference or as the unit a repeat starts with.
///
/// # Errors
///
/// [`HgvsError::ReferenceMismatch`] for a stated residue the protein does
/// not hold, a position past its end included; a `ValidationError` for an
/// interval that runs backwards or names no position.
pub(crate) fn residue_range(
    vp: &PVariant,
    reference: &crate::reference::Reference<'_, '_>,
) -> Result<(usize, usize), HgvsError> {
    use crate::edits::AaEdit;
    let pos = vp
        .posedit
        .pos
        .as_ref()
        .ok_or_else(|| HgvsError::ValidationError("Missing protein position".into()))?;
    let (start, end) = aa_interval_range(pos)?;
    let check = |stated: &str, at: usize, to: usize| -> Result<(), HgvsError> {
        if stated.is_empty() {
            return Ok(());
        }
        let stated = crate::utils::residues_1(stated)?;
        let found = reference.slice(at, to)?;
        // A protein sequence may end before its stop: Ter is then the
        // position just past the last residue.
        let stop_past_the_end = stated == "*"
            && found.is_empty()
            && at > 0
            && reference.base(at - 1)?.is_some_and(|b| b != b'*');
        if stated.eq_ignore_ascii_case(&found) || stop_past_the_end {
            return Ok(());
        }
        Err(HgvsError::ReferenceMismatch {
            accession: reference.accession().to_string(),
            expected: stated,
            found,
            start: at,
            end: to,
        })
    };
    check(&pos.start.aa, start, start + 1)?;
    if let Some(e) = &pos.end {
        check(&e.aa, end - 1, end)?;
    }
    let stated = match &vp.posedit.edit {
        AaEdit::Subst { ref_, .. } | AaEdit::DelIns { ref_, .. } | AaEdit::Del { ref_, .. } => {
            Some(ref_.as_str())
        }
        AaEdit::Dup { ref_, .. } | AaEdit::RefAlt { ref_, .. } => ref_.as_deref(),
        _ => None,
    };
    if let Some(stated) = stated {
        check(stated, start, end)?;
    }
    if let AaEdit::Repeat {
        ref_: Some(unit), ..
    } = &vp.posedit.edit
    {
        if !crate::edits::is_length(unit) {
            let n = crate::utils::residues_1(unit)?.len();
            check(unit, start, start + n)?;
        }
    }
    Ok((start, end))
}

/// The VRS bounds of a g. interval whose breakpoints are uncertain, HGVS
/// `(a_b)_(c_d)` with `?` for an unknown side; `None` when every position is
/// exact. A deletion `(a_b)_(c_d)` starts at an interbase coordinate in
/// `[a-1, b-1]` and ends at one in `[c, d]`; a single `(a_b)` removes one
/// base somewhere in `a..=b`.
fn uncertain_bounds(pos: &SimpleInterval) -> Option<(VrsBound, VrsBound)> {
    let ranged = |p: &SimplePosition| p.end.is_some() || p.base.known().is_none();
    let last = pos.end.as_ref().unwrap_or(&pos.start);
    if !ranged(&pos.start) && !ranged(last) {
        return None;
    }
    let known = |b: crate::coords::GenomicBound| b.known().map(|p| p.0 as usize);
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
            interval(s, e)?
        }
        (start, end) => bounded_interval(start, end)?,
    })
}

/// The g. base the interbase coordinate `n` names, `n + shift`, if g. can
/// number it.
fn g_base(n: usize, shift: i32) -> Result<crate::coords::HgvsGenomicPos, HgvsError> {
    i32::try_from(n)
        .ok()
        .and_then(|n| n.checked_add(shift))
        .and_then(|h| crate::coords::HgvsGenomicPos::try_from(h).ok())
        .ok_or_else(|| HgvsError::ValidationError(format!("coordinate {n} names no g. position")))
}

/// The g. interval over the 0-based half-open range `[s, e)`.
fn interval(s: usize, e: usize) -> Result<SimpleInterval, HgvsError> {
    let position = |i: usize| -> Result<SimplePosition, HgvsError> {
        Ok(SimplePosition {
            base: g_base(i, 1)?.into(),
            end: None,
            uncertain: false,
        })
    };
    let last = e.checked_sub(1).filter(|&last| last > s);
    Ok(SimpleInterval {
        start: position(s)?,
        end: last.map(position).transpose()?,
        uncertain: false,
    })
}

/// The g. interval `(a_b)_(c_d)` from VRS bounds: the inverse of
/// `uncertain_bounds`. Exact bounds come out as parenthesis-free positions.
fn bounded_interval(start: VrsBound, end: VrsBound) -> Result<SimpleInterval, HgvsError> {
    use crate::coords::GenomicBound;
    let position = |b: VrsBound, shift: i32| -> Result<SimplePosition, HgvsError> {
        let bound = |v: Option<usize>| -> Result<GenomicBound, HgvsError> {
            Ok(v.map(|n| g_base(n, shift)).transpose()?.into())
        };
        Ok(match b {
            VrsBound::Exact(n) => SimplePosition {
                base: g_base(n, shift)?.into(),
                end: None,
                uncertain: false,
            },
            VrsBound::Range(lo, hi) => SimplePosition {
                base: bound(lo)?,
                end: Some(bound(hi)?),
                uncertain: true,
            },
        })
    };
    // A single base somewhere in `a..=b` came out as `[a-1, b-1]`, `[a, b]`.
    let single = matches!((start, end), (VrsBound::Range(lo, hi), VrsBound::Range(lo2, hi2))
        if lo.map(|n| n.checked_add(1)) == lo2.map(Some)
            && hi.map(|n| n.checked_add(1)) == hi2.map(Some));
    Ok(SimpleInterval {
        start: position(start, 1)?,
        end: (!single).then(|| position(end, 0)).transpose()?,
        uncertain: false,
    })
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
    /// A statement rather than a change: `p.?`, `p.(Met1?)`, `p.0?`, `p.(=)`.
    Statement(PVariant),
    Change(crate::protein::Reading),
}

/// The protein statement `statement` predicted from `var_c`, bracketed as
/// HGVS writes each prediction: `p.(=)` and `p.(Met1?)`, but `p.?` and `p.0?`
/// bare.
fn predicted_statement(
    var_c: &CVariant,
    protein_ac: &str,
    pos: Option<crate::structs::AaInterval>,
    statement: crate::structs::ProteinStatement,
) -> PVariant {
    use crate::structs::ProteinStatement;
    PVariant {
        ac: protein_ac.to_string(),
        gene: var_c.gene.clone(),
        posedit: crate::structs::PosEdit {
            predicted: match statement {
                ProteinStatement::Unchanged => true,
                ProteinStatement::Unknown => pos.is_some(),
                ProteinStatement::NoProtein { .. } => false,
            },
            pos,
            edit: crate::edits::AaEdit::Statement(statement),
            uncertain: false,
        },
    }
}

/// `p.?`, `p.(Met1?)`, `p.0?` or `p.(=)` for an edit whose codons cannot be
/// read: an intronic position, a start in the 5'UTR (deleting the whole CDS
/// predicts no protein, reaching into it disrupts the start codon, staying
/// upstream says nothing, since an upstream start may be created), or an
/// edit entirely in the 3'UTR, which cannot change the protein. `None`
/// when the edit is to be read on the codons.
fn statement_about_transcript(
    var_c: &CVariant,
    transcript: &TranscriptData,
    protein_ac: &str,
) -> Option<PVariant> {
    use crate::structs::{AAPosition, AaInterval, ProteinPos, ProteinStatement};
    let pos = var_c.posedit.pos.as_ref()?;
    let statement = |pos: Option<AaInterval>, statement: ProteinStatement| {
        predicted_statement(var_c, protein_ac, pos, statement)
    };
    // A start after the stop codon is an edit entirely in the 3'UTR.
    let in_3_utr = pos.start.anchor == Anchor::CdsEnd;
    let in_5_utr = pos.start.anchor == Anchor::CdsStart && pos.start.base.0 < 0;
    let intronic = has_intronic_offset(pos);
    if !(intronic || in_3_utr || in_5_utr) {
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
    if intronic {
        return Some(statement(None, ProteinStatement::Unknown));
    }
    if in_3_utr {
        return Some(statement(None, ProteinStatement::Unchanged));
    }
    let deletes = matches!(var_c.posedit.edit, crate::edits::NaEdit::Del { .. });
    Some(if covers_cds && deletes {
        statement(None, ProteinStatement::NoProtein { uncertain: true })
    } else if reaches_cds {
        let met1 = AaInterval {
            start: AAPosition {
                base: ProteinPos(0).to_hgvs(),
                aa: "Met".to_string(),
            },
            end: None,
            uncertain: false,
        };
        statement(Some(met1), ProteinStatement::Unknown)
    } else {
        statement(None, ProteinStatement::Unknown)
    })
}

/// The first (`from_start`) or last coding base a range end covers, as a
/// 0-based index from the first base of the CDS: negative in the 5'UTR, at
/// or past `cds_len` in the 3'UTR. An end in an intron covers up to the
/// exonic base on its inner side.
fn covered_cds_index(
    pos: &crate::structs::BaseOffsetPosition,
    cds_len: i64,
    from_start: bool,
) -> Result<i64, HgvsError> {
    let base = pos.base.0 as i64;
    let index = match pos.anchor {
        Anchor::CdsStart if base > 0 => base - 1,
        Anchor::CdsStart => base,
        Anchor::CdsEnd => cds_len - 1 + base,
        Anchor::TranscriptStart => {
            return Err(HgvsError::ValidationError(format!(
                "a c. position is numbered from the CDS, not the transcript start: {base}"
            )))
        }
    };
    Ok(match pos.offset.map(|o| o.0) {
        Some(o) if o > 0 && from_start => index + 1,
        Some(o) if o < 0 && !from_start => index - 1,
        _ => index,
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

/// On a transcript whose CDS runs off the 3' end, a change has a protein
/// consequence only if it is in frame and lies wholly within the whole codons
/// the record carries, short of the last: anything else reads into coding
/// sequence the record does not have, or shifts the frame into it.
fn check_within_open_cds(
    ac: &str,
    resolved: &crate::edits::ResolvedEdit,
    cds_start: usize,
    cds_end: usize,
) -> Result<(), HgvsError> {
    let whole_codons = (cds_end + 1 - cds_start) / 3;
    let last_whole_codon_start = cds_start + whole_codons.saturating_sub(1) * 3;
    let in_frame = (resolved.alt.len() as i64 - (resolved.end - resolved.start) as i64) % 3 == 0;
    if in_frame && resolved.end <= last_whole_codon_start {
        return Ok(());
    }
    Err(HgvsError::ValidationError(format!(
        "{ac}: its CDS is open at the 3' end, so the record carries no stop codon; a protein consequence is \
         predicted only for an in-frame change wholly inside the coding bases it carries"
    )))
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
    var_c
        .posedit
        .edit
        .resolve_with(&var_c.ac, start, end, |s, e| Ok(window(s, e).to_string()))
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
/// r. and c./n. Statements about the transcript (r.0, r.spl, r.=, c.=) have
/// no other spelling.
fn relettered_posedit(
    var: &dyn std::fmt::Display,
    posedit: &crate::structs::PosEdit<BaseOffsetInterval, crate::edits::NaEdit>,
    letters: Letters,
    anchor: impl Fn(Anchor) -> Anchor,
) -> Result<crate::structs::PosEdit<BaseOffsetInterval, crate::edits::NaEdit>, HgvsError> {
    if posedit.edit.is_whole_transcript() {
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
pub(crate) fn simple_interval_range(pos: &SimpleInterval) -> Result<(usize, usize), HgvsError> {
    let start_i = pos.start.index()?.0;
    if start_i < 0 {
        return Err(HgvsError::ValidationError(format!(
            "Genomic start position {} is negative or uncertain",
            start_i
        )));
    }
    let start = start_i as usize;
    let last = match &pos.end {
        Some(e) => {
            let end_i = e.index()?.0;
            if end_i < 0 {
                return Err(HgvsError::ValidationError(format!(
                    "Genomic end position {} is negative or uncertain",
                    end_i
                )));
            }
            if (end_i as usize) < start {
                return Err(HgvsError::ValidationError(format!(
                    "Genomic range runs backwards: {} is after {}",
                    pos.start, e
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

/// The certain g. position of the 0-based index `i`.
pub(crate) fn simple_position(i: usize) -> SimplePosition {
    SimplePosition {
        base: GenomicPos(i as i32).to_hgvs().into(),
        end: None,
        uncertain: false,
    }
}

/// Whether a bound of `pos` is known only to lie in a range, `(3_4)`, or is
/// marked uncertain.
pub(crate) fn is_uncertain_linear(pos: &SimpleInterval) -> bool {
    let uncertain = |p: &SimplePosition| p.uncertain || p.end.is_some();
    pos.uncertain || uncertain(&pos.start) || pos.end.as_ref().is_some_and(uncertain)
}

/// The HGVS start index and optional end index to write for `after`: an end
/// for an insertion, which names its two flanking bases, and for an edit of
/// more than one base.
pub(crate) fn hgvs_positions(after: &PlacedEdit) -> (usize, Option<usize>) {
    let (s, e) = after.hgvs_range();
    let last = e - 1;
    (s, (after.is_insertion() || last > s).then_some(last))
}

pub(crate) fn has_intronic_offset(pos: &BaseOffsetInterval) -> bool {
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
        let run = self.repeat_run_g(var_g, pos)?;
        let run = match run {
            Some(run) => Some(self.transcript_run_g(&am, var_g, run)?),
            None => None,
        };
        let (pos, edit) = match &run {
            Some((interval, edit)) => (interval, edit),
            None => (pos, &var_g.posedit.edit),
        };
        if let Some(end) = &pos.end {
            if let (Some(s), Some(e)) = (pos.start.base.known(), end.base.known()) {
                if e.0 < s.0 {
                    return Err(HgvsError::ValidationError(format!(
                        "Genomic range runs backwards: {} is after {}",
                        s.0, e.0
                    )));
                }
            }
        }
        // A g. variant states its bases on the genome, intronic or not.
        if var_g.posedit.edit.states_bases() {
            let (s, e) = simple_interval_range(pos)?;
            let genome = self
                .refs
                .reference(&var_g.ac, IdentifierType::GenomicAccession);
            var_g
                .posedit
                .edit
                .check_stated(&var_g.ac, s, e, |a, b| genome.slice(a, b))?;
        }
        // Each end occupies a half-open transcript range: one base, or none
        // for a base the transcript lacks.
        let fp_start = am.g_to_n_footprint(pos.start.index()?)?;
        let fp_end = match &pos.end {
            Some(end) => am.g_to_n_footprint(end.index()?)?,
            None => fp_start,
        };
        let (mut n_lo, mut off_lo) = (fp_start.range.0, fp_start.offset);
        let (mut n_hi, mut off_hi) = (fp_end.range.0, fp_end.offset);
        // Two intronic ends anchored on one exon boundary share n and differ
        // only in offset, so the order is by (base, offset).
        if (n_lo.0, off_lo.0) > (n_hi.0, off_hi.0) {
            std::mem::swap(&mut n_lo, &mut n_hi);
            std::mem::swap(&mut off_lo, &mut off_hi);
        }
        let (t_lo, t_hi) = (
            fp_start.range.0.min(fp_end.range.0).0 as usize,
            fp_start.range.1.max(fp_end.range.1).0 as usize,
        );
        let position = |n: TranscriptPos, off: crate::structs::IntronicOffset| {
            let (c, c_off, anchor) = am.n_to_c(n)?;
            Ok::<_, HgvsError>(make_base_offset_position(
                c.to_hgvs(),
                c_off.0 + off.0,
                anchor,
            ))
        };

        // An intronic position has no transcript base: the edit is carried as
        // given. An exonic one is read on the genome and written for the record.
        let exonic = off_lo.0 == 0 && off_hi.0 == 0 && n_lo.0 >= 0;
        let (edit, start, end) = if exonic {
            let named = simple_interval_range(pos)?;
            let source = self
                .refs
                .reference(&var_g.ac, IdentifierType::GenomicAccession);
            let target = self
                .refs
                .reference(transcript_ac, IdentifierType::TranscriptAccession);
            // A base the transcript lacks has no index of its own: the
            // transcript can describe only the whole run it belongs to, so an
            // edit touching a run is written over all of it.
            let touches_run = fp_start.gap_run.is_some() || fp_end.gap_run.is_some();
            let runs = [fp_start.gap_run, fp_end.gap_run].into_iter().flatten();
            let wide = runs.fold(named, |(lo, hi), (r_lo, r_hi)| {
                (lo.min(r_lo.0 as usize), hi.max(r_hi.0 as usize))
            });
            let (edit, s, e) = if touches_run {
                project_edit_over_run(
                    edit,
                    (&source, named, wide),
                    &target,
                    (t_lo, t_hi),
                    am.transcript.strand,
                )?
            } else {
                project_edit(
                    edit,
                    Some((&source, named)),
                    &target,
                    (t_lo, t_hi),
                    am.transcript.strand,
                )?
            };
            let repeat = matches!(edit, crate::edits::NaEdit::Repeat { .. });
            if repeat && !within_one_exon(&am, s as i32, e as i32) {
                return Err(HgvsError::UnsupportedOperation(format!(
                    "{var_g} counts a run that the transcript continues across an exon \
                     junction: the transcript's count has no single place on the genome"
                )));
            }
            // An edit over a run is always rewritten for the transcript.
            let unchanged_range = !touches_run && (s, e) == (t_lo, t_hi);
            let end = if unchanged_range {
                pos.end
                    .as_ref()
                    .map(|_| position(n_hi, off_hi))
                    .transpose()?
            } else {
                (e > s + 1)
                    .then(|| position(TranscriptPos(e as i32 - 1), off_hi))
                    .transpose()?
            };
            (edit, position(TranscriptPos(s as i32), off_lo)?, end)
        } else {
            let edit = apply_strand_complement(edit.clone(), am.transcript.strand);
            let end = pos
                .end
                .as_ref()
                .map(|_| position(n_hi, off_hi))
                .transpose()?;
            (edit, position(n_lo, off_lo)?, end)
        };

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
    /// reads as a run from its lower end on either strand. The edit returned
    /// states its unit, which an unstated one takes from the range as written.
    fn repeat_run_tx<V: TranscriptVariant>(
        &self,
        am: &TranscriptMapper,
        var: &V,
        pos: &BaseOffsetInterval,
    ) -> Result<Option<(BaseOffsetInterval, crate::edits::NaEdit)>, HgvsError> {
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
        let run = repeat_over_run(edit, &reference, n_start.0 as usize, n_end.0 as usize)?;
        if let Some((end, _)) = run {
            if !within_one_exon(am, n_start.0, end as i32) {
                return Err(HgvsError::UnsupportedOperation(format!(
                    "{}:{}.{} counts a run that crosses an exon junction: the genome has \
                     no single place for the change",
                    var.ac(),
                    var.coordinate_type(),
                    var.posedit()
                )));
            }
        }
        run.map(|(end, edit)| {
            let interval = BaseOffsetInterval {
                start: pos.start,
                end: (end - 1 > n_start.0 as usize)
                    .then(|| V::position_from_index(am, (end - 1) as i32))
                    .transpose()?,
                uncertain: pos.uncertain,
            };
            Ok((interval, edit))
        })
        .transpose()
    }

    /// The genomic counterpart of [`repeat_run_tx`](Self::repeat_run_tx).
    fn repeat_run_g(
        &self,
        var_g: &GVariant,
        pos: &SimpleInterval,
    ) -> Result<Option<(SimpleInterval, crate::edits::NaEdit)>, HgvsError> {
        if !matches!(var_g.posedit.edit, crate::edits::NaEdit::Repeat { .. }) {
            return Ok(None);
        }
        let (start, end) = simple_interval_range(pos)?;
        let reference = self
            .refs
            .reference(&var_g.ac, IdentifierType::GenomicAccession);
        let run = repeat_over_run(&var_g.posedit.edit, &reference, start, end)?;
        Ok(run.map(|(end, edit)| {
            let interval = SimpleInterval {
                start: pos.start,
                end: (end - 1 > start).then(|| SimplePosition {
                    base: GenomicPos((end - 1) as i32).to_hgvs().into(),
                    end: None,
                    uncertain: false,
                }),
                uncertain: pos.uncertain,
            };
            (interval, edit)
        }))
    }

    /// A genomic repeat over its run, `run`, extended to the run the transcript
    /// counts where the run is intronic on a minus-strand transcript: the
    /// transcript reads it from its upper end downwards, past the genome's run
    /// where the bases below continue it, and those copies join the count. An
    /// exonic run is read on the record instead, by [`project_repeat`].
    fn transcript_run_g(
        &self,
        am: &TranscriptMapper,
        var_g: &GVariant,
        run: (SimpleInterval, crate::edits::NaEdit),
    ) -> Result<(SimpleInterval, crate::edits::NaEdit), HgvsError> {
        let (interval, edit) = run;
        let Some(repeat) = Repeat::of(&edit) else {
            return Ok((interval, edit));
        };
        let (start, end) = simple_interval_range(&interval)?;
        let exonic = |i: usize| -> Result<bool, HgvsError> {
            let fp = am.g_to_n_footprint(GenomicPos(i as i32))?;
            Ok(fp.offset.0 == 0 && fp.range.0 .0 >= 0)
        };
        let minus = am.transcript.strand == crate::data::Strand::Minus;
        if !minus || (exonic(start)? && exonic(end - 1)?) {
            return Ok((interval, edit));
        }
        let genome = self
            .refs
            .reference(&var_g.ac, IdentifierType::GenomicAccession);
        let unit = crate::edits::repeat_unit(&repeat.unit, start, end, |s, e| genome.slice(s, e))?;
        // Above `start` where the unit is absent from the range altogether.
        let lo = run_down(&genome, &unit, end)?;
        if lo >= start {
            return Ok((interval, edit));
        }
        let more = (start - lo) / unit.len();
        let extended = SimpleInterval {
            start: make_simple_position(GenomicPos(lo as i32).to_hgvs()),
            end: Some(make_simple_position(GenomicPos(end as i32 - 1).to_hgvs())),
            uncertain: interval.uncertain,
        };
        Ok((extended, repeat.grown(unit, more)?))
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
        if matches!(var_c.posedit().edit, crate::edits::NaEdit::WholeIdentity) {
            return Err(HgvsError::UnsupportedOperation(format!(
                "{}:{}.{} describes the whole transcript, which is spliced from its exons, \
                 so it has no genomic form",
                var_c.ac(),
                var_c.coordinate_type(),
                var_c.posedit()
            )));
        }
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
        if !pos.is_exact() {
            return Err(HgvsError::UnsupportedOperation(format!(
                "{}:{}.{} has an uncertain position, which a genomic interval cannot carry",
                var_c.ac(),
                var_c.coordinate_type(),
                var_c.posedit()
            )));
        }
        let run = self.repeat_run_tx(&am, var_c, pos)?;
        let (pos, edit) = match &run {
            Some((interval, edit)) => (interval, edit),
            None => (pos, &var_c.posedit().edit),
        };
        if let Some(end) = &pos.end {
            if end.order_key() < pos.start.order_key() {
                return Err(HgvsError::ValidationError(format!(
                    "Transcript range runs backwards: {} is after {end}",
                    pos.start
                )));
            }
        }
        let zero = crate::structs::IntronicOffset(0);
        // Each end occupies a half-open genomic range: one base, or none for
        // a base the genome lacks.
        let project = |p: &BaseOffsetPosition| -> Result<(TranscriptPos, Footprint), HgvsError> {
            let n = am.c_to_n(p.base.to_index(), p.anchor)?;
            Ok((n, am.position_footprint(p)?))
        };
        let (n_start, fp_start) = project(&pos.start)?;
        let (n_end, fp_end) = match &pos.end {
            Some(end) => project(end)?,
            None => (n_start, fp_start),
        };
        let ((a_lo, a_hi), (b_lo, b_hi)) = (fp_start.range, fp_end.range);
        let (g_lo, g_hi) = (a_lo.0.min(b_lo.0), a_hi.0.max(b_hi.0));
        if g_lo < 0 {
            return Err(HgvsError::ValidationError(format!(
                "{}:{} projects before the start of {target_ac}",
                var_c.ac(),
                var_c.posedit()
            )));
        }

        // An exonic edit is read on the record and written for the genome. An
        // intronic one has no record bases: what it states is checked on the
        // genome, read on the transcript's strand, and it is carried as given.
        let exonic = pos.start.offset.unwrap_or(zero).0 == 0
            && pos
                .end
                .as_ref()
                .is_none_or(|e| e.offset.unwrap_or(zero).0 == 0)
            && n_start.0 >= 0;
        let target = self
            .refs
            .reference(&target_ac, IdentifierType::GenomicAccession);
        let source = self
            .refs
            .reference(var_c.ac(), IdentifierType::TranscriptAccession);
        let named = (
            n_start.0.min(n_end.0) as usize,
            n_start.0.max(n_end.0) as usize + 1,
        );
        // A base the genome lacks has no position of its own: the genome can
        // describe only the whole run it belongs to, so an edit touching a
        // run is written over all of it.
        let touches_run = fp_start.gap_run.is_some() || fp_end.gap_run.is_some();
        let runs = [fp_start.gap_run, fp_end.gap_run].into_iter().flatten();
        let wide = runs.fold(named, |(lo, hi), (r_lo, r_hi)| {
            (lo.min(r_lo.0 as usize), hi.max(r_hi.0 as usize))
        });
        if !exonic {
            let (lo, hi) = (g_lo as usize, g_hi as usize);
            // Offsets from `lo` count along the transcript: on the minus
            // strand they run down from `hi`, complemented.
            var_c
                .posedit()
                .edit
                .check_stated(&target_ac, lo, hi, |s, e| {
                    Ok(match am.transcript.strand {
                        crate::data::Strand::Minus => crate::utils::reverse_complement(
                            &target.slice((lo + hi).saturating_sub(e), lo + hi - s)?,
                        ),
                        crate::data::Strand::Plus => target.slice(s, e)?,
                    })
                })?;
        }
        let minus = am.transcript.strand == crate::data::Strand::Minus;
        let intronic_repeat_on_minus = Repeat::of(edit).filter(|_| minus && !exonic);
        let (edit, s, e) = if touches_run {
            project_edit_over_run(
                edit,
                (&source, named, wide),
                &target,
                (g_lo as usize, g_hi as usize),
                am.transcript.strand,
            )?
        } else if let Some(repeat) = intronic_repeat_on_minus {
            intronic_repeat_on_the_minus_strand(&repeat, &target, (g_lo as usize, g_hi as usize))?
        } else {
            project_edit(
                edit,
                exonic.then_some((&source, named)),
                &target,
                (g_lo as usize, g_hi as usize),
                am.transcript.strand,
            )?
        };
        // An edit over a run is always rewritten for the genome.
        let unchanged_range = !touches_run && (s, e) == (g_lo as usize, g_hi as usize);

        Ok(GVariant {
            ac: target_ac,
            gene: var_c.gene().map(str::to_string),
            posedit: crate::structs::PosEdit {
                pos: Some(crate::structs::SimpleInterval {
                    start: make_simple_position(GenomicPos(s as i32).to_hgvs()),
                    // The writer's form survives an edit carried as given; a
                    // rewritten one names exactly the bases it covers.
                    end: if unchanged_range {
                        pos.end
                            .as_ref()
                            .map(|_| make_simple_position(GenomicPos(e as i32 - 1).to_hgvs()))
                    } else {
                        (e > s + 1)
                            .then(|| make_simple_position(GenomicPos(e as i32 - 1).to_hgvs()))
                    },
                    uncertain: false,
                }),
                edit,
                uncertain: var_c.posedit().uncertain,
                predicted: var_c.posedit().predicted,
            },
        })
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
        let start_0 = pos.start.index()?.0;
        let end_0 = match &pos.end {
            Some(e) => e.index()?.0 + 1,
            None => start_0 + 1,
        };

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
                // A base the genome lacks fails on every transcript alike.
                Err(e @ HgvsError::ReferenceMismatch { .. }) => return Err(e),
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
    /// about the protein (`p.0?`, `p.?`, `p.?`, `p.(=)`); `r.(?)`,
    /// which names no change of its own, is refused; anything else is
    /// predicted from its c. spelling. Either way the transcript must have a
    /// CDS.
    pub fn r_to_p(&self, r: &RVariant, protein_ac: Option<&str>) -> Result<PVariant, HgvsError> {
        use crate::structs::{ProteinStatement, RnaStatement};
        let statement = match &r.posedit.edit {
            crate::edits::NaEdit::WholeIdentity => Some((ProteinStatement::Unchanged, true)),
            crate::edits::NaEdit::Statement(RnaStatement::FollowsDna) => {
                return Err(HgvsError::UnsupportedOperation(format!(
                    "{r} says the RNA follows the DNA change it accompanies, which it does not \
                     carry; predict the protein from that change"
                )))
            }
            // No RNA predicts no protein, p.0? as HGVS writes that prediction.
            crate::edits::NaEdit::Statement(RnaStatement::NoRna { .. }) => {
                Some((ProteinStatement::NoProtein { uncertain: true }, false))
            }
            crate::edits::NaEdit::Statement(
                RnaStatement::SplicingAffected { .. } | RnaStatement::Unknown,
            ) => Some((ProteinStatement::Unknown, false)),
            _ => None,
        };
        if let Some((p, predicted)) = statement {
            if !self.has_cds(&r.ac)? {
                return Err(HgvsError::ValidationError(format!(
                    "{} has no CDS, so {r} has no protein consequence",
                    r.ac
                )));
            }
            return Ok(PVariant {
                ac: self.protein_accession(&r.ac, protein_ac)?,
                gene: r.gene.clone(),
                posedit: crate::structs::PosEdit {
                    pos: None,
                    edit: crate::edits::AaEdit::Statement(p),
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
            CodingOutcome::Change(reading) => {
                let mut var_p = reading.describe()?;
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
                crate::edits::AaEdit::Statement(crate::edits::ProteinStatement::NoProtein {
                    ..
                }) => {
                    let whole = self
                        .refs
                        .reference(&p.ac, IdentifierType::ProteinAccession)
                        .whole()?;
                    (p.ac.clone(), strip_stop(whole), String::new())
                }
                // An identity: the reference over the residues it names.
                crate::edits::AaEdit::Identity { .. } => {
                    let pos = p.posedit.pos.as_ref().ok_or_else(|| {
                        HgvsError::ValidationError(format!("{p} names no residue"))
                    })?;
                    let whole = strip_stop(
                        self.refs
                            .reference(&p.ac, IdentifierType::ProteinAccession)
                            .whole()?,
                    );
                    let index = |a: &crate::structs::AAPosition| a.base.to_index().0 as usize;
                    let first = index(&pos.start).min(whole.len());
                    let end = (pos.end.as_ref().map_or(first, index) + 1).clamp(first, whole.len());
                    touched = Some((first, end));
                    (p.ac.clone(), whole.clone(), whole)
                }
                _ => {
                    return Err(HgvsError::UnsupportedOperation(format!(
                        "{p} names no protein sequence, so it has no allele"
                    )))
                }
            },
            CodingOutcome::Change(reading) => {
                let (r, a) = reading.proteins()?;
                let change = reading.change();
                if r == a {
                    let first = change.edit.start / 3;
                    let last = change.edit.end.max(change.edit.start + 1) - 1;
                    touched = Some((
                        first.min(r.len()),
                        (last / 3 + 1).clamp(first.min(r.len()), r.len()),
                    ));
                }
                (change.protein_ac.clone(), r, a)
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
                crate::edits::AaEdit::Statement(crate::edits::ProteinStatement::NoProtein {
                    ..
                }) => Some(String::new()),
                _ => None,
            },
            CodingOutcome::Change(reading) => Some(reading.proteins()?.1),
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
    /// `p.(Met1?)`, `p.0?`, `p.(=)`) when the edit lies outside or across the CDS
    /// boundary or leaves no codon where the protein first changes, else the
    /// resolved change to describe or read.
    fn coding_outcome(
        &self,
        var_c: &CVariant,
        protein_ac: Option<&str>,
    ) -> Result<CodingOutcome, HgvsError> {
        if matches!(var_c.posedit.edit, crate::edits::NaEdit::WholeIdentity) {
            return Ok(CodingOutcome::Statement(
                self.unchanged_protein(var_c, protein_ac)?,
            ));
        }
        let protein_ac = self.protein_accession(&var_c.ac, protein_ac)?;
        let transcript = self.provider().get_transcript(&var_c.ac, None)?;
        if transcript.cds_start_open {
            return Err(HgvsError::ValidationError(format!(
                "{}: its CDS is open at the 5' end, so the record carries no start codon and no protein \
                 consequence can be predicted",
                var_c.ac
            )));
        }
        if var_c.posedit.edit.is_identity() {
            self.check_transcript_stated(var_c)?;
            return Ok(CodingOutcome::Statement(self.identity_prediction(
                var_c,
                &transcript,
                protein_ac,
            )?));
        }
        if let Some(statement) = statement_about_transcript(var_c, &transcript, &protein_ac) {
            // A statement reads no codon, but the bases the edit states must
            // still be there.
            if var_c.posedit.edit.states_bases() {
                self.check_transcript_stated(var_c)?;
            }
            return Ok(CodingOutcome::Statement(statement));
        }
        let ref_seq = self
            .refs
            .reference(&var_c.ac, IdentifierType::TranscriptAccession)
            .whole()?;
        let (cds_start, cds_end) = cds_bounds(&transcript, ref_seq.len())?;
        let cds_end_open = transcript.cds_end_open;
        let am = TranscriptMapper::new(transcript)?;
        let resolved = resolve_in_transcript(var_c, &am, &ref_seq)?;
        if cds_end_open {
            check_within_open_cds(&var_c.ac, &resolved, cds_start, cds_end)?;
        }
        let rel = |i: usize| {
            i.checked_sub(cds_start).ok_or_else(|| {
                HgvsError::ValidationError(format!("Position {} before the CDS start", i))
            })
        };
        let change = crate::protein::CodingChange {
            coding: ref_seq[cds_start..].to_string(),
            cds_len: cds_end + 1 - cds_start,
            edit: crate::edits::ResolvedEdit {
                start: rel(resolved.start)?,
                end: rel(resolved.end)?,
                ref_: resolved.ref_,
                alt: resolved.alt,
            },
            protein_ac,
        };
        let reading = crate::protein::Reading::new(change);
        if reading.runs_off() {
            return Ok(CodingOutcome::Statement(predicted_statement(
                var_c,
                &reading.change().protein_ac,
                None,
                crate::structs::ProteinStatement::Unknown,
            )));
        }
        Ok(CodingOutcome::Change(reading))
    }

    /// What an identity predicts: the residues whose whole codon it covers,
    /// unchanged, `p.(Met1=)`, or `p.(=)` when it covers every codon of a CDS
    /// that ends in a stop. Where its ends lie, in a UTR or an intron, does
    /// not enter.
    ///
    /// # Errors
    ///
    /// [`HgvsError::NoProteinStatement`] when it covers no whole codon: it
    /// states bases, not residues.
    fn identity_prediction(
        &self,
        var_c: &CVariant,
        transcript: &TranscriptData,
        protein_ac: String,
    ) -> Result<PVariant, HgvsError> {
        use crate::structs::{AAPosition, AaInterval, ProteinPos, ProteinStatement};
        let pos = var_c
            .posedit
            .pos
            .as_ref()
            .ok_or_else(|| HgvsError::ValidationError("Missing position".into()))?;
        let ref_seq = self
            .refs
            .reference(&var_c.ac, IdentifierType::TranscriptAccession)
            .whole()?;
        let (cds_start, cds_end) = cds_bounds(transcript, ref_seq.len())?;
        let cds_len = (cds_end + 1 - cds_start) as i64;
        let codons = cds_len / 3;
        let first = covered_cds_index(&pos.start, cds_len, true)?;
        let last = covered_cds_index(pos.end.as_ref().unwrap_or(&pos.start), cds_len, false)?;
        // Whole codons in [first, last], 0-based and half-open.
        let from = (first.max(0) as u64).div_ceil(3) as i64;
        let to = ((last + 1).max(0) / 3).min(codons);
        let every_codon = from == 0 && to == codons && !transcript.cds_end_open;
        if from >= to {
            return Err(HgvsError::NoProteinStatement(format!(
                "{var_c} covers no whole codon, so it states no residue of {protein_ac}"
            )));
        }
        if every_codon {
            return Ok(predicted_statement(
                var_c,
                &protein_ac,
                None,
                ProteinStatement::Unchanged,
            ));
        }
        let residue = |codon: i64| -> Result<AAPosition, HgvsError> {
            let at = cds_start + 3 * codon as usize;
            let aa = crate::utils::translate_codon(&ref_seq[at..at + 3]).ok_or_else(|| {
                HgvsError::ValidationError(format!("{}: no codon at {at}", var_c.ac))
            })?;
            Ok(AAPosition {
                base: ProteinPos(codon as i32).to_hgvs(),
                aa: crate::utils::aa1_to_aa3(aa).to_string(),
            })
        };
        Ok(PVariant {
            ac: protein_ac,
            gene: var_c.gene.clone(),
            posedit: crate::structs::PosEdit {
                pos: Some(AaInterval {
                    start: residue(from)?,
                    end: (to - 1 > from).then(|| residue(to - 1)).transpose()?,
                    uncertain: false,
                }),
                edit: crate::edits::AaEdit::Identity { uncertain: false },
                uncertain: false,
                predicted: true,
            },
        })
    }

    /// `p.(=)` for `c.=`: the whole transcript unchanged predicts the protein
    /// it encodes unchanged, whether or not the record carries its start
    /// codon. The transcript must have a CDS, and `protein_ac`, if given, must
    /// be the protein the provider pairs it with.
    fn unchanged_protein(
        &self,
        var_c: &CVariant,
        protein_ac: Option<&str>,
    ) -> Result<PVariant, HgvsError> {
        if !self.has_cds(&var_c.ac)? {
            return Err(HgvsError::ValidationError(format!(
                "{} has no CDS, so {var_c} has no protein consequence",
                var_c.ac
            )));
        }
        let encoded = self.protein_accession(&var_c.ac, None)?;
        if let Some(ac) = protein_ac.filter(|&ac| ac != encoded) {
            return Err(HgvsError::ValidationError(format!(
                "{} encodes {encoded}, not {ac}, so {var_c} says nothing about {ac}",
                var_c.ac
            )));
        }
        Ok(PVariant {
            ac: encoded,
            gene: var_c.gene.clone(),
            posedit: crate::structs::PosEdit {
                pos: None,
                edit: crate::edits::AaEdit::Statement(crate::edits::ProteinStatement::Unchanged),
                uncertain: false,
                predicted: true,
            },
        })
    }

    /// Back-converts a single amino acid substitution (p.Xxx###Yyy) to a coding variant (c.).
    ///
    /// Only handles `AaEdit::Subst` currently. For ambiguous back-conversions (multiple
    /// codons with the same minimum nucleotide distance), the first alphabetically is chosen
    /// and `is_unique` is set to false. The residue it names is read on the
    /// protein, so the provider must serve the protein's sequence.
    ///
    /// Returns `(CVariant, is_unique)`.
    pub fn p_to_c(
        &self,
        var_p: &PVariant,
        transcript_ac: Option<&str>,
    ) -> Result<(CVariant, bool), HgvsError> {
        let (aa_pos, alt_aa) = substituted_residue(var_p)?;
        residue_range(
            var_p,
            &self
                .refs
                .reference(&var_p.ac, IdentifierType::ProteinAccession),
        )?;
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

    /// Whether the bases and residues a variant states are what its reference
    /// holds: a `g.` or `m.` base on the genome, an exonic `c.` or `n.` base
    /// on the transcript record, an intronic one on the genome read on the
    /// transcript's strand, a `p.` residue on the protein. Every member of a
    /// cis allele must hold, and an allele of nucleotide members must have
    /// one reading ([`read_cis`](Self::read_cis)).
    ///
    /// Returns `Ok(true)` when they hold, or when the edit states none (or
    /// only a length); `Ok(false)` where an operation reading the variant
    /// refuses it with [`HgvsError::ReferenceMismatch`].
    ///
    /// # Errors
    ///
    /// A contradictory allele is a [`HgvsError::ValidationError`] and an
    /// ambiguous one [`HgvsError::Ambiguous`], raised rather than answered
    /// `false`, as any description that cannot be read is; one whose readings
    /// weaver does not list is [`HgvsError::UnsupportedOperation`].
    pub fn validate(&self, var: &crate::SequenceVariant) -> Result<bool, HgvsError> {
        match self.check_stated(var) {
            Ok(()) => Ok(true),
            Err(HgvsError::ReferenceMismatch { .. }) => Ok(false),
            Err(e) => Err(e),
        }
    }

    /// [`validate`](Self::validate) as an error: `Ok` when every base and
    /// residue the variant states is there.
    pub(crate) fn check_stated(&self, var: &crate::SequenceVariant) -> Result<(), HgvsError> {
        match var {
            crate::SequenceVariant::Genomic(v) => self.check_linear_stated(v),
            crate::SequenceVariant::Mitochondrial(v) => self.check_linear_stated(v),
            crate::SequenceVariant::Coding(v)
                if matches!(v.posedit.edit, crate::edits::NaEdit::WholeIdentity) =>
            {
                self.check_whole_identity(v)
            }
            crate::SequenceVariant::Coding(v) => self.check_transcript_stated(v),
            crate::SequenceVariant::NonCoding(v) => self.check_transcript_stated(v),
            // A statement about the whole protein (`p.=`, `p.?`) names no residue.
            crate::SequenceVariant::Protein(v) if v.posedit.pos.is_none() => Ok(()),
            crate::SequenceVariant::Protein(v) => residue_range(
                v,
                &self.refs.reference(&v.ac, IdentifierType::ProteinAccession),
            )
            .map(|_| ()),
            // r.=, r.0, r.spl and their kind state no base.
            crate::SequenceVariant::Rna(r) if r.posedit.edit.is_whole_transcript() => Ok(()),
            crate::SequenceVariant::Rna(r) => self.check_stated(&self.r_as_transcript(r)?),
            crate::SequenceVariant::CisPhased(cis) if cis.coordinate_type() == "p" => {
                cis.members.iter().try_for_each(|m| self.check_stated(m))
            }
            crate::SequenceVariant::CisPhased(cis) => self.read_cis(cis).map(|_| ()),
        }
    }

    /// `c.=` names no bases, so there is nothing to check once its
    /// transcript is found to be a coding one.
    fn check_whole_identity(&self, v: &CVariant) -> Result<(), HgvsError> {
        if !self.has_cds(&v.ac)? {
            return Err(HgvsError::ValidationError(format!(
                "{v} describes a coding sequence, but {} has no CDS",
                v.ac
            )));
        }
        Ok(())
    }

    fn check_linear_stated<L: LinearVariant>(&self, v: &L) -> Result<(), HgvsError> {
        let edit = &v.posedit().edit;
        let pos = v
            .posedit()
            .pos
            .as_ref()
            .ok_or_else(|| HgvsError::ValidationError("Missing position".into()))?;
        let (start, end) = simple_interval_range(pos)?;
        let reference = self
            .refs
            .reference(v.ac(), IdentifierType::GenomicAccession);
        edit.check_stated(v.ac(), start, end, |s, e| reference.slice(s, e))
    }

    /// The bases a c. or n. variant states, on the record where it is exonic
    /// and on the genome, through the projection, where it is not.
    fn check_transcript_stated<V: TranscriptVariant>(&self, v: &V) -> Result<(), HgvsError> {
        let edit = &v.posedit().edit;
        let pos = v
            .posedit()
            .pos
            .as_ref()
            .ok_or_else(|| HgvsError::ValidationError("Missing position".into()))?;
        if has_intronic_offset(pos) {
            if !edit.states_bases() {
                return Ok(());
            }
            return self.tx_to_g(v, None).map(|_| ());
        }
        let am = TranscriptMapper::new(self.provider().get_transcript(v.ac(), None)?)?;
        let (n_start, n_end) = am.interval_to_n(pos)?;
        let start = checked_usize(n_start.0, "transcript start index")?;
        let end = checked_usize(n_end.0, "transcript end index")?;
        let reference = self
            .refs
            .reference(v.ac(), IdentifierType::TranscriptAccession);
        if reference.slice(start, end)?.len() != end - start {
            return Err(HgvsError::ValidationError(format!(
                "{} is shorter than the range {start}..{end} its variant names",
                v.ac()
            )));
        }
        edit.check_stated(v.ac(), start, end, |s, e| reference.slice(s, e))
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
            // (r.0, r.spl, r.=) has nothing to normalise.
            SV::Rna(r) if r.posedit.edit.is_whole_transcript() => SV::Rna(r),
            SV::Rna(r) => {
                let normalised = self.normalize_variant(self.r_as_transcript(&r)?)?;
                SV::Rna(self.tx_to_r(&normalised)?)
            }
            // A p. allele is not normalised, but what each member names is read.
            SV::CisPhased(cis) if cis.coordinate_type() == "p" => {
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
            // Members that meet are written as the one change they make.
            SV::CisPhased(cis) => self.normalize_cis(&cis)?,
            // A p. description is not normalised, but what it names is read.
            SV::Protein(v) => {
                let p = SV::Protein(v);
                self.check_stated(&p)?;
                p
            }
        })
    }

    fn normalize_linear<L: LinearVariant>(&self, mut v: L) -> Result<L, HgvsError> {
        let ac = v.ac().to_string();
        let Some(pos) = &mut v.posedit_mut().pos else {
            return Ok(v);
        };
        if is_uncertain_linear(pos) {
            // A bound known only to lie in a range names no one placement to
            // shift or reduce from. Leave the variant as written.
            return Ok(v);
        }
        let (start, end) = simple_interval_range(pos)?;
        let before = PlacedEdit::from_hgvs_range(start, end, v.posedit().edit.clone());
        let reference = self.refs.reference(&ac, IdentifierType::GenomicAccession);
        let after = normalize::normalize(&reference, before.clone(), &ExonBounds::unspliced())?;
        let posedit = v.posedit_mut();
        let pos = posedit.pos.as_mut().expect("checked above");
        let (s, e) = hgvs_positions(&after);
        pos.start = simple_position(s);
        pos.end = e.map(simple_position);
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
            // normalise against in transcript space. Leave the variant as
            // written, once what it states is found on the genome.
            self.check_transcript_stated(&v)?;
            return Ok(v);
        }
        if pos.uncertain || pos.start.uncertain || pos.end.as_ref().is_some_and(|e| e.uncertain) {
            return Ok(v);
        }
        let (start, end) = self.get_c_indices(pos, &am)?;
        let before = PlacedEdit::from_hgvs_range(start, end, v.posedit().edit.clone());
        let reference = self
            .refs
            .reference(&ac, IdentifierType::TranscriptAccession);
        let after = normalize::normalize(&reference, before.clone(), &am.exon_bounds()?)?;
        let written =
            matches!(v.posedit().edit, crate::edits::NaEdit::Repeat { .. }).then(|| v.clone());
        let posedit = v.posedit_mut();
        let pos = posedit.pos.as_mut().expect("checked above");
        // Positions from indices, so that for c. the c.0 gap and the CDS
        // anchors come out right when a shift crosses the CDS bounds.
        let (s, e) = hgvs_positions(&after);
        pos.start = V::position_from_index(&am, s as i32)?;
        pos.end = e
            .map(|last| V::position_from_index(&am, last as i32))
            .transpose()?;
        posedit.edit = after.edit;
        match written {
            Some(written)
                if v.posedit().edit.is_identity() && genome_may_count_more(&am, s, e)? =>
            {
                self.repeat_statement_or_as_written(written, v)
            }
            _ => Ok(v),
        }
    }

    /// `statement`, the identity a transcript repeat that makes the reference
    /// normalised to, or the repeat as `written` where the genome counts a
    /// different run from it: on the minus strand, copies before the written
    /// start; on either strand, copies the intron continues an exon's run
    /// with. The repeat then states bases on the genome that no identity over
    /// the transcript's run states. Telling the two apart reads the genome.
    fn repeat_statement_or_as_written<V: TranscriptVariant>(
        &self,
        written: V,
        statement: V,
    ) -> Result<V, HgvsError> {
        let stated = |v: &V| -> Result<(usize, usize), HgvsError> {
            let g = self.tx_to_g(v, None)?;
            let pos = g
                .posedit
                .pos
                .as_ref()
                .ok_or_else(|| HgvsError::ValidationError("Missing position".into()))?;
            simple_interval_range(pos)
        };
        Ok(if stated(&written)? == stated(&statement)? {
            statement
        } else {
            written
        })
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
    /// equivalence all derive from. See [`CanonicalAllele`]. A variant whose
    /// position is uncertain names no one change, so it has none.
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
        if !var.has_exact_position() {
            return Err(HgvsError::UnsupportedOperation(format!(
                "{var} has an uncertain position, so it names no one change to make an allele of"
            )));
        }
        if let crate::SequenceVariant::Protein(vp) = var {
            let reference = self
                .refs
                .reference(&vp.ac, IdentifierType::ProteinAccession);
            let (start, end) = residue_range(vp, &reference)?;
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
                | crate::edits::NaEdit::Statement(_)
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
        if let (VrsBound::Exact(s), VrsBound::Exact(e)) = (start, end) {
            let genome = self.refs.reference(&g.ac, IdentifierType::GenomicAccession);
            g.posedit
                .edit
                .check_stated(&g.ac, s, e, |a, b| genome.slice(a, b))?;
        }
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
    ///
    /// An allele of nucleotide members is rendered as
    /// [`to_vrs_variation`](Self::to_vrs_variation) renders it: projected to
    /// the genome, with the members that meet there as one Allele, in order
    /// along the genome. One that is a single change there is that change's
    /// object, not a block.
    pub fn to_vrs_cis_phased(
        &self,
        var: &crate::SequenceVariant,
    ) -> Result<VrsCisPhasedBlock, HgvsError> {
        if !matches!(var, crate::SequenceVariant::CisPhased(_)) {
            return Err(HgvsError::UnsupportedOperation(format!(
                "Only an allele in cis, ac:c.[a;b], renders as a CisPhasedBlock, not {var}"
            )));
        }
        match self.to_vrs_variation(var)? {
            VrsVariation::CisPhasedBlock(block) => Ok(block),
            _ => Err(HgvsError::UnsupportedOperation(format!(
                "{var} is one change on the genome, rendered as one VRS object rather than a \
                 CisPhasedBlock; use to_vrs_variation"
            ))),
        }
    }

    /// The block of `cis`'s members, with `var`, the allele as written, as
    /// its expression.
    fn cis_phased_block(
        &self,
        cis: &crate::structs::CisPhasedVariant,
        var: &crate::SequenceVariant,
    ) -> Result<VrsCisPhasedBlock, HgvsError> {
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
    ///
    /// An allele of nucleotide members is read where it is written, refused
    /// if contradictory or ambiguous, projected to the genome and normalised
    /// there ([`normalize_variant`](Self::normalize_variant) on the `g.`
    /// allele, so no codon groups members): members that overlap, touch or
    /// would shift into each other are one Allele, each fully justified. One
    /// that is a single change there is that change's object, carrying the
    /// allele as written as its expression.
    pub fn to_vrs_variation(
        &self,
        var: &crate::SequenceVariant,
    ) -> Result<VrsVariation, HgvsError> {
        if let crate::SequenceVariant::CisPhased(cis) = var {
            if cis.coordinate_type() == "p" {
                return Ok(VrsVariation::CisPhasedBlock(
                    self.cis_phased_block(cis, var)?,
                ));
            }
            // Read where written, so that a contradictory or ambiguous allele
            // is refused before a block could conflate its readings; then
            // merged on the genome, where no codon groups members, so every
            // spelling of one genomic allele splits into Alleles alike.
            self.read_cis(cis)?;
            return match self.normalize_cis_on_genome(cis)? {
                crate::SequenceVariant::CisPhased(normalised) => Ok(VrsVariation::CisPhasedBlock(
                    self.cis_phased_block(&normalised, var)?,
                )),
                one => {
                    let expressions = vec![crate::vrs::VrsExpression {
                        syntax: format!("hgvs.{}", var.coordinate_type()),
                        value: var.to_string(),
                    }];
                    let mut rendered = self.to_vrs_variation(&one)?;
                    match &mut rendered {
                        VrsVariation::Allele(a) => a.expressions = expressions,
                        VrsVariation::CopyNumberCount(c) => c.expressions = expressions,
                        VrsVariation::CopyNumberChange(c) => c.expressions = expressions,
                        VrsVariation::CisPhasedBlock(b) => b.expressions = expressions,
                    }
                    Ok(rendered)
                }
            };
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
        let past_any_sequence =
            || HgvsError::ValidationError(format!("SPDI {spdi:?} runs past any sequence"));
        let end = if !del.is_empty() && del.chars().all(|c| c.is_ascii_digit()) {
            let length: usize = del.parse().map_err(|_| past_any_sequence())?;
            start.checked_add(length).ok_or_else(past_any_sequence)?
        } else {
            let end = start.checked_add(del.len()).ok_or_else(past_any_sequence)?;
            let actual = self.refs.reference(ac, kind).slice(start, end)?;
            if !actual.eq_ignore_ascii_case(del) {
                return Err(HgvsError::ReferenceMismatch {
                    accession: ac.to_string(),
                    expected: del.to_string(),
                    found: actual,
                    start,
                    end,
                });
            }
            end
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
                    bounded_interval(start, end)?,
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
            let next = start.checked_add(1);
            let after = match next {
                Some(next) => self.refs.reference(ac, kind).slice(start, next)?,
                None => String::new(),
            };
            if after.is_empty() {
                return Err(HgvsError::ValidationError(format!(
                    "{ac} has no base after interbase {start}"
                )));
            }
            let edit = NaEdit::InsLength {
                min,
                max,
                uncertain,
            };
            return Ok(crate::SequenceVariant::Genomic(g_variant(
                ac,
                interval(start - 1, start + 1)?,
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
            interval(start, end)?,
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
        let past_the_end = start > 0 && reference.base(start - 1)?.is_none();
        if ref_.len() != end - start || past_the_end {
            return Err(HgvsError::ValidationError(format!(
                "{ac} is shorter than position {end}"
            )));
        }
        let (edit, s, e) = normalize::minimal_edit(&reference, start, end, &ref_, &alt)?;
        let placed = normalize::normalize(
            &reference,
            PlacedEdit::from_hgvs_range(s, e, edit),
            &ExonBounds::unspliced(),
        )?;
        let (s, e) = placed.hgvs_range();
        Ok(match kind {
            IdentifierType::ProteinAccession => {
                crate::SequenceVariant::Protein(protein_variant(ac, &reference, s, e, placed.edit)?)
            }
            _ => crate::SequenceVariant::Genomic(g_variant(ac, interval(s, e)?, placed.edit)),
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
    /// Either way the bases it states are checked, so a mismatch is a
    /// [`HgvsError::ReferenceMismatch`]. `None` for protein variants and
    /// cis alleles.
    pub fn as_genomic(&self, var: &crate::SequenceVariant) -> Option<Result<GVariant, HgvsError>> {
        use crate::SequenceVariant as SV;
        Some(match var {
            SV::Genomic(v) => self.check_linear_stated(v).map(|()| v.clone()),
            SV::Mitochondrial(v) => self.check_linear_stated(v).map(|()| v.to_genomic()),
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
                if !within_one_exon(&am, start.0, end.0) {
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
