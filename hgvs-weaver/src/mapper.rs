use crate::allele::CanonicalAllele;
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
use crate::vrs::{VrsAllele, VrsBound, VrsMolecule, VrsState};

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

fn apply_strand_complement(
    edit: crate::edits::NaEdit,
    strand: crate::data::Strand,
) -> crate::edits::NaEdit {
    if strand == crate::data::Strand::Minus {
        edit.map_sequence(|s| crate::utils::reverse_complement(s))
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

/// `g.(a_b)_(c_d)del` from VRS bounds: the inverse of `uncertain_bounds`.
fn imprecise_deletion(ac: &str, start: VrsBound, end: VrsBound) -> GVariant {
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
    let pos = SimpleInterval {
        start: position(start, |n| n as i32 + 1),
        end: (!single).then(|| position(end, |n| n as i32)),
        uncertain: false,
    };
    GVariant::from_parts(
        ac.to_string(),
        None,
        crate::structs::PosEdit {
            pos: Some(pos),
            edit: crate::edits::NaEdit::Del {
                ref_: None,
                uncertain: false,
            },
            uncertain: false,
            predicted: false,
        },
    )
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
    /// Data provider used to retrieve transcript and sequence information.
    pub hdp: &'a dyn DataProvider,
    /// Cached, random-access view of every sequence the provider serves.
    pub refs: ReferenceStore<'a>,
}

impl<'a> VariantMapper<'a> {
    /// Creates a new `VariantMapper` with the given data provider.
    pub fn new(hdp: &'a dyn DataProvider) -> Self {
        VariantMapper {
            hdp,
            refs: ReferenceStore::new(hdp),
        }
    }

    /// Transforms a genomic variant (`g.`) to a coding cDNA variant (`c.`).
    pub fn g_to_c(&self, var_g: &GVariant, transcript_ac: &str) -> Result<CVariant, HgvsError> {
        let transcript = self.hdp.get_transcript(transcript_ac, Some(&var_g.ac))?;
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
        let g_start_0 = pos.start.base.to_index();
        let (mut n_pos, mut offset) = am.g_to_n(g_start_0)?;

        if let Some(end_g_simple) = &pos.end {
            let g_end_0 = end_g_simple.base.to_index();
            let (mut n_pos_e, mut offset_e) = am.g_to_n(g_end_0)?;

            if n_pos.0 > n_pos_e.0 {
                std::mem::swap(&mut n_pos, &mut n_pos_e);
                std::mem::swap(&mut offset, &mut offset_e);
            }

            let (c_pos_index, c_offset, anchor) = am.n_to_c(n_pos)?;
            let (c_pos_e_index, c_offset_e, anchor_e) = am.n_to_c(n_pos_e)?;

            let pos_c =
                make_base_offset_position(c_pos_index.to_hgvs(), c_offset.0 + offset.0, anchor);

            let pos_c_e = make_base_offset_position(
                c_pos_e_index.to_hgvs(),
                c_offset_e.0 + offset_e.0,
                anchor_e,
            );

            let edit = apply_strand_complement(var_g.posedit.edit.clone(), am.transcript.strand);

            return Ok(CVariant {
                ac: transcript_ac.to_string(),
                gene: var_g.gene.clone(),
                posedit: crate::structs::PosEdit {
                    pos: Some(crate::structs::BaseOffsetInterval {
                        start: pos_c,
                        end: Some(pos_c_e),
                        uncertain: false,
                    }),
                    edit,
                    uncertain: var_g.posedit.uncertain,
                    predicted: var_g.posedit.predicted,
                },
            });
        }

        let (c_pos_index, c_offset, anchor) = am.n_to_c(n_pos)?;
        let pos_c = make_base_offset_position(c_pos_index.to_hgvs(), c_offset.0 + offset.0, anchor);

        let edit = apply_strand_complement(var_g.posedit.edit.clone(), am.transcript.strand);

        Ok(CVariant {
            ac: transcript_ac.to_string(),
            gene: var_g.gene.clone(),
            posedit: crate::structs::PosEdit {
                pos: Some(crate::structs::BaseOffsetInterval {
                    start: pos_c,
                    end: None,
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
        let transcript = self.hdp.get_transcript(var_c.ac(), reference_ac)?;
        let am = TranscriptMapper::new(transcript)?;

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
        let n_pos = am.c_to_n(pos.start.base.to_index(), pos.start.anchor)?;
        let g_pos = am.n_to_g(
            n_pos,
            pos.start
                .offset
                .unwrap_or(crate::structs::IntronicOffset(0)),
        )?;

        if let Some(end_c) = &pos.end {
            let n_pos_e = am.c_to_n(end_c.base.to_index(), end_c.anchor)?;
            let g_pos_e = am.n_to_g(
                n_pos_e,
                end_c.offset.unwrap_or(crate::structs::IntronicOffset(0)),
            )?;
            let mut pos_g = make_simple_position(g_pos.to_hgvs());
            let mut pos_g_e = make_simple_position(g_pos_e.to_hgvs());

            if pos_g.base.0 > pos_g_e.base.0 {
                std::mem::swap(&mut pos_g, &mut pos_g_e);
            }

            let edit = apply_strand_complement(var_c.posedit().edit.clone(), am.transcript.strand);

            return Ok(GVariant {
                ac: reference_ac
                    .unwrap_or(am.transcript.reference_accession.as_str())
                    .to_string(),
                gene: var_c.gene().map(str::to_string),
                posedit: crate::structs::PosEdit {
                    pos: Some(crate::structs::SimpleInterval {
                        start: pos_g,
                        end: Some(pos_g_e),
                        uncertain: false,
                    }),
                    edit,
                    uncertain: var_c.posedit().uncertain,
                    predicted: var_c.posedit().predicted,
                },
            });
        }

        let pos_g = make_simple_position(g_pos.to_hgvs());
        let edit = apply_strand_complement(var_c.posedit().edit.clone(), am.transcript.strand);

        Ok(GVariant {
            ac: reference_ac
                .unwrap_or(am.transcript.reference_accession.as_str())
                .to_string(),
            gene: var_c.gene().map(str::to_string),
            posedit: crate::structs::PosEdit {
                pos: Some(crate::structs::SimpleInterval {
                    start: pos_g,
                    end: None,
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
            .hdp
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
        let transcript_ac = &var_c.ac;
        let pro_ac_str = self.protein_accession(transcript_ac, protein_ac)?;

        let transcript = self.hdp.get_transcript(transcript_ac, None)?;
        let unknown = |pos: Option<crate::structs::AaInterval>, value: &str| PVariant {
            ac: pro_ac_str.clone(),
            gene: var_c.gene.clone(),
            posedit: crate::structs::PosEdit {
                pos,
                edit: crate::edits::AaEdit::Special {
                    value: value.to_string(),
                    uncertain: false,
                },
                uncertain: false,
                predicted: false,
            },
        };
        if let Some(pos) = &var_c.posedit.pos {
            // Intronic: the protein consequence cannot be predicted.
            let has_offset = pos.start.offset.is_some_and(|o| o.0 != 0)
                || pos
                    .end
                    .as_ref()
                    .is_some_and(|e| e.offset.is_some_and(|o| o.0 != 0));
            if has_offset {
                return Ok(CodingOutcome::Statement(unknown(None, "?")));
            }

            // An edit that starts in the 5'UTR. Deleting the whole CDS predicts
            // no protein (p.0?); reaching into the CDS disrupts the start codon
            // (p.Met1?); staying upstream says nothing about the protein (p.?).
            use crate::coords::Anchor;
            if pos.start.anchor == Anchor::CdsStart && pos.start.base.0 < 0 {
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
                return Ok(CodingOutcome::Statement(if covers_cds && deletes {
                    unknown(None, "0?")
                } else if reaches_cds {
                    let met1 = crate::structs::AaInterval {
                        start: crate::structs::AAPosition {
                            base: crate::structs::ProteinPos(0).to_hgvs(),
                            aa: "Met".to_string(),
                            uncertain: false,
                        },
                        end: None,
                        uncertain: false,
                    };
                    unknown(Some(met1), "?")
                } else {
                    unknown(None, "?")
                }));
            }
        }

        let ref_seq = self
            .refs
            .reference(transcript_ac, IdentifierType::TranscriptAccession)
            .whole()?;

        let cds_start_tx = transcript
            .cds_start_index
            .ok_or_else(|| HgvsError::ValidationError("Missing CDS start".into()))?;
        let cds_end_tx = transcript
            .cds_end_index
            .ok_or_else(|| HgvsError::ValidationError("Missing CDS end".into()))?;
        let cds_start_idx = checked_usize(cds_start_tx.0, "CDS start")?;
        let cds_end_idx = checked_usize(cds_end_tx.0, "CDS end")?;

        if ref_seq.len() < cds_end_idx {
            return Err(HgvsError::ValidationError(format!(
                "Transcript sequence too short (len={}, expected at least {})",
                ref_seq.len(),
                cds_end_idx
            )));
        }

        if cds_start_idx > ref_seq.len() {
            return Err(HgvsError::ValidationError(format!(
                "CDS start {} out of sequence bounds {}",
                cds_start_idx,
                ref_seq.len()
            )));
        }

        let pos = var_c
            .posedit
            .pos
            .as_ref()
            .ok_or_else(|| HgvsError::ValidationError("Missing position".into()))?;
        let am = TranscriptMapper::new(transcript)?;
        let (n_start, n_end) = am.interval_to_n(pos)?;
        if n_start.0 < 0 {
            return Err(HgvsError::ValidationError(format!(
                "Position {} before transcript start",
                n_start.0
            )));
        }
        let (start_idx, end_idx) = (n_start.0 as usize, n_end.0 as usize);
        if end_idx > ref_seq.len() {
            let first_bad = if start_idx >= ref_seq.len() {
                start_idx
            } else {
                end_idx
            };
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
        // The bases the variant says are there must be there.
        if let Some(stated) = var_c.posedit.edit.stated_ref() {
            let actual = window(start_idx, end_idx);
            if actual != stated {
                return Err(HgvsError::TranscriptMismatch {
                    expected: stated.to_string(),
                    found: actual.to_string(),
                    start: start_idx,
                    end: end_idx,
                });
            }
        }
        let resolved = var_c
            .posedit
            .edit
            .resolve_with(start_idx, end_idx, |s, e| Ok(window(s, e).to_string()))?;
        let rel = |i: usize| {
            i.checked_sub(cds_start_idx).ok_or_else(|| {
                HgvsError::ValidationError(format!("Position {} before the CDS start", i))
            })
        };
        Ok(CodingOutcome::Change(crate::protein::CodingChange {
            coding: ref_seq[cds_start_idx..].to_string(),
            cds_len: cds_end_idx + 1 - cds_start_idx,
            edit: crate::edits::ResolvedEdit {
                start: rel(resolved.start)?,
                end: rel(resolved.end)?,
                ref_: resolved.ref_,
                alt: resolved.alt,
            },
            protein_ac: pro_ac_str,
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
        let ac = &var_p.ac;

        // Extract position and edit
        let pos = var_p
            .posedit
            .pos
            .as_ref()
            .ok_or_else(|| HgvsError::ValidationError("Missing protein position".into()))?;

        let edit = &var_p.posedit.edit;

        // Only handle Subst for now
        let alt_aa_str = match edit {
            crate::edits::AaEdit::Subst { alt, .. } => alt.clone(),
            _ => {
                return Err(HgvsError::UnsupportedOperation(
                    "p_to_c only supports single amino acid substitutions".into(),
                ));
            }
        };

        let raw_pos = pos.start.base.0;
        if raw_pos <= 0 {
            return Err(HgvsError::ValidationError(format!(
                "Protein position {raw_pos} is not valid (must be >= 1)"
            )));
        }
        let aa_pos = checked_usize(raw_pos, "protein position")?; // 1-based
        let ref_aa_str = &pos.start.aa;

        // Convert 3-letter or 1-letter to single char
        let ref_aa_1 = crate::utils::aa3_to_aa1(ref_aa_str);
        let alt_aa_1 = crate::utils::aa3_to_aa1(&alt_aa_str);

        let _ref_aa = ref_aa_1
            .chars()
            .next()
            .ok_or_else(|| HgvsError::ValidationError("Invalid reference AA".into()))?;
        let alt_aa = alt_aa_1
            .chars()
            .next()
            .ok_or_else(|| HgvsError::ValidationError("Invalid alternate AA".into()))?;

        // Resolve transcript accession
        let tx_ac = if let Some(ta) = transcript_ac {
            ta.to_string()
        } else {
            self.hdp
                .get_symbol_accessions(ac, IdentifierKind::Protein, IdentifierKind::Transcript)?
                .first()
                .ok_or_else(|| {
                    HgvsError::ValidationError(format!("No transcript accession found for {}", ac))
                })?
                .1
                .clone()
        };

        // Get transcript
        let transcript = self.hdp.get_transcript(&tx_ac, None)?;
        let cds_start = transcript
            .cds_start_index
            .ok_or_else(|| HgvsError::ValidationError("No CDS start for transcript".into()))?
            .0 as usize;

        // Get transcript sequence
        let tx_seq_str = self
            .refs
            .reference(&tx_ac, IdentifierType::TranscriptAccession)
            .whole()?;

        // Calculate codon position (0-based in transcript)
        let codon_start = cds_start + (aa_pos - 1) * 3;
        let codon_end = codon_start + 3;
        if codon_end > tx_seq_str.len() {
            return Err(HgvsError::ValidationError(format!(
                "Codon position {}-{} out of range for sequence length {}",
                codon_start,
                codon_end,
                tx_seq_str.len()
            )));
        }

        let ref_codon = tx_seq_str[codon_start..codon_end].to_uppercase();
        let ref_codon_bytes = ref_codon.as_bytes();

        // Get all codons for the target AA
        let alt_codons = crate::utils::codons_for_aa(alt_aa);
        if alt_codons.is_empty() {
            return Err(HgvsError::ValidationError(format!(
                "No codons found for amino acid '{}'",
                alt_aa
            )));
        }

        // Score each candidate by nucleotide differences, pick minimum
        let mut scored: Vec<(usize, &str)> = alt_codons
            .iter()
            .map(|codon| {
                let diffs = codon
                    .as_bytes()
                    .iter()
                    .zip(ref_codon_bytes.iter())
                    .filter(|(a, b)| a != b)
                    .count();
                (diffs, *codon)
            })
            .collect();
        scored.sort();

        let best_diffs = scored[0].0;
        let best_codons: Vec<&str> = scored
            .iter()
            .filter(|(d, _)| *d == best_diffs)
            .map(|(_, c)| *c)
            .collect();
        let is_unique = best_codons.len() == 1;
        let best_codon = best_codons[0];

        // Build the c. variant for the changed nucleotides
        let mut changes: Vec<(usize, u8, u8)> = Vec::new();
        for i in 0..3 {
            if ref_codon_bytes[i] != best_codon.as_bytes()[i] {
                // c. position is 1-based from CDS start
                let c_pos = (aa_pos - 1) * 3 + i + 1;
                changes.push((c_pos, ref_codon_bytes[i], best_codon.as_bytes()[i]));
            }
        }

        let na_edit = if changes.len() == 1 {
            crate::edits::NaEdit::RefAlt {
                ref_: Some(String::from(changes[0].1 as char)),
                alt: Some(String::from(changes[0].2 as char)),
                uncertain: false,
            }
        } else {
            // Multiple nucleotide changes → delins
            let ref_nts: String = changes.iter().map(|(_, r, _)| *r as char).collect();
            let alt_nts: String = changes.iter().map(|(_, _, a)| *a as char).collect();
            crate::edits::NaEdit::RefAlt {
                ref_: Some(ref_nts),
                alt: Some(alt_nts),
                uncertain: false,
            }
        };

        let start_c_pos = changes.first().map(|(p, _, _)| *p).unwrap_or(1);
        let end_c_pos = changes.last().map(|(p, _, _)| *p).unwrap_or(start_c_pos);

        use crate::coords::{Anchor, HgvsTranscriptPos};
        use crate::structs::BaseOffsetPosition;

        let c_pos = BaseOffsetInterval {
            start: BaseOffsetPosition {
                base: HgvsTranscriptPos(start_c_pos as i32),
                offset: None,
                anchor: Anchor::CdsStart,
                uncertain: false,
            },
            end: if start_c_pos == end_c_pos {
                None
            } else {
                Some(BaseOffsetPosition {
                    base: HgvsTranscriptPos(end_c_pos as i32),
                    offset: None,
                    anchor: Anchor::CdsStart,
                    uncertain: false,
                })
            },
            uncertain: false,
        };

        let c_variant = CVariant {
            ac: tx_ac.clone(),
            gene: var_p.gene.clone(),
            posedit: crate::structs::PosEdit {
                pos: Some(c_pos),
                edit: na_edit,
                uncertain: false,
                predicted: false,
            },
        };

        Ok((c_variant, is_unique))
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
        let transcript = self.hdp.get_transcript(v.ac(), None)?;
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
        let transcript = self.hdp.get_transcript(&ac, None)?;
        let Some(pos) = &v.posedit().pos else {
            return Ok(v);
        };
        if has_intronic_offset(pos) {
            // An intronic base has no transcript index; there is nothing to
            // normalise against in transcript space. Leave the variant as written.
            return Ok(v);
        }
        let (start, end) = self.get_c_indices(pos, &transcript)?;
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
            let am = TranscriptMapper::new(transcript)?;
            let (s, e) = hgvs_positions(&before, &after, pos.end.is_some());
            pos.start = V::position_from_index(&am, s as i32)?;
            pos.end = e
                .map(|last| V::position_from_index(&am, last as i32))
                .transpose()?;
        }
        posedit.edit = after.edit;
        Ok(v)
    }

    pub fn get_c_indices(
        &self,
        pos: &BaseOffsetInterval,
        transcript: &TranscriptData,
    ) -> Result<(usize, usize), HgvsError> {
        let am = TranscriptMapper::new(transcript.clone())?;
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
            let resolved = vp.posedit.edit.resolve(&reference, start, end)?;
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
    pub fn to_vrs(&self, var: &crate::SequenceVariant) -> Result<VrsAllele, HgvsError> {
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
                if !matches!(g.posedit.edit, crate::edits::NaEdit::Del { .. }) {
                    return Err(HgvsError::UnsupportedOperation(format!(
                        "Only a deletion can have uncertain breakpoints in VRS, not {:?}",
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

    /// The variant a GA4GH VRS 2.0 Allele (as JSON) names, written in HGVS on
    /// its own sequence and 3'-normalised. The sequence is identified by its
    /// refget accession: `accession` names it when given, else the provider's
    /// `get_accession_for_refget` must; the digest is checked against the
    /// sequence either way. Range bounds are accepted for a deletion, which
    /// comes back as `g.(a_b)_(c_d)del`.
    pub fn from_vrs(
        &self,
        json: &str,
        accession: Option<&str>,
    ) -> Result<crate::SequenceVariant, HgvsError> {
        let allele = VrsAllele::from_json(json)?;
        let refget = allele.location.sequence_reference.refget_accession.as_str();
        let ac = match accession {
            Some(a) => a.to_string(),
            None => self.hdp.get_accession_for_refget(refget)?.ok_or_else(|| {
                HgvsError::DataProviderError(format!(
                    "No accession is known for {refget}; pass one, or implement DataProvider::get_accession_for_refget"
                ))
            })?,
        };
        let sr = &allele.location.sequence_reference;
        let kind = if sr.residue_alphabet == "aa" || sr.molecule_type == "protein" {
            IdentifierType::ProteinAccession
        } else if sr.residue_alphabet == "na" {
            IdentifierType::GenomicAccession
        } else {
            self.sequence_kind(&ac)?
        };
        let reference = self.refs.reference(&ac, kind);
        let actual = reference.refget_accession()?;
        if actual != refget {
            return Err(HgvsError::ValidationError(format!(
                "{refget} is not the refget accession of {ac}, which is {actual}"
            )));
        }
        match (allele.location.start, allele.location.end) {
            (VrsBound::Exact(start), VrsBound::Exact(end)) => {
                if end < start {
                    return Err(HgvsError::ValidationError(format!(
                        "Location end {end} is before start {start}"
                    )));
                }
                let alt = match &allele.state {
                    VrsState::Literal { sequence, .. } => sequence.clone(),
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
                Ok(crate::SequenceVariant::Genomic(imprecise_deletion(
                    &ac, start, end,
                )))
            }
        }
    }

    /// Whether `ac` names a protein or a nucleotide sequence.
    fn sequence_kind(&self, ac: &str) -> Result<IdentifierType, HgvsError> {
        Ok(match self.hdp.get_identifier_type(ac)? {
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
        use crate::edits::NaEdit;
        let reference = self.refs.reference(ac, kind);
        let ref_ = reference.slice(start, end)?;
        if ref_.len() != end - start {
            return Err(HgvsError::ValidationError(format!(
                "{ac} is shorter than position {end}"
            )));
        }
        // Trim what the allele leaves unchanged, prefix first: a fully
        // justified insertion or deletion then sits at the 3' end of its run,
        // which is where HGVS writes it.
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
        let r = r[..r.len() - suffix].to_string();
        let a = a[..a.len() - suffix].to_string();
        let at = start + prefix;
        let (edit, s, e) = if r.is_empty() && a.is_empty() {
            let edit = NaEdit::RefAlt {
                ref_: None,
                alt: None,
                uncertain: false,
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
                uncertain: false,
            };
            (edit, 0, 1)
        } else if r.is_empty() {
            let edit = NaEdit::Ins {
                alt: Some(a),
                uncertain: false,
            };
            (edit, at - 1, at + 1)
        } else if a.is_empty() {
            let edit = NaEdit::Del {
                ref_: None,
                uncertain: false,
            };
            (edit, at, at + r.len())
        } else if r.len() == 1 && a.len() == 1 {
            let edit = NaEdit::RefAlt {
                ref_: Some(r),
                alt: Some(a),
                uncertain: false,
            };
            (edit, at, at + 1)
        } else if kind != IdentifierType::ProteinAccession
            && r.len() > 1
            && a == crate::utils::reverse_complement(&r)
        {
            let edit = NaEdit::Inv {
                ref_: None,
                uncertain: false,
            };
            (edit, at, at + r.len())
        } else {
            // A delins, written without the deleted bases.
            let edit = NaEdit::RefAlt {
                ref_: Some(String::new()),
                alt: Some(a),
                uncertain: false,
            };
            (edit, at, at + r.len())
        };
        let placed = normalize::normalize(&reference, PlacedEdit::from_hgvs_range(s, e, edit))?;
        let (s, e) = placed.hgvs_range();
        let edit = placed.edit;
        Ok(match kind {
            IdentifierType::ProteinAccession => {
                crate::SequenceVariant::Protein(protein_variant(ac, &reference, s, e, edit)?)
            }
            _ => crate::SequenceVariant::Genomic(GVariant::from_parts(
                ac.to_string(),
                None,
                crate::structs::PosEdit {
                    pos: Some(interval(s, e)),
                    edit,
                    uncertain: false,
                    predicted: false,
                },
            )),
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
            SV::Protein(_) => return None,
        })
    }

    // --- r.: the transcript in RNA letters ---
    //
    // HGVS numbers r. positions like c. on a coding transcript and like n. on
    // a non-coding one, and writes bases in lowercase with u for T. Every r.
    // operation is therefore a conversion to the c. or n. spelling and back.

    /// Whether `ac` has a CDS, which decides whether r. is numbered like c. or n.
    fn has_cds(&self, ac: &str) -> Result<bool, HgvsError> {
        let t = self.hdp.get_transcript(ac, None)?;
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
            let transcript = self.hdp.get_transcript(&r.ac, None)?;
            let exons = transcript.exons.clone();
            let am = TranscriptMapper::new(transcript)?;
            // Intronic positions do not resolve to a transcript range; they
            // name the genome directly and need no guard.
            if let Ok((start, end)) = am.interval_to_n(pos) {
                let within_one_exon = exons
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
