use crate::allele::CanonicalAllele;
use crate::data::{DataProvider, IdentifierKind, IdentifierType, TranscriptData, TranscriptSearch};
use crate::error::HgvsError;
use crate::normalize::{self, PlacedEdit};
use crate::reference::ReferenceStore;
use crate::structs::Variant;
use crate::structs::{
    BaseOffsetInterval, BaseOffsetPosition, CVariant, GVariant, GenomicPos, LinearVariant,
    NVariant, PVariant, SimpleInterval, SimplePosition, TranscriptVariant,
};
use crate::transcript_mapper::TranscriptMapper;
use crate::vrs::VrsAllele;

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
/// be rewritten. A filled-in deletion or duplication reference alone does not.
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

    /// Transforms a coding cDNA variant (`c.`) to a protein variant (`p.`).
    pub fn c_to_p(
        &self,
        var_c: &CVariant,
        protein_ac: Option<&str>,
    ) -> Result<PVariant, HgvsError> {
        let transcript_ac = &var_c.ac;
        let pro_ac_str = if let Some(ac) = protein_ac {
            ac.to_string()
        } else {
            self.hdp
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
                .clone()
        };

        // Handle intronic variants by returning p.?
        if let Some(pos) = &var_c.posedit.pos {
            let has_offset = pos.start.offset.is_some_and(|o| o.0 != 0)
                || pos
                    .end
                    .as_ref()
                    .is_some_and(|e| e.offset.is_some_and(|o| o.0 != 0));

            if has_offset {
                return Ok(PVariant {
                    ac: pro_ac_str,
                    gene: var_c.gene.clone(),
                    posedit: crate::structs::PosEdit {
                        pos: None,
                        edit: crate::edits::AaEdit::Special {
                            value: "?".to_string(),
                            uncertain: false,
                        },
                        uncertain: false,
                        predicted: false,
                    },
                });
            }

            // 5'UTR variants (negative c. position) → p.?
            let start_base = pos.start.base.0;
            if pos.start.anchor == crate::coords::Anchor::CdsStart && start_base < 0 {
                return Ok(PVariant {
                    ac: pro_ac_str,
                    gene: var_c.gene.clone(),
                    posedit: crate::structs::PosEdit {
                        pos: None,
                        edit: crate::edits::AaEdit::Special {
                            value: "?".to_string(),
                            uncertain: false,
                        },
                        uncertain: false,
                        predicted: false,
                    },
                });
            }
        }

        let transcript = self.hdp.get_transcript(transcript_ac, None)?;
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
        let change = crate::protein::CodingChange {
            coding: &ref_seq[cds_start_idx..],
            cds_len: cds_end_idx + 1 - cds_start_idx,
            edit: crate::edits::ResolvedEdit {
                start: rel(resolved.start)?,
                end: rel(resolved.end)?,
                ref_: resolved.ref_,
                alt: resolved.alt,
            },
            protein_ac: pro_ac_str,
        };
        let mut var_p = crate::protein::describe(&change)?;
        var_p.posedit.predicted = true;
        Ok(var_p)
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
            _ => Err(HgvsError::UnsupportedOperation(
                "Validation not implemented for this variant type".into(),
            )),
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
        let g = self.as_genomic(var).ok_or_else(|| {
            HgvsError::UnsupportedOperation(
                "Canonical alleles exist for genomic, mitochondrial, coding and non-coding variants only".into(),
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
        let allele = self.canonical_allele(var)?;
        let refget = self
            .refs
            .reference(&allele.accession, IdentifierType::GenomicAccession)
            .refget_accession()?;
        let syntax = format!("hgvs.{}", var.coordinate_type());
        Ok(VrsAllele::new(
            &allele,
            &refget,
            Some((&syntax, &var.to_string())),
        ))
    }

    pub fn to_spdi(
        &self,
        var: &crate::SequenceVariant,
        unambiguous: bool,
    ) -> Result<String, HgvsError> {
        if unambiguous {
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
            SV::Protein(_) | SV::Rna(_) => return None,
        })
    }
}
