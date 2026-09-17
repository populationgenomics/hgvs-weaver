use crate::error::HgvsError;
use crate::structs::{CVariant, NaEdit, ProteinPos, TranscriptPos};

/// Represents the data for a transcript with a variant applied.
pub struct AltTranscriptData {
    pub transcript_sequence: String,
    pub aa_sequence: String,
    pub cds_start_index: TranscriptPos,
    pub cds_end_index: TranscriptPos,
    pub protein_accession: String,
    pub is_frameshift: bool,
    /// The index of the first affected amino acid.
    pub variant_start_aa: Option<ProteinPos>,
    pub frameshift_start: Option<ProteinPos>,
    pub is_substitution: bool,
    pub variant_start_idx: usize,
    pub variant_end_idx: usize,
    pub is_ambiguous: bool,
    /// The original cDNA variant.
    pub c_variant: CVariant,
}

pub struct AltSeqBuilder<'a> {
    pub var_c: &'a CVariant,
    /// Resolves the variant's c. positions to transcript indices.
    pub mapper: &'a crate::transcript_mapper::TranscriptMapper,
    pub transcript_sequence: &'a str,
    pub cds_start_index: TranscriptPos,
    pub cds_end_index: TranscriptPos,
    pub protein_accession: String,
}

/// `seq` with `[start, end)` replaced by `insert`. Indices past the end clamp.
fn splice(seq: &str, start: usize, end: usize, insert: &str) -> String {
    let start = start.min(seq.len());
    let end = end.min(seq.len()).max(start);
    let mut out = String::with_capacity(seq.len() + insert.len());
    out.push_str(&seq[..start]);
    out.push_str(insert);
    out.push_str(&seq[end..]);
    out
}

/// `seq[start..end]`, clamped to the sequence.
fn window(seq: &str, start: usize, end: usize) -> &str {
    let start = start.min(seq.len());
    let end = end.min(seq.len()).max(start);
    &seq[start..end]
}

impl<'a> AltSeqBuilder<'a> {
    pub fn build_altseq(&self) -> Result<AltTranscriptData, HgvsError> {
        let (start_idx, end_idx) = self.get_variant_indices()?;
        let len = self.transcript_sequence.len();
        if end_idx > len {
            let first_bad = if start_idx >= len { start_idx } else { end_idx };
            return Err(HgvsError::ValidationError(format!(
                "Coordinate out of bounds: index {} is beyond transcript length {}",
                first_bad, len
            )));
        }

        // --- Validate reference sequence ---
        match &self.var_c.posedit.edit {
            NaEdit::RefAlt { ref_: Some(r), .. }
            | NaEdit::Del { ref_: Some(r), .. }
            | NaEdit::Dup { ref_: Some(r), .. } => {
                if !r.is_empty() && !r.chars().all(|c| c.is_ascii_digit()) {
                    let actual_ref = window(self.transcript_sequence, start_idx, end_idx);
                    if actual_ref != r {
                        return Err(HgvsError::TranscriptMismatch {
                            expected: r.to_string(),
                            found: actual_ref.to_string(),
                            start: start_idx,
                            end: end_idx,
                        });
                    }
                }
            }
            _ => {}
        }
        // --- End validation ---

        let cds_start_idx = self.cds_start_index.0 as usize;
        let variant_start_aa = if start_idx >= cds_start_idx {
            Some(ProteinPos(((start_idx - cds_start_idx) / 3) as i32))
        } else {
            Some(ProteinPos(0)) // 5' UTR variant
        };

        let (is_substitution, is_frameshift, alt_transcript) = match &self.var_c.posedit.edit {
            NaEdit::RefAlt { ref_, alt, .. } => {
                let is_identity = ref_.is_none() && alt.is_none();
                let is_ins = ref_.is_none()
                    && !is_identity
                    && self
                        .var_c
                        .posedit
                        .pos
                        .as_ref()
                        .is_some_and(|p| p.end.is_some());
                let alt_str = alt.as_deref().unwrap_or("");

                let (is_subst, is_fs, res) = if is_identity {
                    (false, false, self.transcript_sequence.to_string())
                } else if is_ins {
                    let ins_pos = start_idx + 1;
                    (
                        false,
                        alt_str.len() % 3 != 0,
                        splice(self.transcript_sequence, ins_pos, ins_pos, alt_str),
                    )
                } else {
                    let r_len = if let Some(r) = ref_ {
                        if r.chars().all(|c| c.is_ascii_digit()) {
                            end_idx - start_idx
                        } else {
                            r.len()
                        }
                    } else {
                        end_idx - start_idx
                    };

                    let is_subst =
                        ref_.is_some() && alt.is_some() && r_len == 1 && alt_str.len() == 1;
                    let is_fs = (alt_str.len() as i32 - r_len as i32) % 3 != 0;
                    (
                        is_subst,
                        is_fs,
                        splice(self.transcript_sequence, start_idx, end_idx, alt_str),
                    )
                };

                (is_subst, is_fs, res)
            }
            NaEdit::Del { ref_, .. } => {
                let r_len = if let Some(r) = ref_ {
                    if r.chars().all(|c| c.is_ascii_digit()) {
                        end_idx - start_idx
                    } else {
                        r.len()
                    }
                } else {
                    end_idx - start_idx
                };

                let res = splice(self.transcript_sequence, start_idx, end_idx, "");
                (false, (r_len as i32) % 3 != 0, res)
            }
            NaEdit::Ins { alt: Some(alt), .. } => {
                let ins_pos = start_idx + 1;
                let res = splice(self.transcript_sequence, ins_pos, ins_pos, alt);
                (false, (alt.len() as i32) % 3 != 0, res)
            }
            NaEdit::Dup { ref_, .. } => {
                let dup_str = match ref_ {
                    Some(r) if !r.chars().all(|c| c.is_ascii_digit()) => r.clone(),
                    _ => window(self.transcript_sequence, start_idx, end_idx).to_string(),
                };
                let res = splice(self.transcript_sequence, end_idx, end_idx, &dup_str);
                (false, (dup_str.len() as i32) % 3 != 0, res)
            }
            NaEdit::Inv { .. } => {
                let inv_str = crate::utils::reverse_complement(window(
                    self.transcript_sequence,
                    start_idx,
                    end_idx,
                ));
                let res = splice(self.transcript_sequence, start_idx, end_idx, &inv_str);
                (false, false, res)
            }
            NaEdit::Repeat { max, ref_, .. } => {
                let unit = match ref_ {
                    Some(r) => r.clone(),
                    None => window(self.transcript_sequence, start_idx, end_idx).to_string(),
                };

                // Extend over every existing copy of the unit.
                let seq = self.transcript_sequence;
                let mut current_idx = start_idx;
                while !unit.is_empty()
                    && current_idx + unit.len() <= seq.len()
                    && seq[current_idx..current_idx + unit.len()] == *unit
                {
                    current_idx += unit.len();
                }

                let total_str = unit.repeat(*max as usize);
                let res = splice(self.transcript_sequence, start_idx, current_idx, &total_str);

                let net_change =
                    (unit.len() as i32 * (*max as i32)) - (current_idx as i32 - start_idx as i32);
                (false, net_change % 3 != 0, res)
            }
            NaEdit::None => (false, false, self.transcript_sequence.to_string()),
            _ => {
                return Err(HgvsError::UnsupportedOperation(
                    "Unsupported edit for altseq".into(),
                ))
            }
        };

        let cds_start = self.cds_start_index.0 as usize;
        let alt_transcript_seq = alt_transcript;
        let aa_sequence = if cds_start < alt_transcript_seq.len() {
            crate::utils::translate(&alt_transcript_seq[cds_start..])
        } else {
            String::new()
        };

        // Find cds_end_i in the new sequence.
        // It's original_cds_end + net_change.
        // Or simply finding the length of the alt_transcript_seq.
        // Actually, we should be careful with 3' UTR.
        let net_change = alt_transcript_seq.len() as i32 - self.transcript_sequence.len() as i32;
        let cds_end_i = self.cds_end_index.0 + net_change;

        let is_fs = is_frameshift;

        Ok(AltTranscriptData {
            transcript_sequence: alt_transcript_seq,
            aa_sequence,
            cds_start_index: self.cds_start_index,
            cds_end_index: TranscriptPos(cds_end_i),
            protein_accession: self.protein_accession.clone(),
            is_frameshift: is_fs,
            variant_start_aa,
            frameshift_start: if is_fs { variant_start_aa } else { None },
            is_substitution,
            is_ambiguous: false,
            variant_start_idx: start_idx,
            variant_end_idx: end_idx,
            c_variant: self.var_c.clone(),
        })
    }

    fn get_variant_indices(&self) -> Result<(usize, usize), HgvsError> {
        let pos = self
            .var_c
            .posedit
            .pos
            .as_ref()
            .ok_or_else(|| HgvsError::ValidationError("Missing position".into()))?;
        let (start, end) = self.mapper.interval_to_n(pos).map_err(|e| match e {
            HgvsError::UnsupportedOperation(_) => HgvsError::UnsupportedOperation(
                "Intronic variants not yet supported in c_to_p".into(),
            ),
            other => other,
        })?;
        if start.0 < 0 {
            return Err(HgvsError::ValidationError(format!(
                "Position {} before transcript start",
                start.0
            )));
        }
        Ok((start.0 as usize, end.0 as usize))
    }
}
