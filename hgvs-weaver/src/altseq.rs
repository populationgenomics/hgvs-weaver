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

        let seq = self.transcript_sequence;
        let edit = &self.var_c.posedit.edit;

        // The bases the variant says are there must be there.
        if let Some(stated) = edit.stated_ref() {
            let actual = window(seq, start_idx, end_idx);
            if actual != stated {
                return Err(HgvsError::TranscriptMismatch {
                    expected: stated.to_string(),
                    found: actual.to_string(),
                    start: start_idx,
                    end: end_idx,
                });
            }
        }

        let cds_start_idx = self.cds_start_index.0 as usize;
        let variant_start_aa = if start_idx >= cds_start_idx {
            Some(ProteinPos(((start_idx - cds_start_idx) / 3) as i32))
        } else {
            Some(ProteinPos(0)) // 5' UTR variant
        };

        let (is_substitution, is_frameshift, alt_transcript) = match edit {
            NaEdit::Con { .. } | NaEdit::NACopy { .. } => {
                return Err(HgvsError::UnsupportedOperation(
                    "Unsupported edit for altseq".into(),
                ))
            }
            _ => {
                let resolved = edit
                    .resolve_with(start_idx, end_idx, |s, e| Ok(window(seq, s, e).to_string()))?;
                let is_subst = matches!(
                    edit,
                    NaEdit::RefAlt {
                        ref_: Some(_),
                        alt: Some(_),
                        ..
                    }
                ) && resolved.ref_.len() == 1
                    && resolved.alt.len() == 1;
                let res = splice(seq, resolved.start, resolved.end, &resolved.alt);
                (is_subst, resolved.is_frameshift(), res)
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
