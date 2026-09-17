//! Normalisation of a nucleotide edit against one reference sequence.
//!
//! HGVS wants an edit at its 3'-most equivalent position, with deleted or
//! duplicated bases spelled out, and an insertion that copies the bases just
//! before it written as a duplication. All of that is decided here, on a
//! 0-based half-open range over a [`Reference`]. The coordinate systems (g.,
//! c., n.) only convert their positions in and out; see
//! [`VariantMapper::normalize_variant`](crate::mapper::VariantMapper::normalize_variant).

use crate::edits::NaEdit;
use crate::error::HgvsError;
use crate::reference::Reference;

/// An edit placed on a reference by a 0-based half-open range.
///
/// For an insertion the range is empty: `start == end` is the index of the
/// base the inserted bases go in front of.
#[derive(Debug, Clone, PartialEq)]
pub struct PlacedEdit {
    pub start: usize,
    pub end: usize,
    pub edit: NaEdit,
}

impl PlacedEdit {
    /// Places an edit written over the HGVS range `[start, end)`.
    ///
    /// HGVS names an insertion by its two flanking bases; the edit itself is
    /// the empty range in front of the second one. Every other edit covers
    /// the bases it names.
    pub fn from_hgvs_range(start: usize, end: usize, edit: NaEdit) -> Self {
        if matches!(edit, NaEdit::Ins { .. }) && end > start {
            let anchor = end - 1;
            PlacedEdit {
                start: anchor,
                end: anchor,
                edit,
            }
        } else {
            PlacedEdit { start, end, edit }
        }
    }

    pub fn is_insertion(&self) -> bool {
        matches!(self.edit, NaEdit::Ins { .. })
    }

    /// The half-open range HGVS writes: the two flanking bases for an
    /// insertion, the edited bases for everything else.
    pub fn hgvs_range(&self) -> (usize, usize) {
        if self.is_insertion() {
            (self.start.saturating_sub(1), self.start + 1)
        } else {
            (self.start, self.end)
        }
    }
}

/// Normalises `placed` against `reference`: shifts it as far 3' as the
/// sequence allows, fills in the reference bases of a deletion or duplication,
/// and rewrites an insertion as a duplication when the inserted bases repeat
/// the bases immediately before it.
pub fn normalize(
    reference: &Reference<'_, '_>,
    placed: PlacedEdit,
) -> Result<PlacedEdit, HgvsError> {
    let PlacedEdit {
        start,
        end,
        mut edit,
    } = placed;
    let k = shift_3(reference, start, end, &edit)?;
    let (start, end) = (start + k, end + k);

    if let NaEdit::Del { ref_, .. } | NaEdit::Dup { ref_, .. } = &mut edit {
        *ref_ = Some(reference.slice(start, end)?);
    }

    if let NaEdit::Ins {
        alt: Some(seq),
        uncertain,
    } = &edit
    {
        let n = seq.len();
        if n > 0 && start >= n && reference.slice(start - n, start)? == *seq {
            return Ok(PlacedEdit {
                start: start - n,
                end: start,
                edit: NaEdit::Dup {
                    ref_: Some(seq.clone()),
                    uncertain: *uncertain,
                },
            });
        }
    }

    Ok(PlacedEdit { start, end, edit })
}

/// How many positions the edit over `[start, end)` can move 3' and still
/// describe the same change to the sequence.
pub fn shift_3(
    reference: &Reference<'_, '_>,
    start: usize,
    end: usize,
    edit: &NaEdit,
) -> Result<usize, HgvsError> {
    match shift_pattern(reference, start, end, edit)? {
        Some(pattern) => reference.run_right(end, pattern.as_bytes()),
        None => Ok(0),
    }
}

/// How many positions the edit over `[start, end)` can move 5' and still
/// describe the same change to the sequence.
pub fn shift_5(
    reference: &Reference<'_, '_>,
    start: usize,
    end: usize,
    edit: &NaEdit,
) -> Result<usize, HgvsError> {
    match shift_pattern(reference, start, end, edit)? {
        Some(pattern) => reference.run_left(start, pattern.as_bytes()),
        None => Ok(0),
    }
}

/// The bases whose repetition lets an edit over `[start, end)` slide along the
/// reference, or `None` when the edit cannot be shifted at all.
///
/// Deletions, duplications and length-changing delins slide while the bases
/// beyond the edit repeat its reference bases; pure insertions slide while
/// they repeat the inserted bases.
fn shift_pattern(
    reference: &Reference<'_, '_>,
    start: usize,
    end: usize,
    edit: &NaEdit,
) -> Result<Option<String>, HgvsError> {
    if matches!(
        edit,
        NaEdit::None | NaEdit::Con { .. } | NaEdit::NACopy { .. }
    ) {
        return Ok(None);
    }
    let resolved = edit.resolve(reference, start, end)?;
    let (ref_str, alt_str) = (resolved.ref_, resolved.alt);
    if ref_str == alt_str && matches!(edit, NaEdit::RefAlt { .. }) {
        return Ok(None);
    }
    // A delins (bases out, different bases in) has no shift rule: sliding it
    // by the deletion rule changes the resulting sequence unless the inserted
    // bases happen to be a run of the repeated base.
    let is_del_or_dup = matches!(edit, NaEdit::Del { .. } | NaEdit::Dup { .. });
    if is_del_or_dup || (!ref_str.is_empty() && alt_str.is_empty()) {
        let pattern = if ref_str.is_empty() {
            reference.slice(start, end)?
        } else {
            ref_str
        };
        Ok((!pattern.is_empty()).then_some(pattern))
    } else if start == end && !alt_str.is_empty() {
        Ok(Some(alt_str))
    } else {
        Ok(None)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::{DataProvider, IdentifierKind, IdentifierType, TranscriptData};
    use crate::reference::ReferenceStore;

    struct Fixed(&'static str);

    impl DataProvider for Fixed {
        fn get_transcript(&self, _: &str, _: Option<&str>) -> Result<TranscriptData, HgvsError> {
            unreachable!()
        }
        fn get_seq(
            &self,
            _: &str,
            start: i32,
            end: Option<i32>,
            _: IdentifierType,
        ) -> Result<String, HgvsError> {
            let s = (start.max(0) as usize).min(self.0.len());
            let e = end.map_or(self.0.len(), |e| (e as usize).min(self.0.len()));
            Ok(self.0[s..e.max(s)].to_string())
        }
        fn get_symbol_accessions(
            &self,
            _: &str,
            _: IdentifierKind,
            _: IdentifierKind,
        ) -> Result<Vec<(IdentifierType, String)>, HgvsError> {
            Ok(vec![])
        }
        fn get_identifier_type(&self, _: &str) -> Result<IdentifierType, HgvsError> {
            Ok(IdentifierType::GenomicAccession)
        }
    }

    fn del(ref_: Option<&str>) -> NaEdit {
        NaEdit::Del {
            ref_: ref_.map(str::to_string),
            uncertain: false,
        }
    }
    fn ins(alt: &str) -> NaEdit {
        NaEdit::Ins {
            alt: Some(alt.to_string()),
            uncertain: false,
        }
    }

    fn run(seq: &'static str, placed: PlacedEdit) -> PlacedEdit {
        let hdp = Fixed(seq);
        let store = ReferenceStore::with_block_size(&hdp, 4);
        normalize(
            &store.reference("X", IdentifierType::GenomicAccession),
            placed,
        )
        .unwrap()
    }

    #[test]
    fn insertion_places_at_the_second_flanking_base() {
        let p = PlacedEdit::from_hgvs_range(4, 6, ins("A"));
        assert_eq!((p.start, p.end), (5, 5));
        assert_eq!(p.hgvs_range(), (4, 6));
        let d = PlacedEdit::from_hgvs_range(4, 6, del(None));
        assert_eq!((d.start, d.end), (4, 6));
    }

    #[test]
    fn deletion_shifts_3_prime_and_fills_reference() {
        //             0123456789
        let out = run("TTCAGCAGTT", PlacedEdit::from_hgvs_range(2, 5, del(None)));
        assert_eq!(
            out,
            PlacedEdit {
                start: 5,
                end: 8,
                edit: del(Some("CAG"))
            }
        );
    }

    #[test]
    fn insertion_that_repeats_preceding_bases_becomes_duplication() {
        // Insert CAG in front of index 5 (between the two CAGs): shifts to the
        // end of the run and duplicates the last CAG.
        let out = run("TTCAGCAGTT", PlacedEdit::from_hgvs_range(4, 6, ins("CAG")));
        assert_eq!(
            out,
            PlacedEdit {
                start: 5,
                end: 8,
                edit: NaEdit::Dup {
                    ref_: Some("CAG".into()),
                    uncertain: false
                }
            }
        );
        // A genuine insertion stays one.
        let out = run("TTCAGCAGTT", PlacedEdit::from_hgvs_range(1, 3, ins("GGG")));
        assert_eq!(
            out,
            PlacedEdit {
                start: 2,
                end: 2,
                edit: ins("GGG")
            }
        );
    }

    #[test]
    fn delins_is_not_shifted() {
        // Deleting CAG at [2,5) and inserting TT: the base after the range (C)
        // equals the first deleted base, but sliding would put TT after that C,
        // which is a different sequence. Only the reference gets resolved.
        let delins = NaEdit::RefAlt {
            ref_: Some("3".into()),
            alt: Some("TT".into()),
            uncertain: false,
        };
        let out = run(
            "TTCAGCAGTT",
            PlacedEdit::from_hgvs_range(2, 5, delins.clone()),
        );
        assert_eq!(
            out,
            PlacedEdit {
                start: 2,
                end: 5,
                edit: delins
            }
        );
    }

    #[test]
    fn substitution_is_left_alone() {
        let sub = NaEdit::RefAlt {
            ref_: Some("C".into()),
            alt: Some("T".into()),
            uncertain: false,
        };
        let out = run("TTCAGCAGTT", PlacedEdit::from_hgvs_range(2, 3, sub.clone()));
        assert_eq!(
            out,
            PlacedEdit {
                start: 2,
                end: 3,
                edit: sub
            }
        );
    }

    #[test]
    fn shift_5_mirrors_shift_3() {
        let hdp = Fixed("TTCAGCAGTT");
        let store = ReferenceStore::with_block_size(&hdp, 4);
        let r = store.reference("X", IdentifierType::GenomicAccession);
        assert_eq!(shift_3(&r, 2, 5, &del(None)).unwrap(), 3);
        assert_eq!(shift_5(&r, 5, 8, &del(None)).unwrap(), 3);
        assert_eq!(shift_3(&r, 5, 5, &ins("CAG")).unwrap(), 3);
        assert_eq!(shift_5(&r, 5, 5, &ins("CAG")).unwrap(), 3);
    }
}
