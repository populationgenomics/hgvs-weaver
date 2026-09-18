//! Normalisation of a nucleotide edit against one reference sequence.
//!
//! HGVS wants an edit at its 3'-most equivalent position, with deleted or
//! duplicated bases spelled out, and an insertion that copies the bases just
//! before it written as a duplication. All of that is decided here, on a
//! 0-based half-open range over a [`Reference`]. The coordinate systems (g.,
//! c., n.) only convert their positions in and out; see
//! [`VariantMapper::normalize_variant`](crate::mapper::VariantMapper::normalize_variant).

use crate::edits::{NaEdit, ResolvedEdit};
use crate::error::HgvsError;
use crate::reference::Reference;
use crate::structs::strip_common_prefix_suffix;

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
    // Only a deletion, duplication or insertion is slid: those are written
    // without stating bases that would go stale, and a repeat is written
    // against its run. A delins is left where it is.
    let k = if matches!(
        edit,
        NaEdit::Del { .. } | NaEdit::Dup { .. } | NaEdit::Ins { .. }
    ) {
        shift_3(reference, &edit.resolve(reference, start, end)?)?
    } else {
        0
    };
    let (start, end) = (start + k, end + k);
    // Sliding an insertion 3' by k rotates what is inserted: the k bases it
    // passed over now precede it, and the same k bases of the insert follow.
    if k > 0 {
        if let NaEdit::Ins { alt: Some(seq), .. } = &mut edit {
            let n = seq.len();
            if n > 0 {
                let r = k % n;
                *seq = format!("{}{}", &seq[r..], &seq[..r]);
            }
        }
    }

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

/// The part of a resolved edit that can slide: after dropping the bases the
/// reference and alternate share, either a pure insertion (`[at, at)` plus the
/// inserted bases) or a pure deletion (`[at, at + n)` plus the deleted bases).
/// A substitution, an inversion or a delins with bases on both sides cannot
/// slide and yields `None`.
fn movable(resolved: &ResolvedEdit) -> Option<(usize, usize, String)> {
    if resolved.ref_ == resolved.alt {
        return None;
    }
    let (at, r, a) =
        strip_common_prefix_suffix(resolved.start as i32, &resolved.ref_, &resolved.alt);
    let at = at as usize;
    if r.is_empty() && !a.is_empty() {
        Some((at, at, a))
    } else if a.is_empty() && !r.is_empty() {
        Some((at, at + r.len(), r))
    } else {
        None
    }
}

/// The widest range over which `resolved` is ambiguous: every reference
/// position the change could equally well be written at, 5' and 3'. This is
/// what an unambiguous SPDI spells out. An edit that cannot slide is its own
/// range.
pub fn ambiguous_range(
    reference: &Reference<'_, '_>,
    resolved: &ResolvedEdit,
) -> Result<(usize, usize), HgvsError> {
    let Some((at_start, at_end, pattern)) = movable(resolved) else {
        return Ok((resolved.start, resolved.end));
    };
    let k5 = reference.run_left(at_start, pattern.as_bytes())?;
    let k3 = reference.run_right(at_end, pattern.as_bytes())?;
    Ok((
        resolved.start.min(at_start - k5),
        resolved.end.max(at_end + k3),
    ))
}

/// How many positions the edit can move 3' and still describe the same change,
/// measured at its movable insertion or deletion point.
pub fn shift_3(reference: &Reference<'_, '_>, resolved: &ResolvedEdit) -> Result<usize, HgvsError> {
    match movable(resolved) {
        Some((_, end, pattern)) => reference.run_right(end, pattern.as_bytes()),
        None => Ok(0),
    }
}

/// How many positions the edit can move 5' and still describe the same change,
/// measured at its movable insertion or deletion point. For a duplication that
/// point is the end of the duplicated bases, so the count includes them.
pub fn shift_5(reference: &Reference<'_, '_>, resolved: &ResolvedEdit) -> Result<usize, HgvsError> {
    match movable(resolved) {
        Some((start, _, pattern)) => reference.run_left(start, pattern.as_bytes()),
        None => Ok(0),
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
    fn a_slid_insertion_is_rotated() {
        // Inserting GA before index 2 of GGGGGGGGG slides one base (the next
        // base is G) and becomes an insertion of AG before index 3: the same
        // molecule, GGGAGGGGGGG.
        let out = run("GGGGGGGGG", PlacedEdit::from_hgvs_range(1, 3, ins("GA")));
        assert_eq!(
            out,
            PlacedEdit {
                start: 3,
                end: 3,
                edit: ins("AG")
            }
        );
        // Sliding a whole number of periods leaves the insert as written; here
        // it then reads as a duplication of the last CAG copy, at [5, 8).
        let out = run("TTCAGCAGTT", PlacedEdit::from_hgvs_range(1, 3, ins("CAG")));
        assert_eq!((out.start, out.end), (5, 8));
        assert_eq!(
            out.edit,
            NaEdit::Dup {
                ref_: Some("CAG".into()),
                uncertain: false
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
        let res = |e: NaEdit, s, t| e.resolve(&r, s, t).unwrap();
        assert_eq!(shift_3(&r, &res(del(None), 2, 5)).unwrap(), 3);
        assert_eq!(shift_5(&r, &res(del(None), 5, 8)).unwrap(), 3);
        assert_eq!(shift_3(&r, &res(ins("CAG"), 5, 5)).unwrap(), 3);
        assert_eq!(shift_5(&r, &res(ins("CAG"), 5, 5)).unwrap(), 3);
        // A duplication slides like an insertion of its bases at its end, so
        // its 5' count includes the duplicated bases themselves.
        let dup = NaEdit::Dup {
            ref_: None,
            uncertain: false,
        };
        assert_eq!(shift_3(&r, &res(dup.clone(), 2, 5)).unwrap(), 3);
        assert_eq!(shift_5(&r, &res(dup.clone(), 5, 8)).unwrap(), 6);
        // The ambiguous range covers the whole CAG run either way.
        assert_eq!(ambiguous_range(&r, &res(dup, 5, 8)).unwrap(), (2, 8));
        assert_eq!(ambiguous_range(&r, &res(del(None), 2, 5)).unwrap(), (2, 8));
        assert_eq!(ambiguous_range(&r, &res(ins("CAG"), 5, 5)).unwrap(), (2, 8));
        // A substitution is not ambiguous.
        let sub = NaEdit::RefAlt {
            ref_: Some("C".into()),
            alt: Some("T".into()),
            uncertain: false,
        };
        assert_eq!(ambiguous_range(&r, &res(sub, 2, 3)).unwrap(), (2, 3));
        // A delins that is really an insertion (delCinsCA at index 2, before A)
        // slides like the insertion of A it is.
        let delins = NaEdit::RefAlt {
            ref_: Some("C".into()),
            alt: Some("CA".into()),
            uncertain: false,
        };
        assert_eq!(shift_3(&r, &res(delins, 2, 3)).unwrap(), 1);
    }
}
