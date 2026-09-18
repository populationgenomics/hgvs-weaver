//! The canonical allele: one unambiguous statement of a nucleotide change on a
//! reference sequence.
//!
//! HGVS admits many spellings of the same change (`g.5del`, `g.6del` and
//! `g.5_6delinsA` can all be one event in a run of As). SPDI and GA4GH VRS
//! resolve this with "fully justified" normalisation, adapted from NCBI's
//! Variant Overprecision Correction Algorithm: trim the bases the reference and
//! alternate share, then widen a pure insertion or deletion over every
//! position it could equally well be written at. The result is a value two
//! variants can be compared on, and that SPDI, VRS JSON and normalised HGVS
//! are all renderings of.

use crate::edits::ResolvedEdit;
use crate::error::HgvsError;
use crate::normalize::ambiguous_range;
use crate::reference::Reference;
use crate::structs::strip_common_prefix_suffix;

/// A change to `[start, end)` of the sequence `accession`, widened over its
/// whole region of ambiguity. Coordinates are 0-based inter-residue, as in
/// SPDI and VRS.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CanonicalAllele {
    pub accession: String,
    pub start: usize,
    pub end: usize,
    /// The reference bases over `[start, end)`.
    pub reference: String,
    /// The bases that replace them.
    pub alternate: String,
    /// `Some(unit length)` when the alternate is the reference repeated: a
    /// deletion, or an insertion of copies of reference bases. This is VRS's
    /// `ReferenceLengthExpression`; `None` means the alternate is literal.
    pub repeat_subunit: Option<usize>,
}

impl CanonicalAllele {
    /// The SPDI form, `accession:start:reference:alternate`.
    pub fn spdi(&self) -> String {
        format!(
            "{}:{}:{}:{}",
            self.accession, self.start, self.reference, self.alternate
        )
    }

    /// Whether the allele leaves the reference unchanged.
    pub fn is_reference(&self) -> bool {
        self.reference == self.alternate
    }

    /// Canonicalises `resolved`, an edit already resolved against `reference`.
    ///
    /// The steps are those of VRS allele normalisation:
    /// 1. trim the shared prefix and suffix;
    /// 2. nothing left on either side: a reference allele over the input range;
    ///    bases left on both sides: a substitution, done;
    /// 3. otherwise roll the remaining insertion or deletion left and right
    ///    while the reference keeps matching it cyclically;
    /// 4. widen both sequences over that region;
    /// 5. decide whether the alternate is reference-derived.
    pub fn canonicalize(
        reference: &Reference<'_, '_>,
        accession: &str,
        resolved: &ResolvedEdit,
    ) -> Result<Self, HgvsError> {
        let accession = accession.to_string();
        // The allele describes the sequence, so its reference is what the
        // sequence holds at the range, whatever bases the variant stated. A
        // transcript-derived variant whose stated base differs from the genome
        // canonicalises to what the genome actually changes, possibly nothing.
        let actual = reference.slice(resolved.start, resolved.end)?;
        let (at, r, a) = strip_common_prefix_suffix(resolved.start as i32, &actual, &resolved.alt);
        let at = at as usize;

        if r.is_empty() && a.is_empty() {
            return Ok(CanonicalAllele {
                accession,
                start: resolved.start,
                end: resolved.end,
                reference: actual.clone(),
                alternate: actual,
                repeat_subunit: Some(resolved.end - resolved.start),
            });
        }
        if !r.is_empty() && !a.is_empty() {
            return Ok(CanonicalAllele {
                accession,
                start: at,
                end: at + r.len(),
                reference: r,
                alternate: a,
                repeat_subunit: None,
            });
        }

        // A pure insertion at `at`, or a pure deletion of `[at, at + r.len())`.
        let trimmed = ResolvedEdit {
            start: at,
            end: at + r.len(),
            ref_: r.clone(),
            alt: a.clone(),
        };
        let (u_start, u_end) = ambiguous_range(reference, &trimmed)?;
        let widened_ref = reference.slice(u_start, u_end)?;
        let widened_alt = format!(
            "{}{}{}",
            &widened_ref[..at - u_start],
            a,
            &widened_ref[trimmed.end - u_start..]
        );
        let seed_len = r.len().max(a.len());
        let repeat_subunit = if r.is_empty() && u_start == at && u_end == at {
            None // an insertion nothing can slide: literal
        } else if a.is_empty() {
            Some(seed_len) // a deletion is always reference-derived
        } else {
            reference_derived_unit(&widened_ref, &widened_alt, seed_len)
        };
        Ok(CanonicalAllele {
            accession,
            start: u_start,
            end: u_end,
            reference: widened_ref,
            alternate: widened_alt,
            repeat_subunit,
        })
    }
}

/// The smallest divisor `d` of `seed_len` such that some `d`-base window of
/// `reference` repeated cyclically reproduces `alternate`, if any.
fn reference_derived_unit(reference: &str, alternate: &str, seed_len: usize) -> Option<usize> {
    let alt = alternate.as_bytes();
    let refb = reference.as_bytes();
    (1..=seed_len)
        .filter(|d| seed_len.is_multiple_of(*d) && *d <= refb.len())
        .find(|&d| {
            (0..=refb.len() - d).any(|i| {
                let unit = &refb[i..i + d];
                alt.iter().enumerate().all(|(k, &c)| c == unit[k % d])
            })
        })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::{DataProvider, IdentifierKind, IdentifierType, TranscriptData};
    use crate::edits::NaEdit;
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

    fn canon(seq: &'static str, edit: NaEdit, start: usize, end: usize) -> CanonicalAllele {
        let hdp = Fixed(seq);
        let store = ReferenceStore::with_block_size(&hdp, 4);
        let r = store.reference("X", IdentifierType::GenomicAccession);
        let resolved = edit.resolve(&r, start, end).unwrap();
        CanonicalAllele::canonicalize(&r, "X", &resolved).unwrap()
    }
    fn del() -> NaEdit {
        NaEdit::Del {
            ref_: None,
            uncertain: false,
        }
    }
    fn ins(a: &str) -> NaEdit {
        NaEdit::Ins {
            alt: Some(a.into()),
            uncertain: false,
        }
    }
    fn sub(r: &str, a: &str) -> NaEdit {
        NaEdit::RefAlt {
            ref_: Some(r.into()),
            alt: Some(a.into()),
            uncertain: false,
        }
    }

    //                 0123456789
    const SEQ: &str = "TTCAGCAGTT";

    #[test]
    fn deletion_anywhere_in_a_run_is_the_same_allele() {
        let a = canon(SEQ, del(), 2, 5);
        let b = canon(SEQ, del(), 5, 8);
        assert_eq!(a, b);
        assert_eq!(a.spdi(), "X:2:CAGCAG:CAG");
        assert_eq!(a.repeat_subunit, Some(3));
    }

    #[test]
    fn duplication_and_insertion_of_the_unit_are_the_same_allele() {
        let dup = canon(
            SEQ,
            NaEdit::Dup {
                ref_: None,
                uncertain: false,
            },
            2,
            5,
        );
        let inserted = canon(SEQ, ins("CAG"), 4, 6);
        assert_eq!(dup, inserted);
        assert_eq!(dup.spdi(), "X:2:CAGCAG:CAGCAGCAG");
        assert_eq!(dup.repeat_subunit, Some(3));
    }

    #[test]
    fn literal_insertion_stays_literal_and_substitution_is_trimmed() {
        let i = canon(SEQ, ins("GG"), 1, 3);
        assert_eq!(i.spdi(), "X:2::GG");
        assert_eq!(i.repeat_subunit, None);
        // delAGinsAT is A>? no: only the last base changes.
        let s = canon(SEQ, sub("AG", "AT"), 3, 5);
        assert_eq!(s.spdi(), "X:4:G:T");
        assert_eq!(s.repeat_subunit, None);
    }

    #[test]
    fn a_stated_reference_that_disagrees_with_the_sequence_is_ignored() {
        // The variant claims G>A at index 4, but the sequence has A there:
        // nothing changes, and the allele says so with the sequence's bases.
        let e = canon(SEQ, sub("G", "A"), 3, 4);
        assert!(e.is_reference());
        assert_eq!(e.spdi(), "X:3:A:A");
    }

    #[test]
    fn identity_is_a_reference_allele_over_its_range() {
        let e = canon(
            SEQ,
            NaEdit::RefAlt {
                ref_: None,
                alt: None,
                uncertain: false,
            },
            2,
            5,
        );
        assert!(e.is_reference());
        assert_eq!(e.spdi(), "X:2:CAG:CAG");
        assert_eq!(e.repeat_subunit, Some(3));
    }

    #[test]
    fn reference_derived_unit_finds_the_smallest_period() {
        assert_eq!(reference_derived_unit("CAGCAG", "CAGCAGCAG", 3), Some(3));
        assert_eq!(reference_derived_unit("AAAA", "AAAAAA", 2), Some(1));
        assert_eq!(reference_derived_unit("CAGCAG", "CAGCAGTTT", 3), None);
    }
}
