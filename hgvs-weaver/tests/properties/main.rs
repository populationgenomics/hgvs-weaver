//! Property tests: relationships that must hold for every input the
//! strategies in `strategies.rs` can build. Each one pins a bug class the
//! ClinVar gate cannot see because ClinVar happens not to contain it.

mod strategies;
#[path = "../support/mod.rs"]
mod support;

use hgvs_weaver::data::{IdentifierType, Strand};
use hgvs_weaver::normalize::{ambiguous_range, normalize, ExonBounds, PlacedEdit};
use hgvs_weaver::reference::ReferenceStore;
use hgvs_weaver::structs::GVariant;
use proptest::prelude::*;
use strategies::*;

const G: IdentifierType = IdentifierType::GenomicAccession;

proptest! {
    /// Normalising an edit never changes the molecule it describes, and doing
    /// it twice is the same as doing it once.
    #[test]
    fn normalising_preserves_the_edited_sequence(
        (seq, p) in seq_and_edit(false),
        block in 1usize..=16,
    ) {
        let hdp = Provider::new().sequence("X", &seq);
        let store = ReferenceStore::with_block_size(&hdp, block);
        let r = store.reference("X", G);

        let before = PlacedEdit::from_hgvs_range(p.start, p.end, p.edit.clone());
        let after = normalize(&r, before.clone(), &ExonBounds::unspliced()).unwrap();
        let after_hgvs = placed_to_hgvs(&after);
        prop_assert_eq!(apply(&seq, &p), apply(&seq, &after_hgvs), "normalisation changed the molecule");

        let again = normalize(&r, after.clone(), &ExonBounds::unspliced()).unwrap();
        prop_assert_eq!(again, after, "normalisation is not idempotent");
    }

    /// A deletion, duplication or delins normalises to the same spelling
    /// whether or not it states the bases it replaces.
    #[test]
    fn stating_the_bases_an_edit_replaces_does_not_change_its_normalised_spelling(
        (seq, p) in seq_and_edit(false),
    ) {
        let bases = Some(seq[p.start..p.end].to_string());
        let stated = match &p.edit {
            NaEdit::Del { uncertain, .. } => NaEdit::Del { ref_: bases, uncertain: *uncertain },
            NaEdit::Dup { uncertain, .. } => NaEdit::Dup { ref_: bases, uncertain: *uncertain },
            NaEdit::RefAlt { alt, uncertain, .. } => NaEdit::RefAlt { ref_: bases, alt: alt.clone(), uncertain: *uncertain },
            _ => return Ok(()),
        };
        let hdp = Provider::new().sequence("X", &seq);
        let store = ReferenceStore::with_block_size(&hdp, 4);
        let r = store.reference("X", G);
        let normalised = |edit: NaEdit| {
            normalize(&r, PlacedEdit::from_hgvs_range(p.start, p.end, edit), &ExonBounds::unspliced()).unwrap()
        };
        prop_assert_eq!(normalised(stated.clone()), normalised(p.edit.clone()), "{:?} on {}", stated, seq);
    }

    /// A change and a delins over a wider range that writes out the
    /// reference's bases either side of it are one change, and normalise to
    /// one spelling.
    #[test]
    fn a_change_and_a_wider_delins_spelling_it_normalise_alike(
        (seq, p) in seq_and_edit(false),
        pad_5 in 0usize..4,
        pad_3 in 0usize..4,
    ) {
        let edited = apply(&seq, &p);
        prop_assume!(edited != seq);
        let (lo, hi) = touched(&p);
        let (mut lo, mut hi) = (lo.saturating_sub(pad_5), (hi + pad_3).min(seq.len()));
        // A delins replaces at least one base with at least one.
        if edited.len() + hi - lo == seq.len() || hi == lo {
            if hi < seq.len() { hi += 1 } else { lo -= 1 }
        }
        let alt = edited[lo..edited.len() - (seq.len() - hi)].to_string();
        let wide = PlacedHgvs {
            start: lo,
            end: hi,
            edit: NaEdit::RefAlt { ref_: Some(String::new()), alt: Some(alt), uncertain: false },
        };
        prop_assert_eq!(apply(&seq, &wide), edited.clone(), "the wider delins is another change");

        let hdp = Provider::new().sequence("X", &seq);
        let store = ReferenceStore::with_block_size(&hdp, 4);
        let r = store.reference("X", G);
        let normalised = |q: &PlacedHgvs| {
            normalize(&r, PlacedEdit::from_hgvs_range(q.start, q.end, q.edit.clone()), &ExonBounds::unspliced()).unwrap()
        };
        prop_assert_eq!(normalised(&wide), normalised(&p), "{:?} and {:?} on {}", wide, p, seq);
    }

    /// VRS and SPDI are lossless: an allele read back from either is the same
    /// allele, whatever spelling the input used.
    #[test]
    fn vrs_and_spdi_round_trip((seq, p) in seq_and_edit(false)) {
        let hdp = Provider::new().sequence("X", &seq);
        let mapper = VariantMapper::new(&hdp);
        let v = SequenceVariant::Genomic(genomic_variant("X", &p));
        let vrs = mapper.to_vrs(&v).unwrap();
        let back = mapper.from_vrs(&vrs.to_json(), Some("X")).unwrap();
        prop_assert_eq!(&mapper.to_vrs(&back).unwrap().id, &vrs.id, "{} read back as {}", v, back);
        let spdi = mapper.to_spdi_unambiguous(&v).unwrap();
        let back = mapper.from_spdi(&spdi).unwrap();
        prop_assert_eq!(&mapper.to_vrs(&back).unwrap().id, &vrs.id, "{} as {} read back as {}", v, spdi, back);
        // Reading back is idempotent: the normalised spelling maps to itself.
        prop_assert_eq!(mapper.from_vrs(&mapper.to_vrs(&back).unwrap().to_json(), Some("X")).unwrap().to_string(), back.to_string());
    }

    /// The ambiguity range is exactly the set of positions the change can be
    /// written at: every placement inside gives the same molecule, and the
    /// placements one past either end do not.
    #[test]
    fn ambiguous_range_is_sound_and_complete((seq, p) in seq_and_edit(true)) {
        let hdp = Provider::new().sequence("X", &seq);
        let store = ReferenceStore::with_block_size(&hdp, 4);
        let r = store.reference("X", G);
        let resolved = p.edit.resolve(&r, p.start, p.end).unwrap();
        let (u_start, u_end) = ambiguous_range(&r, &resolved).unwrap();
        let width = resolved.end - resolved.start;
        let expected = apply(&seq, &p);

        // The same change written at position `s`: a deletion of the bases
        // there, or the insertion rotated by how far it moved.
        let is_ins = matches!(p.edit, hgvs_weaver::edits::NaEdit::Ins { .. });
        let placed_at = |s: usize| -> PlacedHgvs {
            let edit = match &p.edit {
                hgvs_weaver::edits::NaEdit::Ins { alt: Some(a), uncertain } => {
                    let n = a.len();
                    let r = ((s as i64 - resolved.start as i64).rem_euclid(n as i64)) as usize;
                    hgvs_weaver::edits::NaEdit::Ins { alt: Some(format!("{}{}", &a[r..], &a[..r])), uncertain: *uncertain }
                }
                other => other.clone(),
            };
            as_hgvs_range(&PlacedHgvs { start: s, end: s + width, edit })
        };
        // Soundness: every placement inside the range is the same molecule. An
        // insertion cannot be written before base 1, so skip position 0.
        for s in u_start..=u_end.saturating_sub(width) {
            if is_ins && s == 0 { continue; }
            prop_assert_eq!(apply(&seq, &placed_at(s)), expected.clone(), "placement at {} inside the range differs", s);
        }
        // Completeness: one past each end is a different molecule (when in bounds).
        if u_end < seq.len() {
            prop_assert_ne!(apply(&seq, &placed_at(u_end - width + 1)), expected.clone(), "range is not maximal on the right");
        }
        if u_start > 0 && !(is_ins && u_start == 1) {
            prop_assert_ne!(apply(&seq, &placed_at(u_start - 1)), expected.clone(), "range is not maximal on the left");
        }
    }

    /// Paging is invisible: whatever the block size, a Reference agrees with
    /// slicing the string directly, including past the end.
    #[test]
    fn reference_store_agrees_with_direct_slicing(
        seq in dna(1, 120),
        block in 1usize..=17,
        ranges in prop::collection::vec((0usize..130, 0usize..130), 1..=8),
        pattern in dna(1, 4),
    ) {
        let hdp = Provider::new().sequence("X", &seq);
        let store = ReferenceStore::with_block_size(&hdp, block);
        let r = store.reference("X", G);
        for (a, b) in ranges {
            let (s, e) = (a.min(b), a.max(b));
            let direct = &seq[s.min(seq.len())..e.min(seq.len())];
            prop_assert_eq!(r.slice(s, e).unwrap(), direct.to_string());
            prop_assert_eq!(r.base(s).unwrap(), seq.as_bytes().get(s).copied());
            // run_right/run_left against a naive cyclic match
            let pat = pattern.as_bytes();
            let mut k = 0;
            while s + k < seq.len() && seq.as_bytes()[s + k] == pat[k % pat.len()] { k += 1; }
            prop_assert_eq!(r.run_right(s, pat).unwrap(), k);
            let mut k = 0;
            while k < s && s - 1 - k < seq.len() && seq.as_bytes()[s - 1 - k] == pat[pat.len() - 1 - (k % pat.len())] { k += 1; }
            prop_assert_eq!(r.run_left(s, pat).unwrap(), k);
        }
        prop_assert_eq!(r.whole().unwrap(), seq.clone());
    }

    /// A Reference learns its length from short reads alone, whatever the
    /// block size and whatever was read first, a range past the end
    /// included.
    #[test]
    fn reference_length_agrees_with_the_sequence(
        seq in dna(0, 120),
        block in 1usize..=17,
        first in prop::option::of(0usize..300),
    ) {
        let hdp = Provider::new().sequence("X", &seq);
        let store = ReferenceStore::with_block_size(&hdp, block);
        let r = store.reference("X", G);
        if let Some(i) = first {
            r.base(i).unwrap();
        }
        prop_assert_eq!(r.length().unwrap(), seq.len());
    }

    /// A range is refused exactly when it runs past the end of the sequence;
    /// an empty one, an insertion point, needs the base in front of which it
    /// inserts.
    #[test]
    fn a_range_is_refused_exactly_when_it_runs_past_the_end(
        seq in dna(0, 120),
        block in 1usize..=17,
        (a, b) in (0usize..130, 0usize..130),
    ) {
        let hdp = Provider::new().sequence("X", &seq);
        let store = ReferenceStore::with_block_size(&hdp, block);
        let r = store.reference("X", G);
        let (s, e) = (a.min(b), a.max(b));
        let inside = if e > s { e <= seq.len() } else { s < seq.len() };
        prop_assert_eq!(r.check_range(s, e).is_ok(), inside);
    }
}

// ---------------------------------------------------------------------------
// Equivalence of single variants
// ---------------------------------------------------------------------------

const SYNTHETIC: &str = "NC_000099.1";

/// An edit of any kind the comparison distinguishes, with identities and
/// inversions besides: an identity, an inversion of a range that is its own
/// reverse complement and a delins of the bases already there all make the
/// reference.
fn single_edit(len: usize) -> BoxedStrategy<PlacedHgvs> {
    let span = (1usize..=len.saturating_sub(2).max(1))
        .prop_flat_map(move |l| (0..=len.saturating_sub(l + 1)).prop_map(move |s| (s, s + l)));
    let identity = span.clone().prop_map(|(start, end)| PlacedHgvs {
        start,
        end,
        edit: NaEdit::RefAlt {
            ref_: None,
            alt: None,
            uncertain: false,
        },
    });
    let inversion = span
        .prop_filter("an inversion of one base is a substitution", |(s, e)| {
            e - s > 1
        })
        .prop_map(|(start, end)| PlacedHgvs {
            start,
            end,
            edit: NaEdit::Inv {
                ref_: None,
                uncertain: false,
            },
        });
    prop_oneof![2 => edit_on_kinds(len, false), 1 => identity, 1 => inversion].boxed()
}

/// `p` written as the delins of the bases it touches, widened by `widen`
/// bases 5' and 3' within the sequence: the same change, or, for a variant
/// that makes the reference, a statement about the widened range.
fn respelled(seq: &str, p: &PlacedHgvs, widen: (usize, usize)) -> PlacedHgvs {
    let (s, e) = touched(p);
    let lo = s.saturating_sub(widen.0);
    let mut hi = (e + widen.1).min(seq.len());
    if lo == hi {
        hi = lo + 1;
    }
    let edited = apply(seq, p);
    let alt = &edited[lo..hi + edited.len() - seq.len()];
    let edit = if alt.is_empty() {
        NaEdit::Del {
            ref_: None,
            uncertain: false,
        }
    } else {
        NaEdit::RefAlt {
            ref_: Some(String::new()),
            alt: Some(alt.to_string()),
            uncertain: false,
        }
    };
    PlacedHgvs {
        start: lo,
        end: hi,
        edit,
    }
}

/// A short reference and three single variants on it. Each after the first
/// is either a fresh edit or the one before respelled, so that chains of one
/// change in several spellings, and statements about neighbouring bases, are
/// common enough to test.
fn seq_and_three_variants() -> impl Strategy<Value = (String, Vec<PlacedHgvs>)> {
    dna(6, 16).prop_flat_map(|seq| {
        let len = seq.len();
        let next = (prop::bool::ANY, single_edit(len), 0usize..=1, 0usize..=1);
        (single_edit(len), next.clone(), next).prop_map(move |(a, b, c)| {
            let mut variants = vec![a];
            for (fresh, edit, l, r) in [b, c] {
                let previous = variants.last().expect("starts with one");
                let variant = if fresh {
                    edit
                } else {
                    respelled(&seq, previous, (l, r))
                };
                variants.push(variant);
            }
            (seq.clone(), variants)
        })
    })
}

/// Every ordered pair's equivalence, by index.
fn pairwise_equivalence(seq: &str, variants: &[PlacedHgvs]) -> Vec<Vec<bool>> {
    let hdp = Provider::new().sequence(SYNTHETIC, seq);
    let mapper = VariantMapper::new(&hdp);
    let eq = VariantEquivalence::new(&mapper, &hdp);
    let parsed: Vec<SequenceVariant> = variants
        .iter()
        .map(|p| SequenceVariant::Genomic(genomic_variant(SYNTHETIC, p)))
        .collect();
    parsed
        .iter()
        .map(|x| {
            parsed
                .iter()
                .map(|y| {
                    eq.equivalent(x, y)
                        .unwrap_or_else(|e| panic!("{x} vs {y}: {e}"))
                })
                .collect()
        })
        .collect()
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(512))]

    /// Two single variants that are each the same change as a third are the
    /// same change as each other, whether they change the sequence or state
    /// its bases.
    #[test]
    fn equivalence_of_single_variants_is_transitive((seq, variants) in seq_and_three_variants()) {
        let same = pairwise_equivalence(&seq, &variants);
        for (x, y, z) in [(0, 1, 2), (0, 2, 1), (1, 0, 2)] {
            if same[x][y] && same[y][z] {
                prop_assert!(same[x][z], "on {}: {:?} ~ {:?} ~ {:?}", seq, variants[x], variants[y], variants[z]);
            }
        }
        for (x, y) in [(0, 1), (0, 2), (1, 2)] {
            prop_assert_eq!(same[x][y], same[y][x], "on {}: not symmetric", seq);
        }
    }

    /// A variant that makes the reference is never the same as one that
    /// changes it.
    #[test]
    fn a_single_statement_never_equals_a_change((seq, variants) in seq_and_three_variants()) {
        let same = pairwise_equivalence(&seq, &variants);
        let states = |p: &PlacedHgvs| apply(&seq, p) == seq;
        for x in 0..3 {
            for y in 0..3 {
                if states(&variants[x]) != states(&variants[y]) {
                    prop_assert!(!same[x][y], "on {}: {:?} equals {:?}", seq, variants[x], variants[y]);
                }
            }
        }
    }

    /// Two variants that make the reference are the same exactly when they
    /// name the same bases, however each is written.
    #[test]
    fn single_statements_are_the_same_when_they_name_the_same_bases(
        (seq, variants) in seq_and_three_variants(),
    ) {
        let same = pairwise_equivalence(&seq, &variants);
        let states = |p: &PlacedHgvs| apply(&seq, p) == seq;
        for x in 0..3 {
            for y in 0..3 {
                let (a, b) = (&variants[x], &variants[y]);
                if states(a) && states(b) {
                    prop_assert_eq!(same[x][y], (a.start, a.end) == (b.start, b.end), "on {}: {:?} vs {:?}", seq, a, b);
                }
            }
        }
    }
}

/// A sequence holding a single variant that makes the reference without
/// being written as `=`, and the half-open range of the bases it names: an
/// inversion of a palindrome, a repeat at the count its run already has, or
/// a delins of the bases already there.
fn seq_and_statement() -> impl Strategy<Value = (String, PlacedHgvs, (usize, usize))> {
    (dna(0, 12), dna(1, 3), 1usize..=4, dna(0, 12), 0u8..3).prop_map(
        |(prefix, unit, copies, suffix, kind)| {
            let at = prefix.len();
            let middle = match kind {
                0 => format!("{unit}{}", revcomp(&unit)),
                _ => unit.repeat(copies),
            };
            let seq = format!("{prefix}{middle}{suffix}");
            let (placed, named) = match kind {
                0 => {
                    let edit = NaEdit::Inv {
                        ref_: None,
                        uncertain: false,
                    };
                    let end = at + middle.len();
                    (
                        PlacedHgvs {
                            start: at,
                            end,
                            edit,
                        },
                        (at, end),
                    )
                }
                1 => {
                    // The run read 3' from `at`, which the suffix can continue.
                    let mut run = at;
                    while seq[run..].starts_with(unit.as_str()) {
                        run += unit.len();
                    }
                    let count = ((run - at) / unit.len()) as i32;
                    let edit = NaEdit::Repeat {
                        ref_: Some(unit.clone()),
                        min: count,
                        max: count,
                        uncertain: false,
                    };
                    (
                        PlacedHgvs {
                            start: at,
                            end: at + unit.len(),
                            edit,
                        },
                        (at, run),
                    )
                }
                _ => {
                    let edit = NaEdit::RefAlt {
                        ref_: Some(String::new()),
                        alt: Some(middle.clone()),
                        uncertain: false,
                    };
                    let end = at + middle.len();
                    (
                        PlacedHgvs {
                            start: at,
                            end,
                            edit,
                        },
                        (at, end),
                    )
                }
            };
            (seq, placed, named)
        },
    )
}

proptest! {
    /// A single variant that makes the reference normalises to the identity
    /// of the bases it names, and that identity normalises to itself.
    #[test]
    fn a_single_variant_that_makes_the_reference_normalises_to_its_statement(
        (seq, p, named) in seq_and_statement(),
    ) {
        let hdp = Provider::new().sequence("X", &seq);
        let store = ReferenceStore::with_block_size(&hdp, 4);
        let r = store.reference("X", G);
        let identity = PlacedEdit {
            start: named.0,
            end: named.1,
            edit: NaEdit::RefAlt { ref_: None, alt: None, uncertain: false },
        };
        let before = PlacedEdit { start: p.start, end: p.end, edit: p.edit.clone() };
        let after = normalize(&r, before, &ExonBounds::unspliced()).unwrap();
        prop_assert_eq!(&after, &identity, "on {}: {:?}", seq, p);
        let again = normalize(&r, after.clone(), &ExonBounds::unspliced()).unwrap();
        prop_assert_eq!(again, after);
    }
}

/// The g. variant on `ac` for an edit over an HGVS range.
fn genomic_variant(ac: &str, p: &PlacedHgvs) -> GVariant {
    let position = |i: usize| serde_json::json!({"base": i, "end": null, "uncertain": false});
    let end = (p.end > p.start + 1).then(|| position(p.end));
    variant_from(ac, interval(position(p.start + 1), end), &p.edit)
}

/// The interval from `start`, to `end` where there is one.
fn interval<P: serde::Serialize>(start: P, end: Option<P>) -> serde_json::Value {
    serde_json::json!({"start": start, "end": end, "uncertain": false})
}

/// The variant on `ac` that changes `pos` by `edit`, built from its fields
/// through serde: outside the crate, the structs cannot be written literally.
fn variant_from<V: serde::de::DeserializeOwned>(
    ac: &str,
    pos: serde_json::Value,
    edit: &NaEdit,
) -> V {
    serde_json::from_value(serde_json::json!({
        "ac": ac,
        "gene": null,
        "posedit": {"pos": pos, "edit": edit, "uncertain": false, "predicted": false},
    }))
    .unwrap()
}

/// The reference bases an edit written over an HGVS range replaces, as a
/// half-open range: empty, at the point the new bases go, for an insertion or
/// a duplication.
fn touched(p: &PlacedHgvs) -> (usize, usize) {
    match p.edit {
        NaEdit::Ins { .. } => (p.end - 1, p.end - 1),
        NaEdit::Dup { .. } => (p.end, p.end),
        _ => (p.start, p.end),
    }
}

/// The HGVS range a placed (normalised) edit is written over.
fn placed_to_hgvs(p: &PlacedEdit) -> PlacedHgvs {
    let (s, e) = p.hgvs_range();
    PlacedHgvs {
        start: s,
        end: e,
        edit: p.edit.clone(),
    }
}

/// For an edit given by its resolved (placed) range, the HGVS range: an
/// insertion at `a` is written between `a - 1` and `a`.
fn as_hgvs_range(p: &PlacedHgvs) -> PlacedHgvs {
    if matches!(p.edit, hgvs_weaver::edits::NaEdit::Ins { .. }) && p.start == p.end {
        PlacedHgvs {
            start: p.start - 1,
            end: p.start + 1,
            edit: p.edit.clone(),
        }
    } else {
        p.clone()
    }
}

// ---------------------------------------------------------------------------
// Coordinates
// ---------------------------------------------------------------------------

use hgvs_weaver::coords::{IntronicOffset, TranscriptPos};
use hgvs_weaver::edits::{AaEdit, NaEdit};
use hgvs_weaver::equivalence::{EquivalenceLevel, VariantEquivalence};
use hgvs_weaver::mapper::VariantMapper;
use hgvs_weaver::structs::{CVariant, PVariant, TranscriptVariant};
use hgvs_weaver::transcript_mapper::TranscriptMapper;
use hgvs_weaver::utils::{aa1_to_aa3, aa3_to_aa1, translate, translate_codon};
use hgvs_weaver::{parse_hgvs_variant, SequenceVariant};

fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        o => o,
    }
}

proptest! {
    /// Every transcript base maps to a genomic base holding the same residue
    /// (complemented on the minus strand), maps back to itself with no
    /// offset, and survives the trip through c. numbering.
    #[test]
    fn transcript_coordinates_round_trip(g in gene()) {
        let am = TranscriptMapper::new(g.transcript_data()).unwrap();
        for i in 0..g.transcript_seq.len() {
            let gpos = am.n_to_g(TranscriptPos(i as i32), IntronicOffset(0)).unwrap();
            let tb = g.transcript_seq.as_bytes()[i];
            let gb = g.genome.as_bytes()[gpos.0 as usize];
            let expected = if g.strand == hgvs_weaver::data::Strand::Minus { complement(tb) } else { tb };
            prop_assert_eq!(gb, expected, "transcript index {} lands on the wrong genomic base", i);
            let (n, off) = am.g_to_n(gpos).unwrap();
            prop_assert_eq!((n, off), (TranscriptPos(i as i32), IntronicOffset(0)));
            let (c, coff, anchor) = am.n_to_c(TranscriptPos(i as i32)).unwrap();
            prop_assert_eq!(coff, IntronicOffset(0));
            prop_assert_eq!(am.c_to_n(c, anchor).unwrap(), TranscriptPos(i as i32));
        }
    }

    /// Counting into an intron from the 3' end of one exon agrees with counting
    /// back from the 5' end of the next: `n+i` and `(n+1)-(L+1-i)` are one base.
    #[test]
    fn intron_offsets_agree_from_either_exon(g in gene()) {
        let am = TranscriptMapper::new(g.transcript_data()).unwrap();
        for (k, &len) in g.introns.iter().enumerate() {
            let last = g.tx_exons[k].1 - 1; // last base of exon k
            let first = g.tx_exons[k + 1].0; // first base of exon k + 1
            prop_assert_eq!(first, last + 1);
            for i in 1..=len {
                let from_left = am.n_to_g(TranscriptPos(last as i32), IntronicOffset(i as i32)).unwrap();
                let from_right = am
                    .n_to_g(TranscriptPos(first as i32), IntronicOffset(-((len + 1 - i) as i32)))
                    .unwrap();
                prop_assert_eq!(from_left, from_right, "intron {} base {} disagrees", k, i);
                // It is an intron base (the generator fills introns with T).
                prop_assert_eq!(g.genome.as_bytes()[from_left.0 as usize], b'T');
                // And the same through c. numbering with an intronic offset.
                let mut pos = CVariant::position_from_index(&am, last as i32).unwrap();
                pos.offset = Some(IntronicOffset(i as i32));
                prop_assert_eq!(am.position_to_g(&pos).unwrap(), from_left);
                // And back from the genome: the base is counted from the nearer
                // exon boundary, from the earlier exon on a tie, on either strand.
                let expected = if i <= len + 1 - i {
                    (TranscriptPos(last as i32), IntronicOffset(i as i32))
                } else {
                    (TranscriptPos(first as i32), IntronicOffset(-((len + 1 - i) as i32)))
                };
                prop_assert_eq!(am.g_to_n(from_left).unwrap(), expected, "intron {} base {} maps back wrongly", k, i);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Protein
// ---------------------------------------------------------------------------

/// A c. variant over transcript indices `[start, end)` (HGVS range) on `g`.
fn coding_variant(
    g: &Gene,
    am: &TranscriptMapper,
    start: usize,
    end: usize,
    edit: NaEdit,
) -> CVariant {
    let s = CVariant::position_from_index(am, start as i32).unwrap();
    let e = if end > start + 1 {
        Some(CVariant::position_from_index(am, (end - 1) as i32).unwrap())
    } else {
        None
    };
    variant_from(&g.ac, interval(s, e), &edit)
}

/// Applies a p. description to the reference protein, the dumb way. Returns
/// `None` for descriptions the oracle does not model (frameshift, extension).
fn apply_protein(ref_aa: &str, vp: &PVariant) -> Option<String> {
    let aa: Vec<char> = ref_aa.chars().collect();
    let Some(pos) = &vp.posedit.pos else {
        return Some(ref_aa.to_string());
    };
    let s = pos.start.base.to_index().0 as usize;
    let e = pos.end.as_ref().map_or(s, |p| p.base.to_index().0 as usize);
    // Three-letter codes concatenated, e.g. "LeuSer" -> "LS".
    let one = |aa3: &str| -> String {
        aa3.as_bytes()
            .chunks(3)
            .map(|c| aa3_to_aa1(std::str::from_utf8(c).unwrap()))
            .collect()
    };
    let (cut_s, cut_e, insert): (usize, usize, String) = match &vp.posedit.edit {
        AaEdit::Identity { .. } => return Some(ref_aa.to_string()),
        AaEdit::Subst { alt, .. } => (s, s + 1, one(alt)),
        AaEdit::Del { .. } => (s, e + 1, String::new()),
        AaEdit::Dup { .. } => (e + 1, e + 1, aa[s..=e].iter().collect()),
        AaEdit::Ins { alt, .. } => (s + 1, s + 1, one(alt)),
        AaEdit::DelIns { alt, .. } => (s, e + 1, one(alt)),
        _ => return None,
    };
    let mut out: String = aa[..cut_s].iter().collect();
    out.push_str(&insert);
    out.extend(aa[cut_e.min(aa.len())..].iter());
    Some(out)
}

proptest! {
    // Three descriptions must meet at one residue for a violation to show.
    #![proptest_config(ProptestConfig::with_cases(2048))]

    /// Judged on the protein, being the same change is transitive: over a
    /// coding substitution at one codon (silent, missense, nonsense or a
    /// stop loss), its prediction, silent and substituted descriptions of
    /// that codon's residue, and the whole protein unchanged, two that are
    /// each the same as a third are the same as each other. Statements that
    /// say only where a change starts, `p.Met1?`, are not generated: they
    /// match changes that differ from each other. So the start codon, whose
    /// every change predicts `p.(Met1?)`, is not a site.
    #[test]
    fn protein_side_equivalence_is_transitive(
        g in gene(),
        site in (1usize..3, 0usize..3, prop::sample::select(BASES.to_vec())),
        picks in prop::collection::vec((0u8..9, prop::sample::select(RESIDUES.to_vec())), 3),
    ) {
        let am = TranscriptMapper::new(g.transcript_data()).unwrap();
        let hdp = g.provider();
        let mapper = VariantMapper::new(&hdp);
        let eq = VariantEquivalence::new(&mapper, &hdp);
        let mut described = Vec::new();
        for &(kind, residue) in &picks {
            let Some(v) = protein_side(&g, &am, &mapper, site, kind, residue) else { return Ok(()) };
            described.push(v);
        }
        let same = |a: &SequenceVariant, b: &SequenceVariant| {
            eq.equivalent_level(a, b).unwrap_or_else(|e| panic!("{a} vs {b}: {e}")).is_equivalent()
        };
        for [x, y, z] in [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]] {
            let (x, y, z) = (&described[x], &described[y], &described[z]);
            if same(x, y) && same(y, z) {
                prop_assert!(same(x, z), "{x} is {y}, and {y} is {z}, but {x} is not {z}");
            }
        }
    }
}

/// The residues a generated protein substitution may name.
const RESIDUES: [char; 21] = [
    'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W',
    'Y', '*',
];

/// A description of kind `kind` on `g` at `(codon, base in codon, base)`,
/// judged on its protein; codon 2 stands for the stop codon. Kind 0 is the
/// coding substitution, 7 one at the same offset in the 3'UTR; 1 is the
/// stated form of kind 0's prediction; at the codon's residue, 2 is
/// `p.Xxx#=`, 3 the substitution by itself, 4 that by `residue`; 5 and 6 are
/// `p.=` and `p.(=)`, and 8 every residue through the stop unchanged. `None`
/// for a substitution by the base already there.
fn protein_side(
    g: &Gene,
    am: &TranscriptMapper,
    mapper: &VariantMapper,
    (codon, offset, base): (usize, usize, char),
    kind: u8,
    residue: char,
) -> Option<SequenceVariant> {
    let protein = translate(&g.transcript_seq[g.cds_start..=g.cds_end]);
    let codon = if codon == 2 { protein.len() - 1 } else { codon };
    let at = match kind {
        7 => g.cds_end + 1 + offset,
        _ => g.cds_start + 3 * codon + offset,
    };
    let reference = g.transcript_seq.as_bytes()[at] as char;
    let snv = || {
        (base != reference).then(|| {
            coding_variant(
                g,
                am,
                at,
                at + 1,
                NaEdit::RefAlt {
                    ref_: Some(reference.to_string()),
                    alt: Some(base.to_string()),
                    uncertain: false,
                },
            )
        })
    };
    let name = aa1_to_aa3(protein.as_bytes()[codon] as char);
    let n = codon + 1;
    let p = |edit: String| parse_hgvs_variant(&format!("NP_PROP.1:p.{edit}")).unwrap();
    Some(match kind {
        0 | 7 => SequenceVariant::Coding(snv()?),
        1 => {
            let predicted = mapper.c_to_p(&snv()?, Some("NP_PROP.1")).unwrap();
            p(predicted
                .to_string()
                .split_once(":p.")
                .unwrap()
                .1
                .replace(['(', ')'], ""))
        }
        2 => p(format!("{name}{n}=")),
        3 => p(format!("{name}{n}{name}")),
        4 => p(format!("{name}{n}{}", aa1_to_aa3(residue))),
        5 => p("=".to_string()),
        6 => p("(=)".to_string()),
        _ => p(format!("Met1_Ter{}=", protein.len())),
    })
}

fn up_to_stop(s: &str) -> String {
    match s.find('*') {
        Some(i) => s[..=i].to_string(),
        None => s.to_string(),
    }
}

// ---------------------------------------------------------------------------
// A coding variant against a protein description
//
// Judged directly, a coding variant is read by the protein its edited
// transcript encodes; judged through its prediction, by what c_to_p's
// description says of that protein. The direct judgement uses the same rules
// on more information, so it is never the less equivalent of the two.
// ---------------------------------------------------------------------------

/// How equivalent `level` says two descriptions are, from Different up.
fn rank(level: EquivalenceLevel) -> u8 {
    match level {
        EquivalenceLevel::Identity => 3,
        EquivalenceLevel::Analogous => 2,
        EquivalenceLevel::Different => 1,
        EquivalenceLevel::Unknown => 0,
    }
}

/// A prediction's edit without its parentheses: `Leu3ThrfsTer8` of
/// `NP_PROP.1:p.(Leu3ThrfsTer8)`.
fn unbracketed(vp: &PVariant) -> String {
    let mut vp = vp.clone();
    vp.posedit.predicted = false;
    vp.to_string().split_once(":p.").unwrap().1.to_string()
}

/// The spellings a description's open forms take in ClinVar: a frameshift
/// without its length, or without its first new residue as well; an
/// extension without its length, or as a substitution of the stop. `None`
/// for a description with neither.
fn open_spelling(edit: &str, which: u8) -> Option<String> {
    if let Some(i) = edit.find("fs") {
        return Some(match which % 2 {
            0 => edit[..i + 2].to_string(),
            _ => format!("{}fs", &edit[..i.checked_sub(3)?]),
        });
    }
    let i = edit.find("ext")?;
    Some(match which % 2 {
        0 => format!("{}extTer?", &edit[..i]),
        _ => edit[..i].to_string(),
    })
}

/// A protein description on `NP_PROP.1`, whose residues before the stop are
/// `protein`, of kind `kind`: 0 and 1 `other`, the prediction of another
/// change, as written and unbracketed; 2 `own`, the prediction under test,
/// unbracketed; 3 and 4 an open spelling of `own` and of `other`; then, at
/// the residue `at` names, a substitution by `residues`' first (5), a
/// deletion of `count + 1` residues (6), their replacement by `residues`
/// (7), a frameshift to `residues`' first with `length` as its length (8), a
/// stop loss to it (9), a change of unknown kind from there (10), statements
/// about the whole protein (11), a duplication (12) and an insertion of
/// `residues` after it (13). `None` where the kind names nothing.
#[allow(clippy::too_many_arguments)]
fn protein_description(
    protein: &str,
    own: &PVariant,
    other: Option<&PVariant>,
    kind: u8,
    at: usize,
    count: usize,
    residues: &[char],
    length: Option<usize>,
) -> Option<String> {
    let at = at % protein.len();
    let end = (at + count).min(protein.len() - 1);
    let name = |i: usize| format!("{}{}", aa1_to_aa3(protein.as_bytes()[i] as char), i + 1);
    let new: String = residues.iter().map(|&r| aa1_to_aa3(r)).collect();
    let first = aa1_to_aa3(residues[0]);
    let range = if end > at {
        format!("{}_{}", name(at), name(end))
    } else {
        name(at)
    };
    let length = length.map_or("?".to_string(), |n| n.to_string());
    let edit = match kind {
        0 => return Some(other?.to_string()),
        1 => unbracketed(other?),
        2 => unbracketed(own),
        3 => open_spelling(&unbracketed(own), count as u8)?,
        4 => open_spelling(&unbracketed(other?), count as u8)?,
        5 => format!("{}{first}", name(at)),
        6 => format!("{range}del"),
        7 => format!("{range}delins{new}"),
        8 => format!("{}{first}fsTer{length}", name(at)),
        9 => format!("Ter{}{first}extTer{length}", protein.len() + 1),
        10 => format!("{}?", name(at)),
        11 => ["=", "?", "0", "(=)"][count % 4].to_string(),
        12 => format!("{range}dup"),
        _ if end + 1 < protein.len() => format!("{}_{}ins{new}", name(end), name(end + 1)),
        _ => format!("{}_Ter{}ins{new}", name(end), end + 2),
    };
    Some(format!("NP_PROP.1:p.{edit}"))
}

proptest! {
    /// A coding variant judged against a protein description is at least as
    /// equivalent as its prediction judged against it, so the direct
    /// judgement says Different only where the prediction's does; except
    /// that a description naming residues a predicted frameshift, extension or
    /// loss of the stop leaves unnamed may match the prediction and not the
    /// protein.
    #[test]
    fn a_coding_variant_is_at_least_as_equivalent_as_its_prediction(
        (g, e, other) in gene_and_two_edits(),
        kind in 0u8..14,
        at in 0usize..10_000,
        count in 0usize..4,
        residues in prop::collection::vec(prop::sample::select(RESIDUES.to_vec()), 1..=4),
        length in prop::option::of(1usize..12),
    ) {
        let am = TranscriptMapper::new(g.transcript_data()).unwrap();
        let hdp = g.provider();
        let mapper = VariantMapper::new(&hdp);
        let eq = VariantEquivalence::new(&mapper, &hdp);
        let vc = tx_edit_variant(&g, &am, &e);
        let Ok(predicted) = mapper.c_to_p(&vc, Some("NP_PROP.1")) else { return Ok(()) };
        let other = mapper.c_to_p(&tx_edit_variant(&g, &am, &other), Some("NP_PROP.1")).ok();
        let protein = translate(&g.transcript_seq[g.cds_start..g.cds_end - 2]);
        let Some(text) = protein_description(
            &protein, &predicted, other.as_ref(), kind, at, count, &residues, length,
        ) else { return Ok(()) };
        let vp = parse_hgvs_variant(&text).unwrap_or_else(|e| panic!("{text}: {e}"));
        let through_prediction = eq.equivalent_level(&SequenceVariant::Protein(predicted.clone()), &vp);
        let direct = eq.equivalent_level(&SequenceVariant::Coding(vc.clone()), &vp);
        // A prediction that does not read back on its own protein has no
        // judgement to compare with: c_to_p names read-through residues past
        // the stop for an in-frame change over it (issue #94), as residues the
        // protein does not hold or as positions past its end.
        let through_prediction = match through_prediction {
            Ok(level) => level,
            Err(
                hgvs_weaver::HgvsError::ReferenceMismatch { .. }
                | hgvs_weaver::HgvsError::PositionOutsideSequence { .. },
            ) => return Ok(()),
            Err(e) => return Err(TestCaseError::fail(format!("{predicted} against {vp}: {e}"))),
        };
        // A predicted frameshift or extension, or a change that removes the
        // stop, leaves residues it does not name, which match any residue, so
        // a description can match those through the prediction and not the
        // residues the transcript encodes. There only Identity carries over.
        let floor = match through_prediction {
            EquivalenceLevel::Analogous if leaves_unnamed(&predicted) && names_residues(kind) => {
                EquivalenceLevel::Different
            }
            level => level,
        };
        prop_assert!(
            direct.as_ref().is_ok_and(|&d| rank(d) >= rank(floor)),
            "{vc} against {vp}: {direct:?} directly, {through_prediction:?} through {predicted}"
        );
    }

    /// A coding variant is equivalent to its predicted frameshift or
    /// extension written open, as ClinVar writes them: without the length,
    /// without the first new residue, or a stop loss as a substitution.
    #[test]
    fn a_coding_variant_is_its_prediction_written_open(
        (g, e) in gene_and_edit(),
        which in 0u8..2,
    ) {
        let am = TranscriptMapper::new(g.transcript_data()).unwrap();
        let hdp = g.provider();
        let mapper = VariantMapper::new(&hdp);
        let eq = VariantEquivalence::new(&mapper, &hdp);
        let vc = tx_edit_variant(&g, &am, &e);
        let Ok(predicted) = mapper.c_to_p(&vc, Some("NP_PROP.1")) else { return Ok(()) };
        let Some(open) = open_spelling(&unbracketed(&predicted), which) else { return Ok(()) };
        let vp = parse_hgvs_variant(&format!("NP_PROP.1:p.{open}")).unwrap_or_else(|e| panic!("{open}: {e}"));
        let direct = eq.equivalent_level(&SequenceVariant::Coding(vc.clone()), &vp);
        prop_assert!(
            direct.as_ref().is_ok_and(EquivalenceLevel::is_equivalent),
            "{vc} against {vp}: {direct:?}, predicted {predicted}"
        );
    }
}

/// Whether a [`protein_description`] of `kind` can name residues another
/// change's prediction leaves unnamed: all but the prediction's own spellings.
fn names_residues(kind: u8) -> bool {
    !matches!(kind, 2 | 3)
}

/// Whether a prediction leaves residues it does not name: a frameshift or an
/// extension, or a change over the stop.
fn leaves_unnamed(predicted: &PVariant) -> bool {
    let over_the_stop = predicted
        .posedit
        .pos
        .as_ref()
        .is_some_and(|p| p.start.aa == "Ter" || p.end.as_ref().is_some_and(|e| e.aa == "Ter"));
    over_the_stop
        || matches!(
            predicted.posedit.edit,
            AaEdit::Fs { .. } | AaEdit::Ext { .. }
        )
}

proptest! {
    /// A single-base substitution in the CDS is described exactly as the codon
    /// table says: silent, missense, nonsense or a stop-loss extension; in the
    /// start codon, an unknown change from Met1.
    #[test]
    fn snv_effect_follows_the_codon_table(g in gene(), pick in 0usize..10_000, alt in prop::sample::select(BASES.to_vec())) {
        let cds_len = g.cds_end + 1 - g.cds_start;
        let pos = g.cds_start + pick % cds_len;
        let ref_base = g.transcript_seq.as_bytes()[pos] as char;
        prop_assume!(alt != ref_base);
        let am = TranscriptMapper::new(g.transcript_data()).unwrap();
        let hdp = g.provider();
        let mapper = VariantMapper::new(&hdp);
        let vc = coding_variant(&g, &am, pos, pos + 1,
            NaEdit::RefAlt { ref_: Some(ref_base.to_string()), alt: Some(alt.to_string()), uncertain: false });
        let vp = mapper.c_to_p(&vc, Some("NP_PROP.1")).unwrap().to_string();

        let ci = (pos - g.cds_start) / 3;
        let codon = &g.transcript_seq[g.cds_start + 3 * ci..g.cds_start + 3 * ci + 3];
        let mut alt_codon: Vec<char> = codon.chars().collect();
        alt_codon[(pos - g.cds_start) % 3] = alt;
        let alt_codon: String = alt_codon.into_iter().collect();
        let (r, a) = (translate_codon(codon).unwrap(), translate_codon(&alt_codon).unwrap());
        let n = ci + 1;
        let expected = if ci == 0 {
            "NP_PROP.1:p.(Met1?)".to_string()
        } else if r == a {
            format!("NP_PROP.1:p.({}{}=)", aa1_to_aa3(r), n)
        } else if a == '*' {
            format!("NP_PROP.1:p.({}{}Ter)", aa1_to_aa3(r), n)
        } else if r == '*' {
            format!("NP_PROP.1:p.(Ter{}{}ext", n, aa1_to_aa3(a))
        } else {
            format!("NP_PROP.1:p.({}{}{})", aa1_to_aa3(r), n, aa1_to_aa3(a))
        };
        prop_assert!(vp.starts_with(&expected), "got {vp}, expected {expected}");
    }

    /// An in-frame deletion or insertion of whole codons inside the CDS is
    /// described as a change that, applied to the reference protein,
    /// reproduces the translated alternate. Never as a frameshift.
    #[test]
    fn in_frame_indels_describe_the_translated_protein(
        g in gene(),
        pick in 0usize..10_000,
        codons in 1usize..=3,
        insert in prop::collection::vec(prop::sample::select(BASES.to_vec()), 3..=9),
        kind in 0u8..3,
    ) {
        let n_codons = (g.cds_end + 1 - g.cds_start) / 3;
        prop_assume!(n_codons > codons + 2);
        // an interior codon index: not the start codon, not the stop
        let ci = 1 + pick % (n_codons - 1 - codons);
        let start = g.cds_start + 3 * ci;
        let am = TranscriptMapper::new(g.transcript_data()).unwrap();
        let hdp = g.provider();
        let mapper = VariantMapper::new(&hdp);

        let (vc, alt_tx) = if kind == 0 {
            let end = start + 3 * codons;
            (coding_variant(&g, &am, start, end, NaEdit::Del { ref_: None, uncertain: false }),
             format!("{}{}", &g.transcript_seq[..start], &g.transcript_seq[end..]))
        } else if kind == 1 {
            // duplication of whole codons
            let end = start + 3 * codons;
            (coding_variant(&g, &am, start, end, NaEdit::Dup { ref_: None, uncertain: false }),
             format!("{}{}{}", &g.transcript_seq[..end], &g.transcript_seq[start..end], &g.transcript_seq[end..]))
        } else {
            let ins: String = insert.iter().take(insert.len() / 3 * 3).collect();
            prop_assume!(!ins.is_empty());
            // between the last base of codon ci-1 and the first of codon ci
            (coding_variant(&g, &am, start - 1, start + 1, NaEdit::Ins { alt: Some(ins.clone()), uncertain: false }),
             format!("{}{}{}", &g.transcript_seq[..start], ins, &g.transcript_seq[start..]))
        };
        let vp = mapper.c_to_p(&vc, Some("NP_PROP.1")).unwrap();
        let described = vp.to_string();
        prop_assert!(!described.contains("fs"), "in-frame change described as frameshift: {described}");

        let ref_aa = translate(&g.transcript_seq[g.cds_start..]);
        let alt_aa = translate(&alt_tx[g.cds_start..]);
        let Some(applied) = apply_protein(&ref_aa, &vp) else {
            prop_assert!(false, "unmodelled description {described}");
            unreachable!()
        };
        prop_assert_eq!(up_to_stop(&applied), up_to_stop(&alt_aa), "{} does not describe the translated protein", described);
    }
}

proptest! {
    /// HGVS writes an in-frame insertion or duplication at its 3'-most
    /// equivalent residues: moving the described residues one position further
    /// towards the C terminus must give a different protein, or the
    /// description was not fully shifted.
    #[test]
    fn protein_insertions_are_written_3_prime_most(
        g in gene(),
        pick in 0usize..10_000,
        codons in 1usize..=3,
        insert in prop::collection::vec(prop::sample::select(BASES.to_vec()), 3..=9),
        dup in prop::bool::ANY,
    ) {
        let n_codons = (g.cds_end + 1 - g.cds_start) / 3;
        prop_assume!(n_codons > codons + 2);
        let ci = 1 + pick % (n_codons - 1 - codons);
        let start = g.cds_start + 3 * ci;
        let am = TranscriptMapper::new(g.transcript_data()).unwrap();
        let hdp = g.provider();
        let mapper = VariantMapper::new(&hdp);
        let vc = if dup {
            coding_variant(&g, &am, start, start + 3 * codons, NaEdit::Dup { ref_: None, uncertain: false })
        } else {
            let ins: String = insert.iter().take(insert.len() / 3 * 3).collect();
            prop_assume!(!ins.is_empty());
            coding_variant(&g, &am, start - 1, start + 1, NaEdit::Ins { alt: Some(ins), uncertain: false })
        };
        let vp = mapper.c_to_p(&vc, Some("NP_PROP.1")).unwrap();
        let ref_aa = translate(&g.transcript_seq[g.cds_start..]);
        let Some(pos) = &vp.posedit.pos else { return Ok(()) };
        let s = pos.start.base.to_index().0 as usize;
        let e = pos.end.as_ref().map_or(s, |p| p.base.to_index().0 as usize);
        let shifted = match &vp.posedit.edit {
            // insertion between s and s+1 -> between s+1 and s+2
            AaEdit::Ins { alt, uncertain } => Some((s + 1, s + 2, AaEdit::Ins { alt: alt.clone(), uncertain: *uncertain })),
            // duplication of s..=e -> of s+1..=e+1
            AaEdit::Dup { ref_, uncertain } => Some((s + 1, e + 1, AaEdit::Dup { ref_: ref_.clone(), uncertain: *uncertain })),
            _ => None,
        };
        let Some((ns, ne, edit)) = shifted else { return Ok(()) };
        if ne >= ref_aa.len() { return Ok(()); }
        let mut moved = vp.clone();
        let p = moved.posedit.pos.as_mut().unwrap();
        p.start.base = hgvs_weaver::coords::ProteinPos(ns as i32).to_hgvs();
        if let Some(end) = p.end.as_mut() { end.base = hgvs_weaver::coords::ProteinPos(ne as i32).to_hgvs(); }
        moved.posedit.edit = edit;
        let (here, there) = (apply_protein(&ref_aa, &vp), apply_protein(&ref_aa, &moved));
        prop_assert!(here.is_some());
        prop_assert_ne!(
            here.clone().map(|x| up_to_stop(&x)),
            there.map(|x| up_to_stop(&x)),
            "{} is not 3'-most: the same protein results one residue further on",
            vp
        );
    }
}

// ---------------------------------------------------------------------------
// Parser and cross-system agreement
// ---------------------------------------------------------------------------

/// Canonically spelled HGVS strings the parser must accept and print back.
fn hgvs_string() -> impl Strategy<Value = String> {
    let bases = |n: usize| {
        prop::collection::vec(prop::sample::select(BASES.to_vec()), 1..=n)
            .prop_map(|v| v.into_iter().collect::<String>())
    };
    let base_pos = 1i32..=5000;
    let c_pos = prop_oneof![
        base_pos.clone().prop_map(|p| p.to_string()),
        (1i32..=300).prop_map(|p| format!("-{p}")),
        (1i32..=300).prop_map(|p| format!("*{p}")),
        (base_pos.clone(), 1i32..=200).prop_map(|(p, o)| format!("{p}+{o}")),
        (base_pos.clone(), 1i32..=200).prop_map(|(p, o)| format!("{p}-{o}")),
    ];
    let range = (1i32..=5000, 1i32..=50).prop_map(|(a, l)| (a, a + l));
    let sys = prop::sample::select(vec![
        ("NC_000001.11", "g"),
        ("NC_012920.1", "m"),
        ("NM_000001.1", "c"),
        ("NR_000001.1", "n"),
    ]);
    let single_pos = move |s: &'static str| -> BoxedStrategy<String> {
        if s == "c" {
            c_pos.clone().boxed()
        } else {
            base_pos.clone().prop_map(|p| p.to_string()).boxed()
        }
    };
    sys.prop_flat_map(move |(ac, s)| {
        let single = single_pos(s);
        prop_oneof![
            (
                single.clone(),
                prop::sample::select(BASES.to_vec()),
                prop::sample::select(BASES.to_vec())
            )
                .prop_filter("sub must change", |(_, r, a)| r != a)
                .prop_map(move |(p, r, a)| format!("{ac}:{s}.{p}{r}>{a}")),
            single.clone().prop_map(move |p| format!("{ac}:{s}.{p}del")),
            single.clone().prop_map(move |p| format!("{ac}:{s}.{p}dup")),
            range
                .clone()
                .prop_map(move |(a, b)| format!("{ac}:{s}.{a}_{b}del")),
            range
                .clone()
                .prop_map(move |(a, b)| format!("{ac}:{s}.{a}_{b}dup")),
            range
                .clone()
                .prop_map(move |(a, b)| format!("{ac}:{s}.{a}_{b}inv")),
            (1i32..=5000, bases(8))
                .prop_map(move |(a, ins)| format!("{ac}:{s}.{a}_{}ins{ins}", a + 1)),
            (range.clone(), bases(8))
                .prop_map(move |((a, b), ins)| format!("{ac}:{s}.{a}_{b}delins{ins}")),
        ]
    })
}

proptest! {
    #[test]
    fn parser_round_trips_canonical_hgvs(s in hgvs_string()) {
        let v = parse_hgvs_variant(&s).unwrap_or_else(|e| panic!("failed to parse {s}: {e}"));
        prop_assert_eq!(v.to_string(), s);
    }

    #[test]
    fn parser_never_panics(s in "\\PC{0,40}") {
        let _ = parse_hgvs_variant(&s);
    }

    /// The protein allele of a coding change, applied to the protein, is the
    /// protein the edited transcript encodes, to its first stop, for every
    /// kind of change: silent, missense, nonsense, in-frame, frameshift,
    /// stop loss. A change to the start codon predicts no protein sequence,
    /// so it has no allele.
    #[test]
    fn the_protein_allele_is_the_edited_translation(g in gene(), pick in 0usize..10_000, alt in prop::sample::select(BASES.to_vec()), kind in 0u8..3) {
        let am = TranscriptMapper::new(g.transcript_data()).unwrap();
        let hdp = g.provider();
        let mapper = VariantMapper::new(&hdp);
        let pos = g.cds_start + pick % (g.cds_end + 1 - g.cds_start);
        let ref_base = g.transcript_seq.as_bytes()[pos] as char;
        let (end, edit) = match kind {
            0 => { prop_assume!(alt != ref_base); (pos + 1, NaEdit::RefAlt { ref_: Some(ref_base.to_string()), alt: Some(alt.to_string()), uncertain: false }) }
            1 => (pos + 1, NaEdit::Del { ref_: None, uncertain: false }),
            _ => (pos + 1, NaEdit::Ins { alt: Some(alt.to_string()), uncertain: false }),
        };
        // An insertion is written between two bases; keep it inside the CDS.
        let (start, end) = if kind == 2 { if pos + 1 > g.cds_end { return Ok(()) } (pos, pos + 2) } else { (pos, end) };
        let vc = coding_variant(&g, &am, start, end, edit.clone());
        let edited = apply(&g.transcript_seq, &PlacedHgvs { start, end, edit });
        let first_codon = g.cds_start..g.cds_start + 3;
        if edited.get(first_codon.clone()) != g.transcript_seq.get(first_codon) {
            let allele = mapper.protein_allele(&vc, Some("NP_PROP.1"));
            prop_assert!(matches!(allele, Err(hgvs_weaver::HgvsError::UnsupportedOperation(_))), "{}: {:?}", vc, allele);
            return Ok(());
        }
        let allele = mapper.protein_allele(&vc, Some("NP_PROP.1")).unwrap();
        let protein = hgvs_weaver::utils::translate(&g.transcript_seq[g.cds_start..=g.cds_end]);
        let protein = protein.trim_end_matches('*');
        prop_assert_eq!(&allele.accession, "NP_PROP.1");
        let applied = format!("{}{}{}", &protein[..allele.start], allele.alternate, &protein[allele.end..]);
        let expected: String = hgvs_weaver::utils::translate(&edited[g.cds_start..]).split('*').next().unwrap().to_string();
        prop_assert_eq!(applied, expected, "{} ({})", vc, mapper.c_to_p(&vc, Some("NP_PROP.1")).map(|p| p.to_string()).unwrap_or_else(|e| e.to_string()));
    }

    /// r. is the c. (or n.) spelling in RNA letters: converting there and back
    /// is the identity, the letters are lowercase with u, and a change inside
    /// one exon has the same allele and the same protein either way.
    #[test]
    fn rna_is_the_transcript_in_other_letters(g in gene(), pick in 0usize..10_000, alt in prop::sample::select(BASES.to_vec()), del in prop::bool::ANY) {
        let am = TranscriptMapper::new(g.transcript_data()).unwrap();
        let hdp = g.provider();
        let mapper = VariantMapper::new(&hdp);
        let pos = g.cds_start + pick % (g.cds_end + 1 - g.cds_start);
        let ref_base = g.transcript_seq.as_bytes()[pos] as char;
        let edit = if del {
            NaEdit::Del { ref_: None, uncertain: false }
        } else {
            prop_assume!(alt != ref_base);
            NaEdit::RefAlt { ref_: Some(ref_base.to_string()), alt: Some(alt.to_string()), uncertain: false }
        };
        let vc = coding_variant(&g, &am, pos, pos + 1, edit);
        let vr = mapper.c_to_r(&vc).unwrap();
        prop_assert_eq!(&mapper.r_to_c(&vr).unwrap(), &vc, "{} did not come back from {}", vc, vr);
        let r_text = vr.to_string();
        let edit_text = r_text.split(":r.").nth(1).unwrap();
        prop_assert!(!edit_text.contains(|ch: char| ch.is_ascii_uppercase() || ch == 't'), "{} is not in RNA letters", vr);

        let sr = SequenceVariant::Rna(vr.clone());
        let sc = SequenceVariant::Coding(vc.clone());
        prop_assert_eq!(mapper.canonical_allele(&sr).unwrap(), mapper.canonical_allele(&sc).unwrap());
        prop_assert_eq!(mapper.to_vrs(&sr).unwrap().id, mapper.to_vrs(&sc).unwrap().id);
        prop_assert_eq!(
            mapper.r_to_p(&vr, Some("NP_PROP.1")).unwrap().to_string(),
            mapper.c_to_p(&vc, Some("NP_PROP.1")).unwrap().to_string()
        );
        prop_assert_eq!(mapper.validate(&sr).unwrap(), mapper.validate(&sc).unwrap());
        let nr = mapper.normalize_variant(sr).unwrap();
        let nc = mapper.normalize_variant(sc).unwrap();
        prop_assert_eq!(nr.to_string(), mapper.tx_to_r(&nc).unwrap().to_string(), "normalising r. differs from normalising c.");
    }

    /// A deletion spanning a splice junction is an RNA finding: it has a
    /// protein prediction but no genomic projection, while its c. spelling
    /// projects to a deletion that includes the intron.
    #[test]
    fn rna_across_a_junction_has_a_protein_but_no_genomic_form(g in gene(), pick in 0usize..10_000) {
        if g.tx_exons.len() < 2 {
            return Ok(());
        }
        let am = TranscriptMapper::new(g.transcript_data()).unwrap();
        let hdp = g.provider();
        let mapper = VariantMapper::new(&hdp);
        let k = pick % (g.tx_exons.len() - 1);
        let last = g.tx_exons[k].1 - 1; // last base of exon k; the next is exon k + 1's first
        let vc = coding_variant(&g, &am, last, last + 2, NaEdit::Del { ref_: None, uncertain: false });
        let vr = mapper.c_to_r(&vc).unwrap();
        prop_assert!(
            matches!(mapper.r_to_g(&vr, None), Err(hgvs_weaver::error::HgvsError::UnsupportedOperation(_))),
            "{} projected to the genome", vr
        );
        prop_assert!(mapper.tx_to_g(&vc, None).is_ok(), "{} did not project", vc);
        if g.cds_start <= last && last + 2 <= g.cds_end + 1 {
            prop_assert_eq!(
                mapper.r_to_p(&vr, Some("NP_PROP.1")).unwrap().to_string(),
                mapper.c_to_p(&vc, Some("NP_PROP.1")).unwrap().to_string()
            );
        }
    }

    /// A coding variant and its genomic projection are the same variant, and
    /// the mitochondrial spelling of a genomic variant is the same allele.
    #[test]
    fn a_variant_agrees_with_its_projection(g in gene(), pick in 0usize..10_000, alt in prop::sample::select(BASES.to_vec()), del in prop::bool::ANY) {
        let am = TranscriptMapper::new(g.transcript_data()).unwrap();
        let hdp = g.provider();
        let mapper = VariantMapper::new(&hdp);
        let eq = VariantEquivalence::new(&mapper, &hdp);
        let pos = g.cds_start + pick % (g.cds_end + 1 - g.cds_start);
        let ref_base = g.transcript_seq.as_bytes()[pos] as char;
        let edit = if del {
            NaEdit::Del { ref_: None, uncertain: false }
        } else {
            prop_assume!(alt != ref_base);
            NaEdit::RefAlt { ref_: Some(ref_base.to_string()), alt: Some(alt.to_string()), uncertain: false }
        };
        let vc = SequenceVariant::Coding(coding_variant(&g, &am, pos, pos + 1, edit));
        let vg = SequenceVariant::Genomic(mapper.as_genomic(&vc).unwrap().unwrap());
        let lvl = eq.equivalent_level(&vc, &vg).unwrap();
        prop_assert!(lvl.is_equivalent(), "c. and its g. projection are not equivalent: {vc} vs {vg} ({lvl:?})");
        prop_assert_eq!(mapper.canonical_allele(&vc).unwrap(), mapper.canonical_allele(&vg).unwrap());

        let m_str = vg.to_string().replacen(":g.", ":m.", 1);
        let vm = parse_hgvs_variant(&m_str).unwrap();
        prop_assert_eq!(mapper.canonical_allele(&vm).unwrap(), mapper.canonical_allele(&vg).unwrap());
        prop_assert_eq!(eq.equivalent_level(&vm, &vg).unwrap(), EquivalenceLevel::Identity);
    }

    /// Normalising a transcript variant moves where it is written, never what
    /// it does to the genome: where the transcript and the genome agree, the
    /// canonical allele of its projection is the same before and after,
    /// identical bases either side of an exon/exon junction included.
    #[test]
    fn normalising_on_the_transcript_keeps_the_genomic_allele((g, e) in gene_and_edit()) {
        if !g.genome_agrees_with_transcript() {
            return Ok(());
        }
        let am = TranscriptMapper::new(g.transcript_data()).unwrap();
        let hdp = g.provider();
        let mapper = VariantMapper::new(&hdp);
        let vc = SequenceVariant::Coding(tx_edit_variant(&g, &am, &e));
        let Ok(before) = mapper.canonical_allele(&vc) else { return Ok(()) };
        let normalised = mapper.normalize_variant(vc.clone()).unwrap();
        let after = mapper
            .canonical_allele(&normalised)
            .map_err(|err| TestCaseError::fail(format!("{vc} normalised to {normalised}: {err}")))?;
        prop_assert_eq!(
            after,
            before,
            "{} normalised to {}",
            vc,
            normalised
        );
    }

    /// On a transcript, normalising twice is normalising once, beside a
    /// junction included.
    #[test]
    fn normalising_a_transcript_variant_twice_is_normalising_it_once((g, e) in gene_and_edit()) {
        let am = TranscriptMapper::new(g.transcript_data()).unwrap();
        let hdp = g.provider();
        let mapper = VariantMapper::new(&hdp);
        let vc = SequenceVariant::Coding(tx_edit_variant(&g, &am, &e));
        let Ok(once) = mapper.normalize_variant(vc.clone()) else { return Ok(()) };
        let twice = mapper
            .normalize_variant(once.clone())
            .map_err(|err| TestCaseError::fail(format!("{vc} normalised to {once}: {err}")))?;
        prop_assert_eq!(twice.to_string(), once.to_string(), "{} normalised to {}", vc, once);
    }

    /// The same, for a deletion, duplication or insertion written just before
    /// an exon/exon junction, where a slide on the transcript could carry it
    /// into the next exon, and for a delins there, up to or across the
    /// junction, that inserts or deletes inside the bases it keeps, which
    /// trimming could carry off the exon.
    #[test]
    fn normalising_beside_a_junction_keeps_the_genomic_allele(
        g in gene(),
        pick in 0usize..10_000,
        back in 0usize..3,
        len in 1usize..=3,
        kind in 0u8..5,
        alt in dna(1, 3),
        over in 0usize..3,
        cut in 0usize..8,
    ) {
        if g.tx_exons.len() < 2 {
            return Ok(());
        }
        let am = TranscriptMapper::new(g.transcript_data()).unwrap();
        let hdp = g.provider();
        let mapper = VariantMapper::new(&hdp);
        // the first index of an exon after the first
        let junction = g.tx_exons[1 + pick % (g.tx_exons.len() - 1)].0;
        let last = junction - 1 - back.min(junction - 1);
        let (start, end, edit) = match kind {
            0 => (last + 1 - len.min(last + 1), last + 1, NaEdit::Del { ref_: None, uncertain: false }),
            1 => (last + 1 - len.min(last + 1), last + 1, NaEdit::Dup { ref_: None, uncertain: false }),
            2 => (last, last + 2, NaEdit::Ins { alt: Some(alt), uncertain: false }),
            _ => {
                let (start, end) = (last + 1 - len.min(last + 1), (last + 1 + over).min(g.transcript_seq.len()));
                let kept = &g.transcript_seq[start..end];
                let at = cut % (kept.len() + 1);
                let inserted = if kind == 3 {
                    format!("{}{alt}{}", &kept[..at], &kept[at..])
                } else {
                    let at = at.min(kept.len() - 1);
                    format!("{}{}", &kept[..at], &kept[at + 1..])
                };
                if inserted.is_empty() {
                    return Ok(());
                }
                (start, end, NaEdit::RefAlt { ref_: Some(String::new()), alt: Some(inserted), uncertain: false })
            }
        };
        let vc = SequenceVariant::Coding(coding_variant(&g, &am, start, end, edit));
        let Ok(before) = mapper.canonical_allele(&vc) else { return Ok(()) };
        let normalised = mapper.normalize_variant(vc.clone()).unwrap();
        let after = mapper
            .canonical_allele(&normalised)
            .map_err(|err| TestCaseError::fail(format!("{vc} normalised to {normalised}: {err}")))?;
        // Compared by the genome each makes: a delins across a junction can
        // change nothing on the genome, and an unchanged range has no one
        // canonical extent.
        prop_assert_eq!(made(&g.genome, &after), made(&g.genome, &before), "{} normalised to {}", vc, normalised);
    }
}

/// The genome with the allele applied.
fn made(genome: &str, allele: &hgvs_weaver::allele::CanonicalAllele) -> String {
    format!(
        "{}{}{}",
        &genome[..allele.start],
        allele.alternate,
        &genome[allele.end..]
    )
}

// ---------------------------------------------------------------------------
// A repeat projected between strands
//
// A repeat's count is of the run read 3' from its first base on the sequence
// it is written on. The oracle reads that run itself and writes the count of
// copies over the bases the run occupies on the other sequence; weaver's
// canonical allele, applied to the genome, must make the same genome.
// ---------------------------------------------------------------------------

/// The exon holding transcript index `i`.
fn exon_holding(g: &Gene, i: usize) -> &hgvs_weaver::data::ExonData {
    g.exons
        .iter()
        .find(|e| e.transcript_start.0 as usize <= i && i < e.transcript_end.0 as usize)
        .unwrap()
}

/// The end of the run of `unit` in `seq` from `start`, half-open.
fn run_end(seq: &str, start: usize, unit: &str) -> usize {
    let mut end = start;
    while seq.get(end..end + unit.len()) == Some(unit) {
        end += unit.len();
    }
    end
}

/// `genome` with `[start, end)` replaced by `alt`.
fn spliced_in(genome: &str, start: usize, end: usize, alt: &str) -> String {
    format!("{}{alt}{}", &genome[..start], &genome[end..])
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(512))]

    /// A transcript repeat, within one exon, makes the genome the transcript's
    /// run with its count of copies would.
    #[test]
    fn a_coding_repeat_does_to_the_genome_what_it_does_to_the_transcript(
        g in gene(),
        pick in 0usize..10_000,
        unit_len in 1usize..=3,
        copies in 0i32..=5,
    ) {
        let tx = &g.transcript_seq;
        let start = pick % tx.len();
        prop_assume!(start + unit_len <= tx.len());
        let unit = &tx[start..start + unit_len];
        let end = run_end(tx, start, unit);
        let exon = exon_holding(&g, start);
        prop_assume!(end <= exon.transcript_end.0 as usize);
        let into = |i: usize| i - exon.transcript_start.0 as usize;
        let (lo, hi) = (exon.reference_start.0 as usize, exon.reference_end.0 as usize + 1);
        let expected = match g.strand {
            Strand::Plus => spliced_in(&g.genome, lo + into(start), lo + into(end), &unit.repeat(copies as usize)),
            Strand::Minus => spliced_in(&g.genome, hi - into(end), hi - into(start), &revcomp(&unit.repeat(copies as usize))),
        };

        let am = TranscriptMapper::new(g.transcript_data()).unwrap();
        let hdp = g.provider();
        let mapper = VariantMapper::new(&hdp);
        let edit = NaEdit::Repeat { ref_: Some(unit.to_string()), min: copies, max: copies, uncertain: false };
        let vc = SequenceVariant::Coding(coding_variant(&g, &am, start, start + 1, edit));
        let allele = mapper
            .canonical_allele(&vc)
            .map_err(|err| TestCaseError::fail(format!("{vc}: {err}")))?;
        let actual = spliced_in(&g.genome, allele.start, allele.end, &allele.alternate);
        prop_assert_eq!(actual, expected, "{} on the {:?} strand", vc, g.strand);
    }

    /// A transcript repeat within one intron makes the genome the run read
    /// in the transcript's sense, with its count of copies, would.
    #[test]
    fn an_intronic_repeat_does_to_the_genome_what_it_does_to_the_pre_mrna(
        g in gene(),
        pick in 0usize..10_000,
        unit_len in 1usize..=3,
        copies in 0i32..=5,
    ) {
        let points = g.points();
        let pre_mrna = g.bases_at(&points);
        let flank = g.flank();
        let intronic: Vec<usize> = (flank..points.len() - flank)
            .filter(|&i| points[i].offset != 0)
            .collect();
        prop_assume!(!intronic.is_empty());
        let start = intronic[pick % intronic.len()];
        prop_assume!(start + unit_len <= points.len() - flank);
        let unit = &pre_mrna[start..start + unit_len];
        let end = run_end(&pre_mrna, start, unit);
        prop_assume!(end <= points.len() - flank && points[start..end].iter().all(|p| p.offset != 0));
        let genomic: Vec<usize> = points[start..end].iter().map(|&p| g.genome_index(p)).collect();
        let (lo, hi) = (*genomic.iter().min().unwrap(), *genomic.iter().max().unwrap() + 1);
        let alt = unit.repeat(copies as usize);
        let alt = match g.strand {
            Strand::Plus => alt,
            Strand::Minus => revcomp(&alt),
        };
        let expected = spliced_in(&g.genome, lo, hi, &alt);

        let am = TranscriptMapper::new(g.transcript_data()).unwrap();
        let hdp = g.provider();
        let mapper = VariantMapper::new(&hdp);
        let edit = NaEdit::Repeat { ref_: Some(unit.to_string()), min: copies, max: copies, uncertain: false };
        let e = TxEdit { start: points[start], end: points[start], edit };
        let vc = SequenceVariant::Coding(tx_edit_variant(&g, &am, &e));
        let allele = mapper
            .canonical_allele(&vc)
            .map_err(|err| TestCaseError::fail(format!("{vc}: {err}")))?;
        let actual = spliced_in(&g.genome, allele.start, allele.end, &allele.alternate);
        prop_assert_eq!(actual, expected, "{} on the {:?} strand", vc, g.strand);
    }

    /// A genomic repeat, within one exon, projected to the transcript makes
    /// the genome the genome's run with its count of copies would.
    #[test]
    fn a_genomic_repeat_projected_to_the_transcript_is_the_same_change(
        g in gene(),
        pick in 0usize..10_000,
        unit_len in 1usize..=3,
        copies in 0i32..=5,
    ) {
        let exon = &g.exons[pick % g.exons.len()];
        let (lo, hi) = (exon.reference_start.0 as usize, exon.reference_end.0 as usize + 1);
        let start = lo + (pick / g.exons.len()) % (hi - lo);
        prop_assume!(start + unit_len <= hi);
        let unit = &g.genome[start..start + unit_len];
        let end = run_end(&g.genome, start, unit);
        prop_assume!(end <= hi);
        // The transcript counts its own run from the projected lower end; one
        // that crosses an exon junction is refused, leaving no c. to compare.
        let tx_start = exon.transcript_start.0 as usize;
        let tx_lo = match g.strand {
            Strand::Plus => tx_start + (start - lo),
            Strand::Minus => tx_start + (hi - end),
        };
        let tx_unit = &g.transcript_seq[tx_lo..tx_lo + unit_len];
        prop_assume!(run_end(&g.transcript_seq, tx_lo, tx_unit) <= exon.transcript_end.0 as usize);
        let expected = spliced_in(&g.genome, start, end, &unit.repeat(copies as usize));

        let hdp = g.provider();
        let mapper = VariantMapper::new(&hdp);
        let vg = format!("{}:g.{}{unit}[{copies}]", g.reference_ac, start + 1);
        let SequenceVariant::Genomic(vg) = parse_hgvs_variant(&vg).unwrap() else { unreachable!() };
        let vc = mapper
            .g_to_c(&vg, &g.ac)
            .map_err(|err| TestCaseError::fail(format!("{vg}: {err}")))?;
        let allele = mapper
            .canonical_allele(&SequenceVariant::Coding(vc.clone()))
            .map_err(|err| TestCaseError::fail(format!("{vg} projected to {vc}: {err}")))?;
        let actual = spliced_in(&g.genome, allele.start, allele.end, &allele.alternate);
        prop_assert_eq!(actual, expected, "{} projected to {} on the {:?} strand", vg, vc, g.strand);
    }
}

// ---------------------------------------------------------------------------
// What a projection writes is HGVS the parser reads back
//
// A caller stores what weaver writes and parses it later, so every projection's
// result must parse, to the value it was written from. One property per
// projection; a projection that refuses its input is not this property's
// concern.
// ---------------------------------------------------------------------------

/// The c. variant a [`TxEdit`] describes on `g`.
fn tx_edit_variant(g: &Gene, am: &TranscriptMapper, e: &TxEdit) -> CVariant {
    let position = |p: TxPoint| {
        let mut pos = CVariant::position_from_index(am, p.index).unwrap();
        pos.offset = (p.offset != 0).then_some(IntronicOffset(p.offset));
        pos
    };
    let end = (e.end != e.start).then(|| position(e.end));
    variant_from(&g.ac, interval(position(e.start), end), &e.edit)
}

/// The variant's HGVS parses, to the variant it was written from.
fn reparses_to_itself(v: SequenceVariant) -> Result<(), TestCaseError> {
    let text = v.to_string();
    let back = parse_hgvs_variant(&text)
        .map_err(|e| TestCaseError::fail(format!("{text} does not parse: {e}")))?;
    prop_assert_eq!(back, v, "{} parses to a different variant", text);
    Ok(())
}

/// The p. variant without the residues its edit replaces: p. text names only
/// the interval's ends, so a parsed consequence never carries them.
fn as_written(mut vp: PVariant) -> PVariant {
    vp.posedit.edit = match vp.posedit.edit {
        AaEdit::Subst { alt, uncertain, .. } => AaEdit::Subst {
            ref_: String::new(),
            alt,
            uncertain,
        },
        AaEdit::Del { uncertain, .. } => AaEdit::Del {
            ref_: String::new(),
            uncertain,
        },
        AaEdit::DelIns { alt, uncertain, .. } => AaEdit::DelIns {
            ref_: String::new(),
            alt,
            uncertain,
        },
        AaEdit::Dup { uncertain, .. } => AaEdit::Dup {
            ref_: None,
            uncertain,
        },
        other => other,
    };
    vp
}

/// The consequence's HGVS parses, to the consequence it was written from, as
/// far as p. text can say ([`as_written`]).
fn protein_reparses_to_itself(vp: PVariant) -> Result<(), TestCaseError> {
    let text = vp.to_string();
    let back = match parse_hgvs_variant(&text) {
        Ok(SequenceVariant::Protein(back)) => back,
        other => {
            return Err(TestCaseError::fail(format!(
                "{text} does not parse to a p. variant: {other:?}"
            )))
        }
    };
    prop_assert_eq!(
        as_written(back),
        as_written(vp),
        "{} parses to a different variant",
        text
    );
    Ok(())
}

proptest! {
    /// The generator's own check: its points name every genome base once, in
    /// transcript order, and the exonic ones hold the transcript's bases.
    #[test]
    fn gene_points_name_every_genome_base_once_in_transcript_order(g in gene()) {
        let named: Vec<usize> = g.points().into_iter().map(|p| g.genome_index(p)).collect();
        let mut expected: Vec<usize> = (0..g.genome.len()).collect();
        if g.strand == hgvs_weaver::data::Strand::Minus {
            expected.reverse();
        }
        prop_assert_eq!(named, expected);
        let exonic: Vec<TxPoint> = g.points().into_iter().filter(|p| p.offset == 0).collect();
        prop_assert_eq!(g.bases_at(&exonic), g.transcript_seq.clone());
    }

    #[test]
    fn a_genomic_projection_reparses_to_itself((g, e) in gene_and_edit()) {
        let am = TranscriptMapper::new(g.transcript_data()).unwrap();
        let hdp = g.provider();
        let mapper = VariantMapper::new(&hdp);
        let vc = tx_edit_variant(&g, &am, &e);
        let Ok(vg) = mapper.c_to_g(&vc, None) else { return Ok(()) };
        reparses_to_itself(SequenceVariant::Genomic(vg))?;
    }

    #[test]
    fn a_coding_projection_reparses_to_itself((g, e) in gene_and_edit()) {
        let am = TranscriptMapper::new(g.transcript_data()).unwrap();
        let hdp = g.provider();
        let mapper = VariantMapper::new(&hdp);
        let vc = tx_edit_variant(&g, &am, &e);
        let Ok(vg) = mapper.c_to_g(&vc, None) else { return Ok(()) };
        let Ok(back) = mapper.g_to_c(&vg, &g.ac) else { return Ok(()) };
        reparses_to_itself(SequenceVariant::Coding(back))?;
    }

    #[test]
    fn an_rna_projection_reparses_to_itself((g, e) in gene_and_edit()) {
        let am = TranscriptMapper::new(g.transcript_data()).unwrap();
        let hdp = g.provider();
        let mapper = VariantMapper::new(&hdp);
        let vc = tx_edit_variant(&g, &am, &e);
        let Ok(vr) = mapper.c_to_r(&vc) else { return Ok(()) };
        reparses_to_itself(SequenceVariant::Rna(vr))?;
    }

    #[test]
    fn a_coding_spelling_of_rna_reparses_to_itself((g, e) in gene_and_edit()) {
        let am = TranscriptMapper::new(g.transcript_data()).unwrap();
        let hdp = g.provider();
        let mapper = VariantMapper::new(&hdp);
        let vc = tx_edit_variant(&g, &am, &e);
        let Ok(vr) = mapper.c_to_r(&vc) else { return Ok(()) };
        let Ok(back) = mapper.r_to_c(&vr) else { return Ok(()) };
        reparses_to_itself(SequenceVariant::Coding(back))?;
    }

    #[test]
    fn a_genomic_projection_of_rna_reparses_to_itself((g, e) in gene_and_edit()) {
        let am = TranscriptMapper::new(g.transcript_data()).unwrap();
        let hdp = g.provider();
        let mapper = VariantMapper::new(&hdp);
        let vc = tx_edit_variant(&g, &am, &e);
        let Ok(vr) = mapper.c_to_r(&vc) else { return Ok(()) };
        let Ok(vg) = mapper.r_to_g(&vr, None) else { return Ok(()) };
        reparses_to_itself(SequenceVariant::Genomic(vg))?;
    }

    #[test]
    fn a_normalised_rna_variant_reparses_to_itself((g, e) in gene_and_edit()) {
        let am = TranscriptMapper::new(g.transcript_data()).unwrap();
        let hdp = g.provider();
        let mapper = VariantMapper::new(&hdp);
        let vc = tx_edit_variant(&g, &am, &e);
        let Ok(vr) = mapper.c_to_r(&vc) else { return Ok(()) };
        let Ok(normalised) = mapper.normalize_variant(SequenceVariant::Rna(vr)) else { return Ok(()) };
        reparses_to_itself(normalised)?;
    }

    #[test]
    fn a_protein_consequence_of_rna_reparses_to_itself((g, e) in gene_and_edit()) {
        let am = TranscriptMapper::new(g.transcript_data()).unwrap();
        let hdp = g.provider();
        let mapper = VariantMapper::new(&hdp);
        let vc = tx_edit_variant(&g, &am, &e);
        let Ok(vr) = mapper.c_to_r(&vc) else { return Ok(()) };
        let Ok(vp) = mapper.r_to_p(&vr, Some("NP_PROP.1")) else { return Ok(()) };
        protein_reparses_to_itself(vp)?;
    }

    #[test]
    fn a_protein_consequence_reparses_to_itself((g, e) in gene_and_edit()) {
        let am = TranscriptMapper::new(g.transcript_data()).unwrap();
        let hdp = g.provider();
        let mapper = VariantMapper::new(&hdp);
        let vc = tx_edit_variant(&g, &am, &e);
        let Ok(vp) = mapper.c_to_p(&vc, Some("NP_PROP.1")) else { return Ok(()) };
        protein_reparses_to_itself(vp)?;
    }

    #[test]
    fn a_normalised_coding_variant_reparses_to_itself((g, e) in gene_and_edit()) {
        let am = TranscriptMapper::new(g.transcript_data()).unwrap();
        let hdp = g.provider();
        let mapper = VariantMapper::new(&hdp);
        let vc = tx_edit_variant(&g, &am, &e);
        let Ok(normalised) = mapper.normalize_variant(SequenceVariant::Coding(vc)) else { return Ok(()) };
        reparses_to_itself(normalised)?;
    }

    #[test]
    fn a_normalised_genomic_variant_reparses_to_itself((g, e) in gene_and_edit()) {
        let am = TranscriptMapper::new(g.transcript_data()).unwrap();
        let hdp = g.provider();
        let mapper = VariantMapper::new(&hdp);
        let vc = tx_edit_variant(&g, &am, &e);
        let Ok(vg) = mapper.c_to_g(&vc, None) else { return Ok(()) };
        let Ok(normalised) = mapper.normalize_variant(SequenceVariant::Genomic(vg)) else { return Ok(()) };
        reparses_to_itself(normalised)?;
    }

    #[test]
    fn a_variant_read_from_spdi_reparses_to_itself((g, e) in gene_and_edit()) {
        let am = TranscriptMapper::new(g.transcript_data()).unwrap();
        let hdp = g.provider();
        let mapper = VariantMapper::new(&hdp);
        let vc = SequenceVariant::Coding(tx_edit_variant(&g, &am, &e));
        let Ok(spdi) = mapper.to_spdi(&vc, false) else { return Ok(()) };
        let Ok(back) = mapper.from_spdi(&spdi) else { return Ok(()) };
        reparses_to_itself(back)?;
    }

    #[test]
    fn a_variant_read_from_vrs_reparses_to_itself((g, e) in gene_and_edit()) {
        let am = TranscriptMapper::new(g.transcript_data()).unwrap();
        let hdp = g.provider();
        let mapper = VariantMapper::new(&hdp);
        let vc = SequenceVariant::Coding(tx_edit_variant(&g, &am, &e));
        let Ok(vrs) = mapper.to_vrs(&vc) else { return Ok(()) };
        let Ok(back) = mapper.from_vrs(&vrs.to_json(), Some(&g.reference_ac)) else { return Ok(()) };
        reparses_to_itself(back)?;
    }
}

proptest! {
    /// An allele in cis projected to the genome and back is the allele it
    /// was, member for member and in the order written, where no member meets
    /// another and every one projects.
    #[test]
    fn a_cis_allele_projected_to_the_genome_and_back_is_itself((g, edits) in gene_and_cis_edits()) {
        let am = TranscriptMapper::new(g.transcript_data()).unwrap();
        let hdp = g.provider();
        let mapper = VariantMapper::new(&hdp);
        let members = edits
            .iter()
            .map(|e| SequenceVariant::Coding(tx_edit_variant(&g, &am, e)))
            .collect();
        let cis = hgvs_weaver::CisPhasedVariant::new(g.ac.clone(), None, members).unwrap();
        let Ok(vg) = mapper.c_to_g_cis(&cis, None) else { return Ok(()) };
        let back = mapper
            .g_to_c_cis(&vg, &g.ac)
            .map_err(|err| TestCaseError::fail(format!("{cis} projected to {vg}: {err}")))?;
        prop_assert_eq!(back.to_string(), cis.to_string(), "through {}", vg);
        prop_assert_eq!(back, cis);
    }

    /// An allele written in transcript order is in genomic order on the
    /// genome, on either strand.
    #[test]
    fn a_cis_allele_in_transcript_order_projects_in_genomic_order((g, edits) in gene_and_cis_edits()) {
        let am = TranscriptMapper::new(g.transcript_data()).unwrap();
        let hdp = g.provider();
        let mapper = VariantMapper::new(&hdp);
        let members = edits
            .iter()
            .map(|e| SequenceVariant::Coding(tx_edit_variant(&g, &am, e)))
            .collect();
        let cis = hgvs_weaver::CisPhasedVariant::new(g.ac.clone(), None, members).unwrap();
        let Ok(vg) = mapper.c_to_g_cis(&cis, None) else { return Ok(()) };
        let starts: Vec<i32> = vg
            .members
            .iter()
            .map(|m| match m {
                SequenceVariant::Genomic(v) => v.posedit.pos.as_ref().unwrap().start.index().unwrap().0,
                other => panic!("{other} is not on the genome"),
            })
            .collect();
        prop_assert!(starts.windows(2).all(|w| w[0] < w[1]), "{} projected to {}", cis, vg);
    }
}

proptest! {
    // Half the genes are on the plus strand, so the filter rejects as many
    // as it keeps.
    #![proptest_config(ProptestConfig { max_local_rejects: u32::MAX, ..ProptestConfig::default() })]

    /// On a minus-strand transcript, an insertion before what one or two
    /// duplications copy makes on the genome what it makes on the transcript,
    /// in any order written: the exon, as the projected allele edits the
    /// genome, read on the minus strand, is the exon as the allele edits the
    /// transcript.
    #[test]
    fn a_duplication_beside_an_insertion_makes_the_same_exon_on_the_genome(
        g in gene().prop_filter("a minus-strand gene with an exon of eight", |g| {
            g.strand == hgvs_weaver::data::Strand::Minus && g.tx_exons.iter().any(|(a, b)| b - a >= 8)
        }),
        at in 0usize..10_000,
        lens in (1usize..=3, 1usize..=3),
        chained in prop::bool::ANY,
        alt in dna(1, 3),
        order in 0usize..6,
    ) {
        let (ts, te) = *g.tx_exons.iter().max_by_key(|(a, b)| b - a).unwrap();
        let reach = lens.0 + if chained { lens.1 } else { 0 };
        // the insertion goes between s - 1 and s, inside the exon
        let s = ts + 1 + at % (te - ts - reach - 1);
        let e = s + lens.0 - 1;
        let am = TranscriptMapper::new(g.transcript_data()).unwrap();
        let mut edits = vec![
            (s, s, Some(alt.clone())),
            (s, e, None),
        ];
        if chained {
            edits.push((e + 1, e + lens.1, None));
        }
        let rotate = order % edits.len();
        edits.rotate_left(rotate);
        if order >= 3 {
            edits.reverse();
        }
        let members = edits
            .iter()
            .map(|(a, b, ins)| {
                let v = match ins {
                    Some(bases) => coding_variant(&g, &am, a - 1, a + 1, NaEdit::Ins { alt: Some(bases.clone()), uncertain: false }),
                    None => coding_variant(&g, &am, *a, b + 1, NaEdit::Dup { ref_: None, uncertain: false }),
                };
                SequenceVariant::Coding(v)
            })
            .collect();
        let cis = hgvs_weaver::CisPhasedVariant::new(g.ac.clone(), None, members).unwrap();
        let hdp = g.provider();
        let mapper = VariantMapper::new(&hdp);
        let vg = mapper.c_to_g_cis(&cis, None).unwrap();

        // On the transcript: each insertion at its 0-based offset, last first.
        let mut tx_inserts: Vec<(usize, String)> = edits
            .iter()
            .map(|(a, b, ins)| match ins {
                Some(bases) => (*a, bases.clone()),
                None => (b + 1, g.transcript_seq[*a..=*b].to_string()),
            })
            .collect();
        let mut on_genome: Vec<(usize, String)> = vg
            .members
            .iter()
            .map(|m| {
                let SequenceVariant::Genomic(v) = m else { panic!("{m} is not on the genome") };
                let pos = v.posedit.pos.as_ref().unwrap();
                let first = pos.start.index().unwrap().0 as usize;
                let last = pos.end.as_ref().map_or(first, |p| p.index().unwrap().0 as usize);
                match &v.posedit.edit {
                    NaEdit::Ins { alt: Some(bases), .. } => (first + 1, bases.clone()),
                    NaEdit::Dup { .. } => (last + 1, g.genome[first..=last].to_string()),
                    other => panic!("{other:?} in {vg}"),
                }
            })
            .collect();
        let apply = |seq: &str, inserts: &mut Vec<(usize, String)>| {
            inserts.sort_by_key(|x| std::cmp::Reverse(x.0));
            let mut out = seq.to_string();
            for (offset, bases) in inserts.iter() {
                out.insert_str(*offset, bases);
            }
            out
        };
        let delta: usize = tx_inserts.iter().map(|i| i.1.len()).sum();
        let edited_tx = apply(&g.transcript_seq, &mut tx_inserts);
        let edited_genome = apply(&g.genome, &mut on_genome);
        let exon = &g.exons.iter().find(|x| x.transcript_start.0 as usize == ts).unwrap();
        let (rs, re) = (exon.reference_start.0 as usize, exon.reference_end.0 as usize);
        prop_assert_eq!(
            revcomp(&edited_genome[rs..=re + delta]),
            edited_tx[ts..te + delta].to_string(),
            "{} projected to {}", cis, vg
        );
    }
}

/// The properties above pass vacuously on a case a projection refuses, so on
/// a fixed sample each chain a property follows must answer most generated
/// edits, and some genes must differ from their genome.
#[test]
fn the_generated_edits_reach_every_projection() {
    use proptest::strategy::ValueTree;
    let mut runner = proptest::test_runner::TestRunner::deterministic();
    let strategy = gene_and_edit();
    let cases = 1000;
    // Each chain, and the tenths of the sample it must answer. Most edits
    // fall outside the CDS, where the protein is a statement.
    let mut chains: [(&str, usize, usize); 13] = [
        ("c_to_g", 7, 0),
        ("g_to_c", 7, 0),
        ("c_to_r", 7, 0),
        ("r_to_c", 7, 0),
        ("r_to_g", 7, 0),
        ("r_to_p", 2, 0),
        ("c_to_p", 2, 0),
        ("normalize c.", 7, 0),
        ("normalize g.", 7, 0),
        ("normalize r.", 7, 0),
        ("canonical allele of c.", 7, 0),
        ("from_spdi", 7, 0),
        ("from_vrs", 7, 0),
    ];
    let mut discordant = 0;
    // a statement (`p.?`, `p.(Met1?)`, `p.(=)`) describes no change
    let describes_a_change =
        |vp: &hgvs_weaver::structs::PVariant| !matches!(vp.posedit.edit, AaEdit::Statement(_));
    for _ in 0..cases {
        let (g, e) = strategy.new_tree(&mut runner).unwrap().current();
        let am = TranscriptMapper::new(g.transcript_data()).unwrap();
        let hdp = g.provider();
        let mapper = VariantMapper::new(&hdp);
        let vc = tx_edit_variant(&g, &am, &e);
        let vg = mapper.c_to_g(&vc, None).ok();
        let vr = mapper.c_to_r(&vc).ok();
        let coding = SequenceVariant::Coding(vc.clone());
        let answered = [
            vg.is_some(),
            vg.as_ref().is_some_and(|vg| mapper.g_to_c(vg, &g.ac).is_ok()),
            vr.is_some(),
            vr.as_ref().is_some_and(|vr| mapper.r_to_c(vr).is_ok()),
            vr.as_ref().is_some_and(|vr| mapper.r_to_g(vr, None).is_ok()),
            vr.as_ref().is_some_and(|vr| {
                matches!(mapper.r_to_p(vr, Some("NP_PROP.1")), Ok(vp) if describes_a_change(&vp))
            }),
            matches!(mapper.c_to_p(&vc, Some("NP_PROP.1")), Ok(vp) if describes_a_change(&vp)),
            mapper.normalize_variant(coding.clone()).is_ok(),
            vg.clone().is_some_and(|vg| mapper.normalize_variant(SequenceVariant::Genomic(vg)).is_ok()),
            vr.clone().is_some_and(|vr| mapper.normalize_variant(SequenceVariant::Rna(vr)).is_ok()),
            mapper.canonical_allele(&coding).is_ok(),
            mapper
                .to_spdi(&coding, false)
                .is_ok_and(|spdi| mapper.from_spdi(&spdi).is_ok()),
            mapper.to_vrs(&coding).is_ok_and(|vrs| {
                mapper
                    .from_vrs(&vrs.to_json(), Some(&g.reference_ac))
                    .is_ok()
            }),
        ];
        for (chain, ok) in chains.iter_mut().zip(answered) {
            chain.2 += usize::from(ok);
        }
        discordant += usize::from(!g.genome_agrees_with_transcript());
    }
    for (name, floor, ok) in chains {
        assert!(
            ok * 10 >= cases * floor,
            "{name} answered {ok} of {cases} generated edits"
        );
    }
    assert!(
        discordant * 10 >= cases,
        "{discordant} of {cases} generated genes differ from their genome"
    );
}

/// The cis round trip passes vacuously on an allele that does not project, so
/// on a fixed sample most generated alleles must project, many with several
/// members and on the minus strand.
#[test]
fn the_generated_cis_alleles_project() {
    use proptest::strategy::ValueTree;
    let mut runner = proptest::test_runner::TestRunner::deterministic();
    let strategy = gene_and_cis_edits();
    let cases = 1000;
    let (mut projected, mut several, mut minus) = (0, 0, 0);
    for _ in 0..cases {
        let (g, edits) = strategy.new_tree(&mut runner).unwrap().current();
        let am = TranscriptMapper::new(g.transcript_data()).unwrap();
        let hdp = g.provider();
        let mapper = VariantMapper::new(&hdp);
        let members = edits
            .iter()
            .map(|e| SequenceVariant::Coding(tx_edit_variant(&g, &am, e)))
            .collect();
        let cis = hgvs_weaver::CisPhasedVariant::new(g.ac.clone(), None, members).unwrap();
        if mapper.c_to_g_cis(&cis, None).is_ok() {
            projected += 1;
            several += usize::from(edits.len() > 1);
            minus += usize::from(g.strand == hgvs_weaver::data::Strand::Minus);
        }
    }
    assert!(
        projected * 10 >= cases * 8,
        "{projected} of {cases} alleles projected"
    );
    assert!(
        several * 10 >= cases * 5,
        "{several} of {cases} projected with several members"
    );
    assert!(
        minus * 10 >= cases * 3,
        "{minus} of {cases} projected on the minus strand"
    );
}

// Reading an allele in cis into one sequence.

/// A member of an allele on a sequence of `len` bases: any of the edits the
/// single-variant properties use, or a substitution, an inversion or a
/// positional identity, so that members overlap in every way the rule reads.
fn cis_member(len: usize) -> impl Strategy<Value = PlacedHgvs> {
    let base = prop::sample::select(BASES.to_vec());
    let sub = (0..len, base).prop_map(|(s, b)| PlacedHgvs {
        start: s,
        end: s + 1,
        edit: NaEdit::RefAlt {
            ref_: None,
            alt: Some(b.to_string()),
            uncertain: false,
        },
    });
    let span = (0..len - 1).prop_flat_map(move |s| (Just(s), s + 2..=len.min(s + 5)));
    let inv = span.clone().prop_map(|(s, e)| PlacedHgvs {
        start: s,
        end: e,
        edit: NaEdit::Inv {
            ref_: None,
            uncertain: false,
        },
    });
    let identity = span.prop_map(|(s, e)| PlacedHgvs {
        start: s,
        end: e,
        edit: NaEdit::RefAlt {
            ref_: None,
            alt: None,
            uncertain: false,
        },
    });
    prop_oneof![4 => edit_on_kinds(len, false), 2 => sub, 1 => inv, 1 => identity]
}

/// A member as [`cis_member`] builds one, or one no sequence can be built
/// from: a copy number, or an insertion of bases it does not state.
fn cis_member_or_unapplied(len: usize) -> impl Strategy<Value = PlacedHgvs> {
    let copy = (0..len - 1)
        .prop_flat_map(move |s| (Just(s), s + 1..=len))
        .prop_map(|(s, e)| PlacedHgvs {
            start: s,
            end: e,
            edit: NaEdit::NACopy {
                copy: 3,
                uncertain: false,
            },
        });
    let unstated = (1..len - 1).prop_map(|b| PlacedHgvs {
        start: b - 1,
        end: b + 1,
        edit: NaEdit::InsLength {
            min: 2,
            max: 2,
            uncertain: false,
        },
    });
    prop_oneof![6 => cis_member(len), 1 => copy, 1 => unstated]
}

/// A short sequence, an allele on it, and either the same members in another
/// order or another allele, so that both equal and different sequences occur.
fn seq_and_two_alleles() -> impl Strategy<Value = (String, Vec<PlacedHgvs>, Vec<PlacedHgvs>)> {
    dna(8, 24).prop_flat_map(|seq| {
        let len = seq.len();
        let allele = || prop::collection::vec(cis_member(len), 1..=3);
        (Just(seq), allele(), allele(), any::<bool>()).prop_flat_map(|(seq, a, b, reorder)| {
            let b = if reorder {
                Just(a.clone()).prop_shuffle().boxed()
            } else {
                Just(b).boxed()
            };
            (Just(seq), Just(a), b)
        })
    })
}

fn genomic_cis(members: &[PlacedHgvs]) -> hgvs_weaver::CisPhasedVariant {
    hgvs_weaver::CisPhasedVariant::new(
        "X".into(),
        None,
        members
            .iter()
            .map(|p| SequenceVariant::Genomic(genomic_variant("X", p)))
            .collect(),
    )
    .unwrap()
}

/// What reading an allele answers, in a form two answers can be compared by:
/// the error's kind, and an ambiguity's readings.
#[derive(Debug, PartialEq)]
enum Outcome {
    Read,
    Contradictory,
    Ambiguous(Vec<String>),
    Unsupported,
    Mismatch,
}

fn outcome(result: &Result<hgvs_weaver::reading::CisReading, hgvs_weaver::HgvsError>) -> Outcome {
    use hgvs_weaver::HgvsError as E;
    match result {
        Ok(_) => Outcome::Read,
        Err(E::ValidationError(m)) if m.contains("has no reading") => Outcome::Contradictory,
        Err(E::Ambiguous { readings, .. }) => {
            Outcome::Ambiguous(readings.iter().map(ToString::to_string).collect())
        }
        Err(E::UnsupportedOperation(_)) => Outcome::Unsupported,
        Err(E::ReferenceMismatch { .. }) => Outcome::Mismatch,
        Err(other) => panic!("reading failed: {other}"),
    }
}

/// The delins over the one stretch where `edited` differs from `seq`, or
/// the deletion where it only loses bases: the same sequence as an allele
/// makes, in one member. `None` when `edited` is `seq`.
fn as_one_member(seq: &str, edited: &str) -> Option<PlacedHgvs> {
    if seq == edited {
        return None;
    }
    let prefix = seq
        .bytes()
        .zip(edited.bytes())
        .take_while(|(a, b)| a == b)
        .count();
    let room = seq.len().min(edited.len()) - prefix;
    let suffix = seq
        .bytes()
        .rev()
        .zip(edited.bytes().rev())
        .take(room)
        .take_while(|(a, b)| a == b)
        .count();
    let (mut start, mut end) = (prefix, seq.len() - suffix);
    // A delins replaces at least one base.
    if start == end {
        if end < seq.len() {
            end += 1;
        } else {
            start -= 1;
        }
    }
    let alt = &edited[start..edited.len() - (seq.len() - end)];
    Some(PlacedHgvs {
        start,
        end,
        edit: if alt.is_empty() {
            NaEdit::Del {
                ref_: None,
                uncertain: false,
            }
        } else {
            NaEdit::RefAlt {
                ref_: Some(String::new()),
                alt: Some(alt.to_string()),
                uncertain: false,
            }
        },
    })
}

/// The bases a reading of an allele on `X` makes of all of `seq`.
fn made_by(mapper: &VariantMapper, seq: &str, cis: &hgvs_weaver::CisPhasedVariant) -> String {
    let reading = mapper
        .read_cis(cis)
        .unwrap_or_else(|e| panic!("{cis}: {e}"));
    let reference = mapper.refs.reference("X", G);
    reading.sequence(&reference, 0, seq.len()).unwrap()
}

/// The property tests over generated alleles are not vacuous: on a fixed
/// sample, a good share of the alleles have one reading.
#[test]
fn the_generated_alleles_are_mostly_read() {
    use proptest::strategy::ValueTree;
    let mut runner = proptest::test_runner::TestRunner::deterministic();
    let strategy = seq_and_two_alleles();
    let cases = 1000;
    let mut read = 0;
    for _ in 0..cases {
        let (seq, a, _) = strategy.new_tree(&mut runner).unwrap().current();
        let hdp = Provider::new().sequence("X", &seq);
        read += usize::from(VariantMapper::new(&hdp).read_cis(&genomic_cis(&a)).is_ok());
    }
    assert!(read * 10 >= cases * 4, "{read} of {cases} alleles read");
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(512))]

    /// Comparing two readings by walking their segments, without building
    /// either sequence, agrees with building both and comparing the strings,
    /// whether the two are split alike or one is the allele's whole change
    /// in one member.
    #[test]
    fn walking_two_readings_agrees_with_building_them((seq, a, b) in seq_and_two_alleles()) {
        let hdp = Provider::new().sequence("X", &seq);
        let mapper = VariantMapper::new(&hdp);
        let Ok(ra) = mapper.read_cis(&genomic_cis(&a)) else {
            return Ok(());
        };
        let reference = mapper.refs.reference("X", G);
        let built_a = ra.sequence(&reference, 0, seq.len()).unwrap();
        prop_assert_eq!(ra.makes_reference(&reference).unwrap(), built_a == seq);
        if let Some(whole) = as_one_member(&seq, &built_a) {
            let one = genomic_cis(std::slice::from_ref(&whole));
            let r_one = mapper.read_cis(&one).unwrap_or_else(|e| panic!("{one}: {e}"));
            prop_assert!(ra.same_sequence(&r_one, &reference).unwrap(), "{} and {}", genomic_cis(&a), one);
        }
        if let Ok(rb) = mapper.read_cis(&genomic_cis(&b)) {
            let built_b = rb.sequence(&reference, 0, seq.len()).unwrap();
            prop_assert_eq!(
                ra.same_sequence(&rb, &reference).unwrap(),
                built_a == built_b,
                "{} makes {} and {} makes {}", genomic_cis(&a), built_a, genomic_cis(&b), built_b
            );
        }
    }

    /// Members that share no base and do not touch make what applying each
    /// on its own makes: the allele's one reading, or, where members that
    /// cancel can be paired more than one way, every reading listed.
    #[test]
    fn members_apart_make_what_each_makes_alone(
        (seq, members) in (dna(12, 40), 1usize..=3).prop_flat_map(|(seq, k)| {
            // each member in a stretch of its own, one base short of the next
            let chunk = seq.len() / k;
            let members: Vec<_> = (0..k)
                .map(|i| {
                    cis_member(chunk - 1).prop_map(move |mut p| {
                        p.start += i * chunk;
                        p.end += i * chunk;
                        p
                    })
                })
                .collect();
            (Just(seq), members).prop_flat_map(|(seq, m)| (Just(seq), Just(m.clone()).prop_shuffle()))
        })
    ) {
        let mut sorted = members.clone();
        sorted.sort_by_key(|p| (p.start, p.end));
        let hdp = Provider::new().sequence("X", &seq);
        let mapper = VariantMapper::new(&hdp);
        let expected = sorted.iter().rev().fold(seq.clone(), |s, p| apply(&s, p));
        let cis = genomic_cis(&members);
        match mapper.read_cis(&cis) {
            Ok(_) => prop_assert_eq!(made_by(&mapper, &seq, &cis), expected),
            Err(hgvs_weaver::HgvsError::Ambiguous { readings, .. }) => {
                for reading in readings {
                    let as_allele = match reading {
                        SequenceVariant::CisPhased(c) => c,
                        one => hgvs_weaver::CisPhasedVariant::new("X".into(), None, vec![one]).unwrap(),
                    };
                    prop_assert_eq!(made_by(&mapper, &seq, &as_allele), expected.clone(), "{}", cis);
                }
            }
            Err(e) => prop_assert!(false, "{cis}: {e}"),
        }
    }

    /// The order members are written in changes nothing: the same sequence and
    /// statements, or the same error with the same readings.
    #[test]
    fn member_order_does_not_change_the_reading(
        (seq, members, shuffled) in dna(8, 24).prop_flat_map(|seq| {
            let len = seq.len();
            prop::collection::vec(cis_member_or_unapplied(len), 1..=4).prop_flat_map(move |m| {
                (Just(seq.clone()), Just(m.clone()), Just(m).prop_shuffle())
            })
        })
    ) {
        let hdp = Provider::new().sequence("X", &seq);
        let mapper = VariantMapper::new(&hdp);
        let (a, b) = (genomic_cis(&members), genomic_cis(&shuffled));
        let (ra, rb) = (mapper.read_cis(&a), mapper.read_cis(&b));
        prop_assert_eq!(outcome(&ra), outcome(&rb), "{} and {}", a, b);
        if let (Ok(ra), Ok(rb)) = (ra, rb) {
            let reference = mapper.refs.reference("X", G);
            prop_assert!(ra.same_sequence(&rb, &reference).unwrap(), "{} and {}", a, b);
            prop_assert_eq!(ra.stated(), rb.stated(), "{} and {}", a, b);
        }
    }
}

// Normalising an allele in cis: the members that meet written as one change.

/// A gene and a `c.` allele of up to four members packed into twelve bases
/// of its transcript, so that members meet, share codons and sit either side
/// of exon junctions.
fn gene_and_dense_cis() -> impl Strategy<Value = (Gene, Vec<PlacedHgvs>)> {
    const WINDOW: usize = 12;
    gene().prop_flat_map(|g| {
        let len = g.transcript_seq.len();
        let members = prop::collection::vec(cis_member(WINDOW), 1..=4);
        (Just(g), 0..len - WINDOW, members).prop_map(|(g, offset, members)| {
            let members = members
                .into_iter()
                .map(|mut p| {
                    p.start += offset;
                    p.end += offset;
                    p
                })
                .collect();
            (g, members)
        })
    })
}

fn coding_cis(g: &Gene, members: &[PlacedHgvs]) -> hgvs_weaver::CisPhasedVariant {
    let am = TranscriptMapper::new(g.transcript_data()).unwrap();
    let position = |i: usize| CVariant::position_from_index(&am, i as i32).unwrap();
    let members = members
        .iter()
        .map(|p| {
            let end = (p.end > p.start + 1).then(|| position(p.end - 1));
            SequenceVariant::Coding(variant_from(
                &g.ac,
                interval(position(p.start), end),
                &p.edit,
            ))
        })
        .collect();
    hgvs_weaver::CisPhasedVariant::new(g.ac.clone(), None, members).unwrap()
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(512))]

    /// Normalising a transcript allele keeps the sequence it makes on the
    /// genome, since no change moves or is merged across an exon junction,
    /// and the bases it states on its transcript; and normalising it again
    /// changes nothing. (The bases an allele states are compared on the
    /// transcript: a member written across a junction names the intron on the
    /// genome, so the two readings of one allele can state different bases.)
    #[test]
    fn normalising_a_transcript_cis_allele_keeps_its_genomic_sequence((g, members) in gene_and_dense_cis()) {
        let hdp = g.provider();
        let mapper = VariantMapper::new(&hdp);
        let cis = coding_cis(&g, &members);
        let (Ok(on_transcript), Ok(before)) = (mapper.read_cis(&cis), mapper.read_cis_on_genome(&cis)) else {
            return Ok(());
        };
        let normalised = match mapper.normalize_variant(SequenceVariant::CisPhased(cis.clone())) {
            Ok(n) => as_allele_on(&g.ac, n),
            // Refused as not yet normalised: an error, not a wrong answer.
            Err(hgvs_weaver::HgvsError::UnsupportedOperation(_)) => return Ok(()),
            Err(e) => {
                return Err(TestCaseError::fail(format!("{cis} reads but does not normalise: {e}")))
            }
        };
        let after = mapper
            .read_cis_on_genome(&normalised)
            .map_err(|e| TestCaseError::fail(format!("{cis} normalises to {normalised}, which does not read on the genome: {e}")))?;
        let reference = mapper.refs.reference(before.accession(), before.kind());
        prop_assert!(before.same_sequence(&after, &reference).unwrap(), "{} and {}", cis, normalised);
        let after_on_transcript = mapper.read_cis(&normalised).unwrap();
        prop_assert_eq!(on_transcript.stated(), after_on_transcript.stated(), "{} and {}", cis, normalised);
        let twice = mapper
            .normalize_variant(SequenceVariant::CisPhased(normalised.clone()))
            .map_err(|e| TestCaseError::fail(format!("{normalised} does not normalise again: {e}")))?;
        prop_assert_eq!(as_allele_on(&g.ac, twice).to_string(), normalised.to_string(), "{}", cis);
    }
}

/// The transcript property is not vacuous: on a fixed sample, a good share
/// of the dense alleles read on both the transcript and the genome, all but
/// one in fifty of those normalise, and some normalise to fewer members than
/// they are written with.
#[test]
fn the_generated_transcript_alleles_are_mostly_read() {
    use proptest::strategy::ValueTree;
    let mut runner = proptest::test_runner::TestRunner::deterministic();
    let strategy = gene_and_dense_cis();
    let cases = 1000;
    let (mut read, mut normalised, mut merged) = (0, 0, 0);
    for _ in 0..cases {
        let (g, members) = strategy.new_tree(&mut runner).unwrap().current();
        let hdp = g.provider();
        let mapper = VariantMapper::new(&hdp);
        let cis = coding_cis(&g, &members);
        if mapper.read_cis(&cis).is_ok() && mapper.read_cis_on_genome(&cis).is_ok() {
            read += 1;
            let n = mapper
                .normalize_variant(SequenceVariant::CisPhased(cis.clone()))
                .map(|n| as_allele_on(&g.ac, n).members.len());
            normalised += usize::from(n.is_ok());
            merged += usize::from(n.is_ok_and(|n| n < cis.members.len()));
        }
    }
    assert!(read * 10 >= cases * 3, "{read} of {cases} alleles read");
    assert!(merged * 10 >= cases, "{merged} of {cases} alleles merged");
    assert!(
        normalised * 50 >= read * 49,
        "{normalised} of {read} alleles that read normalise"
    );
}

fn as_allele_on(ac: &str, v: SequenceVariant) -> hgvs_weaver::CisPhasedVariant {
    match v {
        SequenceVariant::CisPhased(c) => c,
        one => hgvs_weaver::CisPhasedVariant::new(ac.into(), None, vec![one]).unwrap(),
    }
}

/// A short sequence and an allele of up to four members on it, some of which
/// no sequence can be built from.
fn seq_and_allele() -> impl Strategy<Value = (String, Vec<PlacedHgvs>)> {
    dna(8, 24).prop_flat_map(|seq| {
        let len = seq.len();
        (
            Just(seq),
            prop::collection::vec(cis_member_or_unapplied(len), 1..=4),
        )
    })
}

/// A description an allele normalises to, as an allele: a plain variant is
/// the allele of its one member.
fn as_allele(v: SequenceVariant) -> hgvs_weaver::CisPhasedVariant {
    match v {
        SequenceVariant::CisPhased(c) => c,
        one => hgvs_weaver::CisPhasedVariant::new("X".into(), None, vec![one]).unwrap(),
    }
}

/// The allele `members` normalised, where it reads; `None` where reading it
/// refuses it, which normalising it must then do too, or where normalising
/// refuses it as not yet normalised, which
/// [`the_generated_alleles_mostly_normalise`] bounds.
fn normalised_if_read(
    mapper: &VariantMapper,
    members: &[PlacedHgvs],
) -> Result<Option<hgvs_weaver::CisPhasedVariant>, TestCaseError> {
    let cis = genomic_cis(members);
    let read = mapper.read_cis(&cis);
    let normalised = mapper.normalize_variant(SequenceVariant::CisPhased(cis.clone()));
    match (read, normalised) {
        (Ok(_), Ok(n)) => Ok(Some(as_allele(n))),
        (Ok(_), Err(hgvs_weaver::HgvsError::UnsupportedOperation(_))) => Ok(None),
        (Ok(_), Err(e)) => Err(TestCaseError::fail(format!(
            "{cis} reads but does not normalise: {e}"
        ))),
        (Err(_), Err(_)) => Ok(None),
        (Err(e), Ok(n)) => Err(TestCaseError::fail(format!(
            "{cis} is refused ({e}) but normalises to {n}"
        ))),
    }
}

/// Normalising refuses few of the alleles that read: on a fixed sample, all
/// but one in fifty normalise.
#[test]
fn the_generated_alleles_mostly_normalise() {
    use proptest::strategy::ValueTree;
    let mut runner = proptest::test_runner::TestRunner::deterministic();
    let strategy = seq_and_allele();
    let (mut read, mut normalised) = (0, 0);
    for _ in 0..2000 {
        let (seq, members) = strategy.new_tree(&mut runner).unwrap().current();
        let hdp = Provider::new().sequence("X", &seq);
        let mapper = VariantMapper::new(&hdp);
        let cis = genomic_cis(&members);
        if mapper.read_cis(&cis).is_ok() {
            read += 1;
            normalised += usize::from(
                mapper
                    .normalize_variant(SequenceVariant::CisPhased(cis))
                    .is_ok(),
            );
        }
    }
    assert!(
        normalised * 50 >= read * 49,
        "{normalised} of {read} alleles that read normalise"
    );
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(512))]

    /// What an allele normalises to reads, and makes the sequence the allele
    /// makes. `normalize` checks this itself and refuses where it fails, so
    /// the property guards that check; how rarely it refuses is
    /// [`the_generated_alleles_mostly_normalise`]'s.
    #[test]
    fn normalising_a_cis_allele_keeps_its_sequence((seq, members) in seq_and_allele()) {
        let hdp = Provider::new().sequence("X", &seq);
        let mapper = VariantMapper::new(&hdp);
        let Some(normalised) = normalised_if_read(&mapper, &members)? else {
            return Ok(());
        };
        let cis = genomic_cis(&members);
        let reference = mapper.refs.reference("X", G);
        let before = mapper.read_cis(&cis).unwrap();
        let after = mapper
            .read_cis(&normalised)
            .map_err(|e| TestCaseError::fail(format!("{cis} normalises to {normalised}, which does not read: {e}")))?;
        prop_assert!(before.same_sequence(&after, &reference).unwrap(), "{} and {}", cis, normalised);
    }

    /// What an allele normalises to states the bases the allele states, as
    /// `normalize` checks before it answers (see
    /// `normalising_a_cis_allele_keeps_its_sequence`).
    #[test]
    fn normalising_a_cis_allele_keeps_its_stated_bases((seq, members) in seq_and_allele()) {
        let hdp = Provider::new().sequence("X", &seq);
        let mapper = VariantMapper::new(&hdp);
        let Some(normalised) = normalised_if_read(&mapper, &members)? else {
            return Ok(());
        };
        let cis = genomic_cis(&members);
        let before = mapper.read_cis(&cis).unwrap();
        let after = mapper
            .read_cis(&normalised)
            .map_err(|e| TestCaseError::fail(format!("{cis} normalises to {normalised}, which does not read: {e}")))?;
        prop_assert_eq!(before.stated(), after.stated(), "{} and {}", cis, normalised);
    }

    /// Normalising what an allele normalises to changes nothing.
    #[test]
    fn normalising_a_cis_allele_twice_is_normalising_it_once((seq, members) in seq_and_allele()) {
        let hdp = Provider::new().sequence("X", &seq);
        let mapper = VariantMapper::new(&hdp);
        let Some(once) = normalised_if_read(&mapper, &members)? else {
            return Ok(());
        };
        let twice = mapper
            .normalize_variant(SequenceVariant::CisPhased(once.clone()))
            .map_err(|e| TestCaseError::fail(format!("{once} does not normalise again: {e}")))?;
        prop_assert_eq!(as_allele(twice).to_string(), once.to_string(), "{}", genomic_cis(&members));
    }

    /// The order members are written in does not change what the allele
    /// normalises to: its statements, its changes, or the error refusing it.
    #[test]
    fn normalising_a_cis_allele_does_not_depend_on_member_order(
        (seq, members, shuffled) in seq_and_allele().prop_flat_map(|(seq, m)| {
            (Just(seq), Just(m.clone()), Just(m).prop_shuffle())
        })
    ) {
        let hdp = Provider::new().sequence("X", &seq);
        let mapper = VariantMapper::new(&hdp);
        let normalise = |m: &[PlacedHgvs]| {
            mapper
                .normalize_variant(SequenceVariant::CisPhased(genomic_cis(m)))
                .map(|n| n.to_string())
                .map_err(|e| std::mem::discriminant(&e))
        };
        prop_assert_eq!(
            normalise(&members),
            normalise(&shuffled),
            "{} and {}", genomic_cis(&members), genomic_cis(&shuffled)
        );
    }
}

// The protein an allele predicts, and what equivalence makes of alleles that
// state bases beside a change.

/// The changes of a dense transcript allele inside the codons before the
/// stop, or `None` when it does not read, states bases or changes anything
/// outside them: as one delins, the allele that
/// [`an_allele_leaves_the_protein_its_one_delins_leaves`] compares it with.
fn delins_inside_the_codons(
    mapper: &VariantMapper,
    g: &Gene,
    cis: &hgvs_weaver::CisPhasedVariant,
) -> Option<hgvs_weaver::CisPhasedVariant> {
    let reading = mapper.read_cis(cis).ok()?;
    let reference = mapper
        .refs
        .reference(&g.ac, IdentifierType::TranscriptAccession);
    let made = reading
        .sequence(&reference, 0, g.transcript_seq.len())
        .unwrap();
    let one = as_one_member(&g.transcript_seq, &made)?;
    let inside = one.start >= g.cds_start && one.end + 3 <= g.cds_end + 1;
    (reading.stated().is_empty() && inside).then(|| coding_cis(g, std::slice::from_ref(&one)))
}

/// The delins property is not vacuous: on a fixed sample, a good share of
/// the dense alleles change only codons before the stop and state nothing.
#[test]
fn the_generated_alleles_inside_the_codons_are_mostly_compared() {
    use proptest::strategy::ValueTree;
    let mut runner = proptest::test_runner::TestRunner::deterministic();
    let strategy = gene_and_dense_cis();
    let cases = 1000;
    let mut compared = 0;
    for _ in 0..cases {
        let (g, members) = strategy.new_tree(&mut runner).unwrap().current();
        let hdp = g.provider();
        let mapper = VariantMapper::new(&hdp);
        compared +=
            usize::from(delins_inside_the_codons(&mapper, &g, &coding_cis(&g, &members)).is_some());
    }
    assert!(
        compared * 10 >= cases,
        "{compared} of {cases} alleles compared"
    );
}

/// The bases `members` make on `X`, as one delins, followed by an identity
/// of each stretch they state: another spelling of what they make and state.
/// `None` when they do not read or make the reference.
fn respelled_allele(
    mapper: &VariantMapper,
    seq: &str,
    members: &[PlacedHgvs],
) -> Option<Vec<PlacedHgvs>> {
    let reading = mapper.read_cis(&genomic_cis(members)).ok()?;
    let reference = mapper.refs.reference("X", G);
    let made = reading.sequence(&reference, 0, seq.len()).unwrap();
    let mut respelled = vec![as_one_member(seq, &made)?];
    respelled.extend(reading.stated().iter().map(|&(start, end)| PlacedHgvs {
        start,
        end,
        edit: NaEdit::RefAlt {
            ref_: None,
            alt: None,
            uncertain: false,
        },
    }));
    Some(respelled)
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(512))]

    /// Equivalence between nucleotide alleles is the equality of one value,
    /// the sequence made and the bases stated, so it is symmetric and
    /// transitive however the members split the change between them and the
    /// bases they state unchanged.
    #[test]
    fn equivalence_is_transitive_over_statements_and_changes((seq, a, b) in seq_and_two_alleles()) {
        let hdp = Provider::new().sequence("X", &seq);
        let mapper = VariantMapper::new(&hdp);
        let eq = VariantEquivalence::new(&mapper, &hdp);
        let mut alleles: Vec<Vec<PlacedHgvs>> = vec![a.clone(), b];
        alleles.extend(respelled_allele(&mapper, &seq, &a));
        alleles.extend(respelled_allele(&mapper, &seq, &a).map(|mut r| {
            r.truncate(1);
            r
        }));
        let alleles: Vec<SequenceVariant> = alleles
            .iter()
            .map(|m| genomic_cis(m))
            .filter(|cis| mapper.read_cis(cis).is_ok())
            .map(SequenceVariant::CisPhased)
            .collect();
        let same = |x: &SequenceVariant, y: &SequenceVariant| {
            eq.equivalent(x, y).unwrap_or_else(|e| panic!("{x} vs {y}: {e}"))
        };
        for x in &alleles {
            for y in &alleles {
                prop_assert_eq!(same(x, y), same(y, x), "{} and {}: not symmetric", x, y);
                for z in &alleles {
                    if same(x, y) && same(y, z) {
                        prop_assert!(same(x, z), "{x} is {y}, and {y} is {z}, but {x} is not {z}");
                    }
                }
            }
        }
    }

    /// A transcript allele whose changes lie inside the codons before the
    /// stop, and which states nothing, leaves the protein the one delins it
    /// makes leaves: the members are read together, not one at a time.
    #[test]
    fn an_allele_leaves_the_protein_its_one_delins_leaves((g, members) in gene_and_dense_cis()) {
        let hdp = g.provider();
        let mapper = VariantMapper::new(&hdp);
        let cis = coding_cis(&g, &members);
        let Some(delins) = delins_inside_the_codons(&mapper, &g, &cis) else {
            return Ok(());
        };
        let left = |c: &hgvs_weaver::CisPhasedVariant| {
            mapper.predicted_protein_cis(c, None).map_err(|e| std::mem::discriminant(&e))
        };
        prop_assert_eq!(left(&cis), left(&delins), "{} and {}", cis, delins);
    }
}
