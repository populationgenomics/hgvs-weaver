//! Property tests: relationships that must hold for every input the
//! strategies in `strategies.rs` can build. Each one pins a bug class the
//! ClinVar gate cannot see because ClinVar happens not to contain it.

mod strategies;

use hgvs_weaver::data::IdentifierType;
use hgvs_weaver::normalize::{ambiguous_range, normalize, PlacedEdit};
use hgvs_weaver::reference::ReferenceStore;
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
        let hdp = Provider::single("X", &seq);
        let store = ReferenceStore::with_block_size(&hdp, block);
        let r = store.reference("X", G);

        let before = PlacedEdit::from_hgvs_range(p.start, p.end, p.edit.clone());
        let after = normalize(&r, before.clone()).unwrap();
        let after_hgvs = placed_to_hgvs(&after);
        prop_assert_eq!(apply(&seq, &p), apply(&seq, &after_hgvs), "normalisation changed the molecule");

        let again = normalize(&r, after.clone()).unwrap();
        prop_assert_eq!(again, after, "normalisation is not idempotent");
    }

    /// The ambiguity range is exactly the set of positions the change can be
    /// written at: every placement inside gives the same molecule, and the
    /// placements one past either end do not.
    #[test]
    fn ambiguous_range_is_sound_and_complete((seq, p) in seq_and_edit(true)) {
        let hdp = Provider::single("X", &seq);
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
        if u_end + 1 <= seq.len() {
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
        let hdp = Provider::single("X", &seq);
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
use hgvs_weaver::structs::{
    BaseOffsetInterval, BaseOffsetPosition, CVariant, PVariant, PosEdit, TranscriptVariant,
};
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
                let (c, _, anchor) = am.n_to_c(TranscriptPos(last as i32)).unwrap();
                let pos = BaseOffsetPosition { base: c.to_hgvs(), offset: Some(IntronicOffset(i as i32)), anchor, uncertain: false };
                prop_assert_eq!(am.position_to_g(&pos).unwrap(), from_left);
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
    CVariant {
        ac: g.ac.clone(),
        gene: None,
        posedit: PosEdit {
            pos: Some(BaseOffsetInterval {
                start: s,
                end: e,
                uncertain: false,
            }),
            edit,
            uncertain: false,
            predicted: false,
        },
    }
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

fn up_to_stop(s: &str) -> String {
    match s.find('*') {
        Some(i) => s[..=i].to_string(),
        None => s.to_string(),
    }
}

proptest! {
    /// A single-base substitution in the CDS is described exactly as the codon
    /// table says: silent, missense, nonsense or a stop-loss extension.
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
        let expected = if r == a {
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

    /// A coding variant and its genomic projection are the same variant, and
    /// the mitochondrial spelling of a genomic variant is the same allele.
    #[test]
    fn a_variant_agrees_with_its_projection(g in gene(), pick in 0usize..10_000, alt in prop::sample::select(BASES.to_vec()), del in prop::bool::ANY) {
        let am = TranscriptMapper::new(g.transcript_data()).unwrap();
        let hdp = g.provider();
        let mapper = VariantMapper::new(&hdp);
        let eq = VariantEquivalence::new(&hdp, &hdp);
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
}
