//! Property tests: relationships that must hold for every input the
//! strategies in `strategies.rs` can build. Each one pins a bug class the
//! ClinVar gate cannot see because ClinVar happens not to contain it.

mod strategies;
#[path = "../support/mod.rs"]
mod support;

use hgvs_weaver::coords::HgvsGenomicPos;
use hgvs_weaver::data::IdentifierType;
use hgvs_weaver::normalize::{ambiguous_range, normalize, ExonBounds, PlacedEdit};
use hgvs_weaver::reference::ReferenceStore;
use hgvs_weaver::structs::{GVariant, LinearVariant, SimpleInterval, SimplePosition};
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
}

/// The g. variant on `ac` for an edit over an HGVS range.
fn genomic_variant(ac: &str, p: &PlacedHgvs) -> GVariant {
    let position = |i: usize| SimplePosition {
        base: HgvsGenomicPos(i as i32).into(),
        end: None,
        uncertain: false,
    };
    GVariant::from_parts(
        ac.to_string(),
        None,
        PosEdit {
            pos: Some(SimpleInterval {
                start: position(p.start + 1),
                end: (p.end > p.start + 1).then(|| position(p.end)),
                uncertain: false,
            }),
            edit: p.edit.clone(),
            uncertain: false,
            predicted: false,
        },
    )
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

    /// The protein allele of a coding change, applied to the protein, is the
    /// protein the edited transcript encodes, to its first stop, for every
    /// kind of change: silent, missense, nonsense, in-frame, frameshift,
    /// stop loss.
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
        let allele = mapper.protein_allele(&vc, Some("NP_PROP.1")).unwrap();
        let protein = hgvs_weaver::utils::translate(&g.transcript_seq[g.cds_start..=g.cds_end]);
        let protein = protein.trim_end_matches('*');
        prop_assert_eq!(&allele.accession, "NP_PROP.1");
        let applied = format!("{}{}{}", &protein[..allele.start], allele.alternate, &protein[allele.end..]);
        let edited = apply(&g.transcript_seq, &PlacedHgvs { start, end, edit });
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

    /// The same, for a deletion, duplication or insertion written just before
    /// an exon/exon junction, where a slide on the transcript could carry it
    /// into the next exon.
    #[test]
    fn normalising_beside_a_junction_keeps_the_genomic_allele(
        g in gene(),
        pick in 0usize..10_000,
        back in 0usize..3,
        len in 1usize..=3,
        kind in 0u8..3,
        alt in dna(1, 3),
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
            _ => (last, last + 2, NaEdit::Ins { alt: Some(alt), uncertain: false }),
        };
        let vc = SequenceVariant::Coding(coding_variant(&g, &am, start, end, edit));
        let Ok(before) = mapper.canonical_allele(&vc) else { return Ok(()) };
        let normalised = mapper.normalize_variant(vc.clone()).unwrap();
        let after = mapper
            .canonical_allele(&normalised)
            .map_err(|err| TestCaseError::fail(format!("{vc} normalised to {normalised}: {err}")))?;
        prop_assert_eq!(after, before, "{} normalised to {}", vc, normalised);
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
    let position = |p: TxPoint| BaseOffsetPosition {
        offset: (p.offset != 0).then_some(IntronicOffset(p.offset)),
        ..CVariant::position_from_index(am, p.index).unwrap()
    };
    CVariant {
        ac: g.ac.clone(),
        gene: None,
        posedit: PosEdit {
            pos: Some(BaseOffsetInterval {
                start: position(e.start),
                end: (e.end != e.start).then(|| position(e.end)),
                uncertain: false,
            }),
            edit: e.edit.clone(),
            uncertain: false,
            predicted: false,
        },
    }
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
