mod support;

use hgvs_weaver::data::Strand;
use hgvs_weaver::*;
use support::{exon, transcript, Provider};

/// `(transcript, sequence, CDS start index, CDS end index)`. Every transcript
/// is one plus-strand exon whose alignment is marked minus, as the original
/// fixture had it.
const TRANSCRIPTS: &[(&str, &str, i32, i32)] = &[
    // Ten A's, then ATG AAA TAG (Met Lys *), then ATGC repeated.
    ("NM_0001.1", "AAAAAAAAAAATGAAATAG", 10, 19),
    ("NM_SHIFT_BUG", "CCATTTTTTT", 0, 30),
    ("NM_PREMATURE_STOP", "ATGCAACAAGATGATTAA", 0, 18), // M Q Q D D * (18 bases)
    ("NM_INFRAME_DEL", "ATGGCTGCATGCGATTAA", 0, 18),    // M A B C D * (18 bases)
    ("NM_CTERM_SUBST", "ATGGCTGCATGCGATTAA", 0, 18),    // M A B C D *
    ("NM_REPEAT_EXP", "ATGGCTGCTGCTTTTTAA", 0, 18),     // M A A A F *
    ("NM_REPEAT_CON", "ATGGCAGCAGCAGCATTTTAA", 0, 21),  // M A A A A F *
];

fn provider() -> Provider {
    let mut provider = Provider::new().protein_for("NM_0001.1", "NP_0001.1");
    for &(ac, seq, cds_start, cds_end) in TRANSCRIPTS {
        let seq = if ac == "NM_0001.1" {
            format!("{seq}{}", "ATGC".repeat(20))
        } else {
            seq.to_string()
        };
        provider = provider.sequence(ac, &seq).transcript(transcript(
            ac,
            "NC_0001.10",
            Strand::Plus,
            Some((cds_start, cds_end)),
            vec![exon((0, 100), (1000, 1100), Strand::Minus)],
        ));
    }
    provider
}

#[test]
fn test_nonsense_normalization() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);

    // c.4A>T changes AAA (Lys) to TAA (Stop).
    let var_c = parse_hgvs_variant("NM_0001.1:c.4A>T").unwrap();
    if let SequenceVariant::Coding(v) = var_c {
        let var_p = mapper.c_to_p(&v, Some("NP_0001.1")).unwrap();
        // Should be Lys2Ter
        assert_eq!(var_p.to_string(), "NP_0001.1:p.(Lys2Ter)");
    }
}

#[test]
fn test_extension_normalization() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);

    // c.7T>G changes TAG (Stop) to GAG (Glu).
    let var_c = parse_hgvs_variant("NM_0001.1:c.7T>G").unwrap();
    if let SequenceVariant::Coding(v) = var_c {
        let var_p = mapper.c_to_p(&v, Some("NP_0001.1")).unwrap();
        // Should be {alt}extTer{length}
        assert!(var_p.to_string().contains("extTer"));
    }
}

#[test]
fn test_normalization_shift_bug() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);

    // NM_SHIFT_BUG: Sequence "CCAT...". UTR=0.
    // c.1_2delinsAT. Ref=CC. Alt=AT.
    // Result Sequence: ATAT...
    // If shifted -> c.3_4delinsAT. Ref=AT. Alt=AT.
    // Result Sequence: CCAT...
    // The outputs are DIFFERENT. So it MUST NOT shift.

    let var_c = parse_hgvs_variant("NM_SHIFT_BUG:c.1_2delinsAT").unwrap();
    if let SequenceVariant::Coding(v) = var_c {
        let var_norm = mapper
            .normalize_variant(SequenceVariant::Coding(v))
            .unwrap();
        if let SequenceVariant::Coding(v_norm) = var_norm {
            // Should remain c.1_2
            // Because if it shifts to 3_4, it implies AT -> AT, which is Identity.
            let pos = v_norm.posedit.pos.unwrap();
            let start = pos.start.base.to_index().0;
            // 0-based index: 0
            assert_eq!(start, 0, "Variant should not have shifted!");
        }
    }
}

#[test]
fn test_premature_stop_formatting() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);

    // NM_PREMATURE_STOP: M Q Q D D *.
    // c.4_9delinsCATTAA.
    // Replace CAA CAA (Q Q) with CAT TAA (H *).
    // Result: M H *.
    // Without fix: p.Gln2_Asp5delinsHis. (Deletion of QQD...).
    // With fix: p.Gln2_Gln3delinsHisTer.

    let var_c = parse_hgvs_variant("NM_PREMATURE_STOP:c.4_9delinsCATTAA").unwrap();
    if let SequenceVariant::Coding(v) = var_c {
        let var_p = mapper.c_to_p(&v, Some("NP_MOCK")).unwrap();
        assert_eq!(var_p.to_string(), "NP_MOCK:p.(Gln2_Gln3delinsHisTer)");
    }
}

#[test]
fn test_inframe_deletion_tail() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);

    // NM_INFRAME_DEL: M A A C D *.
    // c.4_6del. Del A (second A).
    // Result: M A C D *.
    // Tail match: C D * (3 chars).
    // This is > 1 char. So it should be detected as Original Stop.
    // Result should be p.Ala2del. NOT delins...Ter.

    let var_c = parse_hgvs_variant("NM_INFRAME_DEL:c.4_6del").unwrap();
    if let SequenceVariant::Coding(v) = var_c {
        let var_p = mapper.c_to_p(&v, Some("NP_MOCK")).unwrap();
        assert_eq!(var_p.to_string(), "NP_MOCK:p.(Ala3del)");
    }
}

#[test]
fn test_cterm_substitution() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);

    // NM_CTERM_SUBST: M A B C D *.
    // c.15T>G. D (GAT) -> E (GAG).
    // Tail match: *. (1 char).
    // Ref mismatch length: 1 (D).
    // Alt mismatch length: 1 (E).
    // Should be p.Asp5Glu. NOT delins...Ter.

    let var_c = parse_hgvs_variant("NM_CTERM_SUBST:c.15T>G").unwrap();
    if let SequenceVariant::Coding(v) = var_c {
        let var_p = mapper.c_to_p(&v, Some("NP_MOCK")).unwrap();
        assert_eq!(var_p.to_string(), "NP_MOCK:p.(Asp5Glu)");
    }
}

#[test]
fn test_repeat_expansion() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);

    // NM_REPEAT_EXP: M A A A *. (ATG GCT GCT GCT TAA).
    // c.4GCT[5] -> p.Ala3_Ala4dup? Or p.Ala2_Ala3dup?
    // Ref has 3 copies. Var has 5. Net +2 copies (+6 bases).
    // Expecting duplication notation.

    let var_c = parse_hgvs_variant("NM_REPEAT_EXP:c.4GCT[5]").unwrap();
    if let SequenceVariant::Coding(v) = var_c {
        let var_p = mapper.c_to_p(&v, Some("NP_MOCK")).unwrap();
        assert!(
            var_p.to_string().contains("dup"),
            "Expected dup, got {}",
            var_p
        );
    }
}

#[test]
fn test_repeat_contraction() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);

    // NM_REPEAT_CON: M A A A A *. (ATG GCA GCA GCA GCA TAA).
    // c.4GCA[2].
    // Ref has 4 copies. Var has 2. Net -2 copies (-6 bases).
    // Expecting deletion notation.

    let var_c = parse_hgvs_variant("NM_REPEAT_CON:c.4GCA[2]").unwrap();
    if let SequenceVariant::Coding(v) = var_c {
        let var_p = mapper.c_to_p(&v, Some("NP_MOCK")).unwrap();
        assert!(
            var_p.to_string().contains("del"),
            "Expected del, got {}",
            var_p
        );
    }
}
