mod support;

use hgvs_weaver::data::Strand;
use hgvs_weaver::*;
use support::{exon, transcript, Provider};

/// Ten A's, then ATG (n.11 is c.1), then ATGC repeated; the genome carries the
/// same string at index 1000, where the exon puts it. The CDS is transcript
/// indices 10..=50.
fn provider() -> Provider {
    let seq = format!("AAAAAAAAAAATG{}", "ATGC".repeat(25));
    Provider::new()
        .sequence("NM_0001.3", &seq)
        .sequence("NC_0001.10", &format!("{}{seq}", "N".repeat(1000)))
        .transcript(transcript(
            "NM_0001.3",
            "NC_0001.10",
            Strand::Plus,
            Some((10, 50)),
            vec![exon((0, 100), (1000, 1100), Strand::Plus)],
        ))
        .protein_for("NM_0001.3", "NP_0001.1")
}

/// The prediction `c_to_p` writes for `c` on `NM_0001.3`.
fn predicted(c: &str) -> String {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let SequenceVariant::Coding(v) = parse_hgvs_variant(c).unwrap() else {
        panic!("{c} is not c.");
    };
    mapper.c_to_p(&v, Some("NP_0001.1")).unwrap().to_string()
}

#[test]
fn a_substitution_in_the_atg_start_codon_is_an_unknown_change_from_met1() {
    // ATG becomes GTG.
    assert_eq!(predicted("NM_0001.3:c.1A>G"), "NP_0001.1:p.(Met1?)");
}

#[test]
fn a_frameshift_in_the_start_codon_is_an_unknown_change_from_met1() {
    // ATG ATG CAT becomes AGA TGC ATG.
    assert_eq!(predicted("NM_0001.3:c.2del"), "NP_0001.1:p.(Met1?)");
}

#[test]
fn a_two_base_delins_inside_the_start_codon_is_an_unknown_change_from_met1() {
    // ATG becomes CTG.
    assert_eq!(predicted("NM_0001.3:c.1_2delinsCT"), "NP_0001.1:p.(Met1?)");
}

#[test]
fn an_insertion_inside_the_start_codon_is_an_unknown_change_from_met1() {
    // ATG becomes ATT G..
    assert_eq!(predicted("NM_0001.3:c.2_3insT"), "NP_0001.1:p.(Met1?)");
}

#[test]
fn a_delins_that_rewrites_the_start_codon_as_itself_leaves_met1_unchanged() {
    assert_eq!(predicted("NM_0001.3:c.1_3delinsATG"), "NP_0001.1:p.(Met1=)");
}

#[test]
fn a_deletion_from_the_start_codon_to_the_end_of_the_transcript_is_an_unknown_effect() {
    // Nothing is left to translate.
    assert_eq!(predicted("NM_0001.3:c.1_*62del"), "NP_0001.1:p.?");
}

#[test]
fn an_insertion_inside_the_start_codon_that_leaves_it_is_a_frameshift_from_met2() {
    // c.2_3insG is c.3dup: ATG GAT GCA, the start codon as it was.
    assert_eq!(
        predicted("NM_0001.3:c.2_3insG"),
        "NP_0001.1:p.(Met2AspfsTer?)"
    );
}

#[test]
fn an_insertion_after_the_start_codon_is_predicted_as_an_insertion() {
    // ATG ATG becomes ATG GGG ATG.
    assert_eq!(
        predicted("NM_0001.3:c.3_4insGGG"),
        "NP_0001.1:p.(Met1_Met2insGly)"
    );
}

#[test]
fn test_mapper_c_to_p_missense() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);

    // Mock CDS: ATG(Met1) ATG(Met2) CAT(His3) GCA(Ala4)...
    // c.7C>T changes codon 3 CAT(His) → TAT(Tyr): p.(His3Tyr)
    let var_c = parse_hgvs_variant("NM_0001.3:c.7C>T").unwrap();
    if let SequenceVariant::Coding(v) = var_c {
        let var_p = mapper.c_to_p(&v, Some("NP_0001.1")).unwrap();
        assert_eq!(var_p.to_string(), "NP_0001.1:p.(His3Tyr)");
    } else {
        panic!("Expected coding variant");
    }
}

#[test]
fn test_mapper_c_to_p_frameshift() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);

    // c.7del removes C from codon 3 (CAT=His), causing a frameshift: p.(His3...fsTer...)
    let var_c = parse_hgvs_variant("NM_0001.3:c.7del").unwrap();
    if let SequenceVariant::Coding(v) = var_c {
        let var_p = mapper.c_to_p(&v, Some("NP_0001.1")).unwrap();
        let p_str = var_p.to_string();
        assert!(
            p_str.contains("His3") && p_str.contains("fsTer"),
            "Expected frameshift at His3, got: {}",
            p_str
        );
    } else {
        panic!("Expected coding variant");
    }
}

#[test]
fn test_mapper_g_to_c_3utr() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);

    // Genomic 1052 (index 1051) -> n.52 -> c.*1, a G in both sequences.
    let var_g = parse_hgvs_variant("NC_0001.10:g.1052G>T").unwrap();
    if let SequenceVariant::Genomic(v) = var_g {
        let var_c = mapper.g_to_c(&v, "NM_0001.3").unwrap();
        assert_eq!(var_c.to_string(), "NM_0001.3:c.*1G>T");
    }
}

#[test]
fn test_mapper_c_to_g_3utr() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);

    let var_c = parse_hgvs_variant("NM_0001.3:c.*1G>T").unwrap();
    if let SequenceVariant::Coding(v) = var_c {
        let var_g = mapper.c_to_g(&v, Some("NC_0001.10")).unwrap();
        assert_eq!(var_g.to_string(), "NC_0001.10:g.1052G>T");
    }
}
