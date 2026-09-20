mod support;

use hgvs_weaver::data::Strand;
use hgvs_weaver::*;
use support::{exon, transcript, Provider};

/// Ten A's, then ATG (n.11 is c.1), then ATGC repeated; the same string is
/// the transcript and the genome. The CDS is transcript indices 10..=50.
fn provider() -> Provider {
    let seq = format!("AAAAAAAAAAATG{}", "ATGC".repeat(25));
    Provider::new()
        .sequence("NM_0001.3", &seq)
        .sequence("NC_0001.10", &seq)
        .transcript(transcript(
            "NM_0001.3",
            "NC_0001.10",
            Strand::Plus,
            Some((10, 50)),
            vec![exon((0, 100), (1000, 1100), Strand::Plus)],
        ))
        .protein_for("NM_0001.3", "NP_0001.1")
}

#[test]
fn test_mapper_c_to_p_start_codon_subst() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);

    // c.1A>T changes ATG(Met) → TTG(Leu): predicts the specific amino acid change p.(Met1Leu)
    let var_c = parse_hgvs_variant("NM_0001.3:c.1A>T").unwrap();
    if let SequenceVariant::Coding(v) = var_c {
        let var_p = mapper.c_to_p(&v, Some("NP_0001.1")).unwrap();
        assert_eq!(var_p.to_string(), "NP_0001.1:p.(Met1Leu)");
    }
}

#[test]
fn test_mapper_c_to_p_start_codon_del() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);

    // c.2del removes 'T' from ATG start codon → frameshift from position 1
    let var_c = parse_hgvs_variant("NM_0001.3:c.2del").unwrap();
    if let SequenceVariant::Coding(v) = var_c {
        let var_p = mapper.c_to_p(&v, Some("NP_0001.1")).unwrap();
        assert!(
            var_p.to_string().contains("fsTer"),
            "Expected frameshift annotation, got: {}",
            var_p
        );
    }
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

    // Genomic 1052 (index 1051) -> n.52 -> c.*1
    let var_g = parse_hgvs_variant("NC_0001.10:g.1052A>T").unwrap();
    if let SequenceVariant::Genomic(v) = var_g {
        let var_c = mapper.g_to_c(&v, "NM_0001.3").unwrap();
        assert_eq!(var_c.to_string(), "NM_0001.3:c.*1A>T");
    }
}

#[test]
fn test_mapper_c_to_g_3utr() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);

    let var_c = parse_hgvs_variant("NM_0001.3:c.*1A>T").unwrap();
    if let SequenceVariant::Coding(v) = var_c {
        let var_g = mapper.c_to_g(&v, Some("NC_0001.10")).unwrap();
        assert_eq!(var_g.to_string(), "NC_0001.10:g.1052A>T");
    }
}
