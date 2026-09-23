mod support;

use hgvs_weaver::*;
use support::Provider;

fn provider() -> Provider {
    Provider::from_json_file("../tests/data/toy_data.json")
        .protein_for("NM_PLUS.1", "NP_PLUS.1")
        .protein_for("NM_MINUS.1", "NP_MINUS.1")
}

#[test]
fn test_toy_plus_strand_mapping() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);

    // Genomic variant NC_TOY.1:g.25A>T → c.1A>T (g-to-c positional mapping)
    let var_g = parse_hgvs_variant("NC_TOY.1:g.25A>T").unwrap();
    if let SequenceVariant::Genomic(v) = var_g {
        let var_c = mapper.g_to_c(&v, "NM_PLUS.1").unwrap();
        assert_eq!(var_c.to_string(), "NM_PLUS.1:c.1A>T");
    }

    // c.1A>T changes ATG(Met) → TTG(Leu): predicts p.(Met1Leu)
    let var_c = parse_hgvs_variant("NM_PLUS.1:c.1A>T").unwrap();
    if let SequenceVariant::Coding(v) = var_c {
        let var_p = mapper.c_to_p(&v, Some("NP_PLUS.1")).unwrap();
        assert_eq!(var_p.to_string(), "NP_PLUS.1:p.(Met1Leu)");
    }
}

#[test]
fn test_toy_plus_strand_missense() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);

    // NM_PLUS.1 CDS: ATG(Met1) CGT(Arg2) ACG(Thr3) ...
    // c.7A>G changes codon 3 ACG(Thr) → GCG(Ala): p.(Thr3Ala)
    let var_c = parse_hgvs_variant("NM_PLUS.1:c.7A>G").unwrap();
    if let SequenceVariant::Coding(v) = var_c {
        let var_p = mapper.c_to_p(&v, Some("NP_PLUS.1")).unwrap();
        assert_eq!(var_p.to_string(), "NP_PLUS.1:p.(Thr3Ala)");
    } else {
        panic!("Expected coding variant");
    }
}

#[test]
fn test_toy_minus_strand_mapping() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    // NM_MINUS.1's first exon is genome 200..=235 read on the minus strand, so
    // c.32 (transcript index 35) is genome index 200, g.201: an A on the
    // genome, a T on the transcript.
    let var_g = parse_hgvs_variant("NC_TOY.1:g.201A>G").unwrap();
    if let SequenceVariant::Genomic(v) = var_g {
        let var_c = mapper.g_to_c(&v, "NM_MINUS.1").unwrap();
        assert_eq!(var_c.to_string(), "NM_MINUS.1:c.32T>C");
    }
    let var_c = parse_hgvs_variant("NM_MINUS.1:c.32T>C").unwrap();
    if let SequenceVariant::Coding(v) = var_c {
        let var_g = mapper.c_to_g(&v, Some("NC_TOY.1")).unwrap();
        assert_eq!(var_g.to_string(), "NC_TOY.1:g.201A>G");
    }
}
