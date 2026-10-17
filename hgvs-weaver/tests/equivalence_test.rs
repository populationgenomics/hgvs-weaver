mod support;

use hgvs_weaver::data::Strand;
use hgvs_weaver::equivalence::{EquivalenceLevel, VariantEquivalence};
use hgvs_weaver::error::HgvsError;
use hgvs_weaver::mapper::VariantMapper;
use support::Provider;

/// ACGT repeated as the genome, and a protein NP_0001.1 with Trp2.
fn provider() -> Provider {
    Provider::new()
        .sequence("NC_TEST.1", &"ACGT".repeat(1000))
        .sequence("NP_0001.1", &format!("MW{}", "A".repeat(500)))
}

#[test]
fn test_equivalence_levels() -> Result<(), HgvsError> {
    let hdp = provider();
    let eq_mapper = VariantMapper::new(&hdp);
    let eq = VariantEquivalence::new(&eq_mapper, &hdp);

    // 1. Identity
    let v1 = hgvs_weaver::parse_hgvs_variant("NC_TEST.1:g.1001A>C")?;
    let v2 = hgvs_weaver::parse_hgvs_variant("NC_TEST.1:g.1001A>C")?;
    assert_eq!(eq.equivalent_level(&v1, &v2)?, EquivalenceLevel::Identity);

    // 2. Analogous (Normalization) - ins vs dup. The reference is ACGT repeated,
    // so g.1006 is a C; duplicating it and inserting a C after it are one change.
    let v3 = hgvs_weaver::parse_hgvs_variant("NC_TEST.1:g.1006dupC")?;
    let v4 = hgvs_weaver::parse_hgvs_variant("NC_TEST.1:g.1006_1007insC")?;
    let lvl = eq.equivalent_level(&v3, &v4)?;
    assert!(lvl == EquivalenceLevel::Identity || lvl == EquivalenceLevel::Analogous);

    // 3. Parity (Functional - same protein)
    let p1 = hgvs_weaver::parse_hgvs_variant("NP_0001.1:p.Trp2Ter")?;
    let p2 = hgvs_weaver::parse_hgvs_variant("NP_0001.1:p.Trp2*")?;
    assert_eq!(eq.equivalent_level(&p1, &p2)?, EquivalenceLevel::Identity);

    Ok(())
}

#[test]
fn equivalence_needs_the_sequence() -> Result<(), HgvsError> {
    let hdp = Provider::new().failing_sequences();
    let eq_mapper = VariantMapper::new(&hdp);
    let eq = VariantEquivalence::new(&eq_mapper, &hdp);

    // p.Ala201_Val202insGlyProGlyAla vs p.Gly198_Ala201dup may well be the same
    // change, but only the protein sequence can say: residues 199 and 200 are
    // named by neither. Without a sequence the judgement is an error, not a guess.
    let v1 = hgvs_weaver::parse_hgvs_variant("NP_0001.1:p.Ala201_Val202insGlyProGlyAla")?;
    let v2 = hgvs_weaver::parse_hgvs_variant("NP_0001.1:p.Gly198_Ala201dup")?;
    assert!(eq.equivalent_level(&v1, &v2).is_err());
    Ok(())
}

#[test]
fn a_lone_copy_number_compares_by_projection() -> Result<(), HgvsError> {
    // NM_X.1 is genomic indices 10..50 of NC_X.1, c.1 at index 15, so c.10 is
    // g.25; NR_XR.1 reads the same forty bases on the minus strand, so n.k is
    // g.(51 - k) and g.25_27 is n.24_26.
    let transcript = "GGGGGATGAAACTGGCCTATCGCTAACCGTATAAGTAAGG";
    let minus = "CCTTACTTATACGGTTAGCGATAGGCCAGTTTCATCCCCC";
    let hdp = Provider::new()
        .sequence("NC_X.1", &format!("TTTTTTTTTT{transcript}CCCCCCCCCC"))
        .sequence("NM_X.1", transcript)
        .sequence("NR_XR.1", minus)
        .transcript(support::single_exon_transcript(
            "NM_X.1",
            "NC_X.1",
            10,
            Strand::Plus,
            5,
            25,
            40,
        ))
        .transcript(support::transcript(
            "NR_XR.1",
            "NC_X.1",
            Strand::Minus,
            None,
            vec![support::exon((0, 40), (10, 49), Strand::Minus)],
        ));
    let mapper = VariantMapper::new(&hdp);
    let eq = VariantEquivalence::new(&mapper, &hdp);
    let level = |a: &str, b: &str| -> Result<EquivalenceLevel, HgvsError> {
        eq.equivalent_level(
            &hgvs_weaver::parse_hgvs_variant(a)?,
            &hgvs_weaver::parse_hgvs_variant(b)?,
        )
    };

    assert_eq!(
        level("NM_X.1:c.10_12copy3", "NC_X.1:g.25_27copy3")?,
        EquivalenceLevel::Identity
    );
    // Neither is the other's text; both project to g.25_27copy3.
    assert_eq!(
        level("NM_X.1:c.10_12copy3", "NR_XR.1:n.24_26copy3")?,
        EquivalenceLevel::Analogous
    );
    assert_eq!(
        level("NM_X.1:c.10_12copy3", "NR_XR.1:n.23_25copy3")?,
        EquivalenceLevel::Different
    );
    assert_eq!(
        level("NC_X.1:g.25_27copy3", "NC_X.1:g.25_27copy4")?,
        EquivalenceLevel::Different
    );
    Ok(())
}
