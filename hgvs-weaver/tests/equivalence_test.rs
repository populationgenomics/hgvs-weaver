mod support;

use hgvs_weaver::equivalence::{EquivalenceLevel, VariantEquivalence};
use hgvs_weaver::error::HgvsError;
use hgvs_weaver::mapper::VariantMapper;
use support::Provider;

/// ACGT repeated, served as the genome and as the "protein" of NP_0001.1 alike.
fn provider() -> Provider {
    let seq = "ACGT".repeat(1000);
    Provider::new()
        .sequence("NC_TEST.1", &seq)
        .sequence("NP_0001.1", &seq)
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
