//! Equivalence judgements on the cases that shaped them: ClinVar spellings
//! of truncations, repeats and frameshifts against weaver's predictions.

mod support;

use hgvs_weaver::data::Strand;
use hgvs_weaver::equivalence::{EquivalenceLevel, VariantEquivalence};
use hgvs_weaver::error::HgvsError;
use hgvs_weaver::mapper::VariantMapper;
use hgvs_weaver::parse_hgvs_variant;
use support::{transcript, Provider};

/// A protein of X's with `residues` in place from `at` (0-based), `total` long.
fn protein_with(at: usize, residues: &str, total: usize) -> String {
    format!(
        "{}{residues}{}",
        "X".repeat(at),
        "X".repeat(total - at - residues.len())
    )
}

#[test]
fn test_clinvar_regression_tyr165ter() -> Result<(), HgvsError> {
    // Regression test for NM_001350334.2:c.495_498del (Frameshift at codon 165)
    // Weaver: p.(Tyr165Ter) -> p.Tyr165Ter
    // ClinVar: p.Ala164_Tyr165insTer
    // These should now be Analogous thanks to offset fix.

    // 163: Leu, 164: Ala, 165: Tyr, 166: Arg
    let hdp = Provider::new().sequence("NP_001337263.1", &protein_with(163, "LAYR", 167));
    let eq_mapper = VariantMapper::new(&hdp);
    let eq = VariantEquivalence::new(&eq_mapper, &hdp);

    let v1 = parse_hgvs_variant("NP_001337263.1:p.Tyr165Ter")?;
    let v2 = parse_hgvs_variant("NP_001337263.1:p.Ala164_Tyr165insTer")?;

    let level = eq.equivalent_level(&v1, &v2)?;
    assert_eq!(level, EquivalenceLevel::Analogous);
    Ok(())
}

#[test]
fn test_analogous_protein_truncation() -> Result<(), HgvsError> {
    // User Case: p.Tyr1433_Lys1434delinsTer vs p.(Tyr1433_Val3056del)
    // Both result in truncation at 1433.
    // delinsTer -> ...Tyr1433*
    // del -> ...Tyr1433 (end of sequence)

    // ATM is 3056 residues long; the deletion runs to its end.
    let hdp = Provider::new().sequence("NP_000042.3", &"M".repeat(3056));
    let eq_mapper = VariantMapper::new(&hdp);
    let eq = VariantEquivalence::new(&eq_mapper, &hdp);

    // GT: NP_000042.3:p.Tyr1433_Lys1434delinsTer
    let v_gt = parse_hgvs_variant("NP_000042.3:p.Tyr1433_Lys1434delinsTer")?;
    // W: p.(Tyr1433_Val3056del)
    let v_w = parse_hgvs_variant("NP_000042.3:p.(Tyr1433_Val3056del)")?;

    let lvl = eq.equivalent_level(&v_gt, &v_w)?;
    assert!(matches!(
        lvl,
        EquivalenceLevel::Analogous | EquivalenceLevel::Identity
    ));
    Ok(())
}

#[test]
fn test_analogous_clinvar_tyr165ter_mismatch() -> Result<(), HgvsError> {
    let hdp = Provider::new().sequence("NP_001337263.1", &protein_with(163, "LAYR", 167));
    let eq_mapper = VariantMapper::new(&hdp);
    let eq = VariantEquivalence::new(&eq_mapper, &hdp);

    let v1 = parse_hgvs_variant("NP_001337263.1:p.Tyr165Ter")?;
    let v2 = parse_hgvs_variant("NP_001337263.1:p.Ala164_Tyr165insTer")?;

    let level = eq.equivalent_level(&v1, &v2)?;
    assert_eq!(level, EquivalenceLevel::Analogous);
    Ok(())
}

#[test]
fn test_analogous_repeat_equivalence() -> Result<(), HgvsError> {
    let hdp = Provider::new().sequence("NP_001365049.1", &protein_with(489, "PRS", 592));
    let eq_mapper = VariantMapper::new(&hdp);
    let eq = VariantEquivalence::new(&eq_mapper, &hdp);

    let v1 = parse_hgvs_variant("NP_001365049.1:p.490PRS[1]")?;
    let v2 = parse_hgvs_variant("NP_001365049.1:p.=")?;
    let res = eq.equivalent_level(&v1, &v2)?;
    assert!(res.is_equivalent());

    let v3 = parse_hgvs_variant("NP_001365049.1:p.490PRS[2]")?;
    let v4 = parse_hgvs_variant("NP_001365049.1:p.490_492dup")?;
    let res2 = eq.equivalent_level(&v3, &v4)?;
    assert!(res2.is_equivalent());
    Ok(())
}

#[test]
fn test_analogous_fs_wildcard_unification() -> Result<(), HgvsError> {
    let hdp = Provider::new().sequence("NP_001.1", &protein_with(96, "R", 197));
    let eq_mapper = VariantMapper::new(&hdp);
    let eq = VariantEquivalence::new(&eq_mapper, &hdp);

    let v1 = parse_hgvs_variant("NP_001.1:p.Arg97ProfsTer4")?;
    let v2 = parse_hgvs_variant("NP_001.1:p.Arg97_Arg97delinsProAlaValTer")?;
    let res = eq.equivalent_level(&v1, &v2)?;
    assert!(res.is_equivalent());

    let v3 = parse_hgvs_variant("NP_001.1:p.Arg97_Arg97delinsProAlaValLeuTer")?;
    let res3 = eq.equivalent_level(&v1, &v3)?;
    assert!(!res3.is_equivalent());

    let v4 = parse_hgvs_variant("NP_001.1:p.Arg97_Arg97delinsLeuAlaValTer")?;
    let res4 = eq.equivalent_level(&v1, &v4)?;
    assert!(!res4.is_equivalent());

    let v5 = parse_hgvs_variant("NP_001.1:p.Arg97_Arg97delinsProAlaValLys")?;
    let res5 = eq.equivalent_level(&v1, &v5)?;
    assert!(!res5.is_equivalent());
    Ok(())
}

#[test]
fn test_multi_unit_repeat_equivalence() -> Result<(), HgvsError> {
    let hdp = Provider::new()
        // ...pppsvsatg pgpgpgpgpg pgpgpappny s...
        // Residue 229 starts the GP repeat.
        // 228 X's + 8 GP units (16 chars) + 100 X's
        .sequence(
            "NP_001122316.1",
            &protein_with(228, "GPGPGPGPGPGPGPGP", 344),
        )
        // ...lap apapapap apapvaapap apapapapap apapapdaap...
        // Residue 179 starts the AP repeat.
        // 178 X's + 8 AP units + 100 X's: AP[5] then deletes three units.
        .sequence("NP_000067.1", &protein_with(178, "APAPAPAPAPAPAPAP", 294));
    let eq_mapper = VariantMapper::new(&hdp);
    let eq = VariantEquivalence::new(&eq_mapper, &hdp);

    // Case 1: p.229GP[2] vs p.Gly233_Pro244del
    // 8 units -> 2 units = delete 6 units (12 residues)
    // 229GP refers to 229-230.
    // Gly233_Pro244del is 233 to 244 (12 residues).
    let v1 = parse_hgvs_variant("NP_001122316.1:p.229GP[2]")?;
    let v2 = parse_hgvs_variant("NP_001122316.1:p.Gly233_Pro244del")?;
    let res1 = eq.equivalent_level(&v1, &v2)?;
    assert!(
        res1.is_equivalent(),
        "GP repeat: {} vs {} -> {:?}",
        v1,
        v2,
        res1
    );

    // Case 2: p.179_180AP[5] vs p.Ala189_Pro194del
    // If W is 6 residues (3 units), and GT is AP[5], then initial was 8 units.
    // Ala189_Pro194del (residues 189, 190, 191, 192, 193, 194) is 6 residues.
    let v3 = parse_hgvs_variant("NP_000067.1:p.179_180AP[5]")?;
    let v4 = parse_hgvs_variant("NP_000067.1:p.Ala189_Pro194del")?;
    let res2 = eq.equivalent_level(&v3, &v4)?;
    assert!(
        res2.is_equivalent(),
        "AP repeat: {} vs {} -> {:?}",
        v3,
        v4,
        res2
    );

    Ok(())
}

#[test]
fn test_immediate_stop_normalization() -> Result<(), HgvsError> {
    // AAA (Lys) GGG (Gly)
    let hdp = Provider::new()
        .sequence("NM_1.1", "AAAGGG")
        .transcript(transcript(
            "NM_1.1",
            "NC_1.1",
            Strand::Plus,
            Some((0, 6)),
            vec![],
        ))
        .protein_for("NM_1.1", "NP_1.1");
    let mapper = VariantMapper::new(&hdp);

    // c.1_2delinsTA -> TAG (Stop) instead of ATG (Met)
    let var_c = parse_hgvs_variant("NM_1.1:c.1_2delinsTA")?;
    if let hgvs_weaver::coords::SequenceVariant::Coding(c) = var_c {
        let var_p = mapper.c_to_p(&c, None)?;
        let p_str = format!("{}", var_p);
        // Should be p.Lys1Ter (nonsense) not p.Lys1fsTer1
        assert_eq!(p_str, "NP_1.1:p.(Lys1Ter)");
    } else {
        panic!("Not a coding variant");
    }

    Ok(())
}
