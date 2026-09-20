mod support;

use hgvs_weaver::data::Strand;
use hgvs_weaver::equivalence::{EquivalenceLevel, VariantEquivalence};
use hgvs_weaver::error::HgvsError;
use hgvs_weaver::mapper::VariantMapper;
use hgvs_weaver::parse_hgvs_variant;
use support::{exon, transcript, Provider};

/// 4000 N's with `bases` set at the given indices.
fn placed(bases: &[(usize, u8)]) -> String {
    let mut seq = vec![b'N'; 4000];
    for &(index, base) in bases {
        seq[index] = base;
    }
    String::from_utf8(seq).unwrap()
}

fn provider() -> Provider {
    // Case 9: c.35 is transcript index 34, genomic 3966 (A). The base after it
    // must differ so the insertion cannot shift further.
    // Case 15: c.2673 is transcript index 2672, genomic 1328 (T).
    let cases = placed(&[(3966, b'A'), (34, b'A'), (1328, b'T'), (2672, b'T')]);
    // BRAF Val600 is GTG at transcript indices 1797..=1799.
    let braf = placed(&[(1797, b'G'), (1798, b'T'), (1799, b'G')]);
    let mut provider = Provider::new()
        .sequence("NC_000001.1", &cases)
        .sequence("NP_BRAF", &cases)
        .sequence("NM_BRAF", &braf)
        .sequence("NC_BRAF", &braf)
        .transcript(transcript(
            "NM_BRAF",
            "NC_BRAF",
            Strand::Plus,
            Some((0, 3000)),
            vec![],
        ));
    for ac in ["NM_001166478.1", "NM_005813.3"] {
        // One minus-strand exon: transcript index i is genomic index 4000 - i.
        provider = provider.sequence(ac, &cases).transcript(transcript(
            ac,
            "NC_000001.1",
            Strand::Minus,
            Some((0, 3000)),
            vec![exon((0, 3001), (1000, 4000), Strand::Minus)],
        ));
    }
    provider
}

#[test]
fn test_repro_case9() -> Result<(), HgvsError> {
    let hdp = provider();
    let eq_mapper = VariantMapper::new(&hdp);
    let eq = VariantEquivalence::new(&eq_mapper, &hdp);

    let v1 = parse_hgvs_variant("NM_001166478.1:c.35_36insT")?;
    let v2 = parse_hgvs_variant("NM_001166478.1:c.35dup")?;

    assert_eq!(eq.equivalent_level(&v1, &v2)?, EquivalenceLevel::Analogous);
    Ok(())
}

#[test]
fn test_repro_case15() -> Result<(), HgvsError> {
    let hdp = provider();
    let eq_mapper = VariantMapper::new(&hdp);
    let eq = VariantEquivalence::new(&eq_mapper, &hdp);

    let v1 = parse_hgvs_variant("NM_005813.3:c.2673insA")?;
    let v2 = parse_hgvs_variant("NM_005813.3:c.2673dup")?;

    assert_eq!(eq.equivalent_level(&v1, &v2)?, EquivalenceLevel::Analogous);
    Ok(())
}

#[test]
fn test_braf_identity() -> Result<(), HgvsError> {
    let hdp = provider();
    let eq_mapper = VariantMapper::new(&hdp);
    let eq = VariantEquivalence::new(&eq_mapper, &hdp);

    // c.1799T>A -> p.Val600Glu (predicted)
    let v1 = parse_hgvs_variant("NM_BRAF:c.1799T>A")?;

    // Target 1: p.Val600Glu (Observed/Experimental)
    // Comparisons of c. to p. (observed) are Analogous because c. implies p.(predicted).
    // The notation differs (parens vs no parens), so it's not Strict Identity.
    let v2_observed = parse_hgvs_variant("NP_BRAF:p.Val600Glu")?;
    assert_eq!(
        eq.equivalent_level(&v1, &v2_observed)?,
        EquivalenceLevel::Analogous
    );

    // Target 2: p.(Val600Glu) (Predicted)
    // Comparisons of c. to p.(predicted) should be Identity because c. implies p.(predicted) exactly.
    let v2_predicted = parse_hgvs_variant("NP_BRAF:p.(Val600Glu)")?;

    assert_eq!(
        eq.equivalent_level(&v1, &v2_predicted)?,
        EquivalenceLevel::Identity
    );

    Ok(())
}
