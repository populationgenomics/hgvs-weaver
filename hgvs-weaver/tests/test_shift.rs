mod support;

use hgvs_weaver::error::HgvsError;
use hgvs_weaver::mapper::VariantMapper;
use hgvs_weaver::SequenceVariant;
use support::Provider;

/// A genome of 2000 A's.
fn homopolymer() -> Provider {
    Provider::new().sequence("NC_TEST.1", &"A".repeat(2000))
}

/// A genome of CAG repeated.
fn repeat() -> Provider {
    Provider::new().sequence("NC_TEST.1", &"CAG".repeat(1000))
}

#[test]
fn test_ins_3_prime_shifting() -> Result<(), HgvsError> {
    let hdp = homopolymer();
    let mapper = VariantMapper::new(&hdp);

    // NC_TEST.1:g.1005_1006insA
    let v1 = hgvs_weaver::parse_hgvs_variant("NC_TEST.1:g.1005_1006insA")?;
    let SequenceVariant::Genomic(v1_g) = v1 else {
        panic!()
    };
    let nv1 = mapper.normalize_variant(SequenceVariant::Genomic(v1_g))?;

    // NC_TEST.1:g.1006_1007insA
    let v2 = hgvs_weaver::parse_hgvs_variant("NC_TEST.1:g.1006_1007insA")?;
    let SequenceVariant::Genomic(v2_g) = v2 else {
        panic!()
    };
    let nv2 = mapper.normalize_variant(SequenceVariant::Genomic(v2_g))?;

    assert_eq!(nv1.to_string(), nv2.to_string());
    Ok(())
}

#[test]
fn test_multi_base_ins_3_prime_shifting() -> Result<(), HgvsError> {
    let hdp = repeat();
    let mapper = VariantMapper::new(&hdp);

    // Reference is CAGCAGCAG... (at 1000, 1003, 1006, ...)
    // Insertion g.1002_1003insCAG should shift to 3' end of the repeat.
    let v1 = hgvs_weaver::parse_hgvs_variant("NC_TEST.1:g.1002_1003insCAG")?;
    let SequenceVariant::Genomic(v1_g) = v1 else {
        panic!()
    };
    let nv1 = mapper.normalize_variant(SequenceVariant::Genomic(v1_g))?;

    // If it doesn't shift, nv1.to_string() will be "...g.1002_1003insCAG"
    // If it shifts 1 block, it might be "...g.1005_1006insCAG"
    // Since we have many repeats, it should shift far.
    println!("Normalized: {}", nv1);
    assert!(
        !nv1.to_string().contains("1002_1003"),
        "Should have shifted from original position 1002_1003. Got: {}",
        nv1
    );

    Ok(())
}

#[test]
fn test_ins_5_prime_shifting_no_panic() -> Result<(), HgvsError> {
    let hdp = homopolymer();
    let mapper = VariantMapper::new(&hdp);

    // NC_TEST.1:g.1005_1006insA
    // This will attempt to shift left (5') because the prefix is 'A's.
    let v1 = hgvs_weaver::parse_hgvs_variant("NC_TEST.1:g.1005_1006insA")?;

    // This calls expand_unambiguous_range -> shift_5_prime, which previously panicked.
    let spdi = mapper.to_spdi_unambiguous(&v1)?;
    assert!(spdi.contains("NC_TEST.1:0:")); // Should expand to start of sequence

    Ok(())
}
