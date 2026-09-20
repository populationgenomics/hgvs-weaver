mod support;

use hgvs_weaver::coords::SequenceVariant;
use hgvs_weaver::data::Strand;
use hgvs_weaver::error::HgvsError;
use hgvs_weaver::mapper::VariantMapper;
use support::{exon, transcript, Provider};

/// Two transcripts over one ACGC repeat, with a T at index 690.
fn regression_provider() -> Provider {
    let mut seq = "ACGC".repeat(1000).into_bytes();
    // For NM_058216.3:c.692_694delinsAA
    // Ser231: 691, 692, 693
    // If 691 is T, 692-693 replaced by AA -> TAA (Stop)
    seq[690] = b'T';
    let seq = String::from_utf8(seq).unwrap();
    let mut provider = Provider::new();
    for (ac, np) in [
        ("NM_153046.3", "NP_694591.2"),
        ("NM_058216.3", "NP_478123.1"),
    ] {
        provider = provider
            .sequence(ac, &seq)
            .transcript(transcript(
                ac,
                "NC_TEST",
                Strand::Plus,
                Some((0, 2000)),
                vec![exon((0, 2000), (1000, 3000), Strand::Plus)],
            ))
            .protein_for(ac, np);
    }
    provider
}

#[test]
fn test_regression_c_360_eq() -> Result<(), HgvsError> {
    let hdp = regression_provider();
    let mapper = VariantMapper::new(&hdp);

    // NM_153046.3:c.360=
    let v_raw = hgvs_weaver::parse_hgvs_variant("NM_153046.3:c.360=")?;
    let SequenceVariant::Coding(var_c) = v_raw else {
        panic!()
    };
    let p_var = mapper.c_to_p(&var_c, None)?;

    // Should be synonymous (Thr120=)
    // Current bug makes it a frameshift deletion of c.360
    assert!(
        !p_var.to_string().contains("fs"),
        "Should not be a frameshift, got {}",
        p_var
    );
    Ok(())
}

#[test]
fn test_regression_delins_stop() -> Result<(), HgvsError> {
    let hdp = regression_provider();
    let mapper = VariantMapper::new(&hdp);

    // NM_058216.3:c.692_694delinsAA
    let v_raw = hgvs_weaver::parse_hgvs_variant("NM_058216.3:c.692_694delinsAA")?;
    let SequenceVariant::Coding(var_c) = v_raw else {
        panic!()
    };
    let p_var = mapper.c_to_p(&var_c, None)?;

    // We want it to be a stop codon if it created one, not a frameshift.
    // In our mock: 691=T. Insert AA -> TAA (Stop) at 231.
    assert!(
        p_var.to_string().contains("Ter") || p_var.to_string().contains("*"),
        "Should contain stop codon, got {}",
        p_var
    );
    assert!(
        !p_var.to_string().contains("fs"),
        "Should not be a frameshift if stop is earlier, got {}",
        p_var
    );

    Ok(())
}

/// ATG (1-3) CAG (4-6) CAG (7-9) CAG (10-12) TAG (13-15): M Q Q Q *
fn repeat_provider() -> Provider {
    Provider::new()
        .sequence("NM_001.1", "ATGCAGCAGCAGTAG")
        .transcript(transcript(
            "NM_001.1",
            "NC_001.1",
            Strand::Plus,
            Some((0, 14)),
            vec![],
        ))
        .protein_for("NM_001.1", "NP_001.1")
}

#[test]
fn test_regression_gln4del_vs_ter() -> Result<(), HgvsError> {
    let hdp = repeat_provider();
    let mapper = VariantMapper::new(&hdp);

    // c.4_6del
    let v_c = hgvs_weaver::parse_hgvs_variant("NM_001.1:c.4_6del")?;

    if let SequenceVariant::Coding(vc) = v_c {
        // Generate p.
        let vp = mapper.c_to_p(&vc, Some("NP_001.1"))?;
        let vp_str = vp.to_string();
        println!("Generated: {}", vp_str);

        // Assert we get p.Gln4del (or equivalent del) AND NOT p.Gln4Ter
        assert!(
            vp_str.contains("del"),
            "Expected 'del' in {}, got {}",
            vp_str,
            vp_str
        );
        assert!(
            !vp_str.contains("Ter"),
            "Did not expect 'Ter' in {}, got {}",
            vp_str,
            vp_str
        );
    } else {
        panic!("Parsed wrong type");
    }

    Ok(())
}

/// 5000 G's with CCA (Pro1500) at indices 4497..=4499.
fn delins_mismatch_provider() -> Provider {
    let seq = format!("{}CCA{}", "G".repeat(4497), "G".repeat(500));
    Provider::new()
        .sequence("NM_001008844.3", &seq)
        .transcript(transcript(
            "NM_001008844.3",
            "NC_000001.11",
            Strand::Plus,
            Some((0, 4500)),
            vec![exon((0, 5000), (0, 5000), Strand::Plus)],
        ))
        .protein_for("NM_001008844.3", "NP_001008844.1")
}

#[test]
fn test_regression_pro_ile_mismatch() -> Result<(), HgvsError> {
    let provider = delins_mismatch_provider();
    let mapper = VariantMapper::new(&provider);

    let v_nuc = hgvs_weaver::parse_hgvs_variant("NM_001008844.3:c.4498_4499delinsAT")?;

    if let SequenceVariant::Coding(c_var) = v_nuc {
        let v_p = mapper.c_to_p(&c_var, None)?;
        assert!(v_p.to_string().contains("Pro1500Ile"));
    } else {
        panic!("Parsed variant is not coding");
    }
    Ok(())
}

#[test]
fn test_regression_parse_clinvar_repeat() {
    use hgvs_weaver::parse_hgvs_variant;
    let v = parse_hgvs_variant("NP_001365049.1:p.490PRS[1]");
    assert!(v.is_ok(), "Failed to parse p.490PRS[1]: {:?}", v.err());
    let v_inner = v.unwrap();

    if let SequenceVariant::Protein(p) = v_inner {
        if let hgvs_weaver::edits::AaEdit::Repeat { ref_, min, max, .. } = p.posedit.edit {
            assert_eq!(ref_, Some("PRS".to_string()));
            assert_eq!(min, 1);
            assert_eq!(max, 1);
        } else {
            panic!(
                "Expected Repeat edit, got something else: {:?}",
                p.posedit.edit
            );
        }
    }
}
