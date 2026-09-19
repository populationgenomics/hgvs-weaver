//! Equivalence judgements on the cases that shaped them: ClinVar spellings
//! of truncations, repeats and frameshifts against weaver's predictions.

use hgvs_weaver::mapper::VariantMapper;

#[test]
fn test_clinvar_regression_tyr165ter() -> Result<(), hgvs_weaver::error::HgvsError> {
    use hgvs_weaver::data::{DataProvider, IdentifierKind, IdentifierType, TranscriptData};
    use hgvs_weaver::equivalence::{EquivalenceLevel, VariantEquivalence};

    // Regression test for NM_001350334.2:c.495_498del (Frameshift at codon 165)
    // Weaver: p.(Tyr165Ter) -> p.Tyr165Ter
    // ClinVar: p.Ala164_Tyr165insTer
    // These should now be Analogous thanks to offset fix.

    struct LocalMockProvider;
    impl DataProvider for LocalMockProvider {
        fn get_transcript(
            &self,
            _ac: &str,
            _ref: Option<&str>,
        ) -> Result<TranscriptData, hgvs_weaver::error::HgvsError> {
            Err(hgvs_weaver::error::HgvsError::ValidationError(
                "Not implemented".into(),
            ))
        }
        fn get_seq(
            &self,
            _ac: &str,
            start: i32,
            end: Option<i32>,
            _kind: IdentifierType,
        ) -> Result<String, hgvs_weaver::error::HgvsError> {
            let mut seq = String::new();
            // 163: Leu, 164: Ala, 165: Tyr, 166: Arg
            for i in start..end.unwrap_or(167) {
                match i {
                    163 => seq.push('L'),
                    164 => seq.push('A'),
                    165 => seq.push('Y'),
                    166 => seq.push('R'),
                    _ => seq.push('X'),
                }
            }
            Ok(seq)
        }
        fn get_symbol_accessions(
            &self,
            _: &str,
            _: IdentifierKind,
            _: IdentifierKind,
        ) -> Result<Vec<(IdentifierType, String)>, hgvs_weaver::error::HgvsError> {
            Ok(vec![])
        }
        fn get_identifier_type(
            &self,
            _: &str,
        ) -> Result<IdentifierType, hgvs_weaver::error::HgvsError> {
            Ok(IdentifierType::ProteinAccession)
        }
    }

    struct LocalMockSearch;
    impl hgvs_weaver::data::TranscriptSearch for LocalMockSearch {
        fn get_transcripts_for_region(
            &self,
            _: &str,
            _: i32,
            _: i32,
        ) -> Result<Vec<String>, hgvs_weaver::error::HgvsError> {
            Ok(vec![])
        }
    }

    let hdp = LocalMockProvider;
    let search = LocalMockSearch;
    let eq_mapper = VariantMapper::new(&hdp);
    let eq = VariantEquivalence::new(&eq_mapper, &search);

    let v1 = hgvs_weaver::parse_hgvs_variant("NP_001337263.1:p.Tyr165Ter")?;
    let v2 = hgvs_weaver::parse_hgvs_variant("NP_001337263.1:p.Ala164_Tyr165insTer")?;

    let level = eq.equivalent_level(&v1, &v2)?;
    assert_eq!(level, EquivalenceLevel::Analogous);
    Ok(())
}

#[test]
fn test_analogous_protein_truncation() -> Result<(), hgvs_weaver::error::HgvsError> {
    use hgvs_weaver::data::{DataProvider, IdentifierKind, IdentifierType, TranscriptData};
    use hgvs_weaver::equivalence::{EquivalenceLevel, VariantEquivalence};

    // User Case: p.Tyr1433_Lys1434delinsTer vs p.(Tyr1433_Val3056del)
    // Both result in truncation at 1433.
    // delinsTer -> ...Tyr1433*
    // del -> ...Tyr1433 (end of sequence)

    struct TruncationMockProvider;
    impl DataProvider for TruncationMockProvider {
        fn get_transcript(
            &self,
            _ac: &str,
            _ref: Option<&str>,
        ) -> Result<TranscriptData, hgvs_weaver::error::HgvsError> {
            Err(hgvs_weaver::error::HgvsError::UnsupportedOperation(
                "Not needed".into(),
            ))
        }
        fn get_seq(
            &self,
            _ac: &str,
            _start: i32,
            _end: Option<i32>,
            kind: IdentifierType,
        ) -> Result<String, hgvs_weaver::error::HgvsError> {
            if kind == IdentifierType::ProteinAccession {
                // Mock a long protein sequence
                // ATM is 3056 residues long; the deletion runs to its end.
                return Ok("M".repeat(3056));
            }
            Ok("".to_string())
        }
        fn get_symbol_accessions(
            &self,
            _: &str,
            _: IdentifierKind,
            _: IdentifierKind,
        ) -> Result<Vec<(IdentifierType, String)>, hgvs_weaver::error::HgvsError> {
            Ok(vec![])
        }
        fn get_identifier_type(
            &self,
            _: &str,
        ) -> Result<IdentifierType, hgvs_weaver::error::HgvsError> {
            Ok(IdentifierType::ProteinAccession)
        }
    }

    struct TruncationMockSearch;
    impl hgvs_weaver::data::TranscriptSearch for TruncationMockSearch {
        fn get_transcripts_for_region(
            &self,
            _: &str,
            _: i32,
            _: i32,
        ) -> Result<Vec<String>, hgvs_weaver::error::HgvsError> {
            Ok(vec![])
        }
    }

    let hdp = TruncationMockProvider;
    let search = TruncationMockSearch;
    let eq_mapper = VariantMapper::new(&hdp);
    let eq = VariantEquivalence::new(&eq_mapper, &search);

    // GT: NP_000042.3:p.Tyr1433_Lys1434delinsTer
    let v_gt = hgvs_weaver::parse_hgvs_variant("NP_000042.3:p.Tyr1433_Lys1434delinsTer")?;
    // W: p.(Tyr1433_Val3056del)
    let v_w = hgvs_weaver::parse_hgvs_variant("NP_000042.3:p.(Tyr1433_Val3056del)")?;

    let lvl = eq.equivalent_level(&v_gt, &v_w)?;
    assert!(matches!(
        lvl,
        EquivalenceLevel::Analogous | EquivalenceLevel::Identity
    ));
    Ok(())
}

#[test]
fn test_analogous_clinvar_tyr165ter_mismatch() -> Result<(), hgvs_weaver::error::HgvsError> {
    use hgvs_weaver::data::{DataProvider, IdentifierKind, IdentifierType, TranscriptData};
    use hgvs_weaver::equivalence::{EquivalenceLevel, VariantEquivalence};

    struct MockReproProvider;
    impl DataProvider for MockReproProvider {
        fn get_transcript(
            &self,
            _ac: &str,
            _ref: Option<&str>,
        ) -> Result<TranscriptData, hgvs_weaver::error::HgvsError> {
            Err(hgvs_weaver::error::HgvsError::ValidationError(
                "Not implemented".into(),
            ))
        }
        fn get_seq(
            &self,
            _ac: &str,
            start: i32,
            end: Option<i32>,
            _kind: IdentifierType,
        ) -> Result<String, hgvs_weaver::error::HgvsError> {
            let mut seq = String::new();
            for i in start..end.unwrap_or(167) {
                match i {
                    163 => seq.push('L'),
                    164 => seq.push('A'),
                    165 => seq.push('Y'),
                    166 => seq.push('R'),
                    _ => seq.push('X'),
                }
            }
            Ok(seq)
        }
        fn get_symbol_accessions(
            &self,
            _: &str,
            _: IdentifierKind,
            _: IdentifierKind,
        ) -> Result<Vec<(IdentifierType, String)>, hgvs_weaver::error::HgvsError> {
            Ok(vec![])
        }
        fn get_identifier_type(
            &self,
            _: &str,
        ) -> Result<IdentifierType, hgvs_weaver::error::HgvsError> {
            Ok(IdentifierType::ProteinAccession)
        }
    }

    struct MockSearch;
    impl hgvs_weaver::data::TranscriptSearch for MockSearch {
        fn get_transcripts_for_region(
            &self,
            _: &str,
            _: i32,
            _: i32,
        ) -> Result<Vec<String>, hgvs_weaver::error::HgvsError> {
            Ok(vec![])
        }
    }

    let hdp = MockReproProvider;
    let search = MockSearch;
    let eq_mapper = VariantMapper::new(&hdp);
    let eq = VariantEquivalence::new(&eq_mapper, &search);

    let v1 = hgvs_weaver::parse_hgvs_variant("NP_001337263.1:p.Tyr165Ter")?;
    let v2 = hgvs_weaver::parse_hgvs_variant("NP_001337263.1:p.Ala164_Tyr165insTer")?;

    let level = eq.equivalent_level(&v1, &v2)?;
    assert_eq!(level, EquivalenceLevel::Analogous);
    Ok(())
}

#[test]
fn test_analogous_repeat_equivalence() -> Result<(), hgvs_weaver::error::HgvsError> {
    use hgvs_weaver::data::TranscriptData;
    use hgvs_weaver::equivalence::VariantEquivalence;

    use hgvs_weaver::{parse_hgvs_variant, DataProvider, IdentifierKind};

    struct MockRepeatEqProvider;
    impl DataProvider for MockRepeatEqProvider {
        fn get_seq(
            &self,
            _ac: &str,
            _start: i32,
            _end: Option<i32>,
            _kind: hgvs_weaver::data::IdentifierType,
        ) -> Result<String, hgvs_weaver::error::HgvsError> {
            Ok("X".repeat(489) + "PRS" + &"X".repeat(100))
        }
        fn get_transcript(
            &self,
            _ac: &str,
            _gene: Option<&str>,
        ) -> Result<TranscriptData, hgvs_weaver::error::HgvsError> {
            panic!("Not implemented")
        }
        fn get_identifier_type(
            &self,
            _ac: &str,
        ) -> Result<hgvs_weaver::data::IdentifierType, hgvs_weaver::error::HgvsError> {
            Ok(hgvs_weaver::data::IdentifierType::ProteinAccession)
        }
        fn get_symbol_accessions(
            &self,
            _symbol: &str,
            _source: IdentifierKind,
            _target: IdentifierKind,
        ) -> Result<Vec<(hgvs_weaver::data::IdentifierType, String)>, hgvs_weaver::error::HgvsError>
        {
            Ok(vec![])
        }
    }

    struct MockSearcher;
    impl hgvs_weaver::data::TranscriptSearch for MockSearcher {
        fn get_transcripts_for_region(
            &self,
            _chrom: &str,
            _start: i32,
            _end: i32,
        ) -> Result<Vec<String>, hgvs_weaver::error::HgvsError> {
            Ok(vec![])
        }
    }

    let provider = MockRepeatEqProvider;
    let searcher = MockSearcher;
    let eq_mapper = VariantMapper::new(&provider);
    let eq = VariantEquivalence::new(&eq_mapper, &searcher);

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
fn test_analogous_fs_wildcard_unification() -> Result<(), hgvs_weaver::error::HgvsError> {
    use hgvs_weaver::data::TranscriptData;
    use hgvs_weaver::equivalence::VariantEquivalence;

    use hgvs_weaver::{parse_hgvs_variant, DataProvider, IdentifierKind};

    struct MockWildcardProvider;
    impl DataProvider for MockWildcardProvider {
        fn get_seq(
            &self,
            _ac: &str,
            _start: i32,
            _end: Option<i32>,
            _kind: hgvs_weaver::data::IdentifierType,
        ) -> Result<String, hgvs_weaver::error::HgvsError> {
            Ok("X".repeat(96) + "R" + &"X".repeat(100))
        }
        fn get_transcript(
            &self,
            _ac: &str,
            _gene: Option<&str>,
        ) -> Result<TranscriptData, hgvs_weaver::error::HgvsError> {
            panic!("Not implemented")
        }
        fn get_identifier_type(
            &self,
            _ac: &str,
        ) -> Result<hgvs_weaver::data::IdentifierType, hgvs_weaver::error::HgvsError> {
            Ok(hgvs_weaver::data::IdentifierType::ProteinAccession)
        }
        fn get_symbol_accessions(
            &self,
            _symbol: &str,
            _source: IdentifierKind,
            _target: IdentifierKind,
        ) -> Result<Vec<(hgvs_weaver::data::IdentifierType, String)>, hgvs_weaver::error::HgvsError>
        {
            Ok(vec![])
        }
    }

    struct MockSearcher;
    impl hgvs_weaver::data::TranscriptSearch for MockSearcher {
        fn get_transcripts_for_region(
            &self,
            _chrom: &str,
            _start: i32,
            _end: i32,
        ) -> Result<Vec<String>, hgvs_weaver::error::HgvsError> {
            Ok(vec![])
        }
    }

    let provider = MockWildcardProvider;
    let searcher = MockSearcher;
    let eq_mapper = VariantMapper::new(&provider);
    let eq = VariantEquivalence::new(&eq_mapper, &searcher);

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
fn test_multi_unit_repeat_equivalence() -> Result<(), hgvs_weaver::error::HgvsError> {
    use hgvs_weaver::data::TranscriptData;
    use hgvs_weaver::equivalence::VariantEquivalence;

    use hgvs_weaver::{parse_hgvs_variant, DataProvider, IdentifierKind};

    struct MockMultiRepeatProvider;
    impl DataProvider for MockMultiRepeatProvider {
        fn get_seq(
            &self,
            ac: &str,
            _start: i32,
            _end: Option<i32>,
            _kind: hgvs_weaver::data::IdentifierType,
        ) -> Result<String, hgvs_weaver::error::HgvsError> {
            if ac == "NP_001122316.1" {
                // ...pppsvsatg pgpgpgpgpg pgpgpappny s...
                // Residue 229 starts the GP repeat.
                // 228 X's + 8 GP units (16 chars) + 100 X's
                Ok("X".repeat(228) + "GPGPGPGPGPGPGPGP" + &"X".repeat(100))
            } else if ac == "NP_000067.1" {
                // ...lap apapapap apapvaapap apapapapap apapapdaap...
                // Residue 179 starts the AP repeat.
                // 178 X's + 8 AP units + 100 X's: AP[5] then deletes three units.
                Ok("X".repeat(178) + "APAPAPAPAPAPAPAP" + &"X".repeat(100))
            } else {
                Ok("".to_string())
            }
        }
        fn get_transcript(
            &self,
            _ac: &str,
            _gene: Option<&str>,
        ) -> Result<TranscriptData, hgvs_weaver::error::HgvsError> {
            panic!("Not implemented")
        }
        fn get_identifier_type(
            &self,
            _ac: &str,
        ) -> Result<hgvs_weaver::data::IdentifierType, hgvs_weaver::error::HgvsError> {
            Ok(hgvs_weaver::data::IdentifierType::ProteinAccession)
        }
        fn get_symbol_accessions(
            &self,
            _symbol: &str,
            _source: IdentifierKind,
            _target: IdentifierKind,
        ) -> Result<Vec<(hgvs_weaver::data::IdentifierType, String)>, hgvs_weaver::error::HgvsError>
        {
            Ok(vec![])
        }
    }

    struct MockSearcher;
    impl hgvs_weaver::data::TranscriptSearch for MockSearcher {
        fn get_transcripts_for_region(
            &self,
            _chrom: &str,
            _start: i32,
            _end: i32,
        ) -> Result<Vec<String>, hgvs_weaver::error::HgvsError> {
            Ok(vec![])
        }
    }

    let provider = MockMultiRepeatProvider;
    let searcher = MockSearcher;
    let eq_mapper = VariantMapper::new(&provider);
    let eq = VariantEquivalence::new(&eq_mapper, &searcher);

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
fn test_immediate_stop_normalization() -> Result<(), hgvs_weaver::error::HgvsError> {
    use hgvs_weaver::data::{DataProvider, IdentifierKind, IdentifierType, TranscriptData};
    use hgvs_weaver::mapper::VariantMapper;
    use hgvs_weaver::structs::TranscriptPos;

    struct MockFsProvider;
    impl DataProvider for MockFsProvider {
        fn get_transcript(
            &self,
            _ac: &str,
            _ref: Option<&str>,
        ) -> Result<TranscriptData, hgvs_weaver::error::HgvsError> {
            Ok(TranscriptData {
                ac: "NM_1.1".to_string(),
                gene: "TEST".to_string(),
                cds_start_index: Some(TranscriptPos(0)),
                cds_end_index: Some(TranscriptPos(6)),
                strand: hgvs_weaver::data::Strand::Plus,
                reference_accession: "NC_1.1".to_string(),
                exons: vec![],
            })
        }
        fn get_seq(
            &self,
            _ac: &str,
            _start: i32,
            _end: Option<i32>,
            _kind: IdentifierType,
        ) -> Result<String, hgvs_weaver::error::HgvsError> {
            // AAA (Lys) GGG (Gly)
            Ok("AAAGGG".to_string())
        }
        fn get_symbol_accessions(
            &self,
            _symbol: &str,
            _source: IdentifierKind,
            _target: IdentifierKind,
        ) -> Result<Vec<(IdentifierType, String)>, hgvs_weaver::error::HgvsError> {
            Ok(vec![(
                IdentifierType::ProteinAccession,
                "NP_1.1".to_string(),
            )])
        }
        fn get_identifier_type(
            &self,
            _ac: &str,
        ) -> Result<IdentifierType, hgvs_weaver::error::HgvsError> {
            Ok(IdentifierType::TranscriptAccession)
        }
    }

    let hdp = MockFsProvider;
    let mapper = VariantMapper::new(&hdp);

    // c.1_2delinsTA -> TAG (Stop) instead of ATG (Met)
    let var_c = hgvs_weaver::parse_hgvs_variant("NM_1.1:c.1_2delinsTA")?;
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
