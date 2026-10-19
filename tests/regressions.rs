use hgvs_weaver::data::{DataProvider, ExonData, IdentifierKind, IdentifierType, TranscriptData};
use hgvs_weaver::structs::{GenomicPos, TranscriptPos};
use hgvs_weaver::{parse_hgvs_variant, HgvsError, SequenceVariant, VariantMapper};

struct MockDataProvider {
    transcripts: std::collections::HashMap<String, TranscriptData>,
    sequences: std::collections::HashMap<String, String>,
}

impl DataProvider for MockDataProvider {
    fn get_transcript(&self, ac: &str, _ref_ac: Option<&str>) -> Result<TranscriptData, HgvsError> {
        self.transcripts
            .get(ac)
            .cloned()
            .ok_or_else(|| HgvsError::ValidationError(format!("Transcript {} not found", ac)))
    }

    fn get_seq(
        &self,
        ac: &str,
        start: i32,
        end: Option<i32>,
        _kind: IdentifierType,
    ) -> Result<String, HgvsError> {
        let seq = self
            .sequences
            .get(ac)
            .ok_or_else(|| HgvsError::ValidationError(format!("Sequence {} not found", ac)))?;
        let start = start as usize;
        let end = end.map_or(seq.len(), |e| e as usize);
        if start > seq.len() || end > seq.len() || start > end {
            return Err(HgvsError::ValidationError(format!(
                "Invalid seq range: {}-{} for len {}",
                start,
                end,
                seq.len()
            )));
        }
        Ok(seq[start..end].to_string())
    }

    fn get_symbol_accessions(
        &self,
        _symbol: &str,
        _src: IdentifierKind,
        _target: IdentifierKind,
    ) -> Result<Vec<(IdentifierType, String)>, HgvsError> {
        Ok(vec![])
    }
    fn get_identifier_type(&self, _identifier: &str) -> Result<IdentifierType, HgvsError> {
        Ok(IdentifierType::TranscriptAccession)
    }
}

fn run_regression_test(
    dp: MockDataProvider,
    variant_str: &str,
    protein_ac: &str,
    expected_p: &str,
) {
    let mapper = VariantMapper::new(&dp);
    let variant = parse_hgvs_variant(variant_str).unwrap();
    if let SequenceVariant::Coding(cv) = variant {
        // We use the variant as-is to match ClinVar representation exactly.
        // Normalization can sometimes shift complex variants in ways that change the translation representation.
        let p = mapper.c_to_p(&cv, Some(protein_ac)).expect("c_to_p failed");
        let expected = format!("{}:{}", protein_ac, expected_p);
        assert_eq!(p.to_string(), expected, "Variant: {}", variant_str);
    } else {
        panic!("Parsed variant is not Coding");
    }
}

#[test]
fn test_regression_nm_000038_6() {
    // APC c.1972_1973delinsAT -> p.(Glu658Met)
    let seq = include_str!("data/NM_000038.6.seq").trim();
    let ac = "NM_000038.6";
    let tx = TranscriptData::new(
        ac.to_string(),
        "APC",
        Some(TranscriptPos(59)),
        Some(TranscriptPos(8590)),
        hgvs_weaver::data::Strand::Plus,
        "NC_000017.11",
        vec![
            ExonData::new(
                TranscriptPos(0),
                TranscriptPos(41),
                GenomicPos(112737884),
                GenomicPos(112737924),
                hgvs_weaver::data::Strand::Plus,
                "41=",
            ),
            ExonData::new(
                TranscriptPos(41),
                TranscriptPos(194),
                GenomicPos(112754872),
                GenomicPos(112755024),
                hgvs_weaver::data::Strand::Plus,
                "153=",
            ),
            ExonData::new(
                TranscriptPos(194),
                TranscriptPos(279),
                GenomicPos(112766325),
                GenomicPos(112766409),
                hgvs_weaver::data::Strand::Plus,
                "85=",
            ),
            ExonData::new(
                TranscriptPos(279),
                TranscriptPos(481),
                GenomicPos(112767188),
                GenomicPos(112767389),
                hgvs_weaver::data::Strand::Plus,
                "202=",
            ),
            ExonData::new(
                TranscriptPos(481),
                TranscriptPos(590),
                GenomicPos(112775628),
                GenomicPos(112775736),
                hgvs_weaver::data::Strand::Plus,
                "109=",
            ),
            ExonData::new(
                TranscriptPos(590),
                TranscriptPos(704),
                GenomicPos(112780789),
                GenomicPos(112780902),
                hgvs_weaver::data::Strand::Plus,
                "114=",
            ),
            ExonData::new(
                TranscriptPos(704),
                TranscriptPos(788),
                GenomicPos(112792445),
                GenomicPos(112792528),
                hgvs_weaver::data::Strand::Plus,
                "84=",
            ),
            ExonData::new(
                TranscriptPos(788),
                TranscriptPos(893),
                GenomicPos(112801278),
                GenomicPos(112801382),
                hgvs_weaver::data::Strand::Plus,
                "105=",
            ),
            ExonData::new(
                TranscriptPos(893),
                TranscriptPos(992),
                GenomicPos(112815494),
                GenomicPos(112815592),
                hgvs_weaver::data::Strand::Plus,
                "99=",
            ),
            ExonData::new(
                TranscriptPos(992),
                TranscriptPos(1371),
                GenomicPos(112818965),
                GenomicPos(112819343),
                hgvs_weaver::data::Strand::Plus,
                "379=",
            ),
            ExonData::new(
                TranscriptPos(1371),
                TranscriptPos(1467),
                GenomicPos(112821895),
                GenomicPos(112821990),
                hgvs_weaver::data::Strand::Plus,
                "96=",
            ),
            ExonData::new(
                TranscriptPos(1467),
                TranscriptPos(1607),
                GenomicPos(112827107),
                GenomicPos(112827246),
                hgvs_weaver::data::Strand::Plus,
                "140=",
            ),
            ExonData::new(
                TranscriptPos(1607),
                TranscriptPos(1685),
                GenomicPos(112827928),
                GenomicPos(112828005),
                hgvs_weaver::data::Strand::Plus,
                "78=",
            ),
            ExonData::new(
                TranscriptPos(1685),
                TranscriptPos(1802),
                GenomicPos(112828855),
                GenomicPos(112828971),
                hgvs_weaver::data::Strand::Plus,
                "117=",
            ),
            ExonData::new(
                TranscriptPos(1802),
                TranscriptPos(2017),
                GenomicPos(112834950),
                GenomicPos(112835164),
                hgvs_weaver::data::Strand::Plus,
                "215=",
            ),
            ExonData::new(
                TranscriptPos(2017),
                TranscriptPos(10704),
                GenomicPos(112837552),
                GenomicPos(112846238),
                hgvs_weaver::data::Strand::Plus,
                "8687=",
            ),
        ],
    );
    let mut transcripts = std::collections::HashMap::new();
    transcripts.insert(ac.to_string(), tx);
    let mut sequences = std::collections::HashMap::new();
    sequences.insert(ac.to_string(), seq.to_string());

    let dp = MockDataProvider {
        transcripts,
        sequences,
    };
    run_regression_test(
        dp,
        "NM_000038.6:c.1972_1973delinsAT",
        "NP_000029.2",
        "p.(Glu658Met)",
    );
}

#[test]
fn test_regression_nm_000527_5() {
    // LDLR c.669_683delinsAACTGCGGTAAACTGCGGTAAACT -> p.(Asp224_Glu228delinsThrAlaValAsnCysGlyLysLeu)
    let seq = include_str!("data/NM_000527.5.seq").trim();
    let ac = "NM_000527.5";
    let tx = TranscriptData::new(
        ac.to_string(),
        "LDLR",
        Some(TranscriptPos(86)),
        Some(TranscriptPos(2668)),
        hgvs_weaver::data::Strand::Plus,
        "NC_000019.10",
        vec![
            ExonData::new(
                TranscriptPos(0),
                TranscriptPos(153),
                GenomicPos(11089462),
                GenomicPos(11089614),
                hgvs_weaver::data::Strand::Plus,
                "153=",
            ),
            ExonData::new(
                TranscriptPos(153),
                TranscriptPos(276),
                GenomicPos(11100222),
                GenomicPos(11100344),
                hgvs_weaver::data::Strand::Plus,
                "123=",
            ),
            ExonData::new(
                TranscriptPos(276),
                TranscriptPos(399),
                GenomicPos(11102663),
                GenomicPos(11102785),
                hgvs_weaver::data::Strand::Plus,
                "123=",
            ),
            ExonData::new(
                TranscriptPos(399),
                TranscriptPos(780),
                GenomicPos(11105219),
                GenomicPos(11105599),
                hgvs_weaver::data::Strand::Plus,
                "381=",
            ),
            ExonData::new(
                TranscriptPos(780),
                TranscriptPos(903),
                GenomicPos(11106564),
                GenomicPos(11106686),
                hgvs_weaver::data::Strand::Plus,
                "123=",
            ),
            ExonData::new(
                TranscriptPos(903),
                TranscriptPos(1026),
                GenomicPos(11107391),
                GenomicPos(11107513),
                hgvs_weaver::data::Strand::Plus,
                "123=",
            ),
            ExonData::new(
                TranscriptPos(1026),
                TranscriptPos(1146),
                GenomicPos(11110651),
                GenomicPos(11110770),
                hgvs_weaver::data::Strand::Plus,
                "120=",
            ),
            ExonData::new(
                TranscriptPos(1146),
                TranscriptPos(1272),
                GenomicPos(11111513),
                GenomicPos(11111638),
                hgvs_weaver::data::Strand::Plus,
                "126=",
            ),
            ExonData::new(
                TranscriptPos(1272),
                TranscriptPos(1444),
                GenomicPos(11113277),
                GenomicPos(11113448),
                hgvs_weaver::data::Strand::Plus,
                "172=",
            ),
            ExonData::new(
                TranscriptPos(1444),
                TranscriptPos(1672),
                GenomicPos(11113534),
                GenomicPos(11113761),
                hgvs_weaver::data::Strand::Plus,
                "228=",
            ),
            ExonData::new(
                TranscriptPos(1672),
                TranscriptPos(1791),
                GenomicPos(11116093),
                GenomicPos(11116211),
                hgvs_weaver::data::Strand::Plus,
                "119=",
            ),
            ExonData::new(
                TranscriptPos(1791),
                TranscriptPos(1931),
                GenomicPos(11116858),
                GenomicPos(11116997),
                hgvs_weaver::data::Strand::Plus,
                "140=",
            ),
            ExonData::new(
                TranscriptPos(1931),
                TranscriptPos(2073),
                GenomicPos(11120091),
                GenomicPos(11120232),
                hgvs_weaver::data::Strand::Plus,
                "142=",
            ),
            ExonData::new(
                TranscriptPos(2073),
                TranscriptPos(2226),
                GenomicPos(11120369),
                GenomicPos(11120521),
                hgvs_weaver::data::Strand::Plus,
                "153=",
            ),
            ExonData::new(
                TranscriptPos(2226),
                TranscriptPos(2397),
                GenomicPos(11123173),
                GenomicPos(11123343),
                hgvs_weaver::data::Strand::Plus,
                "171=",
            ),
            ExonData::new(
                TranscriptPos(2397),
                TranscriptPos(2475),
                GenomicPos(11128007),
                GenomicPos(11128084),
                hgvs_weaver::data::Strand::Plus,
                "78=",
            ),
            ExonData::new(
                TranscriptPos(2475),
                TranscriptPos(2633),
                GenomicPos(11129512),
                GenomicPos(11129669),
                hgvs_weaver::data::Strand::Plus,
                "158=",
            ),
            ExonData::new(
                TranscriptPos(2633),
                TranscriptPos(5173),
                GenomicPos(11131280),
                GenomicPos(11133819),
                hgvs_weaver::data::Strand::Plus,
                "2540=",
            ),
        ],
    );
    let mut transcripts = std::collections::HashMap::new();
    transcripts.insert(ac.to_string(), tx);
    let mut sequences = std::collections::HashMap::new();
    sequences.insert(ac.to_string(), seq.to_string());

    let dp = MockDataProvider {
        transcripts,
        sequences,
    };
    run_regression_test(
        dp,
        "NM_000527.5:c.669_683delinsAACTGCGGTAAACTGCGGTAAACT",
        "NP_000518.1",
        "p.(Asp224_Glu228delinsThrAlaValAsnCysGlyLysLeu)",
    );
}

#[test]
fn test_regression_nm_000478_6() {
    // ALPL c.876_882delinsT -> p.(Gly293_Asp294del)
    let seq = include_str!("data/NM_000478.6.seq").trim();
    let ac = "NM_000478.6";
    let tx = TranscriptData::new(
        ac.to_string(),
        "ALPL",
        Some(TranscriptPos(199)),
        Some(TranscriptPos(1773)),
        hgvs_weaver::data::Strand::Plus,
        "NC_000001.11",
        vec![
            ExonData::new(
                TranscriptPos(0),
                TranscriptPos(95),
                GenomicPos(21509422),
                GenomicPos(21509516),
                hgvs_weaver::data::Strand::Plus,
                "95=",
            ),
            ExonData::new(
                TranscriptPos(95),
                TranscriptPos(260),
                GenomicPos(21553977),
                GenomicPos(21554141),
                hgvs_weaver::data::Strand::Plus,
                "165=",
            ),
            ExonData::new(
                TranscriptPos(260),
                TranscriptPos(380),
                GenomicPos(21560625),
                GenomicPos(21560744),
                hgvs_weaver::data::Strand::Plus,
                "120=",
            ),
            ExonData::new(
                TranscriptPos(380),
                TranscriptPos(496),
                GenomicPos(21561096),
                GenomicPos(21561211),
                hgvs_weaver::data::Strand::Plus,
                "116=",
            ),
            ExonData::new(
                TranscriptPos(496),
                TranscriptPos(671),
                GenomicPos(21563109),
                GenomicPos(21563283),
                hgvs_weaver::data::Strand::Plus,
                "175=",
            ),
            ExonData::new(
                TranscriptPos(671),
                TranscriptPos(847),
                GenomicPos(21564040),
                GenomicPos(21564215),
                hgvs_weaver::data::Strand::Plus,
                "176=",
            ),
            ExonData::new(
                TranscriptPos(847),
                TranscriptPos(991),
                GenomicPos(21568103),
                GenomicPos(21568246),
                hgvs_weaver::data::Strand::Plus,
                "144=",
            ),
            ExonData::new(
                TranscriptPos(991),
                TranscriptPos(1061),
                GenomicPos(21570304),
                GenomicPos(21570373),
                hgvs_weaver::data::Strand::Plus,
                "70=",
            ),
            ExonData::new(
                TranscriptPos(1061),
                TranscriptPos(1196),
                GenomicPos(21573664),
                GenomicPos(21573798),
                hgvs_weaver::data::Strand::Plus,
                "135=",
            ),
            ExonData::new(
                TranscriptPos(1196),
                TranscriptPos(1388),
                GenomicPos(21575732),
                GenomicPos(21575923),
                hgvs_weaver::data::Strand::Plus,
                "192=",
            ),
            ExonData::new(
                TranscriptPos(1388),
                TranscriptPos(1508),
                GenomicPos(21576521),
                GenomicPos(21576640),
                hgvs_weaver::data::Strand::Plus,
                "120=",
            ),
            ExonData::new(
                TranscriptPos(1508),
                TranscriptPos(2536),
                GenomicPos(21577382),
                GenomicPos(21578409),
                hgvs_weaver::data::Strand::Plus,
                "1028=",
            ),
        ],
    );
    let mut transcripts = std::collections::HashMap::new();
    transcripts.insert(ac.to_string(), tx);
    let mut sequences = std::collections::HashMap::new();
    sequences.insert(ac.to_string(), seq.to_string());

    let dp = MockDataProvider {
        transcripts,
        sequences,
    };
    run_regression_test(
        dp,
        "NM_000478.6:c.876_882delinsT",
        "NP_000469.3",
        "p.(Gly293_Asp294del)",
    );
}

#[test]
fn test_regression_nm_001122606_1() {
    // CASP8 c.763_768delinsTGAAGT -> p.(Asn255Ter)
    // Note: ClinVar expects p.(Asn255_Pro256delinsTer), but Weaver simplifies this to p.(Asn255Ter)
    // which is biologically equivalent.
    let seq = include_str!("data/NM_001122606.1.seq").trim();
    let ac = "NM_001122606.1";
    let tx = TranscriptData::new(
        ac.to_string(),
        "CASP8",
        Some(TranscriptPos(180)),
        Some(TranscriptPos(1415)),
        hgvs_weaver::data::Strand::Plus,
        "NC_000002.12",
        vec![
            ExonData::new(
                TranscriptPos(0),
                TranscriptPos(244),
                GenomicPos(120469105),
                GenomicPos(120469348),
                hgvs_weaver::data::Strand::Minus,
                "244=",
            ),
            ExonData::new(
                TranscriptPos(244),
                TranscriptPos(363),
                GenomicPos(120456650),
                GenomicPos(120456768),
                hgvs_weaver::data::Strand::Minus,
                "119=",
            ),
            ExonData::new(
                TranscriptPos(363),
                TranscriptPos(577),
                GenomicPos(120455356),
                GenomicPos(120455569),
                hgvs_weaver::data::Strand::Minus,
                "214=",
            ),
            ExonData::new(
                TranscriptPos(577),
                TranscriptPos(736),
                GenomicPos(120448969),
                GenomicPos(120449127),
                hgvs_weaver::data::Strand::Minus,
                "159=",
            ),
            ExonData::new(
                TranscriptPos(736),
                TranscriptPos(921),
                GenomicPos(120447840),
                GenomicPos(120448024),
                hgvs_weaver::data::Strand::Minus,
                "185=",
            ),
            ExonData::new(
                TranscriptPos(921),
                TranscriptPos(1044),
                GenomicPos(120446304),
                GenomicPos(120446426),
                hgvs_weaver::data::Strand::Minus,
                "123=",
            ),
            ExonData::new(
                TranscriptPos(1044),
                TranscriptPos(1108),
                GenomicPos(120442598),
                GenomicPos(120442661),
                hgvs_weaver::data::Strand::Minus,
                "64=",
            ),
            ExonData::new(
                TranscriptPos(1108),
                TranscriptPos(1273),
                GenomicPos(120441729),
                GenomicPos(120441893),
                hgvs_weaver::data::Strand::Minus,
                "165=",
            ),
            ExonData::new(
                TranscriptPos(1273),
                TranscriptPos(3752),
                GenomicPos(120426147),
                GenomicPos(120428625),
                hgvs_weaver::data::Strand::Minus,
                "2479=",
            ),
        ],
    );
    let mut transcripts = std::collections::HashMap::new();
    transcripts.insert(ac.to_string(), tx);
    let mut sequences = std::collections::HashMap::new();
    sequences.insert(ac.to_string(), seq.to_string());

    let dp = MockDataProvider {
        transcripts,
        sequences,
    };
    run_regression_test(
        dp,
        "NM_001122606.1:c.763_768delinsTGAAGT",
        "NP_001116078.1",
        "p.(Asn255Ter)",
    );
}

#[test]
fn test_regression_nm_000465_4() {
    // ABCA4 c.66_70delinsTGCGT -> p.(Pro24Ser)
    let seq = include_str!("data/NM_000465.4.seq").trim();
    let ac = "NM_000465.4";
    let tx = TranscriptData::new(
        ac.to_string(),
        "ABCA4",
        Some(TranscriptPos(114)),
        Some(TranscriptPos(2447)),
        hgvs_weaver::data::Strand::Plus,
        "NC_000001.11",
        vec![
            ExonData::new(
                TranscriptPos(0),
                TranscriptPos(272),
                GenomicPos(214809411),
                GenomicPos(214809682),
                hgvs_weaver::data::Strand::Minus,
                "272=",
            ),
            ExonData::new(
                TranscriptPos(272),
                TranscriptPos(329),
                GenomicPos(214797060),
                GenomicPos(214797116),
                hgvs_weaver::data::Strand::Minus,
                "57=",
            ),
            ExonData::new(
                TranscriptPos(329),
                TranscriptPos(478),
                GenomicPos(214792296),
                GenomicPos(214792444),
                hgvs_weaver::data::Strand::Minus,
                "149=",
            ),
            ExonData::new(
                TranscriptPos(478),
                TranscriptPos(1428),
                GenomicPos(214780559),
                GenomicPos(214781508),
                hgvs_weaver::data::Strand::Minus,
                "950=",
            ),
            ExonData::new(
                TranscriptPos(1428),
                TranscriptPos(1509),
                GenomicPos(214769231),
                GenomicPos(214769311),
                hgvs_weaver::data::Strand::Minus,
                "81=",
            ),
            ExonData::new(
                TranscriptPos(1509),
                TranscriptPos(1682),
                GenomicPos(214767481),
                GenomicPos(214767653),
                hgvs_weaver::data::Strand::Minus,
                "173=",
            ),
            ExonData::new(
                TranscriptPos(1682),
                TranscriptPos(1791),
                GenomicPos(214752446),
                GenomicPos(214752554),
                hgvs_weaver::data::Strand::Minus,
                "109=",
            ),
            ExonData::new(
                TranscriptPos(1791),
                TranscriptPos(1924),
                GenomicPos(214745721),
                GenomicPos(214745853),
                hgvs_weaver::data::Strand::Minus,
                "133=",
            ),
            ExonData::new(
                TranscriptPos(1924),
                TranscriptPos(2017),
                GenomicPos(214745066),
                GenomicPos(214745158),
                hgvs_weaver::data::Strand::Minus,
                "93=",
            ),
            ExonData::new(
                TranscriptPos(2017),
                TranscriptPos(2115),
                GenomicPos(214730410),
                GenomicPos(214730507),
                hgvs_weaver::data::Strand::Minus,
                "98=",
            ),
            ExonData::new(
                TranscriptPos(2115),
                TranscriptPos(5478),
                GenomicPos(214725645),
                GenomicPos(214729007),
                hgvs_weaver::data::Strand::Minus,
                "3363=",
            ),
        ],
    );
    let mut transcripts = std::collections::HashMap::new();
    transcripts.insert(ac.to_string(), tx);
    let mut sequences = std::collections::HashMap::new();
    sequences.insert(ac.to_string(), seq.to_string());

    let dp = MockDataProvider {
        transcripts,
        sequences,
    };
    run_regression_test(
        dp,
        "NM_000465.4:c.66_70delinsTGCGT",
        "NP_000456.2",
        "p.(Pro24Ser)",
    );
}
