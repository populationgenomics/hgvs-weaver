use hgvs_weaver::coords::{GenomicPos, TranscriptPos};
use hgvs_weaver::data::{
    DataProvider, ExonData, IdentifierKind, IdentifierType, TranscriptData, TranscriptSearch,
};
use hgvs_weaver::equivalence::{EquivalenceLevel, VariantEquivalence};
use hgvs_weaver::error::HgvsError;
use hgvs_weaver::parse_hgvs_variant;

struct MockDataProvider;
impl DataProvider for MockDataProvider {
    fn get_transcript(&self, ac: &str, _: Option<&str>) -> Result<TranscriptData, HgvsError> {
        if ac == "NM_001166478.1" || ac == "NM_005813.3" {
            // One minus-strand exon: transcript index i is genomic index 4000 - i.
            Ok(TranscriptData {
                ac: ac.to_string(),
                gene: "TEST".to_string(),
                cds_start_index: Some(TranscriptPos(0)),
                cds_end_index: Some(TranscriptPos(3000)),
                strand: hgvs_weaver::data::Strand::Minus,
                reference_accession: "NC_000001.1".to_string(),
                exons: vec![ExonData {
                    transcript_start: TranscriptPos(0),
                    transcript_end: TranscriptPos(3001),
                    reference_start: GenomicPos(1000),
                    reference_end: GenomicPos(4000),
                    alt_strand: hgvs_weaver::data::Strand::Minus,
                    cigar: "3001M".to_string(),
                }],
            })
        } else if ac == "NM_BRAF" {
            Ok(TranscriptData {
                ac: ac.to_string(),
                gene: "BRAF".to_string(),
                cds_start_index: Some(TranscriptPos(0)),
                cds_end_index: Some(TranscriptPos(3000)),
                strand: hgvs_weaver::data::Strand::Plus, // Plus strand
                reference_accession: "NC_BRAF".to_string(),
                exons: vec![],
            })
        } else {
            Err(HgvsError::ValidationError("Not found".into()))
        }
    }
    fn get_seq(
        &self,
        ac: &str,
        s: i32,
        e: Option<i32>,
        _k: IdentifierType,
    ) -> Result<String, HgvsError> {
        // 4000 N's with the few bases the cases below depend on.
        let mut seq = vec![b'N'; 4000];
        if ac == "NM_BRAF" || ac == "NC_BRAF" {
            // BRAF Val600 is GTG at transcript indices 1797..=1799.
            seq[1797] = b'G';
            seq[1798] = b'T';
            seq[1799] = b'G';
        } else {
            // Case 9: c.35 is transcript index 34, genomic 3966 (A). The base
            // after it must differ so the insertion cannot shift further.
            seq[3966] = b'A';
            seq[34] = b'A';
            // Case 15: c.2673 is transcript index 2672, genomic 1328 (T).
            seq[1328] = b'T';
            seq[2672] = b'T';
        }
        let start = (s.max(0) as usize).min(seq.len());
        let end = e.map_or(seq.len(), |e| (e as usize).min(seq.len()));
        Ok(String::from_utf8(seq[start..end.max(start)].to_vec()).unwrap())
    }
    fn get_symbol_accessions(
        &self,
        _s: &str,
        _f: IdentifierKind,
        _t: IdentifierKind,
    ) -> Result<Vec<(IdentifierType, String)>, HgvsError> {
        Ok(vec![])
    }
    fn get_identifier_type(&self, _id: &str) -> Result<IdentifierType, HgvsError> {
        Ok(IdentifierType::TranscriptAccession)
    }
}

struct MockSearch;
impl TranscriptSearch for MockSearch {
    fn get_transcripts_for_region(
        &self,
        _: &str,
        _: i32,
        _: i32,
    ) -> Result<Vec<String>, HgvsError> {
        Ok(vec![])
    }
}

#[test]
fn test_repro_case9() -> Result<(), HgvsError> {
    let hdp = MockDataProvider;
    let search = MockSearch;
    let eq = VariantEquivalence::new(&hdp, &search);

    let v1 = parse_hgvs_variant("NM_001166478.1:c.35_36insT")?;
    let v2 = parse_hgvs_variant("NM_001166478.1:c.35dup")?;

    assert_eq!(eq.equivalent_level(&v1, &v2)?, EquivalenceLevel::Analogous);
    Ok(())
}

#[test]
fn test_repro_case15() -> Result<(), HgvsError> {
    let hdp = MockDataProvider;
    let search = MockSearch;
    let eq = VariantEquivalence::new(&hdp, &search);

    let v1 = parse_hgvs_variant("NM_005813.3:c.2673insA")?;
    let v2 = parse_hgvs_variant("NM_005813.3:c.2673dup")?;

    assert_eq!(eq.equivalent_level(&v1, &v2)?, EquivalenceLevel::Analogous);
    Ok(())
}

#[test]
fn test_braf_identity() -> Result<(), HgvsError> {
    let hdp = MockDataProvider;
    let search = MockSearch;
    let eq = VariantEquivalence::new(&hdp, &search);

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
