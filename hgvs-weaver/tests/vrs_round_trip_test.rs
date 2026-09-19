//! VRS and SPDI read back: `from_vrs` and `from_spdi` invert `to_vrs` and
//! `to_spdi`, on nucleotide and protein sequences.

use hgvs_weaver::data::{DataProvider, IdentifierKind, IdentifierType, TranscriptData};
use hgvs_weaver::error::HgvsError;
use hgvs_weaver::mapper::VariantMapper;
use hgvs_weaver::parse_hgvs_variant;
use hgvs_weaver::refget::Refget;
use hgvs_weaver::vrs::refget_accession;

const GENOME: &str = "ACGTTTGCAAGGCTAGCTAGCTTTTAACGGGATCGATCGA";
const OTHER: &str = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
const PROTEIN: &str = "MKLAAAYRQ";

/// Serves three sequences and can look an accession up by refget digest.
struct Provider {
    lookup: bool,
}

fn sequence(ac: &str) -> Result<&'static str, HgvsError> {
    match ac {
        "NC_TEST.1" => Ok(GENOME),
        "NC_OTHER.1" => Ok(OTHER),
        "NP_TEST.1" => Ok(PROTEIN),
        _ => Err(HgvsError::DataProviderError(format!("no sequence {ac}"))),
    }
}

impl DataProvider for Provider {
    fn get_transcript(&self, ac: &str, _: Option<&str>) -> Result<TranscriptData, HgvsError> {
        Err(HgvsError::DataProviderError(format!("no transcript {ac}")))
    }

    fn get_seq(
        &self,
        ac: &str,
        start: i32,
        end: Option<i32>,
        _kind: IdentifierType,
    ) -> Result<String, HgvsError> {
        let seq = sequence(ac)?;
        let start = (start.max(0) as usize).min(seq.len());
        let end = end.map_or(seq.len(), |e| (e.max(0) as usize).min(seq.len()));
        Ok(seq[start..end.max(start)].to_string())
    }

    fn get_symbol_accessions(
        &self,
        _: &str,
        _: IdentifierKind,
        _: IdentifierKind,
    ) -> Result<Vec<(IdentifierType, String)>, HgvsError> {
        Ok(vec![])
    }

    fn get_identifier_type(&self, id: &str) -> Result<IdentifierType, HgvsError> {
        Ok(if id.starts_with("NP_") {
            IdentifierType::ProteinAccession
        } else {
            IdentifierType::GenomicAccession
        })
    }
}

impl Refget for Provider {
    fn refget_accession(&self, _ac: &str) -> Result<Option<String>, HgvsError> {
        Ok(None) // let the store compute it
    }

    fn accession_for_refget(&self, refget: &str) -> Result<Option<String>, HgvsError> {
        if !self.lookup {
            return Ok(None);
        }
        Ok(["NC_TEST.1", "NC_OTHER.1", "NP_TEST.1"]
            .into_iter()
            .find(|ac| refget_accession(sequence(ac).unwrap()) == refget)
            .map(str::to_string))
    }
}

/// `hgvs` through VRS and back, and through SPDI and back; both must give an
/// allele with the same identifier and the expected HGVS spelling.
fn round_trip(mapper: &VariantMapper, hgvs: &str, expected: &str) {
    let var = parse_hgvs_variant(hgvs).unwrap();
    let vrs = mapper.to_vrs(&var).unwrap();
    let back = mapper.from_vrs(&vrs.to_json(), None).unwrap();
    assert_eq!(back.to_string(), expected, "{hgvs} via VRS");
    assert_eq!(mapper.to_vrs(&back).unwrap().id, vrs.id, "{hgvs} via VRS");
    for spdi in [
        mapper.to_spdi_unambiguous(&var).unwrap(),
        mapper.to_spdi(&var, false).unwrap(),
    ] {
        if spdi.ends_with("::") {
            // The plain SPDI of an identity keeps no range to read back.
            continue;
        }
        let back = mapper.from_spdi(&spdi).unwrap();
        assert_eq!(back.to_string(), expected, "{hgvs} via {spdi}");
        assert_eq!(
            mapper.to_vrs(&back).unwrap().id,
            vrs.id,
            "{hgvs} via {spdi}"
        );
    }
}

#[test]
fn nucleotide_variants_round_trip_to_their_normalised_form() {
    let hdp = Provider { lookup: true };
    let mapper = VariantMapper::with_refget(&hdp, &hdp);
    // ACGTTTGCAAGGCTAGCTAGCTTTTAACGGGATCGATCGA
    // 1234567890123456789012345678901234567890
    round_trip(&mapper, "NC_TEST.1:g.7G>C", "NC_TEST.1:g.7G>C");
    round_trip(&mapper, "NC_TEST.1:g.4del", "NC_TEST.1:g.6del");
    round_trip(&mapper, "NC_TEST.1:g.4dup", "NC_TEST.1:g.6dup");
    round_trip(&mapper, "NC_TEST.1:g.3_4insT", "NC_TEST.1:g.6dup");
    round_trip(&mapper, "NC_TEST.1:g.12_13insAG", "NC_TEST.1:g.12_13insAG");
    // CTAG at 13-16 rolls over CTAG CT: the same molecule written 3'-most.
    round_trip(&mapper, "NC_TEST.1:g.13_16del", "NC_TEST.1:g.19_22del");
    round_trip(
        &mapper,
        "NC_TEST.1:g.7_8delinsAA",
        "NC_TEST.1:g.7_8delinsAA",
    );
    round_trip(
        &mapper,
        "NC_TEST.1:g.7_9delinsGAT",
        "NC_TEST.1:g.8_9delinsAT",
    );
    round_trip(&mapper, "NC_TEST.1:g.10=", "NC_TEST.1:g.10=");
    round_trip(&mapper, "NC_TEST.1:g.10_12=", "NC_TEST.1:g.10_12=");
    round_trip(&mapper, "NC_TEST.1:g.6_7inv", "NC_TEST.1:g.6_7inv");
}

#[test]
fn protein_variants_round_trip_to_their_normalised_form() {
    let hdp = Provider { lookup: true };
    let mapper = VariantMapper::with_refget(&hdp, &hdp);
    round_trip(&mapper, "NP_TEST.1:p.Lys2Leu", "NP_TEST.1:p.Lys2Leu");
    round_trip(&mapper, "NP_TEST.1:p.Ala4del", "NP_TEST.1:p.Ala6del");
    round_trip(&mapper, "NP_TEST.1:p.Ala4dup", "NP_TEST.1:p.Ala6dup");
    round_trip(
        &mapper,
        "NP_TEST.1:p.Lys2_Leu3insGlyGly",
        "NP_TEST.1:p.Lys2_Leu3insGlyGly",
    );
    round_trip(
        &mapper,
        "NP_TEST.1:p.Leu3_Ala4delinsTrp",
        "NP_TEST.1:p.Leu3_Ala4delinsTrp",
    );
    round_trip(&mapper, "NP_TEST.1:p.Lys2=", "NP_TEST.1:p.Lys2=");
}

#[test]
fn imprecise_deletions_round_trip_as_written() {
    let hdp = Provider { lookup: true };
    let mapper = VariantMapper::with_refget(&hdp, &hdp);
    for hgvs in [
        "NC_TEST.1:g.(?_5)_(10_?)del",
        "NC_TEST.1:g.(3_5)_(10_12)del",
        "NC_TEST.1:g.5_(10_12)del",
        "NC_TEST.1:g.(3_5)_12del",
        "NC_TEST.1:g.(3_5)del",
    ] {
        let vrs = mapper.to_vrs(&parse_hgvs_variant(hgvs).unwrap()).unwrap();
        let back = mapper.from_vrs(&vrs.to_json(), None).unwrap();
        assert_eq!(back.to_string(), hgvs);
        assert_eq!(mapper.to_vrs(&back).unwrap().id, vrs.id);
    }
}

#[test]
fn the_sequence_is_named_by_the_caller_or_looked_up_and_always_checked() {
    let no_lookup = Provider { lookup: false };
    let mapper = VariantMapper::new(&no_lookup);
    let json = mapper
        .to_vrs(&parse_hgvs_variant("NC_TEST.1:g.7G>C").unwrap())
        .unwrap()
        .to_json();
    assert!(matches!(
        mapper.from_vrs(&json, None),
        Err(HgvsError::DataProviderError(_))
    ));
    assert_eq!(
        mapper
            .from_vrs(&json, Some("NC_TEST.1"))
            .unwrap()
            .to_string(),
        "NC_TEST.1:g.7G>C"
    );
    // The digest must be the named sequence's.
    assert!(matches!(
        mapper.from_vrs(&json, Some("NC_OTHER.1")),
        Err(HgvsError::ValidationError(_))
    ));
}

#[test]
fn alleles_from_other_producers_parse() {
    let hdp = Provider { lookup: true };
    let mapper = VariantMapper::with_refget(&hdp, &hdp);
    let refget = refget_accession(GENOME);
    // No ids or digests, no sequence on the reference-length state, extra
    // properties, a Range location elsewhere: what another tool may emit.
    let json = format!(
        r#"{{"type":"Allele","name":"x","location":{{"type":"SequenceLocation",
            "sequenceReference":{{"type":"SequenceReference","refgetAccession":"{refget}"}},
            "start":3,"end":6}},
            "state":{{"type":"ReferenceLengthExpression","length":2,"repeatSubunitLength":1}}}}"#
    );
    assert_eq!(
        mapper.from_vrs(&json, None).unwrap().to_string(),
        "NC_TEST.1:g.6del"
    );

    let json = format!(
        r#"{{"type":"Allele","location":{{"type":"SequenceLocation",
            "sequenceReference":{{"type":"SequenceReference","refgetAccession":"{refget}"}},
            "start":[null,4],"end":[10,null]}},
            "state":{{"type":"LiteralSequenceExpression","sequence":""}}}}"#
    );
    assert_eq!(
        mapper.from_vrs(&json, None).unwrap().to_string(),
        "NC_TEST.1:g.(?_5)_(10_?)del"
    );

    let json = format!(
        r#"{{"type":"Allele","location":{{"type":"SequenceLocation",
            "sequenceReference":{{"type":"SequenceReference","refgetAccession":"{refget}"}},
            "start":[2,4],"end":[10,12]}},
            "state":{{"type":"LiteralSequenceExpression","sequence":"A"}}}}"#
    );
    assert!(matches!(
        mapper.from_vrs(&json, None),
        Err(HgvsError::UnsupportedOperation(_))
    ));

    assert!(matches!(
        mapper.from_vrs(r#"{"type":"SequenceLocation"}"#, None),
        Err(HgvsError::ValidationError(_))
    ));
    assert!(matches!(
        mapper.from_spdi("NC_TEST.1:7"),
        Err(HgvsError::ValidationError(_))
    ));
    assert!(matches!(
        mapper.from_spdi("NC_TEST.1:6:A:C"),
        Err(HgvsError::ValidationError(_))
    ));
    assert_eq!(
        mapper.from_spdi("NC_TEST.1:6:1:C").unwrap().to_string(),
        "NC_TEST.1:g.7G>C"
    );
    assert_eq!(
        mapper.from_spdi("NP_TEST.1:1:K:L").unwrap().to_string(),
        "NP_TEST.1:p.Lys2Leu"
    );
}
