//! VRS and SPDI read back: `from_vrs` and `from_spdi` invert `to_vrs` and
//! `to_spdi`, on nucleotide and protein sequences.

mod support;

use hgvs_weaver::error::HgvsError;
use hgvs_weaver::mapper::VariantMapper;
use hgvs_weaver::parse_hgvs_variant;
use hgvs_weaver::vrs::refget_accession;
use support::Provider;

const GENOME: &str = "ACGTTTGCAAGGCTAGCTAGCTTTTAACGGGATCGATCGA";
const OTHER: &str = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
const PROTEIN: &str = "MKLAAAYRQ";

/// Serves three sequences; as a `Refget` it looks accessions up by digest.
fn provider() -> Provider {
    Provider::new()
        .sequence("NC_TEST.1", GENOME)
        .sequence("NC_OTHER.1", OTHER)
        .sequence("NP_TEST.1", PROTEIN)
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
    let hdp = provider();
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
    let hdp = provider();
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
    let hdp = provider();
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
    let no_lookup = provider();
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
    let hdp = provider();
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

/// One test per case: the SPDI `$spdi` names a range no sequence can hold
/// and is refused, rather than read as a different allele or overflowing.
macro_rules! an_spdi_past_any_sequence_is_refused {
    ($($name:ident: $spdi:literal;)*) => {$(
        #[test]
        fn $name() {
            let hdp = provider();
            let mapper = VariantMapper::with_refget(&hdp, &hdp);
            let result = mapper.from_spdi($spdi);
            assert!(matches!(result, Err(HgvsError::ValidationError(_))), "{result:?}");
        }
    )*};
}

an_spdi_past_any_sequence_is_refused! {
    an_spdi_deletion_length_past_usize: "NC_TEST.1:6:99999999999999999999999:C";
    an_spdi_deletion_length_running_past_usize: "NC_TEST.1:18446744073709551615:5:C";
    an_spdi_deletion_running_past_usize: "NC_TEST.1:18446744073709551615:A:C";
    an_spdi_insertion_past_the_sequence: "NC_TEST.1:41:0:C";
    an_spdi_insertion_past_i32: "NC_TEST.1:2147483648:0:C";
    an_spdi_insertion_past_u32: "NC_TEST.1:4294967296:0:C";
    an_spdi_insertion_at_the_last_index: "NC_TEST.1:18446744073709551615:0:C";
    an_spdi_deletion_length_past_the_sequence: "NC_TEST.1:0:99999999999999:C";
}

/// A deletion from VRS bounds `start`, `end` on GENOME, as JSON.
fn bounded_deletion_json(start: &str, end: &str) -> String {
    let refget = refget_accession(GENOME);
    format!(
        r#"{{"type":"Allele","location":{{"type":"SequenceLocation",
            "sequenceReference":{{"type":"SequenceReference","refgetAccession":"{refget}"}},
            "start":{start},"end":{end}}},
            "state":{{"type":"LiteralSequenceExpression","sequence":""}}}}"#
    )
}

/// A location on GENOME from `start` to `end`, as the JSON of a VRS object.
fn location_json(start: &str, end: &str) -> String {
    let refget = refget_accession(GENOME);
    format!(
        r#""location":{{"type":"SequenceLocation",
            "sequenceReference":{{"type":"SequenceReference","refgetAccession":"{refget}"}},
            "start":{start},"end":{end}}}"#
    )
}

/// One test per case: the VRS object `$json` locates its change at
/// coordinates no g. position holds, and is refused rather than written at a
/// wrapped or zero position, or panicking.
macro_rules! a_vrs_location_past_any_position_is_refused {
    ($($name:ident: $json:expr;)*) => {$(
        #[test]
        fn $name() {
            let hdp = provider();
            let mapper = VariantMapper::with_refget(&hdp, &hdp);
            let result = mapper.from_vrs(&$json, None);
            assert!(matches!(result, Err(HgvsError::ValidationError(_))), "{result:?}");
        }
    )*};
}

a_vrs_location_past_any_position_is_refused! {
    a_copy_number_at_the_last_u32: format!(
        r#"{{"type":"CopyNumberCount","copies":3,{}}}"#,
        location_json("4294967295", "4294967296")
    );
    a_copy_number_past_i32: format!(
        r#"{{"type":"CopyNumberCount","copies":3,{}}}"#,
        location_json("2147483647", "2147483648")
    );
    a_copy_number_past_u32: format!(
        r#"{{"type":"CopyNumberCount","copies":3,{}}}"#,
        location_json("4294967300", "4294967310")
    );
    a_deletion_range_from_the_last_index: format!(
        r#"{{"type":"Allele",{},"state":{{"type":"LiteralSequenceExpression","sequence":""}}}}"#,
        location_json("[18446744073709551615,null]", "[null,null]")
    );
    a_length_expression_at_the_last_index: format!(
        r#"{{"type":"Allele",{},"state":{{"type":"LengthExpression","length":5}}}}"#,
        location_json("18446744073709551615", "18446744073709551615")
    );
}

#[test]
fn a_copy_number_at_the_last_g_position_keeps_it() {
    let hdp = provider();
    let mapper = VariantMapper::with_refget(&hdp, &hdp);
    let json = format!(
        r#"{{"type":"CopyNumberCount","copies":3,{}}}"#,
        location_json("2147483646", "2147483647")
    );
    assert_eq!(
        mapper.from_vrs(&json, None).unwrap().to_string(),
        "NC_TEST.1:g.2147483647copy3"
    );
}

#[test]
fn a_vrs_deletion_ending_at_interbase_zero_is_refused() {
    let hdp = provider();
    let mapper = VariantMapper::with_refget(&hdp, &hdp);
    let result = mapper.from_vrs(&bounded_deletion_json("[null,4]", "0"), None);
    assert!(
        matches!(result, Err(HgvsError::ValidationError(_))),
        "{result:?}"
    );
}

#[test]
fn a_vrs_deletion_whose_end_range_starts_at_interbase_zero_is_refused() {
    let hdp = provider();
    let mapper = VariantMapper::with_refget(&hdp, &hdp);
    let result = mapper.from_vrs(&bounded_deletion_json("[null,4]", "[0,10]"), None);
    assert!(
        matches!(result, Err(HgvsError::ValidationError(_))),
        "{result:?}"
    );
}

#[test]
fn an_insertion_before_the_first_base_reads_back_as_a_delins_of_that_base() {
    // HGVS has no insertion before base 1. An allele that puts T in front of
    // the sequence (which starts with A) is written as g.1delinsTA.
    let hdp = provider();
    let mapper = VariantMapper::with_refget(&hdp, &hdp);
    let refget = refget_accession(GENOME);
    let json = format!(
        r#"{{"type":"Allele","location":{{"type":"SequenceLocation",
            "sequenceReference":{{"type":"SequenceReference","refgetAccession":"{refget}"}},
            "start":0,"end":0}},
            "state":{{"type":"LiteralSequenceExpression","sequence":"T"}}}}"#
    );
    assert_eq!(
        mapper.from_vrs(&json, None).unwrap().to_string(),
        "NC_TEST.1:g.1delinsTA"
    );
}
