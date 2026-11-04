//! VRS and SPDI read back: `from_vrs` and `from_spdi` invert `to_vrs` and
//! `to_spdi`, on nucleotide and protein sequences.

mod support;

use std::collections::BTreeSet;

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
        Err(HgvsError::ReferenceMismatch { .. })
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

/// A `CisPhasedBlock` of `members`, each an Allele on `sequence` from
/// interbase `start` to `end` whose state is `alt`.
fn block_on(sequence: &str, members: &[(usize, usize, &str)]) -> String {
    let refget = refget_accession(sequence);
    let alleles: Vec<String> = members
        .iter()
        .map(|(start, end, alt)| {
            format!(
                r#"{{"type":"Allele","location":{{"type":"SequenceLocation",
                    "sequenceReference":{{"type":"SequenceReference","refgetAccession":"{refget}"}},
                    "start":{start},"end":{end}}},
                    "state":{{"type":"LiteralSequenceExpression","sequence":"{alt}"}}}}"#
            )
        })
        .collect();
    format!(
        r#"{{"type":"CisPhasedBlock","members":[{}]}}"#,
        alleles.join(",")
    )
}

/// The readings `from_vrs` lists for the ambiguous `block` on `X`, which is
/// `sequence`.
fn ambiguous_readings(sequence: &str, block: &str) -> BTreeSet<String> {
    let hdp = Provider::new().sequence("X", sequence);
    match VariantMapper::with_refget(&hdp, &hdp).from_vrs(block, None) {
        Err(HgvsError::Ambiguous { readings, .. }) => {
            readings.iter().map(|r| r.to_string()).collect()
        }
        other => panic!("expected an ambiguity, got {other:?}"),
    }
}

#[test]
fn a_block_inserting_twice_at_one_point_reads_back_as_ambiguous() {
    // An A and a G between g.4 and g.5 of TTCAGCAGTT, in either order: AG,
    // or GA, a copy of g.4_5.
    let block = block_on("TTCAGCAGTT", &[(4, 4, "A"), (4, 4, "G")]);
    assert_eq!(
        ambiguous_readings("TTCAGCAGTT", &block),
        BTreeSet::from(["X:g.4_5insAG".to_string(), "X:g.4_5dup".to_string()])
    );
}

#[test]
fn a_block_of_two_identical_insertions_reads_back_as_ambiguous() {
    // One T written twice, or two Ts, in the run of TTCAGCAGTT's g.1_2.
    let block = block_on("TTCAGCAGTT", &[(2, 2, "T"), (2, 2, "T")]);
    assert_eq!(
        ambiguous_readings("TTCAGCAGTT", &block),
        BTreeSet::from(["X:g.2dup".to_string(), "X:g.1_2dup".to_string()])
    );
}

#[test]
fn a_block_of_two_identical_insertions_after_a_base_reads_back_as_ambiguous() {
    // One A after g.4 of TTCAGCAGTT written twice, or two As.
    let block = block_on("TTCAGCAGTT", &[(4, 4, "A"), (4, 4, "A")]);
    assert_eq!(
        ambiguous_readings("TTCAGCAGTT", &block),
        BTreeSet::from(["X:g.4dup".to_string(), "X:g.4_5insAA".to_string()])
    );
}

#[test]
fn a_protein_block_reads_back_in_residue_order() {
    // MKLAAAYRQ: Tyr7His listed before Lys2Asn, as a block's digest order may.
    let hdp = provider();
    let mapper = VariantMapper::with_refget(&hdp, &hdp);
    let block = block_on(PROTEIN, &[(6, 7, "H"), (1, 2, "N")]);
    assert_eq!(
        mapper.from_vrs(&block, None).unwrap().to_string(),
        "NP_TEST.1:p.[Lys2Asn;Tyr7His]"
    );
}

/// The protein allele `hgvs` rendered by `to_vrs_variation` and read back:
/// the reading must render to the same identifier.
fn protein_allele_read_back(mapper: &VariantMapper, hgvs: &str) -> String {
    let vrs = mapper
        .to_vrs_variation(&parse_hgvs_variant(hgvs).unwrap())
        .unwrap();
    let back = mapper.from_vrs(&vrs.to_json(), None).unwrap();
    assert_eq!(
        mapper.to_vrs_variation(&back).unwrap().id(),
        vrs.id(),
        "{hgvs} read back as {back}"
    );
    back.to_string()
}

#[test]
fn a_protein_deletion_in_a_run_beside_a_change_to_its_last_residue_reads_back() {
    // MKLAAAYRQ: an Ala deleted from Ala4_Ala6 and Ala6 changed leave MKLAGYRQ.
    let hdp = provider();
    let mapper = VariantMapper::with_refget(&hdp, &hdp);
    assert_eq!(
        protein_allele_read_back(&mapper, "NP_TEST.1:p.[Ala4del;Ala6Gly]"),
        "NP_TEST.1:p.Ala5_Ala6delinsGly"
    );
}

#[test]
fn a_protein_duplication_in_a_run_beside_a_change_to_its_last_residue_reads_back() {
    // MKLAAAYRQ: an Ala added to Ala4_Ala6 and Ala6 changed leave MKLAAAGYRQ.
    let hdp = provider();
    let mapper = VariantMapper::with_refget(&hdp, &hdp);
    assert_eq!(
        protein_allele_read_back(&mapper, "NP_TEST.1:p.[Ala4dup;Ala6Gly]"),
        "NP_TEST.1:p.Ala6_Tyr7insGly"
    );
}

#[test]
fn two_protein_deletions_in_one_run_read_back_as_both() {
    // MKLAAAYRQ: two of the three Ala deleted leave MKLAYRQ.
    let hdp = provider();
    let mapper = VariantMapper::with_refget(&hdp, &hdp);
    assert_eq!(
        protein_allele_read_back(&mapper, "NP_TEST.1:p.[Ala4del;Ala5del]"),
        "NP_TEST.1:p.Ala5_Ala6del"
    );
}

#[test]
fn a_protein_member_meeting_none_stays_its_own_member_of_the_block() {
    let hdp = provider();
    let mapper = VariantMapper::with_refget(&hdp, &hdp);
    assert_eq!(
        protein_allele_read_back(&mapper, "NP_TEST.1:p.[Lys2Asn;Ala4del;Ala6Gly]"),
        "NP_TEST.1:p.[Lys2Asn;Ala5_Ala6delinsGly]"
    );
}

#[test]
fn a_protein_deletion_in_a_run_beside_an_identity_reads_back_short_of_it() {
    let hdp = provider();
    let mapper = VariantMapper::with_refget(&hdp, &hdp);
    assert_eq!(
        protein_allele_read_back(&mapper, "NP_TEST.1:p.[Ala4del;Ala6=]"),
        "NP_TEST.1:p.[Ala5del;Ala6=]"
    );
}

#[test]
fn a_protein_duplication_in_a_run_beside_an_identity_reads_back_short_of_it() {
    let hdp = provider();
    let mapper = VariantMapper::with_refget(&hdp, &hdp);
    assert_eq!(
        protein_allele_read_back(&mapper, "NP_TEST.1:p.[Ala4dup;Ala6=]"),
        "NP_TEST.1:p.[Ala5dup;Ala6=]"
    );
}

#[test]
fn spellings_of_one_protein_allele_beside_an_identity_have_one_vrs_identifier() {
    // MKLAAAYRQ: an Ala added to Ala4_Ala6, written anywhere in the run.
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let id = |hgvs: &str| {
        let var = parse_hgvs_variant(hgvs).unwrap();
        mapper.to_vrs_variation(&var).unwrap().id().to_string()
    };
    assert_eq!(
        id("NP_TEST.1:p.[Ala4dup;Tyr7=]"),
        id("NP_TEST.1:p.[Ala6dup;Tyr7=]")
    );
}

#[test]
fn protein_changes_that_touch_render_as_the_one_change_they_make() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let allele = parse_hgvs_variant("NP_TEST.1:p.[Lys2Asn;Leu3Pro]").unwrap();
    let one = parse_hgvs_variant("NP_TEST.1:p.Lys2_Leu3delinsAsnPro").unwrap();
    assert_eq!(
        mapper.to_vrs_variation(&allele).unwrap().id(),
        mapper.to_vrs(&one).unwrap().id
    );
}

#[test]
fn a_protein_change_written_twice_renders_once() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let allele = parse_hgvs_variant("NP_TEST.1:p.[Lys2Asn;Lys2delinsAsn]").unwrap();
    let one = parse_hgvs_variant("NP_TEST.1:p.Lys2Asn").unwrap();
    assert_eq!(
        mapper.to_vrs_variation(&allele).unwrap().id(),
        mapper.to_vrs(&one).unwrap().id
    );
}

#[test]
fn a_protein_duplication_written_twice_renders_once() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let allele = parse_hgvs_variant("NP_TEST.1:p.[Ala4dup;Ala4dup]").unwrap();
    let one = parse_hgvs_variant("NP_TEST.1:p.Ala4dup").unwrap();
    assert_eq!(
        mapper.to_vrs_variation(&allele).unwrap().id(),
        mapper.to_vrs(&one).unwrap().id
    );
}

#[test]
fn a_protein_insertion_copying_a_stated_residue_reads_back() {
    let hdp = provider();
    let mapper = VariantMapper::with_refget(&hdp, &hdp);
    assert_eq!(
        protein_allele_read_back(&mapper, "NP_TEST.1:p.[Ala5=;Ala5_Ala6insAla;Ala6=]"),
        "NP_TEST.1:p.[Ala5=;Ala5dup;Ala6=]"
    );
}

#[test]
fn a_protein_duplication_of_a_stated_residue_has_a_reading() {
    // It inserts a copy beside Ala5, which it leaves as it was.
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let var = parse_hgvs_variant("NP_TEST.1:p.[Ala5=;Ala5dup]").unwrap();
    assert!(mapper.validate(&var).unwrap());
}

#[test]
fn a_protein_extension_past_a_stated_stop_has_no_reading() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let var = parse_hgvs_variant("NP_TEST.1:p.[Ter10=;Ter10GlnextTer5]").unwrap();
    assert!(has_no_reading(mapper.validate(&var)));
}

#[test]
fn a_protein_extension_before_a_stated_met1_has_no_reading() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let var = parse_hgvs_variant("NP_TEST.1:p.[Met1=;Met1ext-5]").unwrap();
    assert!(has_no_reading(mapper.validate(&var)));
}

#[test]
fn to_vrs_refuses_protein_changes_that_together_leave_the_protein_as_it_is() {
    // An Ala deleted from Ala4_Ala6 and one added back.
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let var = parse_hgvs_variant("NP_TEST.1:p.[Ala4del;Ala5dup]").unwrap();
    assert!(matches!(
        mapper.to_vrs_variation(&var),
        Err(HgvsError::UnsupportedOperation(m)) if m.contains("leave the protein as it is")
    ));
}

#[test]
fn a_protein_duplication_and_an_insertion_at_its_copy_have_no_reading() {
    // Ala4dup puts its copy between Ala4 and Ala5, where the Gly goes too.
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let var = parse_hgvs_variant("NP_TEST.1:p.[Ala4dup;Ala4_Ala5insGly]").unwrap();
    assert!(has_no_reading(mapper.validate(&var)));
}

#[test]
fn a_protein_insertion_beside_an_identity_has_a_reading() {
    // The same molecule as p.[Ala6dup;Tyr7=]: the insertion changes no residue Tyr7= states.
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let var = parse_hgvs_variant("NP_TEST.1:p.[Ala6_Tyr7insAla;Tyr7=]").unwrap();
    assert!(mapper.validate(&var).unwrap());
}

#[test]
fn a_protein_insertion_inside_an_identity_has_no_reading() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let var = parse_hgvs_variant("NP_TEST.1:p.[Ala4_Ala6=;Ala4_Ala5insGly]").unwrap();
    assert!(has_no_reading(mapper.validate(&var)));
}

#[test]
fn a_protein_allele_that_is_one_change_renders_as_its_allele() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let allele = parse_hgvs_variant("NP_TEST.1:p.[Ala4del;Ala6Gly]").unwrap();
    let one = parse_hgvs_variant("NP_TEST.1:p.Ala5_Ala6delinsGly").unwrap();
    assert_eq!(
        mapper.to_vrs_variation(&allele).unwrap().id(),
        mapper.to_vrs(&one).unwrap().id
    );
}

/// Whether `result` is the refusal of an allele with no one reading.
fn has_no_reading<T: std::fmt::Debug>(result: Result<T, HgvsError>) -> bool {
    matches!(&result, Err(HgvsError::ValidationError(m)) if m.contains("no reading"))
}

#[test]
fn to_vrs_refuses_a_protein_allele_inserting_twice_at_one_point() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let var = parse_hgvs_variant("NP_TEST.1:p.[Lys2_Leu3insGln;Lys2_Leu3insGly]").unwrap();
    assert!(has_no_reading(mapper.to_vrs_variation(&var)));
}

#[test]
fn to_vrs_refuses_a_protein_allele_changing_one_residue_two_ways() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let var = parse_hgvs_variant("NP_TEST.1:p.[Lys2Asn;Lys2Thr]").unwrap();
    assert!(has_no_reading(mapper.to_vrs_variation(&var)));
}

#[test]
fn validate_refuses_a_protein_allele_changing_one_residue_two_ways() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let var = parse_hgvs_variant("NP_TEST.1:p.[Lys2Asn;Lys2Thr]").unwrap();
    assert!(has_no_reading(mapper.validate(&var)));
}

#[test]
fn a_protein_block_changing_one_residue_two_ways_does_not_read_back() {
    let hdp = provider();
    let mapper = VariantMapper::with_refget(&hdp, &hdp);
    let block = block_on(PROTEIN, &[(1, 2, "N"), (1, 2, "T")]);
    assert!(has_no_reading(mapper.from_vrs(&block, None)));
}

#[test]
fn a_protein_block_inserting_twice_after_the_last_residue_does_not_read_back() {
    // MKLAAAYRQ is served without its stop, so index 9 is before the stop.
    let hdp = provider();
    let mapper = VariantMapper::with_refget(&hdp, &hdp);
    let block = block_on(PROTEIN, &[(9, 9, "Q"), (9, 9, "G")]);
    assert!(has_no_reading(mapper.from_vrs(&block, None)));
}

#[test]
fn a_protein_insertion_beside_a_change_on_its_right_has_a_reading() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let var = parse_hgvs_variant("NP_TEST.1:p.[Lys2_Leu3insGln;Leu3Pro]").unwrap();
    assert!(mapper.validate(&var).unwrap());
}

#[test]
fn a_protein_insertion_beside_a_change_on_its_left_has_a_reading() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let var = parse_hgvs_variant("NP_TEST.1:p.[Lys2Pro;Lys2_Leu3insGln]").unwrap();
    assert!(mapper.validate(&var).unwrap());
}

#[test]
fn a_protein_change_written_twice_is_read_once() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let var = parse_hgvs_variant("NP_TEST.1:p.[Lys2Asn;Tyr7Asp;Tyr7Asp]").unwrap();
    assert!(mapper.validate(&var).unwrap());
}

#[test]
fn a_protein_member_naming_no_residues_does_not_hide_two_that_clash() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let var = parse_hgvs_variant("NP_TEST.1:p.[Met1?;Lys2Asn;Lys2Thr]").unwrap();
    assert!(has_no_reading(mapper.validate(&var)));
}

#[test]
fn normalize_refuses_a_protein_allele_changing_one_residue_two_ways() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let var = parse_hgvs_variant("NP_TEST.1:p.[Lys2Asn;Lys2Thr]").unwrap();
    assert!(has_no_reading(mapper.normalize_variant(var)));
}
