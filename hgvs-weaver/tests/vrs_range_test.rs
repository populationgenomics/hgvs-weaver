//! VRS `Range` bounds for deletions whose breakpoints are uncertain.

mod support;

use hgvs_weaver::error::HgvsError;
use hgvs_weaver::mapper::VariantMapper;
use hgvs_weaver::parse_hgvs_variant;
use hgvs_weaver::vrs::VrsBound::{Exact, Range};
use support::Provider;

const GENOME: &str = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

/// The same genome under a nuclear and the mitochondrial accession.
fn provider() -> Provider {
    Provider::new()
        .sequence("NC_TEST.1", GENOME)
        .sequence("NC_012920.1", GENOME)
}

fn bounds(
    mapper: &VariantMapper,
    hgvs: &str,
) -> Result<(hgvs_weaver::vrs::VrsBound, hgvs_weaver::vrs::VrsBound), HgvsError> {
    let var = parse_hgvs_variant(hgvs).unwrap();
    let vrs = mapper.to_vrs(&var)?;
    assert!(vrs.id.starts_with("ga4gh:VA."));
    assert_eq!(vrs.expressions[0].value, hgvs);
    Ok((vrs.location.start, vrs.location.end))
}

#[test]
fn uncertain_breakpoints_become_ranges() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let b = |h| bounds(&mapper, h).unwrap();
    assert_eq!(
        b("NC_TEST.1:g.(3_5)_(10_12)del"),
        (Range(Some(2), Some(4)), Range(Some(10), Some(12)))
    );
    assert_eq!(
        b("NC_TEST.1:g.(?_5)_(10_?)del"),
        (Range(None, Some(4)), Range(Some(10), None))
    );
    assert_eq!(
        b("NC_TEST.1:g.5_(10_12)del"),
        (Exact(4), Range(Some(10), Some(12)))
    );
    assert_eq!(
        b("NC_TEST.1:g.(3_5)_12del"),
        (Range(Some(2), Some(4)), Exact(12))
    );
    // One base somewhere in 3..=5.
    assert_eq!(
        b("NC_TEST.1:g.(3_5)del"),
        (Range(Some(2), Some(4)), Range(Some(3), Some(5)))
    );
    // Parenthesised exact positions are exact, so the deletion is normalised
    // like any other: here it rolls over the whole ACGT run.
    assert_eq!(b("NC_TEST.1:g.(5)_(12)del"), b("NC_TEST.1:g.5_12del"));
    assert_eq!(
        b("NC_012920.1:m.(3_5)_(10_12)del"),
        (Range(Some(2), Some(4)), Range(Some(10), Some(12)))
    );
}

#[test]
fn an_imprecise_deletion_has_an_empty_literal_state_and_json_ranges() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let var = parse_hgvs_variant("NC_TEST.1:g.(?_5)_(10_?)del").unwrap();
    let json = mapper.to_vrs(&var).unwrap().to_json();
    assert!(json.contains(r#""start":[null,4]"#), "{json}");
    assert!(json.contains(r#""end":[10,null]"#), "{json}");
    assert!(
        json.contains(r#""state":{"type":"LiteralSequenceExpression","sequence":""}"#),
        "{json}"
    );
}

#[test]
fn only_deletions_may_have_uncertain_breakpoints() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    for hgvs in [
        "NC_TEST.1:g.(3_5)_(10_12)dup",
        "NC_TEST.1:g.(3_5)_(10_12)inv",
        "NC_TEST.1:g.(3_5)_(10_12)delinsA",
    ] {
        let err = bounds(&mapper, hgvs).unwrap_err();
        assert!(
            matches!(err, HgvsError::UnsupportedOperation(_)),
            "{hgvs}: {err}"
        );
    }
    // An exact deletion still goes through normalisation.
    let var = parse_hgvs_variant("NC_TEST.1:g.5_12del").unwrap();
    assert_eq!(mapper.to_vrs(&var).unwrap().location.start, Exact(0));
}

#[test]
fn unknown_bounds_round_trip_as_question_marks() {
    for hgvs in [
        "NC_TEST.1:g.(?_5)_(10_?)del",
        "NC_012920.1:m.(?_5)_(10_12)del",
        "NC_TEST.1:g.?_12del",
    ] {
        assert_eq!(parse_hgvs_variant(hgvs).unwrap().to_string(), hgvs);
    }
}
