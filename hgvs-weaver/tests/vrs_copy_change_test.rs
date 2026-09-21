//! Duplications and deletions, `g.(a_b)_(c_d)dup`, rendered as VRS
//! `CopyNumberChange` objects (a gain or a loss of copies) and read back.

mod support;

use hgvs_weaver::error::HgvsError;
use hgvs_weaver::mapper::VariantMapper;
use hgvs_weaver::parse_hgvs_variant;
use hgvs_weaver::vrs::VrsBound::{Exact, Range};
use hgvs_weaver::vrs::{VrsCopyChange, VrsVariation};
use support::Provider;

//                     1234567890123456789012345678901234567890
const GENOME: &str = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
const MITO: &str = "TTGCAAGGCTAGCTAGCTTTTAACGGGATCGATCGAACGT";

/// The genome and a distinct mitochondrial sequence, so digests tell them
/// apart.
fn provider() -> Provider {
    Provider::new()
        .sequence("NC_TEST.1", GENOME)
        .sequence("NC_012920.1", MITO)
}

#[test]
fn an_imprecise_duplication_is_a_gain_over_range_bounds() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let var = parse_hgvs_variant("NC_TEST.1:g.(3_5)_(10_12)dup").unwrap();
    let change = mapper.to_vrs_copy_number_change(&var).unwrap();
    assert_eq!(change.type_, "CopyNumberChange");
    assert!(change.id.starts_with("ga4gh:CX."), "{}", change.id);
    assert_eq!(change.id, format!("ga4gh:CX.{}", change.digest));
    assert_eq!(
        (change.location.start, change.location.end),
        (Range(Some(2), Some(4)), Range(Some(10), Some(12)))
    );
    assert_eq!(change.copy_change, VrsCopyChange::Gain);
    assert_eq!(change.expressions.len(), 1);
    assert_eq!(change.expressions[0].syntax, "hgvs.g");
    assert_eq!(change.expressions[0].value, "NC_TEST.1:g.(3_5)_(10_12)dup");

    let json = change.to_json();
    for field in [
        r#""type":"CopyNumberChange""#,
        r#""copyChange":"gain""#,
        r#""location":{"id":"ga4gh:SL."#,
        r#""start":[2,4],"end":[10,12]}"#,
        r#""expressions":[{"syntax":"hgvs.g","value":"NC_TEST.1:g.(3_5)_(10_12)dup"}]"#,
    ] {
        assert!(json.contains(field), "{field} missing from {json}");
    }
    assert!(!json.contains("state"), "{json}");
    assert!(!json.contains("copies"), "{json}");
}

#[test]
fn exact_duplications_and_deletions_have_changes_too() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let dup = parse_hgvs_variant("NC_TEST.1:g.5_12dup").unwrap();
    let del = parse_hgvs_variant("NC_TEST.1:g.5_12del").unwrap();
    let gain = mapper.to_vrs_copy_number_change(&dup).unwrap();
    let loss = mapper.to_vrs_copy_number_change(&del).unwrap();
    assert_eq!(
        (gain.location.start, gain.location.end),
        (Exact(4), Exact(12))
    );
    assert_eq!(gain.copy_change, VrsCopyChange::Gain);
    assert_eq!(loss.copy_change, VrsCopyChange::Loss);
    // Same location, different change, different identifier.
    assert_eq!(gain.location, loss.location);
    assert_ne!(gain.id, loss.id);
    // The location is the one a copy-number count over the range has.
    let count = parse_hgvs_variant("NC_TEST.1:g.5_12copy3").unwrap();
    assert_eq!(
        mapper.to_vrs_copy_number(&count).unwrap().location.id,
        gain.location.id
    );
    // Nothing is normalised: a dup that would roll along the ACGT run as an
    // Allele stays where it was written as a change of copies.
    assert_eq!(mapper.to_vrs(&dup).unwrap().location.start, Exact(0));
}

#[test]
fn to_vrs_variation_gives_a_change_for_an_imprecise_dup_only() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let imprecise_dup = parse_hgvs_variant("NC_TEST.1:g.(3_5)_(10_12)dup").unwrap();
    let imprecise_del = parse_hgvs_variant("NC_TEST.1:g.(3_5)_(10_12)del").unwrap();
    let exact_dup = parse_hgvs_variant("NC_TEST.1:g.5_12dup").unwrap();
    let one_side = parse_hgvs_variant("NC_TEST.1:g.5_(10_12)dup").unwrap();
    for var in [&imprecise_dup, &one_side] {
        match mapper.to_vrs_variation(var).unwrap() {
            VrsVariation::CopyNumberChange(c) => {
                assert_eq!(c.id, mapper.to_vrs_copy_number_change(var).unwrap().id)
            }
            other => panic!("{var}: expected a CopyNumberChange, got {other:?}"),
        }
    }
    // An imprecise deletion keeps its Allele form, an exact dup its own.
    for var in [&imprecise_del, &exact_dup] {
        match mapper.to_vrs_variation(var).unwrap() {
            VrsVariation::Allele(a) => assert_eq!(a.id, mapper.to_vrs(var).unwrap().id),
            other => panic!("{var}: expected an Allele, got {other:?}"),
        }
    }
    // to_vrs itself still refuses, and says where to go.
    let err = mapper.to_vrs(&imprecise_dup).unwrap_err();
    assert!(
        matches!(&err, HgvsError::UnsupportedOperation(m) if m.contains("to_vrs_copy_number_change")),
        "{err}"
    );
}

#[test]
fn copy_number_changes_read_back_as_the_same_hgvs() {
    let hdp = provider();
    let mapper = VariantMapper::with_refget(&hdp, &hdp);
    for hgvs in [
        "NC_TEST.1:g.(3_5)_(10_12)dup",
        "NC_TEST.1:g.(?_5)_(10_?)dup",
        "NC_TEST.1:g.5_(10_12)dup",
        "NC_TEST.1:g.(3_5)dup",
        "NC_TEST.1:g.5_12dup",
        "NC_TEST.1:g.5dup",
        "NC_TEST.1:g.(3_5)_(10_12)del",
        "NC_TEST.1:g.5_12del",
    ] {
        let var = parse_hgvs_variant(hgvs).unwrap();
        let change = mapper.to_vrs_copy_number_change(&var).unwrap();
        for accession in [Some("NC_TEST.1"), None] {
            let back = mapper.from_vrs(&change.to_json(), accession).unwrap();
            assert_eq!(back.to_string(), hgvs);
            assert_eq!(
                mapper.to_vrs_copy_number_change(&back).unwrap().id,
                change.id
            );
        }
    }
    // Through the dispatching pair as well.
    let var = parse_hgvs_variant("NC_TEST.1:g.(3_5)_(10_12)dup").unwrap();
    let json = mapper.to_vrs_variation(&var).unwrap().to_json();
    assert_eq!(
        mapper.from_vrs(&json, None).unwrap().to_string(),
        "NC_TEST.1:g.(3_5)_(10_12)dup"
    );
}

#[test]
fn mitochondrial_changes_are_projected_to_genomic() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let m = parse_hgvs_variant("NC_012920.1:m.(3_5)_(10_12)dup").unwrap();
    let g = parse_hgvs_variant("NC_012920.1:g.(3_5)_(10_12)dup").unwrap();
    let from_m = mapper.to_vrs_copy_number_change(&m).unwrap();
    assert_eq!(from_m.id, mapper.to_vrs_copy_number_change(&g).unwrap().id);
    assert_eq!(from_m.expressions[0].syntax, "hgvs.m");
    assert_eq!(
        from_m.expressions[0].value,
        "NC_012920.1:m.(3_5)_(10_12)dup"
    );
    assert!(matches!(
        mapper.to_vrs_variation(&m).unwrap(),
        VrsVariation::CopyNumberChange(_)
    ));
    let back = mapper
        .from_vrs(&from_m.to_json(), Some("NC_012920.1"))
        .unwrap();
    assert_eq!(back.to_string(), "NC_012920.1:g.(3_5)_(10_12)dup");
}

#[test]
fn the_whole_gain_family_is_a_dup_and_the_loss_family_a_del() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let var = parse_hgvs_variant("NC_TEST.1:g.(3_5)_(10_12)dup").unwrap();
    let json = mapper.to_vrs_copy_number_change(&var).unwrap().to_json();
    let with =
        |term: &str| json.replace(r#""copyChange":"gain""#, &format!(r#""copyChange":{term}"#));
    for term in [
        r#""gain""#,
        r#""low-level gain""#,
        r#""high-level gain""#,
        r#""EFO:0030070""#,
        r#""EFO:0030071""#,
        r#""EFO:0030072""#,
        // The VRS 2.0.0 MappableConcept form.
        r#"{"primaryCoding":{"code":"EFO:0030070","system":"https://www.ebi.ac.uk/efo/"}}"#,
    ] {
        let back = mapper.from_vrs(&with(term), Some("NC_TEST.1")).unwrap();
        assert_eq!(back.to_string(), "NC_TEST.1:g.(3_5)_(10_12)dup", "{term}");
    }
    for term in [
        r#""loss""#,
        r#""low-level loss""#,
        r#""high-level loss""#,
        r#""complete genomic loss""#,
        r#""EFO:0030067""#,
        r#""EFO:0030068""#,
        r#""EFO:0030069""#,
        r#""EFO:0020073""#,
        r#"{"primaryCode":"EFO:0030067"}"#,
    ] {
        let back = mapper.from_vrs(&with(term), Some("NC_TEST.1")).unwrap();
        assert_eq!(back.to_string(), "NC_TEST.1:g.(3_5)_(10_12)del", "{term}");
    }
    for term in [
        r#""regional base ploidy""#,
        r#""EFO:0030064""#,
        r#""EFO:0000001""#,
        r#""amplification""#,
    ] {
        let err = mapper.from_vrs(&with(term), Some("NC_TEST.1")).unwrap_err();
        assert!(
            matches!(err, HgvsError::UnsupportedOperation(_)),
            "{term}: {err}"
        );
    }
    // Not a term at all.
    assert!(matches!(
        mapper.from_vrs(&with("7"), Some("NC_TEST.1")).unwrap_err(),
        HgvsError::ValidationError(_)
    ));
    // The refget accession must be the sequence's.
    assert!(matches!(
        mapper
            .from_vrs(&json.replace("SQ.", "SQ.x"), Some("NC_TEST.1"))
            .unwrap_err(),
        HgvsError::ValidationError(_)
    ));
}

#[test]
fn only_genomic_duplications_and_deletions_have_changes() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    for hgvs in [
        "NC_TEST.1:g.5_12inv",
        "NC_TEST.1:g.5A>T",
        "NC_TEST.1:g.5_12copy3",
        "NC_TEST.1:g.(3_5)_(10_12)copy3",
        "NM_TEST.1:c.5_12dup",
        "NP_TEST.1:p.Ala5dup",
    ] {
        let var = parse_hgvs_variant(hgvs).unwrap();
        let err = mapper.to_vrs_copy_number_change(&var).unwrap_err();
        assert!(
            matches!(err, HgvsError::UnsupportedOperation(_)),
            "{hgvs}: {err}"
        );
    }
    // Imprecise breakpoints on anything but a dup or del are still refused
    // everywhere.
    let inv = parse_hgvs_variant("NC_TEST.1:g.(3_5)_(10_12)inv").unwrap();
    assert!(matches!(
        mapper.to_vrs_variation(&inv).unwrap_err(),
        HgvsError::UnsupportedOperation(_)
    ));
}
