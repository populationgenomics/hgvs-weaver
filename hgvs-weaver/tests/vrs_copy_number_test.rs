//! Copy-number edits, `g.1000_2000copy3`, rendered as VRS `CopyNumberCount`
//! objects and read back.

use hgvs_weaver::allele::CanonicalAllele;
use hgvs_weaver::data::{DataProvider, IdentifierKind, IdentifierType, TranscriptData};
use hgvs_weaver::error::HgvsError;
use hgvs_weaver::mapper::VariantMapper;
use hgvs_weaver::parse_hgvs_variant;
use hgvs_weaver::vrs::VrsBound::{Exact, Range};
use hgvs_weaver::vrs::{VrsAllele, VrsMolecule, VrsVariation};

const GENOME: &str = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

struct Provider;

impl DataProvider for Provider {
    fn get_transcript(&self, ac: &str, _: Option<&str>) -> Result<TranscriptData, HgvsError> {
        Err(HgvsError::DataProviderError(format!("no transcript {ac}")))
    }

    fn get_seq(
        &self,
        _ac: &str,
        start: i32,
        end: Option<i32>,
        _kind: IdentifierType,
    ) -> Result<String, HgvsError> {
        let start = (start.max(0) as usize).min(GENOME.len());
        let end = end.map_or(GENOME.len(), |e| (e.max(0) as usize).min(GENOME.len()));
        Ok(GENOME[start..end.max(start)].to_string())
    }

    fn get_symbol_accessions(
        &self,
        _: &str,
        _: IdentifierKind,
        _: IdentifierKind,
    ) -> Result<Vec<(IdentifierType, String)>, HgvsError> {
        Ok(vec![])
    }

    fn get_identifier_type(&self, _: &str) -> Result<IdentifierType, HgvsError> {
        Ok(IdentifierType::GenomicAccession)
    }
}

#[test]
fn a_copy_number_edit_is_a_copy_number_count_over_the_range() {
    let hdp = Provider;
    let mapper = VariantMapper::new(&hdp);
    let var = parse_hgvs_variant("NC_TEST.1:g.5_12copy3").unwrap();
    let count = mapper.to_vrs_copy_number(&var).unwrap();
    assert_eq!(count.type_, "CopyNumberCount");
    assert!(count.id.starts_with("ga4gh:CN."), "{}", count.id);
    assert_eq!(count.id, format!("ga4gh:CN.{}", count.digest));
    assert_eq!(
        (count.location.start, count.location.end),
        (Exact(4), Exact(12))
    );
    assert_eq!(count.copies, Exact(3));
    assert_eq!(count.location.sequence_reference.residue_alphabet, "na");
    assert_eq!(count.expressions.len(), 1);
    assert_eq!(count.expressions[0].syntax, "hgvs.g");
    assert_eq!(count.expressions[0].value, "NC_TEST.1:g.5_12copy3");

    let json = count.to_json();
    for field in [
        r#""type":"CopyNumberCount""#,
        r#""copies":3"#,
        r#""location":{"id":"ga4gh:SL."#,
        r#""start":4,"end":12}"#,
        r#""expressions":[{"syntax":"hgvs.g","value":"NC_TEST.1:g.5_12copy3"}]"#,
    ] {
        assert!(json.contains(field), "{field} missing from {json}");
    }
    assert!(!json.contains("state"), "{json}");
}

#[test]
fn the_location_is_the_one_an_allele_over_the_same_range_has() {
    let hdp = Provider;
    let mapper = VariantMapper::new(&hdp);
    let var = parse_hgvs_variant("NC_TEST.1:g.5_12copy3").unwrap();
    let count = mapper.to_vrs_copy_number(&var).unwrap();
    let allele = VrsAllele::new(
        &CanonicalAllele {
            accession: "NC_TEST.1".into(),
            start: 4,
            end: 12,
            reference: GENOME[4..12].into(),
            alternate: String::new(),
            repeat_subunit: None,
        },
        &count.location.sequence_reference.refget_accession,
        VrsMolecule::Genomic,
        None,
    );
    assert_eq!(count.location.id, allele.location.id);
    assert_eq!(count.location, allele.location);
}

#[test]
fn copy_number_counts_read_back_as_the_same_hgvs() {
    let hdp = Provider;
    let mapper = VariantMapper::new(&hdp);
    for hgvs in [
        "NC_TEST.1:g.5_12copy3",
        "NC_TEST.1:g.5copy2",
        "NC_TEST.1:g.(3_5)_(10_12)copy4",
        "NC_TEST.1:g.(?_5)_(10_?)copy0",
    ] {
        let var = parse_hgvs_variant(hgvs).unwrap();
        let json = mapper.to_vrs_copy_number(&var).unwrap().to_json();
        let back = mapper.from_vrs(&json, Some("NC_TEST.1")).unwrap();
        assert_eq!(back.to_string(), hgvs);
        assert_eq!(
            mapper.to_vrs_copy_number(&back).unwrap().id,
            mapper.to_vrs_copy_number(&var).unwrap().id
        );
    }
}

#[test]
fn uncertain_breakpoints_become_ranges() {
    let hdp = Provider;
    let mapper = VariantMapper::new(&hdp);
    let var = parse_hgvs_variant("NC_TEST.1:g.(3_5)_(10_12)copy4").unwrap();
    let count = mapper.to_vrs_copy_number(&var).unwrap();
    assert_eq!(
        (count.location.start, count.location.end),
        (Range(Some(2), Some(4)), Range(Some(10), Some(12)))
    );
}

#[test]
fn mitochondrial_copy_numbers_are_projected_to_genomic() {
    let hdp = Provider;
    let mapper = VariantMapper::new(&hdp);
    let m = parse_hgvs_variant("NC_012920.1:m.5_12copy3").unwrap();
    let g = parse_hgvs_variant("NC_012920.1:g.5_12copy3").unwrap();
    let from_m = mapper.to_vrs_copy_number(&m).unwrap();
    let from_g = mapper.to_vrs_copy_number(&g).unwrap();
    assert_eq!(from_m.id, from_g.id);
    assert_eq!(from_m.expressions[0].syntax, "hgvs.m");
    assert_eq!(from_m.expressions[0].value, "NC_012920.1:m.5_12copy3");
    let back = mapper
        .from_vrs(&from_m.to_json(), Some("NC_012920.1"))
        .unwrap();
    assert_eq!(back.to_string(), "NC_012920.1:g.5_12copy3");
}

#[test]
fn only_copy_number_edits_on_genomic_sequences_have_counts() {
    let hdp = Provider;
    let mapper = VariantMapper::new(&hdp);
    for hgvs in [
        "NC_TEST.1:g.5_12del",
        "NC_TEST.1:g.5_12dup",
        "NC_TEST.1:g.5A>T",
        "NM_TEST.1:c.5_12copy2",
        "NP_TEST.1:p.Ala5del",
    ] {
        let var = parse_hgvs_variant(hgvs).unwrap();
        let err = mapper.to_vrs_copy_number(&var).unwrap_err();
        assert!(
            matches!(err, HgvsError::UnsupportedOperation(_)),
            "{hgvs}: {err}"
        );
    }
}

#[test]
fn to_vrs_variation_picks_the_object_by_edit() {
    let hdp = Provider;
    let mapper = VariantMapper::new(&hdp);
    let copy = parse_hgvs_variant("NC_TEST.1:g.5_12copy3").unwrap();
    let del = parse_hgvs_variant("NC_TEST.1:g.5_12del").unwrap();
    match mapper.to_vrs_variation(&copy).unwrap() {
        VrsVariation::CopyNumberCount(c) => assert!(c.id.starts_with("ga4gh:CN.")),
        other => panic!("expected a CopyNumberCount, got {other:?}"),
    }
    match mapper.to_vrs_variation(&del).unwrap() {
        VrsVariation::Allele(a) => assert!(a.id.starts_with("ga4gh:VA.")),
        other => panic!("expected an Allele, got {other:?}"),
    }
    // to_vrs itself still refuses: it returns an Allele.
    assert!(matches!(
        mapper.to_vrs(&copy).unwrap_err(),
        HgvsError::UnsupportedOperation(_)
    ));
}

#[test]
fn from_vrs_checks_the_count_it_is_given() {
    let hdp = Provider;
    let mapper = VariantMapper::new(&hdp);
    let var = parse_hgvs_variant("NC_TEST.1:g.5_12copy3").unwrap();
    let json = mapper.to_vrs_copy_number(&var).unwrap().to_json();

    // HGVS has no range of counts.
    let ranged = json.replace(r#""copies":3"#, r#""copies":[3,null]"#);
    assert!(matches!(
        mapper.from_vrs(&ranged, Some("NC_TEST.1")).unwrap_err(),
        HgvsError::UnsupportedOperation(_)
    ));
    // The refget accession must be the sequence's.
    let wrong = json.replace("SQ.", "SQ.x");
    assert!(matches!(
        mapper.from_vrs(&wrong, Some("NC_TEST.1")).unwrap_err(),
        HgvsError::ValidationError(_)
    ));
    // Other types are still refused as before.
    let other = json.replace("CopyNumberCount", "Haplotype");
    let err = mapper.from_vrs(&other, Some("NC_TEST.1")).unwrap_err();
    assert!(
        matches!(&err, HgvsError::ValidationError(m) if m.contains("Allele")),
        "{err}"
    );
}
