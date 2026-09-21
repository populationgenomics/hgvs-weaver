//! Insertions of bases known only by number, `insN[20]`, rendered as VRS
//! Alleles with a `LengthExpression` state and read back.

mod support;

use hgvs_weaver::data::Strand;
use hgvs_weaver::error::HgvsError;
use hgvs_weaver::mapper::VariantMapper;
use hgvs_weaver::parse_hgvs_variant;
use hgvs_weaver::vrs::VrsBound::{Exact, Range};
use hgvs_weaver::vrs::{VrsState, VrsVariation};
use support::{single_exon_transcript, Provider};

//                     1234567890123456789012345678901234567890
const GENOME: &str = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
const MITO: &str = "TTGCAAGGCTAGCTAGCTTTTAACGGGATCGATCGAACGT";

/// The genome, a distinct mitochondrial sequence (so digests tell them
/// apart), and one transcript over all of the genome on the plus strand, so
/// `c.` numbers as `g.`.
fn provider() -> Provider {
    Provider::new()
        .sequence("NC_TEST.1", GENOME)
        .sequence("NC_012920.1", MITO)
        .sequence("NM_TEST.1", GENOME)
        .transcript(single_exon_transcript(
            "NM_TEST.1",
            "NC_TEST.1",
            0,
            Strand::Plus,
            0,
            39,
            40,
        ))
}

#[test]
fn a_length_insertion_is_an_allele_at_the_insertion_point() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let var = parse_hgvs_variant("NC_TEST.1:g.10_11insN[20]").unwrap();
    let vrs = mapper.to_vrs(&var).unwrap();
    assert!(vrs.id.starts_with("ga4gh:VA."), "{}", vrs.id);
    // Between bases 10 and 11: the empty interbase range at 10.
    assert_eq!((vrs.location.start, vrs.location.end), (Exact(10), Exact(10)));
    assert!(matches!(
        vrs.state,
        VrsState::Length {
            length: Exact(20),
            ..
        }
    ));
    assert_eq!(vrs.expressions[0].syntax, "hgvs.g");
    assert_eq!(vrs.expressions[0].value, "NC_TEST.1:g.10_11insN[20]");
    let json = vrs.to_json();
    assert!(
        json.contains(r#""state":{"type":"LengthExpression","length":20}"#),
        "{json}"
    );
    assert!(json.contains(r#""start":10,"end":10}"#), "{json}");
}

#[test]
fn the_older_spelling_gives_the_same_allele() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let recommended = parse_hgvs_variant("NC_TEST.1:g.10_11insN[20]").unwrap();
    let older = parse_hgvs_variant("NC_TEST.1:g.10_11ins(20)").unwrap();
    assert_eq!(
        mapper.to_vrs(&older).unwrap().id,
        mapper.to_vrs(&recommended).unwrap().id
    );
    // The expression carries the recommended spelling either way.
    assert_eq!(
        mapper.to_vrs(&older).unwrap().expressions[0].value,
        "NC_TEST.1:g.10_11insN[20]"
    );
}

#[test]
fn an_uncertain_length_is_a_range() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    for hgvs in [
        "NC_TEST.1:g.10_11insN[(20_30)]",
        "NC_TEST.1:g.10_11insN[20_30]",
        "NC_TEST.1:g.10_11ins(20_30)",
    ] {
        let var = parse_hgvs_variant(hgvs).unwrap();
        let vrs = mapper.to_vrs(&var).unwrap();
        assert!(
            matches!(
                vrs.state,
                VrsState::Length {
                    length: Range(Some(20), Some(30)),
                    ..
                }
            ),
            "{hgvs}: {:?}",
            vrs.state
        );
        assert!(vrs.to_json().contains(r#""length":[20,30]"#), "{hgvs}");
    }
}

#[test]
fn a_delins_of_a_stated_length_covers_the_deleted_range() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let var = parse_hgvs_variant("NC_TEST.1:g.10_12delinsN[5]").unwrap();
    let vrs = mapper.to_vrs(&var).unwrap();
    assert_eq!((vrs.location.start, vrs.location.end), (Exact(9), Exact(12)));
    assert!(matches!(
        vrs.state,
        VrsState::Length {
            length: Exact(5),
            ..
        }
    ));
    // Stated deleted bases or their count do not change the allele.
    for other in [
        "NC_TEST.1:g.10_12del3insN[5]",
        "NC_TEST.1:g.10_12delCGTinsN[5]",
        "NC_TEST.1:g.10_12delins(5)",
    ] {
        let v = parse_hgvs_variant(other).unwrap();
        assert_eq!(mapper.to_vrs(&v).unwrap().id, vrs.id, "{other}");
    }
}

#[test]
fn length_expressions_are_not_normalised() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    // An insertion of known bases here would slide 3' along the ACGT run;
    // unknown bases stay where they were written.
    let var = parse_hgvs_variant("NC_TEST.1:g.4_5insN[4]").unwrap();
    let vrs = mapper.to_vrs(&var).unwrap();
    assert_eq!((vrs.location.start, vrs.location.end), (Exact(4), Exact(4)));
    let known = parse_hgvs_variant("NC_TEST.1:g.4_5insACGT").unwrap();
    assert_ne!(mapper.to_vrs(&known).unwrap().location.start, Exact(4));
}

#[test]
fn transcript_and_mitochondrial_insertions_are_projected_to_genomic() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let g = parse_hgvs_variant("NC_TEST.1:g.10_11insN[20]").unwrap();
    let c = parse_hgvs_variant("NM_TEST.1:c.10_11insN[20]").unwrap();
    let from_c = mapper.to_vrs(&c).unwrap();
    assert_eq!(from_c.id, mapper.to_vrs(&g).unwrap().id);
    assert_eq!(from_c.expressions[0].syntax, "hgvs.c");
    assert_eq!(from_c.expressions[0].value, "NM_TEST.1:c.10_11insN[20]");

    let m = parse_hgvs_variant("NC_012920.1:m.10_11insN[20]").unwrap();
    let from_m = mapper.to_vrs(&m).unwrap();
    assert_eq!(from_m.expressions[0].syntax, "hgvs.m");
    let back = mapper
        .from_vrs(&from_m.to_json(), Some("NC_012920.1"))
        .unwrap();
    assert_eq!(back.to_string(), "NC_012920.1:g.10_11insN[20]");
}

#[test]
fn length_expressions_read_back_as_the_recommended_spelling() {
    let hdp = provider();
    let mapper = VariantMapper::with_refget(&hdp, &hdp);
    for (hgvs, expected) in [
        ("NC_TEST.1:g.10_11insN[20]", "NC_TEST.1:g.10_11insN[20]"),
        ("NC_TEST.1:g.10_11ins(20)", "NC_TEST.1:g.10_11insN[20]"),
        (
            "NC_TEST.1:g.10_11insN[(20_30)]",
            "NC_TEST.1:g.10_11insN[(20_30)]",
        ),
        ("NC_TEST.1:g.10_11ins(20_30)", "NC_TEST.1:g.10_11insN[(20_30)]"),
        ("NC_TEST.1:g.10_12delinsN[5]", "NC_TEST.1:g.10_12delinsN[5]"),
        ("NC_TEST.1:g.10_12del3insN[5]", "NC_TEST.1:g.10_12delinsN[5]"),
        ("NC_TEST.1:g.10delinsN[(2_3)]", "NC_TEST.1:g.10delinsN[(2_3)]"),
        ("NC_TEST.1:g.39_40insN[1]", "NC_TEST.1:g.39_40insN[1]"),
    ] {
        let var = parse_hgvs_variant(hgvs).unwrap();
        let vrs = mapper.to_vrs(&var).unwrap();
        // With the accession given and looked up from the digest alike.
        for accession in [Some("NC_TEST.1"), None] {
            let back = mapper.from_vrs(&vrs.to_json(), accession).unwrap();
            assert_eq!(back.to_string(), expected, "{hgvs}");
            assert_eq!(mapper.to_vrs(&back).unwrap().id, vrs.id, "{hgvs}");
        }
    }
}

#[test]
fn to_vrs_variation_gives_the_allele() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let var = parse_hgvs_variant("NC_TEST.1:g.10_11insN[20]").unwrap();
    match mapper.to_vrs_variation(&var).unwrap() {
        VrsVariation::Allele(a) => assert_eq!(a.id, mapper.to_vrs(&var).unwrap().id),
        other => panic!("expected an Allele, got {other:?}"),
    }
}

#[test]
fn a_length_expression_the_sequence_cannot_hold_is_refused() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let var = parse_hgvs_variant("NC_TEST.1:g.10_11insN[20]").unwrap();
    let json = mapper.to_vrs(&var).unwrap().to_json();
    // HGVS has no insertion before the first base.
    let at_zero = json.replace(r#""start":10,"end":10"#, r#""start":0,"end":0"#);
    assert!(matches!(
        mapper.from_vrs(&at_zero, Some("NC_TEST.1")).unwrap_err(),
        HgvsError::UnsupportedOperation(_)
    ));
    // Nor one after the last.
    let past_end = json.replace(r#""start":10,"end":10"#, r#""start":40,"end":40"#);
    assert!(matches!(
        mapper.from_vrs(&past_end, Some("NC_TEST.1")).unwrap_err(),
        HgvsError::ValidationError(_)
    ));
    // An open-ended range of lengths has no HGVS spelling here.
    let open = json.replace(r#""length":20"#, r#""length":[20,null]"#);
    assert!(matches!(
        mapper.from_vrs(&open, Some("NC_TEST.1")).unwrap_err(),
        HgvsError::UnsupportedOperation(_)
    ));
}

#[test]
fn a_length_insertion_has_no_canonical_allele_or_spdi() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let var = parse_hgvs_variant("NC_TEST.1:g.10_11insN[20]").unwrap();
    assert!(matches!(
        mapper.canonical_allele(&var).unwrap_err(),
        HgvsError::UnsupportedOperation(_)
    ));
    assert!(matches!(
        mapper.to_spdi(&var, false).unwrap_err(),
        HgvsError::UnsupportedOperation(_)
    ));
    // Validation has nothing to check and normalisation nothing to move.
    assert!(mapper.validate(&var).unwrap());
    assert_eq!(
        mapper.normalize_variant(var.clone()).unwrap().to_string(),
        var.to_string()
    );
}
