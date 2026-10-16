//! Positions are numbered from 1: a 0 is refused, and an unknown position is
//! `?`, never a number.

use hgvs_weaver::error::HgvsError;
use hgvs_weaver::{parse_hgvs_variant, SequenceVariant};

fn parse(s: &str) -> SequenceVariant {
    parse_hgvs_variant(s).unwrap_or_else(|e| panic!("{s}: {e}"))
}

fn from_json(json: &str) -> Result<SequenceVariant, serde_json::Error> {
    serde_json::from_str(json)
}

/// A g. deletion on NC_W.1 whose interval is `pos`, as JSON.
fn g_del_json(pos: &str) -> String {
    format!(
        r#"{{"variant_type": "Genomic", "ac": "NC_W.1", "gene": null, "posedit": {{"pos": {pos},
            "edit": {{"type": "Del", "ref_": null, "uncertain": false}}, "uncertain": false, "predicted": false}}}}"#
    )
}

/// One test per case: `$hgvs` names a position 0 and does not parse.
macro_rules! position_zero_is_a_parse_error {
    ($($name:ident: $hgvs:literal;)*) => {$(
        #[test]
        fn $name() {
            let result = parse_hgvs_variant($hgvs);
            assert!(matches!(result, Err(HgvsError::PestError(_))), "{result:?}");
        }
    )*};
}

position_zero_is_a_parse_error! {
    a_genomic_substitution_at_zero: "NC_W.1:g.0A>G";
    a_genomic_range_from_zero: "NC_W.1:g.0_5del";
    an_uncertain_genomic_bound_at_zero: "NC_W.1:g.(0_5)_(10_12)del";
    a_mitochondrial_substitution_at_zero: "NC_W.1:m.0A>G";
    a_coding_substitution_at_zero: "NM_W.1:c.0A>G";
    a_coding_substitution_at_minus_zero: "NM_W.1:c.-0A>G";
    a_coding_substitution_at_star_zero: "NM_W.1:c.*0A>G";
    a_coding_range_to_zero: "NM_W.1:c.-5_0del";
    a_non_coding_substitution_at_zero: "NR_W.1:n.0A>G";
    an_rna_substitution_at_zero: "NM_W.1:r.0a>g";
}

#[test]
fn a_dict_for_a_coding_position_of_zero_is_refused() {
    let json = r#"{"variant_type": "Coding", "ac": "NM_W.1", "gene": null, "posedit": {"pos": {"start":
        {"base": 0, "offset": null, "anchor": "CdsStart", "uncertain": false}, "end": null, "uncertain": false},
        "edit": {"type": "Del", "ref_": null, "uncertain": false}, "uncertain": false, "predicted": false}}"#;
    assert!(from_json(json).is_err(), "{json}");
}

#[test]
fn an_unknown_genomic_bound_reparses_to_itself() {
    let v = parse("NC_W.1:g.(?_5)_(10_?)del");
    assert_eq!(parse(&v.to_string()), v);
}

#[test]
fn an_unknown_genomic_bound_survives_a_dict() {
    let v = parse("NC_W.1:g.(?_5)_(10_?)del");
    let json = serde_json::to_string(&v).unwrap();
    assert_eq!(from_json(&json).unwrap(), v);
}

#[test]
fn a_dict_for_an_unknown_genomic_bound_names_it_with_a_question_mark() {
    let json = g_del_json(
        r#"{"start": {"base": "?", "end": 5, "uncertain": true},
            "end": {"base": 10, "end": "?", "uncertain": true}, "uncertain": false}"#,
    );
    assert_eq!(from_json(&json).unwrap(), parse("NC_W.1:g.(?_5)_(10_?)del"));
}

#[test]
fn a_dict_for_a_genomic_position_of_zero_is_refused() {
    let json = g_del_json(
        r#"{"start": {"base": 0, "end": null, "uncertain": false}, "end": null, "uncertain": false}"#,
    );
    assert!(from_json(&json).is_err(), "{json}");
}
