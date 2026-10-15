//! Statements about a whole transcript — `r.0`, `r.spl`, `r.?` — read as what
//! they mean, print as HGVS writes them, and survive JSON, including the
//! `Special` shape older dictionaries hold.

use hgvs_weaver::edits::{NaEdit, RnaStatement};
use hgvs_weaver::{parse_hgvs_variant, SequenceVariant};

fn parse(s: &str) -> SequenceVariant {
    parse_hgvs_variant(s).unwrap_or_else(|e| panic!("{s}: {e}"))
}

fn na_edit(v: &SequenceVariant) -> &NaEdit {
    match v {
        SequenceVariant::Rna(r) => &r.posedit.edit,
        other => panic!("{other} is not r."),
    }
}

fn from_json(json: &str) -> Result<SequenceVariant, serde_json::Error> {
    serde_json::from_str(json)
}

fn json_round_trip(v: &SequenceVariant) -> SequenceVariant {
    from_json(&serde_json::to_string(v).unwrap()).unwrap()
}

/// An r. variant on NM_W.1 whose edit is `edit`, as JSON.
fn r_json(edit: &str) -> String {
    format!(
        r#"{{"variant_type": "Rna", "ac": "NM_W.1", "gene": null, "posedit": {{"pos": null,
            "edit": {edit}, "uncertain": false, "predicted": false}}}}"#
    )
}

/// One test per case: parsing `$hgvs` gives the edit `$edit`.
macro_rules! reads_as {
    ($edit_of:ident; $($name:ident: $hgvs:literal => $edit:expr;)*) => {$(
        #[test]
        fn $name() {
            assert_eq!($edit_of(&parse($hgvs)), &$edit);
        }
    )*};
}

/// One test per case: `$hgvs` survives JSON unchanged.
macro_rules! survives_json {
    ($($name:ident: $hgvs:literal;)*) => {$(
        #[test]
        fn $name() {
            let v = parse($hgvs);
            assert_eq!(json_round_trip(&v), v);
        }
    )*};
}

/// One test per case: the old `Special` value `$value` loads as `$edit`.
macro_rules! special_loads_as {
    ($edit_of:ident, $json:expr; $($name:ident: $value:literal => $edit:expr;)*) => {$(
        #[test]
        fn $name() {
            let special = format!(r#"{{"type": "Special", "value": "{}", "uncertain": false}}"#, $value);
            let v = from_json(&$json(&special)).unwrap();
            assert_eq!($edit_of(&v), &$edit);
        }
    )*};
}

// --- Parsing reads a statement by what it says ---

reads_as! { na_edit;
    r_0_is_no_rna: "NM_W.1:r.0" => NaEdit::Statement(RnaStatement::NoRna { uncertain: false });
    r_0_uncertain_is_probably_no_rna: "NM_W.1:r.0?" => NaEdit::Statement(RnaStatement::NoRna { uncertain: true });
    r_spl_is_splicing_affected: "NM_W.1:r.spl" => NaEdit::Statement(RnaStatement::SplicingAffected { uncertain: false });
    r_spl_uncertain_is_splicing_probably_affected: "NM_W.1:r.spl?" => NaEdit::Statement(RnaStatement::SplicingAffected { uncertain: true });
    r_unknown_is_an_unknown_effect: "NM_W.1:r.?" => NaEdit::Statement(RnaStatement::Unknown);
    r_spl_with_a_space_before_its_question_mark_is_uncertain: "NM_W.1:r.spl ?" => NaEdit::Statement(RnaStatement::SplicingAffected { uncertain: true });
}

// --- Printing writes each statement as HGVS does ---

#[test]
fn an_r_statement_written_with_a_space_prints_without_it() {
    assert_eq!(parse("NM_W.1:r.spl ?").to_string(), "NM_W.1:r.spl?");
}

// --- JSON in the current shape ---

survives_json! {
    r_0_survives_json: "NM_W.1:r.0";
    r_0_uncertain_survives_json: "NM_W.1:r.0?";
    r_spl_survives_json: "NM_W.1:r.spl";
    r_spl_uncertain_survives_json: "NM_W.1:r.spl?";
    r_unknown_survives_json: "NM_W.1:r.?";
}

#[test]
fn a_json_statement_weaver_does_not_know_is_refused() {
    let json = r_json(r#"{"type": "Statement", "statement": "Maybe"}"#);
    let err = from_json(&json).unwrap_err().to_string();
    assert!(err.contains("unknown variant `Maybe`"), "{err}");
}

#[test]
fn a_json_rna_statement_on_a_coding_variant_is_refused() {
    let json = r#"{"variant_type": "Coding", "ac": "NM_W.1", "gene": null, "posedit": {"pos": null,
        "edit": {"type": "Statement", "statement": "SplicingAffected", "uncertain": false},
        "uncertain": false, "predicted": false}}"#;
    let err = from_json(json).unwrap_err().to_string();
    assert!(err.contains("c. cannot write"), "{err}");
}

// --- JSON in the older Special shape still loads ---

special_loads_as! { na_edit, r_json;
    special_0_loads_as_no_rna: "0" => NaEdit::Statement(RnaStatement::NoRna { uncertain: false });
    special_0_uncertain_loads_as_probably_no_rna: "0?" => NaEdit::Statement(RnaStatement::NoRna { uncertain: true });
    special_spl_loads_as_splicing_affected: "spl" => NaEdit::Statement(RnaStatement::SplicingAffected { uncertain: false });
    special_spl_uncertain_loads_as_splicing_probably_affected: "spl?" => NaEdit::Statement(RnaStatement::SplicingAffected { uncertain: true });
    special_spl_with_a_space_loads_as_splicing_probably_affected: "spl ?" => NaEdit::Statement(RnaStatement::SplicingAffected { uncertain: true });
    special_unknown_loads_as_an_unknown_rna_effect: "?" => NaEdit::Statement(RnaStatement::Unknown);
    special_identity_loads_as_the_rna_unchanged: "=" => NaEdit::WholeIdentity;
}

#[test]
fn a_special_value_the_r_grammar_does_not_write_is_refused() {
    let json = r_json(r#"{"type": "Special", "value": "??", "uncertain": false}"#);
    let err = from_json(&json).unwrap_err().to_string();
    assert!(err.contains("not a statement about RNA"), "{err}");
}

#[test]
fn a_special_with_its_uncertain_flag_set_is_refused() {
    // The `?` was part of the value; no description set the flag.
    let json = r_json(r#"{"type": "Special", "value": "spl", "uncertain": true}"#);
    let err = from_json(&json).unwrap_err().to_string();
    assert!(err.contains("never uncertain"), "{err}");
}

#[test]
fn a_special_identity_on_a_coding_variant_is_refused() {
    let json = r#"{"variant_type": "Coding", "ac": "NM_W.1", "gene": null, "posedit": {"pos": null,
        "edit": {"type": "Special", "value": "=", "uncertain": false},
        "uncertain": false, "predicted": false}}"#;
    let err = from_json(json).unwrap_err().to_string();
    assert!(err.contains("c. cannot write"), "{err}");
}
