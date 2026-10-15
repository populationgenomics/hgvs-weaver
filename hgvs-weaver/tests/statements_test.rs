//! Statements about a whole transcript or protein — `r.0`, `r.spl`, `r.?`,
//! `p.0`, `p.?`, `p.=` — read as what they mean, print as HGVS writes them,
//! and survive JSON, in their own shape and in the `Special` one.

mod support;

use hgvs_weaver::data::Strand;
use hgvs_weaver::edits::{AaEdit, NaEdit, ProteinStatement, RnaStatement};
use hgvs_weaver::equivalence::{EquivalenceLevel, VariantEquivalence};
use hgvs_weaver::mapper::VariantMapper;
use hgvs_weaver::{
    parse_hgvs_variant, transform_variant, SequenceVariant, StartCodonConvention,
    VariantTransformSettings,
};
use support::{exon, transcript, Provider};

/// NM_W.1, CDS at transcript indices 5..=94 encoding NP_W.1, spliced from
/// genome[10, 60) and genome[100, 160) of a CGT repeat with ATG at 15, AAAA
/// at 30 and TAA at 142.
fn provider() -> Provider {
    let mut g: Vec<u8> = b"CGT".iter().copied().cycle().take(200).collect();
    g[15..18].copy_from_slice(b"ATG");
    g[30..34].copy_from_slice(b"AAAA");
    g[142..145].copy_from_slice(b"TAA");
    let g = String::from_utf8(g).unwrap();
    let spliced = format!("{}{}", &g[10..60], &g[100..160]);
    Provider::new()
        .sequence("NC_W.1", &g)
        .sequence("NM_W.1", &spliced)
        .sequence("NP_W.1", "MRRRRKSRRRRRRRRVVVVVVVVVVVVVV")
        .transcript(transcript(
            "NM_W.1",
            "NC_W.1",
            Strand::Plus,
            Some((5, 94)),
            vec![
                exon((0, 50), (10, 59), Strand::Plus),
                exon((50, 110), (100, 159), Strand::Plus),
            ],
        ))
        .protein_for("NM_W.1", "NP_W.1")
}

fn parse(s: &str) -> SequenceVariant {
    parse_hgvs_variant(s).unwrap_or_else(|e| panic!("{s}: {e}"))
}

fn c_to_p(s: &str) -> SequenceVariant {
    let SequenceVariant::Coding(c) = parse(s) else {
        panic!("{s} is not c.");
    };
    let hdp = provider();
    let p = VariantMapper::new(&hdp).c_to_p(&c, None).unwrap();
    SequenceVariant::Protein(p)
}

fn r_to_p(s: &str) -> SequenceVariant {
    let SequenceVariant::Rna(r) = parse(s) else {
        panic!("{s} is not r.");
    };
    let hdp = provider();
    let p = VariantMapper::new(&hdp).r_to_p(&r, None).unwrap();
    SequenceVariant::Protein(p)
}

fn na_edit(v: &SequenceVariant) -> &NaEdit {
    match v {
        SequenceVariant::Rna(r) => &r.posedit.edit,
        other => panic!("{other} is not r."),
    }
}

fn aa_edit(v: &SequenceVariant) -> &AaEdit {
    match v {
        SequenceVariant::Protein(p) => &p.posedit.edit,
        other => panic!("{other} is not p."),
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

/// A p. variant on NP_W.1 whose edit is `edit`, as JSON.
fn p_json(edit: &str, predicted: bool) -> String {
    format!(
        r#"{{"variant_type": "Protein", "ac": "NP_W.1", "gene": null, "posedit": {{"pos": null,
            "edit": {edit}, "uncertain": false, "predicted": {predicted}}}}}"#
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

reads_as! { aa_edit;
    p_0_is_no_protein: "NP_W.1:p.0" => AaEdit::Statement(ProteinStatement::NoProtein { uncertain: false });
    p_unknown_is_an_unknown_effect: "NP_W.1:p.?" => AaEdit::Statement(ProteinStatement::Unknown);
    p_identity_is_the_protein_unchanged: "NP_W.1:p.=" => AaEdit::Statement(ProteinStatement::Unchanged);
    predicted_p_identity_is_the_protein_unchanged: "NP_W.1:p.(=)" => AaEdit::Statement(ProteinStatement::Unchanged);
    p_0_uncertain_is_probably_no_protein: "NP_W.1:p.0?" => AaEdit::Statement(ProteinStatement::NoProtein { uncertain: true });
    an_unknown_change_from_a_residue_is_an_unknown_effect: "NP_W.1:p.Arg12?" => AaEdit::Statement(ProteinStatement::Unknown);
}

#[test]
fn p_met1_unknown_parses_to_what_the_mapper_builds() {
    assert_eq!(parse("NP_W.1:p.Met1?"), c_to_p("NM_W.1:c.-3_4del"));
}

#[test]
fn an_unknown_change_over_a_range_does_not_parse() {
    assert!(parse_hgvs_variant("NP_W.1:p.Met1_Arg3?").is_err());
}

#[test]
fn a_predicted_probably_no_protein_does_not_parse() {
    assert!(parse_hgvs_variant("NP_W.1:p.(0?)").is_err());
}

#[test]
fn a_residue_numbered_0_does_not_parse() {
    assert!(parse_hgvs_variant("NP_W.1:p.Met0Val").is_err());
}

#[test]
fn stated_p_identity_is_not_a_prediction() {
    let SequenceVariant::Protein(p) = parse("NP_W.1:p.=") else {
        panic!("not p.");
    };
    assert!(!p.posedit.predicted);
}

#[test]
fn predicted_p_identity_carries_the_prediction_on_its_posedit() {
    let SequenceVariant::Protein(p) = parse("NP_W.1:p.(=)") else {
        panic!("not p.");
    };
    assert!(p.posedit.predicted);
}

// --- Printing writes each statement as HGVS does ---

#[test]
fn p_0_prints_back_as_parsed() {
    assert_eq!(parse("NP_W.1:p.0").to_string(), "NP_W.1:p.0");
}

#[test]
fn p_unknown_prints_back_as_parsed() {
    assert_eq!(parse("NP_W.1:p.?").to_string(), "NP_W.1:p.?");
}

#[test]
fn p_identity_prints_back_as_parsed() {
    assert_eq!(parse("NP_W.1:p.=").to_string(), "NP_W.1:p.=");
}

#[test]
fn predicted_p_identity_prints_back_as_parsed() {
    assert_eq!(parse("NP_W.1:p.(=)").to_string(), "NP_W.1:p.(=)");
}

#[test]
fn p_0_uncertain_prints_back_as_parsed() {
    assert_eq!(parse("NP_W.1:p.0?").to_string(), "NP_W.1:p.0?");
}

#[test]
fn an_unknown_change_from_a_residue_prints_back_as_parsed() {
    assert_eq!(parse("NP_W.1:p.Arg12?").to_string(), "NP_W.1:p.Arg12?");
}

#[test]
fn a_predicted_unknown_change_from_met1_prints_back_as_parsed() {
    assert_eq!(parse("NP_W.1:p.(Met1?)").to_string(), "NP_W.1:p.(Met1?)");
}

#[test]
fn an_r_statement_written_with_a_space_prints_without_it() {
    assert_eq!(parse("NM_W.1:r.spl ?").to_string(), "NM_W.1:r.spl?");
}

#[test]
fn a_deleted_cds_is_probably_no_protein() {
    let p = c_to_p("NM_W.1:c.-3_*2del");
    assert_eq!(
        aa_edit(&p),
        &AaEdit::Statement(ProteinStatement::NoProtein { uncertain: true })
    );
}

#[test]
fn the_mappers_probably_no_protein_reparses_to_itself() {
    let p = c_to_p("NM_W.1:c.-3_*2del");
    assert_eq!(parse(&p.to_string()), p);
}

#[test]
fn the_mappers_no_protein_reparses_to_itself() {
    let p = r_to_p("NM_W.1:r.0");
    assert_eq!(parse(&p.to_string()), p);
}

#[test]
fn the_mappers_unknown_effect_reparses_to_itself() {
    let p = c_to_p("NM_W.1:c.45+1G>A");
    assert_eq!(parse(&p.to_string()), p);
}

#[test]
fn the_mappers_unknown_change_from_met1_reparses_to_itself() {
    let p = c_to_p("NM_W.1:c.-3_4del");
    assert_eq!(parse(&p.to_string()), p);
}

#[test]
fn the_start_codon_transforms_unknown_change_from_met1_reparses_to_itself() {
    let settings = VariantTransformSettings {
        start_codon: StartCodonConvention::HgvsQuestion,
    };
    let p = transform_variant(&parse("NP_W.1:p.(Met1Val)"), &settings);
    assert_eq!(parse(&p.to_string()), p);
}

#[test]
fn a_start_loss_is_the_parsed_unknown_change_from_met1() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let level = VariantEquivalence::new(&mapper, &hdp)
        .equivalent_level(&parse("NM_W.1:c.-3_4del"), &parse("NP_W.1:p.Met1?"))
        .unwrap();
    assert_eq!(level, EquivalenceLevel::Identity);
}

#[test]
fn a_3_utr_change_predicts_the_protein_unchanged_and_reparses_to_it() {
    let p = c_to_p("NM_W.1:c.*5T>A");
    assert_eq!(p.to_string(), "NP_W.1:p.(=)");
    assert_eq!(parse(&p.to_string()), p);
}

// --- JSON in the current shape ---

survives_json! {
    r_0_survives_json: "NM_W.1:r.0";
    r_0_uncertain_survives_json: "NM_W.1:r.0?";
    r_spl_survives_json: "NM_W.1:r.spl";
    r_spl_uncertain_survives_json: "NM_W.1:r.spl?";
    r_unknown_survives_json: "NM_W.1:r.?";
    p_0_survives_json: "NP_W.1:p.0";
    p_unknown_survives_json: "NP_W.1:p.?";
    p_identity_survives_json: "NP_W.1:p.=";
    predicted_p_identity_survives_json: "NP_W.1:p.(=)";
}

#[test]
fn probably_no_protein_survives_json() {
    let p = c_to_p("NM_W.1:c.-3_*2del");
    assert_eq!(json_round_trip(&p), p);
}

#[test]
fn an_unknown_change_from_met1_survives_json() {
    let p = c_to_p("NM_W.1:c.-3_4del");
    assert_eq!(p.to_string(), "NP_W.1:p.Met1?");
    assert_eq!(json_round_trip(&p), p);
}

#[test]
fn a_json_statement_weaver_does_not_know_is_refused() {
    let json = r_json(r#"{"type": "Statement", "statement": "Maybe"}"#);
    let err = from_json(&json).unwrap_err().to_string();
    assert!(err.contains("Maybe"), "{err}");
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
}

special_loads_as! { aa_edit, |edit| p_json(edit, false);
    special_0_loads_as_no_protein: "0" => AaEdit::Statement(ProteinStatement::NoProtein { uncertain: false });
    special_0_uncertain_loads_as_probably_no_protein: "0?" => AaEdit::Statement(ProteinStatement::NoProtein { uncertain: true });
    special_unknown_loads_as_an_unknown_protein_effect: "?" => AaEdit::Statement(ProteinStatement::Unknown);
    special_identity_loads_as_the_protein_unchanged: "=" => AaEdit::Statement(ProteinStatement::Unchanged);
    special_identity_with_spaces_loads_as_the_protein_unchanged: " = " => AaEdit::Statement(ProteinStatement::Unchanged);
    special_0_with_a_space_loads_as_probably_no_protein: "0 ?" => AaEdit::Statement(ProteinStatement::NoProtein { uncertain: true });
}

#[test]
fn a_predicted_special_identity_loads_as_predicted_p_identity() {
    let json = p_json(
        r#"{"type": "Special", "value": "=", "uncertain": false}"#,
        true,
    );
    assert_eq!(from_json(&json).unwrap(), parse("NP_W.1:p.(=)"));
}

#[test]
fn a_special_value_the_r_grammar_does_not_write_is_refused() {
    let json = r_json(r#"{"type": "Special", "value": "??", "uncertain": false}"#);
    let err = from_json(&json).unwrap_err().to_string();
    assert!(err.contains("not a statement about RNA"), "{err}");
}

#[test]
fn a_special_rna_statement_given_to_a_protein_is_refused() {
    let json = p_json(
        r#"{"type": "Special", "value": "spl", "uncertain": false}"#,
        false,
    );
    let err = from_json(&json).unwrap_err().to_string();
    assert!(err.contains("not a statement about a protein"), "{err}");
}

#[test]
fn a_special_with_its_uncertain_flag_set_is_refused() {
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

#[test]
fn a_special_unknown_at_met1_loads_as_the_mappers_met1_unknown() {
    let json = r#"{"variant_type": "Protein", "ac": "NP_W.1", "gene": null, "posedit": {
        "pos": {"start": {"base": 1, "aa": "Met", "uncertain": false}, "end": null, "uncertain": false},
        "edit": {"type": "Special", "value": "?", "uncertain": false}, "uncertain": false, "predicted": false}}"#;
    assert_eq!(from_json(json).unwrap(), c_to_p("NM_W.1:c.-3_4del"));
}

#[test]
fn a_special_protein_statement_with_its_uncertain_flag_set_is_refused() {
    let json = p_json(
        r#"{"type": "Special", "value": "0", "uncertain": true}"#,
        false,
    );
    let err = from_json(&json).unwrap_err().to_string();
    assert!(err.contains("never uncertain"), "{err}");
}

#[test]
fn an_rna_unknown_naming_an_uncertainty_it_cannot_have_is_refused() {
    let json = r_json(r#"{"type": "Statement", "statement": "Unknown", "uncertain": true}"#);
    let err = from_json(&json).unwrap_err().to_string();
    assert!(err.contains("fields its statement does not have"), "{err}");
}

#[test]
fn a_protein_unchanged_naming_an_uncertainty_it_cannot_have_is_refused() {
    let json = p_json(
        r#"{"type": "Statement", "statement": "Unchanged", "uncertain": false}"#,
        false,
    );
    let err = from_json(&json).unwrap_err().to_string();
    assert!(err.contains("fields its statement does not have"), "{err}");
}
