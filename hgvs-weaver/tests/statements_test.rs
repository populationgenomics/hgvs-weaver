//! Statements about a whole transcript or protein — `r.0`, `r.spl`, `r.?`,
//! `p.0`, `p.?`, `p.=` — read as what they mean, print as HGVS writes them,
//! and survive JSON, in their own shape and in the `Special` one.

mod support;

use hgvs_weaver::data::Strand;
use hgvs_weaver::edits::{AaEdit, NaEdit, ProteinStatement, RnaStatement};
use hgvs_weaver::equivalence::{EquivalenceLevel, VariantEquivalence};
use hgvs_weaver::error::HgvsError;
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
    SequenceVariant::Protein(try_c_to_p(s).unwrap_or_else(|e| panic!("{s}: {e}")))
}

fn try_c_to_p(s: &str) -> Result<hgvs_weaver::PVariant, HgvsError> {
    let SequenceVariant::Coding(c) = parse(s) else {
        panic!("{s} is not c.");
    };
    let hdp = provider();
    VariantMapper::new(&hdp).c_to_p(&c, None)
}

/// Whether `s` is refused as making no statement about the protein.
fn makes_no_protein_statement(s: &str) -> bool {
    matches!(try_c_to_p(s), Err(HgvsError::NoProteinStatement(_)))
}

/// Whether `s` fails to parse, rather than parsing and failing later.
fn is_parse_error(s: &str) -> bool {
    matches!(parse_hgvs_variant(s), Err(HgvsError::PestError(_)))
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
    r_follows_the_dna_is_its_own_statement: "NM_W.1:r.(?)" => NaEdit::Statement(RnaStatement::FollowsDna);
    r_spl_with_a_space_before_its_question_mark_is_uncertain: "NM_W.1:r.spl ?" => NaEdit::Statement(RnaStatement::SplicingAffected { uncertain: true });
}

reads_as! { aa_edit;
    p_0_is_no_protein: "NP_W.1:p.0" => AaEdit::Statement(ProteinStatement::NoProtein { uncertain: false });
    p_unknown_is_an_unknown_effect: "NP_W.1:p.?" => AaEdit::Statement(ProteinStatement::Unknown);
    p_identity_is_the_protein_unchanged: "NP_W.1:p.=" => AaEdit::Statement(ProteinStatement::Unchanged);
    predicted_p_identity_is_the_protein_unchanged: "NP_W.1:p.(=)" => AaEdit::Statement(ProteinStatement::Unchanged);
    p_0_uncertain_is_probably_no_protein: "NP_W.1:p.0?" => AaEdit::Statement(ProteinStatement::NoProtein { uncertain: true });
    predicted_p_unknown_is_an_unknown_effect: "NP_W.1:p.(?)" => AaEdit::Statement(ProteinStatement::Unknown);
    an_unknown_change_from_a_residue_is_an_unknown_effect: "NP_W.1:p.Arg12?" => AaEdit::Statement(ProteinStatement::Unknown);
}

#[test]
fn a_predicted_p_met1_unknown_parses_to_what_the_mapper_builds() {
    assert_eq!(parse("NP_W.1:p.(Met1?)"), c_to_p("NM_W.1:c.-3_4del"));
}

#[test]
fn a_start_loss_is_predicted_bracketed_as_hgvs_writes_it() {
    assert_eq!(c_to_p("NM_W.1:c.-3_4del").to_string(), "NP_W.1:p.(Met1?)");
}

#[test]
fn a_change_upstream_of_the_cds_is_an_unknown_effect_written_bare() {
    assert_eq!(c_to_p("NM_W.1:c.-3C>A").to_string(), "NP_W.1:p.?");
}

#[test]
fn a_predicted_change_is_written_bracketed() {
    assert_eq!(c_to_p("NM_W.1:c.10C>G").to_string(), "NP_W.1:p.(Arg4Gly)");
}

#[test]
fn an_unknown_change_over_a_range_does_not_parse() {
    assert!(is_parse_error("NP_W.1:p.Met1_Arg3?"));
}

#[test]
fn an_unknown_change_from_an_uncertain_residue_does_not_parse() {
    assert!(is_parse_error("NP_W.1:p.(Met1)?"));
}

#[test]
fn a_predicted_probably_no_protein_does_not_parse() {
    assert!(is_parse_error("NP_W.1:p.(0?)"));
}

#[test]
fn a_residue_numbered_0_does_not_parse() {
    assert!(is_parse_error("NP_W.1:p.Met0Val"));
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
fn predicted_p_unknown_prints_back_as_parsed() {
    assert_eq!(parse("NP_W.1:p.(?)").to_string(), "NP_W.1:p.(?)");
}

#[test]
fn predicted_p_unknown_is_the_unknown_effect_in_another_spelling() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let level = VariantEquivalence::new(&mapper, &hdp)
        .equivalent_level(&parse("NP_W.1:p.(?)"), &parse("NP_W.1:p.?"))
        .unwrap();
    assert_eq!(level, EquivalenceLevel::Identity);
}

#[test]
fn predicted_p_unknown_does_not_match_a_coding_change() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let level = VariantEquivalence::new(&mapper, &hdp)
        .equivalent_level(&parse("NP_W.1:p.(?)"), &parse("NM_W.1:c.10C>G"))
        .unwrap();
    assert_eq!(level, EquivalenceLevel::Different);
}

#[test]
fn r_follows_the_dna_prints_back_as_parsed() {
    assert_eq!(parse("NM_W.1:r.(?)").to_string(), "NM_W.1:r.(?)");
}

#[test]
fn a_dict_for_r_follows_the_dna_names_its_statement() {
    let json = r_json(r#"{"type": "Statement", "statement": "FollowsDna"}"#);
    assert_eq!(from_json(&json).unwrap(), parse("NM_W.1:r.(?)"));
}

#[test]
fn a_dict_for_r_unknown_marked_predicted_is_refused() {
    let json = r#"{"variant_type": "Rna", "ac": "NM_W.1", "gene": null, "posedit": {"pos": null,
        "edit": {"type": "Statement", "statement": "Unknown"}, "uncertain": false, "predicted": true}}"#;
    let err = from_json(json).unwrap_err().to_string();
    assert!(err.contains("never predicted"), "{err}");
}

#[test]
fn a_dict_for_a_special_unknown_marked_predicted_is_refused() {
    let json = r#"{"variant_type": "Rna", "ac": "NM_W.1", "gene": null, "posedit": {"pos": null,
        "edit": {"type": "Special", "value": "?", "uncertain": false}, "uncertain": false, "predicted": true}}"#;
    let err = from_json(json).unwrap_err().to_string();
    assert!(err.contains("never predicted"), "{err}");
}

#[test]
fn r_follows_the_dna_has_no_protein_consequence_of_its_own() {
    let SequenceVariant::Rna(r) = parse("NM_W.1:r.(?)") else {
        panic!("not r.");
    };
    let hdp = provider();
    let result = VariantMapper::new(&hdp).r_to_p(&r, None);
    assert!(
        matches!(result, Err(HgvsError::UnsupportedOperation(_))),
        "{result:?}"
    );
}

#[test]
fn r_follows_the_dna_is_not_the_unknown_effect() {
    assert_eq!(
        level("NM_W.1:r.(?)", "NM_W.1:r.?"),
        EquivalenceLevel::Different
    );
}

#[test]
fn r_follows_the_dna_is_its_one_member_cis_allele() {
    assert_eq!(
        level("NM_W.1:r.[(?)]", "NM_W.1:r.(?)"),
        EquivalenceLevel::Identity
    );
}

#[test]
fn a_one_member_cis_allele_of_r_follows_the_dna_is_not_the_unknown_effect() {
    assert_eq!(
        level("NM_W.1:r.[(?)]", "NM_W.1:r.?"),
        EquivalenceLevel::Different
    );
}

#[test]
fn one_member_cis_alleles_of_r_follows_the_dna_and_the_unknown_effect_differ() {
    assert_eq!(
        level("NM_W.1:r.[(?)]", "NM_W.1:r.[?]"),
        EquivalenceLevel::Different
    );
}

#[test]
fn a_cis_allele_holding_r_follows_the_dna_is_not_one_holding_the_unknown_effect() {
    assert_eq!(
        level("NM_W.1:r.[(?);10c>g]", "NM_W.1:r.[?;10c>g]"),
        EquivalenceLevel::Different
    );
}

#[test]
fn r_follows_the_dna_is_itself() {
    assert_eq!(
        level("NM_W.1:r.(?)", "NM_W.1:r.(?)"),
        EquivalenceLevel::Identity
    );
}

/// The level `a` and `b` are judged at over the provider.
fn level(a: &str, b: &str) -> EquivalenceLevel {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    VariantEquivalence::new(&mapper, &hdp)
        .equivalent_level(&parse(a), &parse(b))
        .unwrap()
}

#[test]
fn an_unknown_change_from_a_residue_matches_a_change_starting_there() {
    assert_eq!(
        level("NP_W.1:p.Arg12?", "NP_W.1:p.Arg12Gly"),
        EquivalenceLevel::Analogous
    );
}

#[test]
fn an_unknown_change_from_a_residue_does_not_match_a_change_starting_after_it() {
    assert_eq!(
        level("NP_W.1:p.Arg12?", "NP_W.1:p.Arg13Gly"),
        EquivalenceLevel::Different
    );
}

#[test]
fn an_unknown_change_from_a_residue_matches_a_coding_change_predicted_there() {
    assert_eq!(
        level("NM_W.1:c.34C>G", "NP_W.1:p.Arg12?"),
        EquivalenceLevel::Analogous
    );
}

#[test]
fn a_predicted_unknown_change_from_a_residue_is_the_stated_one() {
    assert_eq!(
        level("NP_W.1:p.(Arg12?)", "NP_W.1:p.Arg12?"),
        EquivalenceLevel::Identity
    );
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
fn an_unknown_change_from_an_unknown_residue_prints_back_as_parsed() {
    assert_eq!(parse("NP_W.1:p.Xaa12?").to_string(), "NP_W.1:p.Xaa12?");
}

#[test]
fn a_substitution_of_an_unknown_residue_prints_back_as_parsed() {
    assert_eq!(parse("NP_W.1:p.Xaa12Val").to_string(), "NP_W.1:p.Xaa12Val");
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
fn the_mappers_unknown_effect_reparses_to_itself() {
    let p = c_to_p("NM_W.1:c.45+1C>A");
    assert_eq!(parse(&p.to_string()), p);
}

#[test]
fn the_start_codon_transform_of_a_one_letter_description_reparses_to_itself() {
    let settings = VariantTransformSettings {
        start_codon: StartCodonConvention::HgvsQuestion,
    };
    let p = transform_variant(&parse("NP_W.1:p.(M1V)"), &settings);
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
fn a_start_loss_is_the_parsed_predicted_unknown_change_from_met1() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let level = VariantEquivalence::new(&mapper, &hdp)
        .equivalent_level(&parse("NM_W.1:c.-3_4del"), &parse("NP_W.1:p.(Met1?)"))
        .unwrap();
    assert_eq!(level, EquivalenceLevel::Identity);
}

#[test]
fn a_start_loss_is_analogous_to_the_unknown_change_from_met1_stated() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let level = VariantEquivalence::new(&mapper, &hdp)
        .equivalent_level(&parse("NM_W.1:c.-3_4del"), &parse("NP_W.1:p.Met1?"))
        .unwrap();
    assert_eq!(level, EquivalenceLevel::Analogous);
}

#[test]
fn the_protein_unchanged_past_the_cds_reparses_to_itself() {
    let p = c_to_p("NM_W.1:c.92T>A");
    assert_eq!(parse(&p.to_string()), p);
}

#[test]
fn the_protein_unchanged_past_the_cds_survives_json() {
    let p = c_to_p("NM_W.1:c.92T>A");
    assert_eq!(json_round_trip(&p), p);
}

#[test]
fn an_unknown_change_from_met1_is_not_the_same_text_as_an_unknown_residue_there() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let level = VariantEquivalence::new(&mapper, &hdp)
        .equivalent_level(&parse("NP_W.1:p.Met1?"), &parse("NP_W.1:p.Met1Xaa"))
        .unwrap();
    assert_eq!(level, EquivalenceLevel::Analogous);
}

#[test]
fn a_3_utr_change_predicts_the_protein_unchanged_and_reparses_to_it() {
    let p = c_to_p("NM_W.1:c.*5T>A");
    assert_eq!(p.to_string(), "NP_W.1:p.(=)");
    assert_eq!(parse(&p.to_string()), p);
}

#[test]
fn an_identity_across_the_whole_transcript_predicts_the_protein_unchanged() {
    assert_eq!(c_to_p("NM_W.1:c.-5_*15="), parse("NP_W.1:p.(=)"));
}

#[test]
fn an_identity_from_the_5_utr_into_the_cds_is_silent_at_the_codons_it_covers() {
    assert_eq!(c_to_p("NM_W.1:c.-2_3="), parse("NP_W.1:p.(Met1=)"));
}

#[test]
fn an_identity_ending_in_an_intron_is_silent_at_the_whole_codons_it_covers() {
    assert_eq!(c_to_p("NM_W.1:c.43_45+1="), parse("NP_W.1:p.(Arg15=)"));
}

#[test]
fn an_identity_starting_in_an_intron_is_silent_at_the_whole_codons_it_covers() {
    assert_eq!(c_to_p("NM_W.1:c.45+2_48="), parse("NP_W.1:p.(Val16=)"));
}

#[test]
fn an_identity_from_an_intron_in_the_5_utr_is_silent_at_the_whole_codons_it_covers() {
    // NM_U.1 is NM_W.1 with a one-base intron, genome base 12, between c.-3
    // and c.-2.
    let mut g: Vec<u8> = b"CGT".iter().copied().cycle().take(200).collect();
    g[15..18].copy_from_slice(b"ATG");
    g[30..34].copy_from_slice(b"AAAA");
    g[142..145].copy_from_slice(b"TAA");
    let g = String::from_utf8(g).unwrap();
    let spliced = format!("{}{}{}", &g[10..12], &g[13..60], &g[100..160]);
    let hdp = Provider::new()
        .sequence("NC_W.1", &g)
        .sequence("NM_U.1", &spliced)
        .sequence("NP_W.1", "MRRRRKSRRRRRRRRVVVVVVVVVVVVVV")
        .transcript(transcript(
            "NM_U.1",
            "NC_W.1",
            Strand::Plus,
            Some((4, 93)),
            vec![
                exon((0, 2), (10, 11), Strand::Plus),
                exon((2, 49), (13, 59), Strand::Plus),
                exon((49, 109), (100, 159), Strand::Plus),
            ],
        ))
        .protein_for("NM_U.1", "NP_W.1");
    let SequenceVariant::Coding(c) = parse("NM_U.1:c.-2-1_3=") else {
        panic!("not c.");
    };
    let p = VariantMapper::new(&hdp).c_to_p(&c, None).unwrap();
    assert_eq!(SequenceVariant::Protein(p), parse("NP_W.1:p.(Met1=)"));
}

#[test]
fn an_identity_covering_no_whole_codon_makes_no_protein_statement() {
    assert!(makes_no_protein_statement("NM_W.1:c.-2_2="));
}

#[test]
fn an_identity_inside_one_codon_makes_no_protein_statement() {
    assert!(makes_no_protein_statement("NM_W.1:c.5_6="));
}

#[test]
fn one_base_unchanged_makes_no_protein_statement() {
    assert!(makes_no_protein_statement("NM_W.1:c.5="));
}

#[test]
fn an_identity_of_every_codon_predicts_the_protein_unchanged() {
    assert_eq!(c_to_p("NM_W.1:c.1_90="), parse("NP_W.1:p.(=)"));
}

#[test]
fn an_identity_of_every_codon_but_the_stop_names_its_residues() {
    assert_eq!(c_to_p("NM_W.1:c.1_87="), parse("NP_W.1:p.(Met1_Val29=)"));
}

#[test]
fn an_identity_in_the_5_utr_makes_no_protein_statement() {
    assert!(makes_no_protein_statement("NM_W.1:c.-5_-1="));
}

#[test]
fn an_identity_in_the_3_utr_makes_no_protein_statement() {
    assert!(makes_no_protein_statement("NM_W.1:c.*3="));
}

#[test]
fn an_intronic_identity_makes_no_protein_statement() {
    assert!(makes_no_protein_statement("NM_W.1:c.45+1="));
}

#[test]
fn an_intronic_conversion_is_an_unknown_effect() {
    assert_eq!(
        c_to_p("NM_W.1:c.45+1_45+3conNM_W.1:c.10_12"),
        parse("NP_W.1:p.?")
    );
}

#[test]
fn an_upstream_conversion_is_an_unknown_effect() {
    assert_eq!(
        c_to_p("NM_W.1:c.-5_-3conNM_W.1:c.10_12"),
        parse("NP_W.1:p.?")
    );
}

#[test]
fn an_identity_inside_the_codons_is_written_silently_at_its_residues() {
    assert_eq!(c_to_p("NM_W.1:c.4_6="), parse("NP_W.1:p.(Arg2=)"));
}

#[test]
fn an_identity_across_the_whole_transcript_is_not_no_protein() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let level = VariantEquivalence::new(&mapper, &hdp)
        .equivalent_level(&parse("NM_W.1:c.-5_*15="), &parse("NP_W.1:p.0"))
        .unwrap();
    assert_eq!(level, EquivalenceLevel::Different);
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
    predicted_p_unknown_survives_json: "NP_W.1:p.(?)";
    r_follows_the_dna_survives_json: "NM_W.1:r.(?)";
    an_unknown_change_from_a_residue_survives_json: "NP_W.1:p.Arg12?";
    a_predicted_unknown_change_from_a_residue_survives_json: "NP_W.1:p.(Arg12?)";
}

#[test]
fn probably_no_protein_survives_json() {
    let p = c_to_p("NM_W.1:c.-3_*2del");
    assert_eq!(json_round_trip(&p), p);
}

#[test]
fn an_unknown_change_from_met1_survives_json() {
    let p = c_to_p("NM_W.1:c.-3_4del");
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
fn a_special_unknown_at_met1_loads_as_the_unknown_change_from_met1() {
    let json = r#"{"variant_type": "Protein", "ac": "NP_W.1", "gene": null, "posedit": {
        "pos": {"start": {"base": 1, "aa": "Met", "uncertain": false}, "end": null, "uncertain": false},
        "edit": {"type": "Special", "value": "?", "uncertain": false}, "uncertain": false, "predicted": false}}"#;
    assert_eq!(from_json(json).unwrap(), parse("NP_W.1:p.Met1?"));
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

// --- A p. dict describes only what p. writes ---

/// A p. variant on NP_W.1 with posedit fields `pos`, `edit`, `uncertain` and
/// `predicted`, as JSON.
fn p_json_at(pos: &str, edit: &str, uncertain: bool, predicted: bool) -> String {
    format!(
        r#"{{"variant_type": "Protein", "ac": "NP_W.1", "gene": null, "posedit": {{"pos": {pos},
            "edit": {edit}, "uncertain": {uncertain}, "predicted": {predicted}}}}}"#
    )
}

/// Whether `json` was refused because p. cannot write what it describes.
fn refused_by_p(json: &str) -> bool {
    from_json(json).is_err_and(|e| e.to_string().contains("p. cannot write"))
}

#[test]
fn a_dict_for_no_protein_at_a_residue_is_refused() {
    let json = p_json_at(
        r#"{"start": {"base": 1, "aa": "Met", "uncertain": false}, "end": null, "uncertain": false}"#,
        r#"{"type": "Statement", "statement": "NoProtein", "uncertain": false}"#,
        false,
        false,
    );
    assert!(refused_by_p(&json), "{json}");
}

#[test]
fn a_dict_for_the_protein_unchanged_at_a_residue_is_refused() {
    let json = p_json_at(
        r#"{"start": {"base": 2, "aa": "Arg", "uncertain": false}, "end": null, "uncertain": false}"#,
        r#"{"type": "Statement", "statement": "Unchanged"}"#,
        false,
        true,
    );
    assert!(refused_by_p(&json), "{json}");
}

#[test]
fn a_dict_for_an_unknown_change_over_a_range_is_refused() {
    let json = p_json_at(
        r#"{"start": {"base": 1, "aa": "Met", "uncertain": false},
            "end": {"base": 3, "aa": "Arg", "uncertain": false}, "uncertain": false}"#,
        r#"{"type": "Statement", "statement": "Unknown"}"#,
        false,
        false,
    );
    assert!(refused_by_p(&json), "{json}");
}

#[test]
fn a_predicted_no_protein_does_not_parse() {
    assert!(is_parse_error("NP_W.1:p.(0)"));
}

#[test]
fn a_dict_for_a_predicted_no_protein_is_refused() {
    let json = p_json_at(
        "null",
        r#"{"type": "Statement", "statement": "NoProtein", "uncertain": false}"#,
        false,
        true,
    );
    assert!(refused_by_p(&json), "{json}");
}

#[test]
fn a_dict_for_a_predicted_unknown_effect_loads_as_it() {
    let json = p_json_at(
        "null",
        r#"{"type": "Statement", "statement": "Unknown"}"#,
        false,
        true,
    );
    assert_eq!(from_json(&json).unwrap(), parse("NP_W.1:p.(?)"));
}

#[test]
fn a_dict_for_a_change_with_no_position_is_refused() {
    let json = p_json_at(
        "null",
        r#"{"type": "Subst", "ref_": "", "alt": "Gly", "uncertain": false}"#,
        false,
        false,
    );
    assert!(refused_by_p(&json), "{json}");
}

#[test]
fn a_dict_for_a_change_at_residue_0_is_refused() {
    let json = p_json_at(
        r#"{"start": {"base": 0, "aa": "", "uncertain": false}, "end": null, "uncertain": false}"#,
        r#"{"type": "Identity", "uncertain": false}"#,
        false,
        false,
    );
    assert!(refused_by_p(&json), "{json}");
}

#[test]
fn a_dict_for_an_uncertain_p_posedit_is_refused() {
    let json = p_json_at(
        r#"{"start": {"base": 2, "aa": "Arg", "uncertain": false}, "end": null, "uncertain": false}"#,
        r#"{"type": "Subst", "ref_": "", "alt": "Gly", "uncertain": false}"#,
        true,
        false,
    );
    assert!(refused_by_p(&json), "{json}");
}

#[test]
fn a_dict_naming_no_edit_is_refused() {
    let json = p_json_at(
        r#"{"start": {"base": 2, "aa": "Arg", "uncertain": false}, "end": null, "uncertain": false}"#,
        r#"{"type": "None"}"#,
        false,
        false,
    );
    assert!(refused_by_p(&json), "{json}");
}

#[test]
fn a_dict_for_a_change_at_a_residue_still_loads() {
    let json = p_json_at(
        r#"{"start": {"base": 2, "aa": "Arg", "uncertain": false}, "end": null, "uncertain": false}"#,
        r#"{"type": "Subst", "ref_": "", "alt": "Gly", "uncertain": false}"#,
        false,
        true,
    );
    assert_eq!(from_json(&json).unwrap(), parse("NP_W.1:p.(Arg2Gly)"));
}

#[test]
fn a_dict_for_a_residue_with_no_uncertain_flag_loads() {
    let json = p_json_at(
        r#"{"start": {"base": 2, "aa": "Arg"}, "end": null, "uncertain": false}"#,
        r#"{"type": "Subst", "ref_": "", "alt": "Gly", "uncertain": false}"#,
        false,
        true,
    );
    assert_eq!(from_json(&json).unwrap(), parse("NP_W.1:p.(Arg2Gly)"));
}

#[test]
fn a_dict_in_one_letter_code_loads_as_its_printed_description() {
    let json = p_json_at(
        r#"{"start": {"base": 2, "aa": "R"}, "end": null, "uncertain": false}"#,
        r#"{"type": "Fs", "alt": "P", "stop": {"kind": "Residues", "count": 5}}"#,
        false,
        false,
    );
    assert_eq!(from_json(&json).unwrap(), parse("NP_W.1:p.Arg2ProfsTer5"));
}

#[test]
fn a_dict_for_a_frameshift_naming_no_stop_loads_as_the_bare_frameshift() {
    let json = p_json_at(
        r#"{"start": {"base": 2, "aa": "Arg"}, "end": null, "uncertain": false}"#,
        r#"{"type": "Fs", "alt": "Pro", "stop": null}"#,
        false,
        false,
    );
    assert_eq!(from_json(&json).unwrap(), parse("NP_W.1:p.Arg2Profs"));
}

#[test]
fn a_dict_for_a_frameshift_is_the_variant_its_printed_text_parses_to() {
    let json = p_json_at(
        r#"{"start": {"base": 2, "aa": "Arg"}, "end": null, "uncertain": false}"#,
        r#"{"type": "Fs", "ref_": "Arg", "alt": "Pro", "stop": null, "uncertain": true}"#,
        false,
        false,
    );
    let loaded = from_json(&json).unwrap();
    assert_eq!(parse(&loaded.to_string()), loaded);
}

#[test]
fn a_dict_for_a_frameshift_with_its_stop_as_text_is_refused() {
    let json = p_json_at(
        r#"{"start": {"base": 2, "aa": "Arg"}, "end": null, "uncertain": false}"#,
        r#"{"type": "Fs", "ref_": "", "alt": "Pro", "term": "Ter", "length": "5", "uncertain": false}"#,
        false,
        false,
    );
    assert!(
        from_json(&json).is_err_and(|e| e.to_string().contains("missing field `stop`")),
        "{json}"
    );
}

#[test]
fn a_dict_for_an_extension_from_a_stop_written_x_loads_at_ter() {
    let json = p_json_at(
        r#"{"start": {"base": 30, "aa": "X"}, "end": null, "uncertain": false}"#,
        r#"{"type": "Ext", "alt": "Q", "end": {"kind": "Stop", "distance": {"kind": "Residues", "count": 5}}}"#,
        false,
        false,
    );
    assert_eq!(from_json(&json).unwrap(), parse("NP_W.1:p.Ter30GlnextTer5"));
}

#[test]
fn a_dict_for_an_extension_by_a_count_after_x_loads_as_the_distance_to_a_stop() {
    let json = p_json_at(
        r#"{"start": {"base": 30, "aa": "Ter"}, "end": null, "uncertain": false}"#,
        r#"{"type": "Ext", "alt": "Gln", "end": {"kind": "Count", "residue": "X", "count": 5}}"#,
        false,
        false,
    );
    assert_eq!(from_json(&json).unwrap(), parse("NP_W.1:p.Ter30GlnextTer5"));
}

#[test]
fn a_dict_for_an_n_terminal_extension_from_x_loads_from_an_unknown_residue() {
    let json = p_json_at(
        r#"{"start": {"base": 1, "aa": "M"}, "end": null, "uncertain": false}"#,
        r#"{"type": "Ext", "alt": "", "end": {"kind": "Upstream", "residue": "X", "distance": 5}}"#,
        false,
        false,
    );
    assert_eq!(from_json(&json).unwrap(), parse("NP_W.1:p.Met1extXaa-5"));
}

#[test]
fn a_dict_for_an_n_terminal_extension_to_a_stop_is_refused() {
    let json = p_json_at(
        r#"{"start": {"base": 1, "aa": "Met"}, "end": null, "uncertain": false}"#,
        r#"{"type": "Ext", "alt": "", "end": {"kind": "Upstream", "residue": "*", "distance": 5}}"#,
        false,
        false,
    );
    assert!(refused_by_p(&json), "{json}");
}

#[test]
fn a_dict_for_an_n_terminal_extension_naming_an_empty_residue_is_refused() {
    let json = p_json_at(
        r#"{"start": {"base": 1, "aa": "Met"}, "end": null, "uncertain": false}"#,
        r#"{"type": "Ext", "alt": "", "end": {"kind": "Upstream", "residue": "", "distance": 5}}"#,
        false,
        false,
    );
    assert!(refused_by_p(&json), "{json}");
}

#[test]
fn a_dict_for_an_n_terminal_extension_naming_a_letter_no_residue_has_is_refused() {
    let json = p_json_at(
        r#"{"start": {"base": 1, "aa": "Met"}, "end": null, "uncertain": false}"#,
        r#"{"type": "Ext", "alt": "", "end": {"kind": "Upstream", "residue": "J", "distance": 5}}"#,
        false,
        false,
    );
    assert!(refused_by_p(&json), "{json}");
}

#[test]
fn a_dict_for_an_extension_by_a_count_after_two_residues_is_refused() {
    let json = p_json_at(
        r#"{"start": {"base": 30, "aa": "Ter"}, "end": null, "uncertain": false}"#,
        r#"{"type": "Ext", "alt": "Gln", "end": {"kind": "Count", "residue": "LeuSer", "count": 5}}"#,
        false,
        false,
    );
    assert!(refused_by_p(&json), "{json}");
}

#[test]
fn a_dict_for_an_extension_with_its_end_as_text_is_refused() {
    let json = p_json_at(
        r#"{"start": {"base": 30, "aa": "Ter"}, "end": null, "uncertain": false}"#,
        r#"{"type": "Ext", "alt": "Gln", "aaterm": "Ter", "length": "5"}"#,
        false,
        false,
    );
    assert!(
        from_json(&json).is_err_and(|e| e.to_string().contains("missing field `end`")),
        "{json}"
    );
}

// --- Older p. dicts give p.Met1? and p.0? as a substitution to ? ---

#[test]
fn an_older_dict_for_p_met1_unknown_loads_as_the_unknown_change_from_met1() {
    let json = p_json_at(
        r#"{"start": {"base": 1, "aa": "Met", "uncertain": false}, "end": null, "uncertain": false}"#,
        r#"{"type": "Subst", "ref_": "", "alt": "?", "uncertain": false}"#,
        false,
        false,
    );
    assert_eq!(from_json(&json).unwrap(), parse("NP_W.1:p.Met1?"));
}

#[test]
fn an_older_dict_for_a_predicted_p_met1_unknown_loads_as_its_prediction() {
    let json = p_json_at(
        r#"{"start": {"base": 1, "aa": "Met", "uncertain": false}, "end": null, "uncertain": false}"#,
        r#"{"type": "Subst", "ref_": "", "alt": "?", "uncertain": false}"#,
        false,
        true,
    );
    assert_eq!(from_json(&json).unwrap(), parse("NP_W.1:p.(Met1?)"));
}

#[test]
fn an_older_dict_for_p_0_uncertain_loads_as_probably_no_protein() {
    let json = p_json_at(
        r#"{"start": {"base": 0, "aa": "", "uncertain": false}, "end": null, "uncertain": false}"#,
        r#"{"type": "Subst", "ref_": "", "alt": "?", "uncertain": false}"#,
        false,
        false,
    );
    assert_eq!(from_json(&json).unwrap(), parse("NP_W.1:p.0?"));
}

#[test]
fn an_older_dict_for_a_predicted_p_0_uncertain_is_refused() {
    let json = p_json_at(
        r#"{"start": {"base": 0, "aa": "", "uncertain": false}, "end": null, "uncertain": false}"#,
        r#"{"type": "Subst", "ref_": "", "alt": "?", "uncertain": false}"#,
        false,
        true,
    );
    assert!(refused_by_p(&json), "{json}");
}

#[test]
fn an_older_dict_for_an_unknown_substitution_over_a_range_is_refused() {
    let json = p_json_at(
        r#"{"start": {"base": 1, "aa": "Met", "uncertain": false},
            "end": {"base": 3, "aa": "Arg", "uncertain": false}, "uncertain": false}"#,
        r#"{"type": "Subst", "ref_": "", "alt": "?", "uncertain": false}"#,
        false,
        false,
    );
    assert!(refused_by_p(&json), "{json}");
}

#[test]
fn a_special_dict_for_the_protein_unchanged_at_a_residue_is_refused() {
    let json = p_json_at(
        r#"{"start": {"base": 2, "aa": "Arg", "uncertain": false}, "end": null, "uncertain": false}"#,
        r#"{"type": "Special", "value": "=", "uncertain": false}"#,
        false,
        false,
    );
    assert!(refused_by_p(&json), "{json}");
}

#[test]
fn an_older_dict_for_an_unknown_substitution_with_no_position_is_refused() {
    let json = p_json_at(
        "null",
        r#"{"type": "Subst", "ref_": "", "alt": "?", "uncertain": false}"#,
        false,
        false,
    );
    assert!(refused_by_p(&json), "{json}");
}

#[test]
fn an_older_dict_for_an_unknown_substitution_at_an_uncertain_residue_0_is_refused() {
    let json = p_json_at(
        r#"{"start": {"base": 0, "aa": "", "uncertain": true}, "end": null, "uncertain": false}"#,
        r#"{"type": "Subst", "ref_": "", "alt": "?", "uncertain": false}"#,
        false,
        false,
    );
    assert!(refused_by_p(&json), "{json}");
}

#[test]
fn a_dict_for_an_unknown_change_from_an_uncertain_residue_is_refused() {
    let json = p_json_at(
        r#"{"start": {"base": 1, "aa": "Met", "uncertain": true}, "end": null, "uncertain": false}"#,
        r#"{"type": "Statement", "statement": "Unknown"}"#,
        false,
        false,
    );
    assert!(refused_by_p(&json), "{json}");
}

#[test]
fn a_dict_for_a_change_at_an_uncertain_residue_is_refused() {
    let json = p_json_at(
        r#"{"start": {"base": 2, "aa": "Arg", "uncertain": true}, "end": null, "uncertain": false}"#,
        r#"{"type": "Subst", "ref_": "", "alt": "Gly", "uncertain": false}"#,
        false,
        true,
    );
    assert!(refused_by_p(&json), "{json}");
}

#[test]
fn an_older_dict_for_the_protein_unchanged_past_the_stop_loads_as_its_statement() {
    let json = p_json_at(
        "null",
        r#"{"type": "Identity", "uncertain": false}"#,
        false,
        true,
    );
    assert_eq!(from_json(&json).unwrap(), parse("NP_W.1:p.(=)"));
}
