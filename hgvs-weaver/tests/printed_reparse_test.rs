//! A description weaver prints parses back to the value it printed from.

use hgvs_weaver::{parse_hgvs_variant, SequenceVariant};

fn parse(s: &str) -> SequenceVariant {
    parse_hgvs_variant(s).unwrap_or_else(|e| panic!("{s}: {e}"))
}

/// The value `s` parses to, printed and parsed again.
fn reparsed(s: &str) -> SequenceVariant {
    let printed = parse(s).to_string();
    parse_hgvs_variant(&printed).unwrap_or_else(|e| panic!("{s} printed {printed}: {e}"))
}

/// One test per case: `$hgvs` parses, prints and parses back to the same value.
macro_rules! reparses_to_itself {
    ($($name:ident: $hgvs:literal;)*) => {$(
        #[test]
        fn $name() {
            assert_eq!(reparsed($hgvs), parse($hgvs));
        }
    )*};
}

reparses_to_itself! {
    a_frameshift_to_a_stop_at_an_unknown_distance: "NP_W.1:p.Arg97fsTer?";
    a_frameshift_to_a_stop_at_a_stated_distance: "NP_W.1:p.Arg97fsTer5";
    a_frameshift_naming_no_stop: "NP_W.1:p.Arg97Profs";
    a_frameshift_to_a_stop_at_no_stated_distance: "NP_W.1:p.Arg97Profs*";
    a_frameshift_count_with_leading_zeros: "NP_W.1:p.Arg97ProfsTer007";
    an_unknown_change_from_a_residue_past_met1: "NP_W.1:p.Arg12?";
    a_predicted_unknown_change_from_a_residue_past_met1: "NP_W.1:p.(Arg12?)";
}

reparses_to_itself! {
    a_protein_repeat_of_a_range_of_copies: "NP_W.1:p.Gln18(70_80)";
    a_coding_repeat_of_a_range_of_copies: "NM_W.1:c.10_12CAG(70_80)";
    a_genomic_repeat_of_a_range_of_copies: "NC_W.1:g.10_11TG(7_9)";
    an_rna_repeat_of_a_range_of_copies: "NM_W.1:r.10_12cag(70_80)";
}

#[test]
fn a_repeat_of_a_range_of_copies_parses_as_hgvs_writes_it() {
    assert_eq!(
        parse("NP_W.1:p.Gln18[(70_80)]"),
        parse("NP_W.1:p.Gln18(70_80)")
    );
}

#[test]
fn a_repeat_of_a_range_of_copies_prints_as_hgvs_writes_it() {
    assert_eq!(
        parse("NM_W.1:c.10_12CAG(70_80)").to_string(),
        "NM_W.1:c.10_12CAG[(70_80)]"
    );
}

reparses_to_itself! {
    a_one_letter_residue_at_a_position: "NP_W.1:p.R2G";
    a_one_letter_unknown_change_from_met1: "NP_W.1:p.M1?";
    a_one_letter_predicted_unknown_change_from_met1: "NP_W.1:p.(M1?)";
    a_one_letter_nonsense_change: "NP_W.1:p.R2*";
    a_one_letter_frameshift: "NP_W.1:p.R97Pfs*23";
    a_one_letter_frameshift_naming_no_new_residue: "NP_W.1:p.R97fs*5";
    a_frameshift_to_a_stop_written_x: "NP_W.1:p.Arg97ProfsX5";
    a_one_letter_extension: "NP_W.1:p.*30Qext*5";
    an_extension_to_a_stop_written_x: "NP_W.1:p.Ter30GlnextX5";
    an_ambiguous_residue_written_b: "NP_W.1:p.B2G";
    a_selenocysteine_written_u: "NP_W.1:p.U2G";
}

#[test]
fn a_one_letter_unknown_change_from_met1_is_the_three_letter_one() {
    assert_eq!(parse("NP_W.1:p.M1?"), parse("NP_W.1:p.Met1?"));
}

#[test]
fn a_frameshift_count_with_leading_zeros_is_the_number() {
    assert_eq!(
        parse("NP_W.1:p.Arg97ProfsTer007"),
        parse("NP_W.1:p.Arg97ProfsTer7")
    );
}

#[test]
fn a_one_letter_frameshift_is_the_three_letter_one() {
    assert_eq!(
        parse("NP_W.1:p.R97Pfs*23"),
        parse("NP_W.1:p.Arg97ProfsTer23")
    );
}

#[test]
fn a_one_letter_ambiguous_residue_prints_as_its_three_letter_code() {
    assert_eq!(parse("NP_W.1:p.B2Z").to_string(), "NP_W.1:p.Asx2Glx");
}

reparses_to_itself! {
    an_extension_from_a_stop_written_x: "NP_W.1:p.X30Qext*5";
    an_n_terminal_extension_from_an_unknown_residue: "NP_W.1:p.Met1extXaa-5";
    an_n_terminal_extension_from_a_named_residue: "NP_W.1:p.M1extM-5";
}

reparses_to_itself! {
    an_extension_with_nothing_after_ext: "NP_W.1:p.Ter30Glnext";
    an_extension_to_a_stop_at_no_stated_distance: "NP_W.1:p.Ter30Glnext*";
    an_extension_to_a_stop_at_an_unknown_distance: "NP_W.1:p.Ter30Glnext*?";
    an_n_terminal_extension: "NP_W.1:p.Met1ext-5";
    an_extension_by_a_bare_count: "NP_W.1:p.Ter30Glnext5";
    an_extension_by_a_count_after_a_residue: "NP_W.1:p.Ter30GlnextLeu5";
    an_extension_by_a_count_written_plus: "NP_W.1:p.Ter30Glnext+5";
    an_extension_by_a_count_written_plus_after_x: "NP_W.1:p.Ter30GlnextX+5";
    an_extension_count_with_leading_zeros: "NP_W.1:p.Ter30GlnextTer007";
}

#[test]
fn an_extension_count_written_plus_is_the_bare_count() {
    assert_eq!(
        parse("NP_W.1:p.Ter30Glnext+5"),
        parse("NP_W.1:p.Ter30Glnext5")
    );
}

#[test]
fn an_extension_count_with_leading_zeros_is_the_number() {
    assert_eq!(
        parse("NP_W.1:p.Ter30GlnextTer007"),
        parse("NP_W.1:p.Ter30GlnextTer7")
    );
}

#[test]
fn an_extension_count_after_x_is_the_distance_to_a_stop() {
    assert_eq!(
        parse("NP_W.1:p.Ter30GlnextX+5"),
        parse("NP_W.1:p.Ter30GlnextTer5")
    );
}

#[test]
fn an_extension_from_a_stop_written_x_is_the_one_written_ter() {
    assert_eq!(
        parse("NP_W.1:p.X30Qext*5"),
        parse("NP_W.1:p.Ter30GlnextTer5")
    );
}
