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
