//! A protein description writes every residue in one amino acid code, the
//! one-letter or the three-letter; `*`, a stop in both, belongs to either.

use hgvs_weaver::error::HgvsError;
use hgvs_weaver::{parse_hgvs_variant, SequenceVariant};

fn parse(s: &str) -> SequenceVariant {
    parse_hgvs_variant(s).unwrap_or_else(|e| panic!("{s}: {e}"))
}

/// One test per case: `$hgvs` is refused for mixing the codes.
macro_rules! refused_as_mixed {
    ($($name:ident: $hgvs:literal;)*) => {$(
        #[test]
        fn $name() {
            let result = parse_hgvs_variant($hgvs);
            assert!(
                matches!(result, Err(HgvsError::MixedResidueCodes { .. })),
                "{}: {result:?}",
                $hgvs
            );
        }
    )*};
}

refused_as_mixed! {
    a_substitution_to_x_after_a_three_letter_residue: "NP_060360.3:p.Leu2781X";
    a_frameshift_to_x_after_a_three_letter_residue: "NP_060360.3:p.Leu2781fsX5";
    a_frameshift_to_a_three_letter_residue_and_x: "NP_060360.3:p.Leu2781ProfsX5";
    a_nonsense_change_written_x_as_clinvar_has_it: "NP_000021.1:p.Arg111X";
    a_frameshift_written_fsx_as_clinvar_has_it: "NP_004619.3:p.Ile223MetfsX45";
    a_one_letter_insertion_between_three_letter_residues: "NP_000081.1:p.Ala356_Gly357insDPA";
    a_one_letter_deletion_insertion_at_a_three_letter_residue: "NP_W.1:p.Arg2delinsKL";
    a_one_letter_repeat_unit_over_three_letter_residues: "NP_W.1:p.Pro490_Ser492PRS[1]";
    a_three_letter_repeat_unit_over_one_letter_residues: "NP_W.1:p.Q18_A20GlnAla[3]";
    a_three_letter_insertion_between_one_letter_residues: "NP_W.1:p.A2_G3insTerGlu";
    a_three_letter_substitution_at_a_one_letter_residue: "NP_W.1:p.R2Lys";
    a_range_whose_ends_use_different_codes: "NP_W.1:p.Lys2_L4del";
    a_one_letter_frameshift_to_a_stop_written_ter: "NP_005206.2:p.V191GfsTer35";
    an_extension_to_a_stop_written_x_from_three_letter_residues: "NP_W.1:p.*30GlnextX5";
    an_n_terminal_extension_from_a_three_letter_residue: "NP_W.1:p.M1extLeu-5";
    a_predicted_change_mixing_the_codes: "NP_W.1:p.(Arg2X)";
    an_allele_whose_members_use_different_codes: "NP_W.1:p.[Arg2Lys;L3P]";
    a_predicted_allele_whose_members_use_different_codes: "NP_W.1:p.[(Arg2Lys;L3P)]";
}

/// The residue or run written in each code that the refusal of `hgvs` names.
fn named(hgvs: &str) -> (String, String) {
    match parse_hgvs_variant(hgvs) {
        Err(HgvsError::MixedResidueCodes {
            one_letter,
            three_letter,
            ..
        }) => (one_letter, three_letter),
        other => panic!("{hgvs}: {other:?}"),
    }
}

#[test]
fn the_refusal_names_a_lone_x_and_the_three_letter_residue() {
    assert_eq!(
        named("NP_000021.1:p.Arg111X"),
        ("X".to_string(), "Arg".to_string())
    );
}

#[test]
fn the_refusal_names_a_one_letter_run_as_written() {
    assert_eq!(
        named("NP_000081.1:p.Ala356_Gly357insDPA"),
        ("DPA".to_string(), "Ala".to_string())
    );
}

#[test]
fn the_refusal_names_an_x_before_an_earlier_one_letter_residue() {
    assert_eq!(
        named("NP_W.1:p.Lys2_L4delinsX"),
        ("X".to_string(), "Lys".to_string())
    );
}

#[test]
fn the_refusal_of_an_x_says_how_a_stop_and_an_unknown_residue_are_written() {
    let message = parse_hgvs_variant("NP_000021.1:p.Arg111X")
        .unwrap_err()
        .to_string();
    assert!(message.contains("a stop is Ter or *"), "{message}");
}

#[test]
fn the_refusal_names_the_description() {
    let err = parse_hgvs_variant("NP_000021.1:p.Arg111X").unwrap_err();
    assert!(
        matches!(&err, HgvsError::MixedResidueCodes { description, .. }
            if description == "NP_000021.1:p.Arg111X"),
        "{err:?}"
    );
}

/// One test per case: `$hgvs`, with a stop written `*`, is the description
/// `$ter` that writes it `Ter`.
macro_rules! a_star_is_ter {
    ($($name:ident: $hgvs:literal => $ter:literal;)*) => {$(
        #[test]
        fn $name() {
            assert_eq!(parse($hgvs), parse($ter));
        }
    )*};
}

a_star_is_ter! {
    a_three_letter_nonsense_change_written_star: "NP_W.1:p.Arg2*" => "NP_W.1:p.Arg2Ter";
    a_one_letter_nonsense_change_written_star: "NP_W.1:p.R2*" => "NP_W.1:p.Arg2Ter";
    a_three_letter_frameshift_to_a_stop_written_star: "NP_W.1:p.Arg97Profs*23" => "NP_W.1:p.Arg97ProfsTer23";
    a_one_letter_frameshift_to_a_stop_written_star: "NP_W.1:p.R97Pfs*23" => "NP_W.1:p.Arg97ProfsTer23";
    a_three_letter_extension_from_a_stop_written_star: "NP_W.1:p.*30Glnext*5" => "NP_W.1:p.Ter30GlnextTer5";
    a_one_letter_extension_from_a_stop_written_star: "NP_W.1:p.*30Qext*5" => "NP_W.1:p.Ter30GlnextTer5";
    a_stop_written_star_inserted_between_three_letter_residues: "NP_W.1:p.Ala2_Gly3delins*" => "NP_W.1:p.Ala2_Gly3delinsTer";
    a_stop_written_star_inserted_between_one_letter_residues: "NP_W.1:p.A2_G3delins*" => "NP_W.1:p.Ala2_Gly3delinsTer";
}

#[test]
fn a_one_letter_frameshift_to_x_is_to_a_stop() {
    assert_eq!(
        parse("NP_004619.3:p.I223MfsX45"),
        parse("NP_004619.3:p.Ile223MetfsTer45")
    );
}

#[test]
fn a_one_letter_substitution_to_x_is_to_an_unknown_residue() {
    assert_eq!(
        parse("NP_060360.3:p.L2781X"),
        parse("NP_060360.3:p.Leu2781Xaa")
    );
}

/// One test per case: `$hgvs`, which writes no residue or writes them in one
/// code, still parses.
macro_rules! parses {
    ($($name:ident: $hgvs:literal;)*) => {$(
        #[test]
        fn $name() {
            let result = parse_hgvs_variant($hgvs);
            assert!(result.is_ok(), "{}: {result:?}", $hgvs);
        }
    )*};
}

parses! {
    an_unknown_effect: "NP_W.1:p.?";
    a_predicted_unknown_effect: "NP_W.1:p.(?)";
    no_protein: "NP_W.1:p.0";
    probably_no_protein: "NP_W.1:p.0?";
    the_protein_unchanged: "NP_W.1:p.=";
    the_protein_predicted_unchanged: "NP_W.1:p.(=)";
    a_one_letter_repeat_at_a_bare_number: "NP_W.1:p.490PRS[1]";
    a_three_letter_insertion_of_xaa: "NP_W.1:p.Ala2_Gly3insXaa";
}

/// One test per case: `$hgvs`, whose residues are in one-letter code,
/// prints `$printed`, in three-letter code throughout.
macro_rules! prints_in_three_letter_code {
    ($($name:ident: $hgvs:literal => $printed:literal;)*) => {$(
        #[test]
        fn $name() {
            assert_eq!(parse($hgvs).to_string(), $printed);
        }
    )*};
}

prints_in_three_letter_code! {
    a_one_letter_insertion: "NP_000081.1:p.A356_G357insDPA" => "NP_000081.1:p.Ala356_Gly357insAspProAla";
    a_one_letter_deletion_insertion: "NP_W.1:p.R2delinsKL" => "NP_W.1:p.Arg2delinsLysLeu";
    a_one_letter_insertion_of_x: "NP_W.1:p.A2_L3insX" => "NP_W.1:p.Ala2_Leu3insXaa";
    a_one_letter_repeat_unit: "NP_W.1:p.490PRS[1]" => "NP_W.1:p.490ProArgSer[1]";
    a_one_letter_repeat_unit_over_a_range: "NP_W.1:p.Q18_A20QA[3]" => "NP_W.1:p.Gln18_Ala20GlnAla[3]";
    a_one_letter_insertion_that_spells_a_three_letter_code: "NP_W.1:p.A2_G3insHIS" => "NP_W.1:p.Ala2_Gly3insHisIleSer";
}
