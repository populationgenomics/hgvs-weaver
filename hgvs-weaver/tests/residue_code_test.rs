//! A protein description writes every residue in one amino acid code, the
//! one-letter or the three-letter; `*`, `Xaa` and `X` belong to either, and
//! an `X` is read where only one meaning fits it.

use hgvs_weaver::error::{HgvsError, ResidueToken};
use hgvs_weaver::{
    parse_hgvs_variant, parse_hgvs_variant_with, ParseOptions, SequenceVariant, XTreatment,
};

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
    a_one_letter_insertion_between_three_letter_residues: "NP_000081.1:p.Ala356_Gly357insDPA";
    a_one_letter_deletion_insertion_at_a_three_letter_residue: "NP_W.1:p.Arg2delinsKL";
    a_one_letter_repeat_unit_over_three_letter_residues: "NP_W.1:p.Pro490_Ser492PRS[1]";
    a_three_letter_repeat_unit_over_one_letter_residues: "NP_W.1:p.Q18_A20GlnAla[3]";
    a_three_letter_insertion_between_one_letter_residues: "NP_W.1:p.A2_G3insTerGlu";
    a_three_letter_substitution_at_a_one_letter_residue: "NP_W.1:p.R2Lys";
    a_range_whose_ends_use_different_codes: "NP_W.1:p.Lys2_L4del";
    a_one_letter_frameshift_to_a_stop_written_ter: "NP_005206.2:p.V191GfsTer35";
    an_n_terminal_extension_from_a_three_letter_residue: "NP_W.1:p.M1extLeu-5";
    a_predicted_change_mixing_the_codes: "NP_W.1:p.(Arg2K)";
    an_allele_whose_members_use_different_codes: "NP_W.1:p.[Arg2Lys;L3P]";
    a_predicted_allele_whose_members_use_different_codes: "NP_W.1:p.[(Arg2Lys;L3P)]";
}

/// A residue token as `(text, start, end)`.
type Token = (String, usize, usize);

/// The one-letter and the three-letter tokens the refusal of `hgvs` carries.
fn tokens(hgvs: &str) -> (Vec<Token>, Vec<Token>) {
    let as_tuples = |ts: Vec<ResidueToken>| -> Vec<Token> {
        ts.into_iter().map(|t| (t.text, t.start, t.end)).collect()
    };
    match parse_hgvs_variant(hgvs) {
        Err(HgvsError::MixedResidueCodes {
            one_letter,
            three_letter,
            ..
        }) => (as_tuples(one_letter), as_tuples(three_letter)),
        other => panic!("{hgvs}: {other:?}"),
    }
}

fn token(text: &str, start: usize, end: usize) -> Token {
    (text.to_string(), start, end)
}

#[test]
fn the_refusal_carries_a_lone_residue_and_the_three_letter_residue_where_written() {
    assert_eq!(
        tokens("NP_W.1:p.Arg111K"),
        (vec![token("K", 15, 16)], vec![token("Arg", 9, 12)])
    );
}

#[test]
fn the_refusal_leaves_a_lone_x_in_neither_code() {
    assert_eq!(
        tokens("NP_W.1:p.[Arg2X;L3P]"),
        (
            vec![token("L", 16, 17), token("P", 18, 19)],
            vec![token("Arg", 10, 13)]
        )
    );
}

#[test]
fn the_refusal_carries_a_one_letter_run_as_one_token() {
    assert_eq!(
        tokens("NP_000081.1:p.Ala356_Gly357insDPA"),
        (
            vec![token("DPA", 30, 33)],
            vec![token("Ala", 14, 17), token("Gly", 21, 24)]
        )
    );
}

#[test]
fn the_refusal_carries_an_upper_case_run_as_one_one_letter_token() {
    assert_eq!(
        tokens("NP_W.1:p.Ala2_Gly3insHIS"),
        (
            vec![token("HIS", 21, 24)],
            vec![token("Ala", 9, 12), token("Gly", 14, 17)]
        )
    );
}

#[test]
fn the_refusal_of_an_allele_carries_the_tokens_of_every_member() {
    assert_eq!(
        tokens("NP_W.1:p.[Arg2Lys;L3P]"),
        (
            vec![token("L", 18, 19), token("P", 20, 21)],
            vec![token("Arg", 10, 13), token("Lys", 14, 17)]
        )
    );
}

#[test]
fn a_one_letter_run_ending_in_a_stop_written_star_is_one_one_letter_token() {
    assert_eq!(
        tokens("NP_W.1:p.Arg2_Lys3delinsL*"),
        (
            vec![token("L*", 24, 26)],
            vec![token("Arg", 9, 12), token("Lys", 14, 17)]
        )
    );
}

#[test]
fn the_refusal_of_a_predicted_allele_carries_offsets_in_the_whole_description() {
    assert_eq!(
        tokens("NP_W.1:p.[(Arg2Lys;L3P)]"),
        (
            vec![token("L", 19, 20), token("P", 21, 22)],
            vec![token("Arg", 11, 14), token("Lys", 15, 18)]
        )
    );
}

#[test]
fn the_refusal_counts_offsets_from_the_start_of_the_input_whitespace_included() {
    let err = parse_hgvs_variant("  NP_W.1:p.Arg2K").unwrap_err();
    let HgvsError::MixedResidueCodes {
        description,
        one_letter,
        ..
    } = err
    else {
        panic!("{err:?}");
    };
    assert_eq!(description, "  NP_W.1:p.Arg2K");
    assert_eq!(
        one_letter,
        vec![ResidueToken {
            text: "K".to_string(),
            start: 15,
            end: 16
        }]
    );
}

#[test]
fn the_message_names_a_run_holding_x_before_an_earlier_one_letter_residue() {
    let message = parse_hgvs_variant("NP_W.1:p.Lys2_L4delinsKX")
        .unwrap_err()
        .to_string();
    assert!(
        message.contains("and KX in the one-letter code"),
        "{message}"
    );
}

#[test]
fn the_refusal_of_a_run_holding_x_says_how_a_stop_and_an_unknown_residue_are_written() {
    let message = parse_hgvs_variant("NP_W.1:p.Arg2_Lys3delinsKX")
        .unwrap_err()
        .to_string();
    assert!(message.contains("a stop is Ter or *"), "{message}");
}

#[test]
fn the_refusal_names_the_description() {
    let err = parse_hgvs_variant("NP_W.1:p.Arg2K").unwrap_err();
    assert!(
        matches!(&err, HgvsError::MixedResidueCodes { description, .. }
            if description == "NP_W.1:p.Arg2K"),
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

/// One test per case: `$hgvs` is the description `$same`, `Xaa` belonging
/// to either code.
macro_rules! xaa_is_either_code {
    ($($name:ident: $hgvs:literal => $same:literal;)*) => {$(
        #[test]
        fn $name() {
            assert_eq!(parse($hgvs), parse($same));
        }
    )*};
}

xaa_is_either_code! {
    a_substitution_to_xaa_at_a_one_letter_residue: "NP_W.1:p.R50Xaa" => "NP_W.1:p.Arg50Xaa";
    a_one_letter_allele_with_a_member_to_xaa: "NP_W.1:p.[R2Xaa;L3P]" => "NP_W.1:p.[Arg2Xaa;Leu3Pro]";
}

#[test]
fn a_mixed_allele_does_not_count_xaa_in_either_code() {
    assert_eq!(
        tokens("NP_W.1:p.[Arg2Xaa;L3P]"),
        (
            vec![token("L", 18, 19), token("P", 20, 21)],
            vec![token("Arg", 10, 13)]
        )
    );
}

/// The `X`s the ambiguity refusal of `hgvs` carries, as `(text, start, end)`.
fn ambiguous(hgvs: &str) -> Vec<Token> {
    match parse_hgvs_variant(hgvs) {
        Err(HgvsError::AmbiguousResidue { ambiguous, .. }) => ambiguous
            .into_iter()
            .map(|t| (t.text, t.start, t.end))
            .collect(),
        other => panic!("{hgvs}: {other:?}"),
    }
}

/// One test per case: `$hgvs`, a one-letter description, is refused for the
/// `X`s at `$spans`, each of which could be a stop or any residue.
macro_rules! an_x_is_ambiguous {
    ($($name:ident: $hgvs:literal => [$(($start:literal, $end:literal)),*];)*) => {$(
        #[test]
        fn $name() {
            assert_eq!(ambiguous($hgvs), vec![$(token("X", $start, $end)),*]);
        }
    )*};
}

an_x_is_ambiguous! {
    a_one_letter_substitution_to_x: "NP_000026.2:p.R60X" => [(17, 18)];
    a_one_letter_insertion_of_x: "NP_W.1:p.A2_G3insX" => [(17, 18)];
    a_one_letter_deletion_insertion_ending_in_x: "NP_W.1:p.R2delinsKX" => [(18, 19)];
    every_x_of_a_one_letter_allele: "NP_W.1:p.[R2X;A4_G5insKXX]" => [(12, 13), (23, 24), (24, 25)];
    an_x_in_a_predicted_one_letter_allele: "NP_W.1:p.[(R2X;L3P)]" => [(13, 14)];
    a_substitution_to_x_at_a_residue_written_xaa: "NP_W.1:p.Xaa5X" => [(13, 14)];
    a_three_letter_nonsense_change_written_x_as_clinvar_has_it: "NP_000021.1:p.Arg111X" => [(20, 21)];
    a_three_letter_substitution_to_x: "NP_060360.3:p.Leu2781X" => [(21, 22)];
    a_predicted_three_letter_substitution_to_x: "NP_W.1:p.(Arg2X)" => [(14, 15)];
    an_x_in_a_three_letter_insertion: "NP_W.1:p.Ala2_Gly3insLeuX" => [(24, 25)];
    every_x_of_a_three_letter_allele: "NP_W.1:p.[Arg2X;Ala4_Gly5insLysXX]" => [(14, 15), (31, 32), (32, 33)];
}

#[test]
fn the_ambiguity_message_offers_a_stop_and_any_residue() {
    let message = parse_hgvs_variant("NP_000026.2:p.R60X")
        .unwrap_err()
        .to_string();
    assert!(
        message.contains("write * for a stop and Xaa for any residue"),
        "{message}"
    );
}

/// One test per case: `$hgvs` reads its one-letter `X` the one way it can,
/// as the description `$read` does.
macro_rules! an_x_with_one_reading {
    ($($name:ident: $hgvs:literal => $read:literal;)*) => {$(
        #[test]
        fn $name() {
            assert_eq!(parse($hgvs), parse($read));
        }
    )*};
}

an_x_with_one_reading! {
    a_frameshift_to_a_stop_written_x: "NP_W.1:p.R210fsX5" => "NP_W.1:p.Arg210fsTer5";
    an_extension_to_a_stop_written_x: "NP_W.1:p.*30QextX5" => "NP_W.1:p.Ter30GlnextTer5";
    a_stop_replaced_by_x_which_cannot_be_a_stop: "NP_W.1:p.*30X" => "NP_W.1:p.Ter30Xaa";
    a_repeat_unit_holding_x_checked_against_the_sequence: "NP_W.1:p.Q18_A20QX[3]" => "NP_W.1:p.Gln18_Ala20GlnXaa[3]";
    a_frameshift_from_x_whose_first_residue_is_no_stop: "NP_W.1:p.R2Xfs*5" => "NP_W.1:p.Arg2XaafsTer5";
    an_n_terminal_extension_from_x_which_starts_no_protein_at_a_stop: "NP_W.1:p.M1extX-5" => "NP_W.1:p.Met1extXaa-5";
    a_reference_residue_written_x_checked_against_the_sequence: "NP_W.1:p.X5_L7del" => "NP_W.1:p.Xaa5_Leu7del";
    a_three_letter_frameshift_written_fsx_as_clinvar_has_it: "NP_004619.3:p.Ile223MetfsX45" => "NP_004619.3:p.Ile223MetfsTer45";
    a_three_letter_frameshift_to_x_with_no_new_residue: "NP_060360.3:p.Leu2781fsX5" => "NP_060360.3:p.Leu2781fsTer5";
    a_three_letter_extension_to_a_stop_written_x: "NP_W.1:p.*30GlnextX5" => "NP_W.1:p.Ter30GlnextTer5";
    a_three_letter_stop_replaced_by_x: "NP_W.1:p.Ter30X" => "NP_W.1:p.Ter30Xaa";
    a_three_letter_frameshift_from_x: "NP_W.1:p.Arg2Xfs*5" => "NP_W.1:p.Arg2XaafsTer5";
    a_three_letter_range_from_a_residue_written_x: "NP_W.1:p.X5_Leu7del" => "NP_W.1:p.Xaa5_Leu7del";
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
    an_insertion_of_xaa_between_one_letter_residues: "NP_W.1:p.A2_L3insXaa" => "NP_W.1:p.Ala2_Leu3insXaa";
    a_one_letter_run_holding_xaa: "NP_W.1:p.A2_G3insKXaaXaa" => "NP_W.1:p.Ala2_Gly3insLysXaaXaa";
    a_one_letter_run_holding_xaa_and_a_stop: "NP_W.1:p.R2delinsXaaL*" => "NP_W.1:p.Arg2delinsXaaLeuTer";
    a_one_letter_repeat_unit: "NP_W.1:p.490PRS[1]" => "NP_W.1:p.490ProArgSer[1]";
    a_one_letter_repeat_unit_over_a_range: "NP_W.1:p.Q18_A20QA[3]" => "NP_W.1:p.Gln18_Ala20GlnAla[3]";
    a_one_letter_insertion_that_spells_a_three_letter_code: "NP_W.1:p.A2_G3insHIS" => "NP_W.1:p.Ala2_Gly3insHisIleSer";
}

fn parse_treating_x(s: &str, x_treatment: XTreatment) -> Result<SequenceVariant, HgvsError> {
    parse_hgvs_variant_with(s, ParseOptions::default().with_x_treatment(x_treatment))
}

/// One test per case: `$hgvs`, read with `$treatment`, is the description
/// `$read`.
macro_rules! an_x_read_by_its_treatment {
    ($($name:ident: $hgvs:literal, $treatment:ident => $read:literal;)*) => {$(
        #[test]
        fn $name() {
            let read = parse_treating_x($hgvs, XTreatment::$treatment)
                .unwrap_or_else(|e| panic!("{}: {e}", $hgvs));
            assert_eq!(read, parse($read));
        }
    )*};
}

an_x_read_by_its_treatment! {
    a_one_letter_substitution_to_x_read_as_a_stop: "NP_000026.2:p.R60X", Stop => "NP_000026.2:p.Arg60Ter";
    a_one_letter_substitution_to_x_read_as_any_residue: "NP_000026.2:p.R60X", Unknown => "NP_000026.2:p.Arg60Xaa";
    a_three_letter_substitution_to_x_read_as_a_stop: "NP_000021.1:p.Arg111X", Stop => "NP_000021.1:p.Arg111Ter";
    a_three_letter_substitution_to_x_read_as_any_residue: "NP_000021.1:p.Arg111X", Unknown => "NP_000021.1:p.Arg111Xaa";
    a_one_letter_run_of_x_read_as_stops: "NP_W.1:p.A2_G3insKXX", Stop => "NP_W.1:p.Ala2_Gly3insLysTerTer";
    a_one_letter_run_of_x_read_as_any_residues: "NP_W.1:p.A2_G3insKXX", Unknown => "NP_W.1:p.Ala2_Gly3insLysXaaXaa";
    a_three_letter_run_ending_in_x_read_as_any_residue: "NP_W.1:p.Ala2_Gly3insLeuX", Unknown => "NP_W.1:p.Ala2_Gly3insLeuXaa";
    every_x_of_an_allele_read_as_a_stop: "NP_W.1:p.[R2X;A4_G5insKXX]", Stop => "NP_W.1:p.[Arg2Ter;Ala4_Gly5insLysTerTer]";
    an_x_after_leading_whitespace_read_as_a_stop: "  NP_000026.2:p.R60X", Stop => "NP_000026.2:p.Arg60Ter";
    an_x_after_fs_stays_a_stop_read_as_any_residue: "NP_W.1:p.R210fsX5", Unknown => "NP_W.1:p.Arg210fsTer5";
    a_stop_replaced_by_x_stays_any_residue_read_as_a_stop: "NP_W.1:p.*30X", Stop => "NP_W.1:p.Ter30Xaa";
    an_x_at_a_position_stays_any_residue_read_as_a_stop: "NP_W.1:p.X5_L7del", Stop => "NP_W.1:p.Xaa5_Leu7del";
}

#[test]
fn the_default_treatment_refuses_an_ambiguous_x() {
    let result = parse_treating_x("NP_000026.2:p.R60X", XTreatment::default());
    assert!(
        matches!(result, Err(HgvsError::AmbiguousResidue { .. })),
        "{result:?}"
    );
}

#[test]
fn a_treatment_of_x_leaves_a_mix_of_codes_refused() {
    let result = parse_treating_x("NP_000081.1:p.Ala356_Gly357insDPA", XTreatment::Stop);
    assert!(
        matches!(result, Err(HgvsError::MixedResidueCodes { .. })),
        "{result:?}"
    );
}
