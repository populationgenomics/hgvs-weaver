//! `normalize` writes one change one way: a delins is reduced to the edit it
//! is before it is shifted, and the bases an edit can state are checked and
//! then not stated.

mod support;

use hgvs_weaver::error::HgvsError;
use hgvs_weaver::mapper::VariantMapper;
use hgvs_weaver::parse_hgvs_variant;
use support::Provider;

/// `NC_000099.1` is the ten bases `TTCAGCAGTT`: g.3 to g.5 and g.6 to g.8 are
/// both `CAG`.
fn provider() -> Provider {
    Provider::new().sequence("NC_000099.1", "TTCAGCAGTT")
}

fn normalized(s: &str) -> String {
    let hdp = provider();
    let v = parse_hgvs_variant(s).unwrap_or_else(|e| panic!("{s}: {e}"));
    VariantMapper::new(&hdp)
        .normalize_variant(v)
        .unwrap_or_else(|e| panic!("normalising {s}: {e}"))
        .to_string()
}

#[test]
fn a_delins_that_keeps_its_first_base_is_an_insertion() {
    assert_eq!(
        normalized("NC_000099.1:g.3delinsCGG"),
        "NC_000099.1:g.3_4insGG"
    );
}

#[test]
fn a_delins_that_is_an_insertion_is_shifted_3_prime() {
    // C becomes CAA: an insertion of AA after g.3, which slides past the A at g.4.
    assert_eq!(
        normalized("NC_000099.1:g.3delinsCAA"),
        "NC_000099.1:g.4_5insAA"
    );
}

#[test]
fn a_delins_that_is_an_insertion_repeating_the_bases_before_it_is_a_duplication() {
    // C becomes CAG: AG inserted after g.6 slides to follow g.8, and repeats g.7_8.
    assert_eq!(
        normalized("NC_000099.1:g.6delinsCAG"),
        "NC_000099.1:g.7_8dup"
    );
}

#[test]
fn a_delins_that_keeps_its_last_base_is_an_insertion() {
    assert_eq!(
        normalized("NC_000099.1:g.5delinsTG"),
        "NC_000099.1:g.4_5insT"
    );
}

#[test]
fn a_delins_that_is_a_deletion_is_shifted_3_prime() {
    // TCAG becomes T: the CAG at g.3_5 deleted, the 3'-most copy g.6_8.
    assert_eq!(
        normalized("NC_000099.1:g.2_5delinsT"),
        "NC_000099.1:g.6_8del"
    );
}

#[test]
fn a_delins_that_changes_one_base_is_a_substitution() {
    assert_eq!(
        normalized("NC_000099.1:g.3_5delinsCTG"),
        "NC_000099.1:g.4A>T"
    );
}

#[test]
fn a_delins_that_is_the_reverse_complement_is_an_inversion() {
    // CA becomes TG.
    assert_eq!(
        normalized("NC_000099.1:g.3_4delinsTG"),
        "NC_000099.1:g.3_4inv"
    );
}

#[test]
fn a_delins_that_changes_a_base_at_each_end_stays_where_it_is() {
    assert_eq!(
        normalized("NC_000099.1:g.3_5delinsTT"),
        "NC_000099.1:g.3_5delinsTT"
    );
}

#[test]
fn a_delins_is_written_without_its_deleted_bases() {
    assert_eq!(
        normalized("NC_000099.1:g.3_5delCAGinsTT"),
        "NC_000099.1:g.3_5delinsTT"
    );
}

#[test]
fn a_delins_that_changes_nothing_is_an_identity() {
    assert_eq!(
        normalized("NC_000099.1:g.3_5delinsCAG"),
        "NC_000099.1:g.3_5="
    );
}

#[test]
fn a_deletion_that_does_not_move_is_written_without_its_bases() {
    assert_eq!(
        normalized("NC_000099.1:g.6_8delCAG"),
        "NC_000099.1:g.6_8del"
    );
}

#[test]
fn a_deletion_written_with_its_length_is_written_without_it() {
    assert_eq!(normalized("NC_000099.1:g.6_8del3"), "NC_000099.1:g.6_8del");
}

#[test]
fn a_duplication_that_does_not_move_is_written_without_its_bases() {
    assert_eq!(
        normalized("NC_000099.1:g.6_8dupCAG"),
        "NC_000099.1:g.6_8dup"
    );
}

#[test]
fn an_inversion_is_written_without_its_bases() {
    assert_eq!(normalized("NC_000099.1:g.3_4invCA"), "NC_000099.1:g.3_4inv");
}

#[test]
fn an_identity_of_one_base_is_written_without_it() {
    assert_eq!(normalized("NC_000099.1:g.3C="), "NC_000099.1:g.3=");
}

#[test]
fn an_identity_of_a_range_is_written_without_its_bases() {
    assert_eq!(normalized("NC_000099.1:g.3_5CAG="), "NC_000099.1:g.3_5=");
}

#[test]
fn a_delins_of_a_number_of_bases_is_written_without_its_deleted_bases() {
    assert_eq!(
        normalized("NC_000099.1:g.3_5delCAGinsN[2]"),
        "NC_000099.1:g.3_5delinsN[2]"
    );
}

/// The bases normalising `s` refused, as `(expected, found)`.
fn refused(s: &str) -> (String, String) {
    let hdp = provider();
    let v = parse_hgvs_variant(s).unwrap_or_else(|e| panic!("{s}: {e}"));
    match VariantMapper::new(&hdp).normalize_variant(v) {
        Err(HgvsError::ReferenceMismatch {
            expected, found, ..
        }) => (expected, found),
        other => panic!("normalising {s} gave {other:?}, not a mismatch"),
    }
}

#[test]
fn a_deletion_stating_bases_the_reference_lacks_is_refused() {
    assert_eq!(
        refused("NC_000099.1:g.3_5delAAA"),
        ("AAA".to_string(), "CAG".to_string())
    );
}

#[test]
fn a_delins_whose_deleted_bases_the_reference_lacks_is_refused() {
    assert_eq!(
        refused("NC_000099.1:g.3_5delAAAinsTT"),
        ("AAA".to_string(), "CAG".to_string())
    );
}

#[test]
fn a_delins_of_a_number_of_bases_whose_deleted_bases_the_reference_lacks_is_refused() {
    assert_eq!(
        refused("NC_000099.1:g.3_5delAAAinsN[2]"),
        ("AAA".to_string(), "CAG".to_string())
    );
}

#[test]
fn an_inversion_stating_bases_the_reference_lacks_is_refused() {
    assert_eq!(
        refused("NC_000099.1:g.3_4invGG"),
        ("GG".to_string(), "CA".to_string())
    );
}

#[test]
fn an_identity_stating_a_base_the_reference_lacks_is_refused() {
    assert_eq!(
        refused("NC_000099.1:g.3G="),
        ("G".to_string(), "C".to_string())
    );
}

#[test]
fn a_substitution_keeps_its_stated_base() {
    assert_eq!(normalized("NC_000099.1:g.3C>T"), "NC_000099.1:g.3C>T");
}

#[test]
fn one_base_stated_over_two_is_refused() {
    assert_eq!(
        refused("NC_000099.1:g.3_4delCinsT"),
        ("C".to_string(), "CA".to_string())
    );
}

#[test]
fn a_range_of_one_base_is_written_as_the_base() {
    assert_eq!(normalized("NC_000099.1:g.6_6dup"), "NC_000099.1:g.6dup");
}

#[test]
fn an_insertion_does_not_slide_past_the_last_base_of_the_sequence() {
    // TA before g.10, a T, could slide to follow it as AT, after a base that
    // does not exist.
    assert_eq!(
        normalized("NC_000099.1:g.9_10insTA"),
        "NC_000099.1:g.9_10insTA"
    );
}

#[test]
fn a_delins_reduced_to_an_insertion_after_the_last_base_keeps_the_last_base() {
    assert_eq!(
        normalized("NC_000099.1:g.10delinsTA"),
        "NC_000099.1:g.10delinsTA"
    );
}

#[test]
fn a_change_with_an_uncertain_position_is_left_as_written() {
    assert_eq!(
        normalized("NC_000099.1:g.(3_4)delinsCGG"),
        "NC_000099.1:g.(3_4)delinsCGG"
    );
}

#[test]
fn a_change_with_uncertain_ends_is_left_as_written() {
    assert_eq!(
        normalized("NC_000099.1:g.(1_3)_(5_6)delinsT"),
        "NC_000099.1:g.(1_3)_(5_6)delinsT"
    );
}
