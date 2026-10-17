mod support;

use hgvs_weaver::data::Strand;
use hgvs_weaver::equivalence::{EquivalenceLevel, VariantEquivalence};
use hgvs_weaver::error::HgvsError;
use hgvs_weaver::mapper::VariantMapper;
use hgvs_weaver::parse_hgvs_variant;
use support::{single_exon_transcript, Provider};

/// ACGT repeated as the genome, and a protein NP_0001.1 with Trp2.
fn provider() -> Provider {
    Provider::new()
        .sequence("NC_TEST.1", &"ACGT".repeat(1000))
        .sequence("NP_0001.1", &format!("MW{}", "A".repeat(500)))
}

#[test]
fn test_equivalence_levels() -> Result<(), HgvsError> {
    let hdp = provider();
    let eq_mapper = VariantMapper::new(&hdp);
    let eq = VariantEquivalence::new(&eq_mapper, &hdp);

    // 1. Identity
    let v1 = hgvs_weaver::parse_hgvs_variant("NC_TEST.1:g.1001A>C")?;
    let v2 = hgvs_weaver::parse_hgvs_variant("NC_TEST.1:g.1001A>C")?;
    assert_eq!(eq.equivalent_level(&v1, &v2)?, EquivalenceLevel::Identity);

    // 2. Analogous (Normalization) - ins vs dup. The reference is ACGT repeated,
    // so g.1006 is a C; duplicating it and inserting a C after it are one change.
    let v3 = hgvs_weaver::parse_hgvs_variant("NC_TEST.1:g.1006dupC")?;
    let v4 = hgvs_weaver::parse_hgvs_variant("NC_TEST.1:g.1006_1007insC")?;
    let lvl = eq.equivalent_level(&v3, &v4)?;
    assert!(lvl == EquivalenceLevel::Identity || lvl == EquivalenceLevel::Analogous);

    // 3. Parity (Functional - same protein)
    let p1 = hgvs_weaver::parse_hgvs_variant("NP_0001.1:p.Trp2Ter")?;
    let p2 = hgvs_weaver::parse_hgvs_variant("NP_0001.1:p.Trp2*")?;
    assert_eq!(eq.equivalent_level(&p1, &p2)?, EquivalenceLevel::Identity);

    Ok(())
}

#[test]
fn equivalence_needs_the_sequence() -> Result<(), HgvsError> {
    let hdp = Provider::new().failing_sequences();
    let eq_mapper = VariantMapper::new(&hdp);
    let eq = VariantEquivalence::new(&eq_mapper, &hdp);

    // p.Ala201_Val202insGlyProGlyAla vs p.Gly198_Ala201dup may well be the same
    // change, but only the protein sequence can say: residues 199 and 200 are
    // named by neither. Without a sequence the judgement is an error, not a guess.
    let v1 = hgvs_weaver::parse_hgvs_variant("NP_0001.1:p.Ala201_Val202insGlyProGlyAla")?;
    let v2 = hgvs_weaver::parse_hgvs_variant("NP_0001.1:p.Gly198_Ala201dup")?;
    assert!(eq.equivalent_level(&v1, &v2).is_err());
    Ok(())
}

#[test]
fn a_lone_copy_number_compares_by_projection() -> Result<(), HgvsError> {
    // NM_X.1 is genomic indices 10..50 of NC_X.1, c.1 at index 15, so c.10 is
    // g.25; NR_XR.1 reads the same forty bases on the minus strand, so n.k is
    // g.(51 - k) and g.25_27 is n.24_26.
    let transcript = "GGGGGATGAAACTGGCCTATCGCTAACCGTATAAGTAAGG";
    let minus = "CCTTACTTATACGGTTAGCGATAGGCCAGTTTCATCCCCC";
    let hdp = Provider::new()
        .sequence("NC_X.1", &format!("TTTTTTTTTT{transcript}CCCCCCCCCC"))
        .sequence("NM_X.1", transcript)
        .sequence("NR_XR.1", minus)
        .transcript(support::single_exon_transcript(
            "NM_X.1",
            "NC_X.1",
            10,
            Strand::Plus,
            5,
            25,
            40,
        ))
        .transcript(support::transcript(
            "NR_XR.1",
            "NC_X.1",
            Strand::Minus,
            None,
            vec![support::exon((0, 40), (10, 49), Strand::Minus)],
        ));
    let mapper = VariantMapper::new(&hdp);
    let eq = VariantEquivalence::new(&mapper, &hdp);
    let level = |a: &str, b: &str| -> Result<EquivalenceLevel, HgvsError> {
        eq.equivalent_level(
            &hgvs_weaver::parse_hgvs_variant(a)?,
            &hgvs_weaver::parse_hgvs_variant(b)?,
        )
    };

    assert_eq!(
        level("NM_X.1:c.10_12copy3", "NC_X.1:g.25_27copy3")?,
        EquivalenceLevel::Identity
    );
    // Neither is the other's text; both project to g.25_27copy3.
    assert_eq!(
        level("NM_X.1:c.10_12copy3", "NR_XR.1:n.24_26copy3")?,
        EquivalenceLevel::Analogous
    );
    assert_eq!(
        level("NM_X.1:c.10_12copy3", "NR_XR.1:n.23_25copy3")?,
        EquivalenceLevel::Different
    );
    assert_eq!(
        level("NC_X.1:g.25_27copy3", "NC_X.1:g.25_27copy4")?,
        EquivalenceLevel::Different
    );
    Ok(())
}

/// `NM_X.1`: a five-base 5'UTR, the CDS `ATG AAA CTG GCC TAT CGC TAA` (Met Lys
/// Leu Ala Tyr Arg Ter), a fourteen-base 3'UTR; it encodes `NP_X.1`, `MKLAYR`.
fn coding_provider() -> Provider {
    let transcript = "GGGGGATGAAACTGGCCTATCGCTAACCGTATAAGTAAGG";
    Provider::new()
        .sequence("NM_X.1", transcript)
        .sequence("NP_X.1", "MKLAYR")
        .transcript(single_exon_transcript(
            "NM_X.1",
            "NC_X.1",
            100,
            Strand::Plus,
            5,
            25,
            transcript.len() as i32,
        ))
        .protein_for("NM_X.1", "NP_X.1")
}

/// The level at which `a` and `b` are judged the same on `NM_X.1`/`NP_X.1`.
fn coding_level(a: &str, b: &str) -> EquivalenceLevel {
    let hdp = coding_provider();
    let mapper = VariantMapper::new(&hdp);
    let parse = |s: &str| parse_hgvs_variant(s).unwrap_or_else(|e| panic!("{s}: {e}"));
    VariantEquivalence::new(&mapper, &hdp)
        .equivalent_level(&parse(a), &parse(b))
        .unwrap()
}

#[test]
fn silent_statements_at_different_residues_differ() {
    assert_eq!(
        coding_level("NP_X.1:p.(Leu3=)", "NP_X.1:p.(Ala4=)"),
        EquivalenceLevel::Different
    );
}

#[test]
fn a_silent_residue_is_not_the_whole_protein_unchanged() {
    assert_eq!(
        coding_level("NP_X.1:p.(Leu3=)", "NP_X.1:p.(=)"),
        EquivalenceLevel::Different
    );
}

#[test]
fn a_silent_residue_is_its_statement_in_any_spelling() {
    assert_eq!(
        coding_level("NP_X.1:p.(Leu3=)", "NP_X.1:p.Leu3Leu"),
        EquivalenceLevel::Analogous
    );
}

#[test]
fn a_stop_inserted_before_the_stop_states_its_two_flanks() {
    assert_eq!(
        coding_level("NP_X.1:p.Arg6_Ter7insTer", "NP_X.1:p.Arg6_Ter7="),
        EquivalenceLevel::Analogous
    );
}

#[test]
fn every_residue_through_the_stop_unchanged_is_the_whole_protein_unchanged() {
    assert_eq!(
        coding_level("NP_X.1:p.Met1_Ter7=", "NP_X.1:p.(=)"),
        EquivalenceLevel::Analogous
    );
}

#[test]
fn every_residue_but_the_stop_unchanged_is_not_the_whole_protein_unchanged() {
    assert_eq!(
        coding_level("NP_X.1:p.Met1_Arg6=", "NP_X.1:p.(=)"),
        EquivalenceLevel::Different
    );
}

#[test]
fn an_identity_of_every_codon_is_the_whole_protein_unchanged() {
    assert_eq!(
        coding_level("NM_X.1:c.1_21=", "NP_X.1:p.(=)"),
        EquivalenceLevel::Identity
    );
}

#[test]
fn an_identity_of_every_codon_from_the_5_utr_is_every_residue_unchanged() {
    assert_eq!(
        coding_level("NM_X.1:c.-2_21=", "NP_X.1:p.Met1_Ter7="),
        EquivalenceLevel::Analogous
    );
}

#[test]
fn a_stop_inserted_before_the_stop_is_not_the_stop_unchanged() {
    assert_eq!(
        coding_level("NP_X.1:p.Arg6_Ter7insTer", "NP_X.1:p.(Ter7=)"),
        EquivalenceLevel::Different
    );
}

#[test]
fn a_silent_change_in_the_stop_codon_is_not_a_stop_inserted_before_it() {
    assert_eq!(
        coding_level("NM_X.1:c.21A>G", "NP_X.1:p.Arg6_Ter7insTer"),
        EquivalenceLevel::Different
    );
}

#[test]
fn one_base_unchanged_is_its_residue_unchanged_read_outward() {
    assert_eq!(
        coding_level("NM_X.1:c.9=", "NP_X.1:p.Leu3="),
        EquivalenceLevel::Analogous
    );
}

#[test]
fn one_base_unchanged_is_not_another_residue_unchanged() {
    assert_eq!(
        coding_level("NM_X.1:c.9=", "NP_X.1:p.Ala4="),
        EquivalenceLevel::Different
    );
}

#[test]
fn an_identity_across_two_codons_is_both_residues_read_outward() {
    assert_eq!(
        coding_level("NM_X.1:c.8_10=", "NP_X.1:p.Leu3_Ala4="),
        EquivalenceLevel::Analogous
    );
}

#[test]
fn an_identity_across_two_codons_is_not_the_first_residue_alone() {
    assert_eq!(
        coding_level("NM_X.1:c.8_10=", "NP_X.1:p.(Leu3=)"),
        EquivalenceLevel::Different
    );
}

#[test]
fn an_identity_is_the_whole_codons_it_covers_as_predicted() {
    assert_eq!(
        coding_level("NM_X.1:c.7_10=", "NP_X.1:p.(Leu3=)"),
        EquivalenceLevel::Identity
    );
}

#[test]
fn an_identity_is_also_every_residue_it_touches() {
    assert_eq!(
        coding_level("NM_X.1:c.7_10=", "NP_X.1:p.Leu3_Ala4="),
        EquivalenceLevel::Analogous
    );
}

#[test]
fn an_identity_reaching_into_the_start_codon_is_it_read_outward() {
    assert_eq!(
        coding_level("NM_X.1:c.-2_2=", "NP_X.1:p.Met1="),
        EquivalenceLevel::Analogous
    );
}

#[test]
fn an_identity_in_the_3_utr_is_no_residue_unchanged() {
    assert_eq!(
        coding_level("NM_X.1:c.*3=", "NP_X.1:p.Ter7="),
        EquivalenceLevel::Different
    );
}

#[test]
fn an_identity_in_the_3_utr_is_not_the_whole_protein_unchanged() {
    assert_eq!(
        coding_level("NM_X.1:c.*3=", "NP_X.1:p.(=)"),
        EquivalenceLevel::Different
    );
}

#[test]
fn a_coding_change_predicted_silent_is_its_statement() {
    assert_eq!(
        coding_level("NM_X.1:c.9G>A", "NP_X.1:p.(Leu3=)"),
        EquivalenceLevel::Identity
    );
}

#[test]
fn a_coding_change_predicted_silent_says_nothing_of_another_residue() {
    assert_eq!(
        coding_level("NM_X.1:c.9G>A", "NP_X.1:p.(Ala4=)"),
        EquivalenceLevel::Different
    );
}

#[test]
fn a_three_prime_utr_change_is_the_whole_protein_unchanged() {
    assert_eq!(
        coding_level("NM_X.1:c.*3G>A", "NP_X.1:p.(=)"),
        EquivalenceLevel::Identity
    );
}

#[test]
fn a_three_prime_utr_change_is_not_the_stop_unchanged() {
    assert_eq!(
        coding_level("NM_X.1:c.*3G>A", "NP_X.1:p.Ter7="),
        EquivalenceLevel::Different
    );
}

#[test]
fn a_silent_change_in_the_stop_codon_is_the_stop_unchanged() {
    assert_eq!(
        coding_level("NM_X.1:c.21A>G", "NP_X.1:p.Ter7="),
        EquivalenceLevel::Analogous
    );
}

#[test]
fn a_silent_change_in_the_stop_codon_is_the_stop_written_as_itself() {
    assert_eq!(
        coding_level("NM_X.1:c.21A>G", "NP_X.1:p.Ter7Ter"),
        EquivalenceLevel::Analogous
    );
}

#[test]
fn a_silent_change_in_the_stop_codon_is_not_a_stop_loss() {
    assert_eq!(
        coding_level("NM_X.1:c.21A>G", "NP_X.1:p.Ter7Trp"),
        EquivalenceLevel::Different
    );
}

#[test]
fn the_stop_unchanged_is_the_stop_written_as_itself() {
    assert_eq!(
        coding_level("NP_X.1:p.Ter7=", "NP_X.1:p.Ter7Ter"),
        EquivalenceLevel::Analogous
    );
}

#[test]
fn replacing_the_stop_over_a_range_is_not_the_protein_unchanged() {
    assert_eq!(
        coding_level("NP_X.1:p.Arg6_Ter7delinsArg", "NP_X.1:p.Arg6_Ter7="),
        EquivalenceLevel::Different
    );
}

#[test]
fn an_insertion_before_the_stop_is_the_residue_added_before_it() {
    assert_eq!(
        coding_level(
            "NP_X.1:p.Arg6_Ter7insGly",
            "NP_X.1:p.Arg6_Ter7delinsArgGlyTer"
        ),
        EquivalenceLevel::Analogous
    );
}

#[test]
fn a_silent_prediction_is_not_read_again_on_the_protein() {
    // NP_X.1 served with Ile3 where the transcript encodes Leu3: the silent
    // prediction, like a missense one, is the transcript's.
    let transcript = "GGGGGATGAAACTGGCCTATCGCTAACCGTATAAGTAAGG";
    let hdp = Provider::new()
        .sequence("NM_X.1", transcript)
        .sequence("NP_X.1", "MKIAYR")
        .transcript(single_exon_transcript(
            "NM_X.1",
            "NC_X.1",
            100,
            Strand::Plus,
            5,
            25,
            transcript.len() as i32,
        ))
        .protein_for("NM_X.1", "NP_X.1");
    let mapper = VariantMapper::new(&hdp);
    let level = VariantEquivalence::new(&mapper, &hdp)
        .equivalent_level(
            &parse_hgvs_variant("NM_X.1:c.9G>A").unwrap(),
            &parse_hgvs_variant("NP_X.1:p.Ala4Gly").unwrap(),
        )
        .unwrap();
    assert_eq!(level, EquivalenceLevel::Different);
}

#[test]
fn an_identity_from_the_5_utr_into_the_cds_is_the_start_codon_unchanged() {
    assert_eq!(
        coding_level("NM_X.1:c.-2_3=", "NP_X.1:p.Met1="),
        EquivalenceLevel::Analogous
    );
}

#[test]
fn an_identity_from_the_5_utr_into_the_cds_is_not_the_whole_protein_unchanged() {
    assert_eq!(
        coding_level("NM_X.1:c.-2_3=", "NP_X.1:p.(=)"),
        EquivalenceLevel::Different
    );
}

#[test]
fn an_unknown_change_from_the_stop_is_not_the_protein_unchanged() {
    assert_eq!(
        coding_level("NP_X.1:p.Ter7?", "NP_X.1:p.(=)"),
        EquivalenceLevel::Different
    );
}
