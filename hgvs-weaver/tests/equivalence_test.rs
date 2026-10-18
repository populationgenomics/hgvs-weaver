mod support;

use hgvs_weaver::data::{Strand, TranscriptData};
use hgvs_weaver::equivalence::{EquivalenceLevel, VariantEquivalence};
use hgvs_weaver::error::HgvsError;
use hgvs_weaver::mapper::VariantMapper;
use hgvs_weaver::{parse_hgvs_variant, SequenceVariant};
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
fn one_base_unchanged_is_not_its_residue_unchanged() {
    assert_eq!(
        coding_level("NM_X.1:c.9=", "NP_X.1:p.Leu3="),
        EquivalenceLevel::Different
    );
}

#[test]
fn one_base_unchanged_is_not_the_whole_protein_unchanged() {
    assert_eq!(
        coding_level("NM_X.1:c.9=", "NP_X.1:p.(=)"),
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
fn an_identity_is_its_prediction_in_another_spelling() {
    assert_eq!(
        coding_level("NM_X.1:c.7_10=", "NP_X.1:p.L3="),
        EquivalenceLevel::Identity
    );
}

#[test]
fn an_identity_is_not_every_residue_it_touches() {
    assert_eq!(
        coding_level("NM_X.1:c.7_10=", "NP_X.1:p.Leu3_Ala4="),
        EquivalenceLevel::Different
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
        EquivalenceLevel::Identity
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
        EquivalenceLevel::Identity
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

// --- Parentheses ---
//
// Only a prediction's parentheses, around a whole p. description, are
// spelling; those of an uncertain position say the change is not known
// exactly, so the description names no one change.

/// The level of `a` against `b` on `NC_TEST.1`, `ACGT` repeated.
fn genomic_level(a: &str, b: &str) -> EquivalenceLevel {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let parse = |s: &str| parse_hgvs_variant(s).unwrap_or_else(|e| panic!("{s}: {e}"));
    VariantEquivalence::new(&mapper, &hdp)
        .equivalent_level(&parse(a), &parse(b))
        .unwrap()
}

/// [`coding_provider`] with `NM_X.1` at genomic index 100 of `NC_X.1`, so
/// c.1 is g.106.
fn placed_provider() -> Provider {
    let transcript = "GGGGGATGAAACTGGCCTATCGCTAACCGTATAAGTAAGG";
    coding_provider().sequence(
        "NC_X.1",
        &format!("{}{transcript}{}", "T".repeat(100), "T".repeat(20)),
    )
}

/// The level of `a` against `b` over [`placed_provider`].
fn placed_level(a: &str, b: &str) -> EquivalenceLevel {
    let hdp = placed_provider();
    let mapper = VariantMapper::new(&hdp);
    let parse = |s: &str| parse_hgvs_variant(s).unwrap_or_else(|e| panic!("{s}: {e}"));
    VariantEquivalence::new(&mapper, &hdp)
        .equivalent_level(&parse(a), &parse(b))
        .unwrap()
}

#[test]
fn one_base_deleted_somewhere_in_a_range_is_not_the_range_deleted() {
    assert_eq!(
        genomic_level("NC_TEST.1:g.(1003_1004)del", "NC_TEST.1:g.1003_1004del"),
        EquivalenceLevel::Different
    );
}

#[test]
fn one_base_deleted_somewhere_in_a_range_is_not_its_first_base_deleted() {
    assert_eq!(
        genomic_level("NC_TEST.1:g.(1003_1004)del", "NC_TEST.1:g.1003del"),
        EquivalenceLevel::Different
    );
}

#[test]
fn one_base_deleted_somewhere_in_a_range_is_itself() {
    assert_eq!(
        genomic_level("NC_TEST.1:g.(1003_1004)del", "NC_TEST.1:g.(1003_1004)del"),
        EquivalenceLevel::Identity
    );
}

#[test]
fn a_deletion_between_uncertain_breakpoints_is_not_one_from_their_first_bases() {
    assert_eq!(
        genomic_level(
            "NC_TEST.1:g.(1001_1003)_(1010_1012)del",
            "NC_TEST.1:g.1001_1010del"
        ),
        EquivalenceLevel::Different
    );
}

#[test]
fn a_copy_number_somewhere_in_a_range_is_not_one_of_the_range() {
    assert_eq!(
        genomic_level("NC_TEST.1:g.(1001_1010)copy3", "NC_TEST.1:g.1001_1010copy3"),
        EquivalenceLevel::Different
    );
}

#[test]
fn a_coding_deletion_over_an_uncertain_range_is_not_the_exact_one() {
    assert_eq!(
        placed_level("NM_X.1:c.(4_6)del", "NM_X.1:c.4_6del"),
        EquivalenceLevel::Different
    );
}

#[test]
fn a_coding_deletion_over_an_uncertain_range_is_not_the_exact_genomic_one() {
    // c.4_6 is g.109_111.
    assert_eq!(
        placed_level("NM_X.1:c.(4_6)del", "NC_X.1:g.109_111del"),
        EquivalenceLevel::Different
    );
}

#[test]
fn one_base_deleted_somewhere_in_a_range_has_no_canonical_allele() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let var = parse_hgvs_variant("NC_TEST.1:g.(1003_1004)del").unwrap();
    let result = mapper.canonical_allele(&var);
    assert!(
        matches!(result, Err(HgvsError::UnsupportedOperation(_))),
        "{result:?}"
    );
}

#[test]
fn a_coding_deletion_over_an_uncertain_range_has_no_genomic_form() {
    let hdp = placed_provider();
    let mapper = VariantMapper::new(&hdp);
    let hgvs_weaver::SequenceVariant::Coding(c) = parse_hgvs_variant("NM_X.1:c.(4_6)del").unwrap()
    else {
        panic!("not c.");
    };
    let result = mapper.c_to_g(&c, None);
    assert!(
        matches!(result, Err(HgvsError::UnsupportedOperation(_))),
        "{result:?}"
    );
}

#[test]
fn a_coding_deletion_over_an_uncertain_range_is_not_a_protein_change() {
    // c.4_6del, exact, would be p.(Lys2del).
    assert_eq!(
        placed_level("NM_X.1:c.(4_6)del", "NP_X.1:p.(Lys2del)"),
        EquivalenceLevel::Different
    );
}

#[test]
fn a_coding_deletion_over_an_uncertain_range_predicts_no_protein() {
    let hdp = placed_provider();
    let mapper = VariantMapper::new(&hdp);
    let hgvs_weaver::SequenceVariant::Coding(c) = parse_hgvs_variant("NM_X.1:c.(4_6)del").unwrap()
    else {
        panic!("not c.");
    };
    let result = mapper.c_to_p(&c, None);
    assert!(
        matches!(result, Err(HgvsError::UnsupportedOperation(_))),
        "{result:?}"
    );
}

#[test]
fn one_base_deleted_somewhere_in_a_range_has_no_spdi() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let var = parse_hgvs_variant("NC_TEST.1:g.(1003_1004)del").unwrap();
    let result = mapper.to_spdi(&var, false);
    assert!(
        matches!(result, Err(HgvsError::UnsupportedOperation(_))),
        "{result:?}"
    );
}

#[test]
fn a_deletion_between_parenthesised_single_bases_has_the_exact_ones_spdi() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let spdi = |s: &str| {
        mapper
            .to_spdi(&parse_hgvs_variant(s).unwrap(), false)
            .unwrap()
    };
    assert_eq!(
        spdi("NC_TEST.1:g.(1005)_(1012)del"),
        spdi("NC_TEST.1:g.1005_1012del")
    );
}

#[test]
fn a_genomic_deletion_somewhere_in_a_range_has_no_transcript_form() {
    let hdp = placed_provider();
    let mapper = VariantMapper::new(&hdp);
    let hgvs_weaver::SequenceVariant::Genomic(g) =
        parse_hgvs_variant("NC_X.1:g.(109_110)del").unwrap()
    else {
        panic!("not g.");
    };
    let result = mapper.g_to_c(&g, "NM_X.1");
    assert!(
        matches!(result, Err(HgvsError::UnsupportedOperation(_))),
        "{result:?}"
    );
}

#[test]
fn a_copy_number_between_parenthesised_single_bases_is_the_exact_one() {
    assert_eq!(
        genomic_level(
            "NC_TEST.1:g.(1005)_(1012)copy3",
            "NC_TEST.1:g.1005_1012copy3"
        ),
        EquivalenceLevel::Analogous
    );
}

#[test]
fn a_deletion_between_parenthesised_single_bases_is_the_exact_one() {
    assert_eq!(
        genomic_level("NC_TEST.1:g.(1005)_(1012)del", "NC_TEST.1:g.1005_1012del"),
        EquivalenceLevel::Analogous
    );
}

#[test]
fn a_protein_deletion_over_an_uncertain_range_is_not_the_exact_one() {
    assert_eq!(
        coding_level("NP_X.1:p.(Lys2_Leu3)del", "NP_X.1:p.Lys2_Leu3del"),
        EquivalenceLevel::Different
    );
}

#[test]
fn a_predicted_protein_change_is_its_confirmation() {
    assert_eq!(
        coding_level("NP_X.1:p.(Lys2Asn)", "NP_X.1:p.Lys2Asn"),
        EquivalenceLevel::Identity
    );
}

#[test]
fn a_cis_allele_of_predicted_protein_changes_is_its_confirmation() {
    assert_eq!(
        coding_level(
            "NP_X.1:p.[(Lys2Asn);(Leu3Pro)]",
            "NP_X.1:p.[Lys2Asn;Leu3Pro]"
        ),
        EquivalenceLevel::Identity
    );
}

#[test]
fn a_coding_change_is_its_prediction_written_unbracketed() {
    // c.7C>G turns CTG into GTG: p.(Leu3Val), confirmed.
    assert_eq!(
        coding_level("NM_X.1:c.7C>G", "NP_X.1:p.Leu3Val"),
        EquivalenceLevel::Identity
    );
}

#[test]
fn a_cis_allele_with_an_uncertain_coding_member_has_no_vrs_on_the_genome() {
    let hdp = placed_provider();
    let mapper = VariantMapper::new(&hdp);
    let var = parse_hgvs_variant("NM_X.1:c.[(4_6)del;10G>A]").unwrap();
    let result = mapper.to_vrs_variation(&var);
    assert!(
        matches!(
            &result,
            Err(HgvsError::CisMember { index: 0, source, .. })
                if matches!(**source, HgvsError::UnsupportedOperation(_))
        ),
        "{result:?}"
    );
}

#[test]
fn a_cis_allele_with_a_member_somewhere_in_a_range_is_not_the_exact_allele() {
    // g.115 is the G of CTG.
    assert_eq!(
        placed_level(
            "NC_X.1:g.[(109_110)del;115G>A]",
            "NC_X.1:g.[109_110del;115G>A]"
        ),
        EquivalenceLevel::Different
    );
}

#[test]
fn a_cis_allele_with_a_member_between_parenthesised_single_bases_normalises_as_the_exact_one() {
    let hdp = placed_provider();
    let mapper = VariantMapper::new(&hdp);
    let normalised = |s: &str| {
        mapper
            .normalize_variant(parse_hgvs_variant(s).unwrap())
            .unwrap()
            .to_string()
    };
    assert_eq!(
        normalised("NC_X.1:g.[(105)_(112)del;115G>A]"),
        normalised("NC_X.1:g.[105_112del;115G>A]")
    );
}

// --- Alleles in cis with an uncertain member ---
//
// The members of an allele are unordered: an uncertain member compares by its
// spelling, an exact one by its canonical allele, and a member written twice
// counts once. On `NC_TEST.1`, g.1001 is A and g.1007 G; on `NM_X.1`, c.10
// is the G of GCC and c.13 the T of TAT; `NP_X.1` is MKLAYR.

#[test]
fn a_cis_allele_with_an_uncertain_member_is_itself_in_another_order() {
    assert_eq!(
        genomic_level(
            "NC_TEST.1:g.[(1003_1004)del;1007G>A]",
            "NC_TEST.1:g.[1007G>A;(1003_1004)del]"
        ),
        EquivalenceLevel::Identity
    );
}

#[test]
fn a_coding_cis_allele_with_an_uncertain_member_is_itself_in_another_order() {
    assert_eq!(
        placed_level("NM_X.1:c.[(4_6)del;10G>A]", "NM_X.1:c.[10G>A;(4_6)del]"),
        EquivalenceLevel::Identity
    );
}

#[test]
fn a_cis_allele_with_an_uncertain_member_in_another_order_is_not_the_exact_allele() {
    assert_eq!(
        genomic_level(
            "NC_TEST.1:g.[(1003_1004)del;1007G>A]",
            "NC_TEST.1:g.[1007G>A;1003del]"
        ),
        EquivalenceLevel::Different
    );
}

#[test]
fn a_coding_cis_allele_with_an_uncertain_member_is_not_the_exact_allele() {
    assert_eq!(
        placed_level("NM_X.1:c.[(4_6)del;10G>A]", "NM_X.1:c.[10G>A;4_6del]"),
        EquivalenceLevel::Different
    );
}

#[test]
fn cis_alleles_with_different_uncertain_members_differ() {
    assert_eq!(
        genomic_level(
            "NC_TEST.1:g.[(1003_1004)del;1007G>A]",
            "NC_TEST.1:g.[1007G>A;(1003_1005)del]"
        ),
        EquivalenceLevel::Different
    );
}

#[test]
fn an_uncertain_member_written_twice_counts_once() {
    assert_eq!(
        genomic_level(
            "NC_TEST.1:g.[(1003_1004)del;1007G>A;(1003_1004)del]",
            "NC_TEST.1:g.[1007G>A;(1003_1004)del]"
        ),
        EquivalenceLevel::Analogous
    );
}

#[test]
fn the_exact_members_beside_an_uncertain_one_compare_by_their_alleles() {
    // g.1001dup and g.1001_1002insA insert one A after g.1001.
    assert_eq!(
        genomic_level(
            "NC_TEST.1:g.[(1003_1004)del;1001_1002insA]",
            "NC_TEST.1:g.[1001dup;(1003_1004)del]"
        ),
        EquivalenceLevel::Analogous
    );
}

/// `NM_C.1` and `NP_C.1`: a five-base 5'UTR, then `cds`, then the fourteen-base
/// 3'UTR of `NM_X.1`; `protein` is the protein the record names.
fn record_provider(cds: &str, protein: &str) -> Provider {
    let transcript = format!("GGGGG{cds}CCGTATAAGTAAGG");
    Provider::new()
        .sequence("NM_C.1", &transcript)
        .sequence("NP_C.1", protein)
        .transcript(single_exon_transcript(
            "NM_C.1",
            "NC_C.1",
            100,
            Strand::Plus,
            5,
            5 + cds.len() as i32 - 1,
            transcript.len() as i32,
        ))
        .protein_for("NM_C.1", "NP_C.1")
}

/// The level at which `a` and `b` are judged the same over `hdp`.
fn level_over(hdp: &Provider, a: &str, b: &str) -> EquivalenceLevel {
    let mapper = VariantMapper::new(hdp);
    let parse = |s: &str| parse_hgvs_variant(s).unwrap_or_else(|e| panic!("{s}: {e}"));
    VariantEquivalence::new(&mapper, hdp)
        .equivalent_level(&parse(a), &parse(b))
        .unwrap()
}

/// A CDS that starts at CTG, as some do; the record reads the initiator as Met,
/// `MKLAYR`.
const CUG_START: &str = "CTGAAACTGGCCTATCGCTAA";

#[test]
fn a_frameshift_from_a_non_aug_start_is_the_frameshift_written_without_its_length() {
    // c.6dup reads CTG AAA ACT GGC CTA TCG CTA ACC GTA TAA: p.(Leu3ThrfsTer8).
    let hdp = record_provider(CUG_START, "MKLAYR");
    assert_eq!(
        level_over(&hdp, "NM_C.1:c.6dup", "NP_C.1:p.Leu3fs"),
        EquivalenceLevel::Analogous
    );
}

#[test]
fn a_frameshift_from_a_non_aug_start_is_not_one_from_a_later_residue() {
    let hdp = record_provider(CUG_START, "MKLAYR");
    assert_eq!(
        level_over(&hdp, "NM_C.1:c.6dup", "NP_C.1:p.Ala4fs"),
        EquivalenceLevel::Different
    );
}

#[test]
fn a_non_aug_start_predicts_the_protein_the_record_names() {
    let hdp = record_provider(CUG_START, "MKLAYR");
    let mapper = VariantMapper::new(&hdp);
    let SequenceVariant::Coding(c) = parse_hgvs_variant("NM_C.1:c.6dup").unwrap() else {
        panic!("not c.");
    };
    assert_eq!(
        mapper
            .predicted_protein(&c, None)
            .unwrap()
            .map(|p| p.residues)
            .as_deref(),
        Some("MKTGLSLTV")
    );
}

/// A selenoprotein: codon 3 is TGA, selenocysteine in the record, `MKUAYR`.
const SEC_AT_3: &str = "ATGAAATGAGCCTATCGCTAA";

#[test]
fn a_frameshift_past_a_selenocysteine_is_the_frameshift_written_without_its_length() {
    // c.12dup reads ATG AAA TGA GCC CTA TCG CTA ACC GTA TAA: p.(Tyr5LeufsTer6).
    let hdp = record_provider(SEC_AT_3, "MKUAYR");
    assert_eq!(
        level_over(&hdp, "NM_C.1:c.12dup", "NP_C.1:p.Tyr5fs"),
        EquivalenceLevel::Analogous
    );
}

#[test]
fn the_exact_coding_members_beside_an_uncertain_one_compare_by_their_alleles() {
    // c.13dup and c.13_14insT insert one T after c.13.
    assert_eq!(
        placed_level("NM_X.1:c.[(4_6)del;13dup]", "NM_X.1:c.[13_14insT;(4_6)del]"),
        EquivalenceLevel::Analogous
    );
}

#[test]
fn a_mitochondrial_uncertain_member_is_its_genomic_spelling() {
    assert_eq!(
        genomic_level(
            "NC_TEST.1:m.[(1003_1004)del;1007G>A]",
            "NC_TEST.1:g.[1007G>A;(1003_1004)del]"
        ),
        EquivalenceLevel::Analogous
    );
}

#[test]
fn an_rna_uncertain_member_is_its_coding_spelling() {
    assert_eq!(
        placed_level("NM_X.1:r.[(4_6)del;10g>a]", "NM_X.1:c.[10G>A;(4_6)del]"),
        EquivalenceLevel::Analogous
    );
}

#[test]
fn a_protein_cis_allele_with_an_uncertain_member_is_itself_in_another_order() {
    assert_eq!(
        coding_level(
            "NP_X.1:p.[(Lys2_Leu3)del;Tyr5Asp]",
            "NP_X.1:p.[Y5D;(K2_L3)del]"
        ),
        EquivalenceLevel::Identity
    );
}

#[test]
fn the_exact_protein_members_beside_an_uncertain_one_compare_by_their_alleles() {
    assert_eq!(
        coding_level(
            "NP_X.1:p.[(Lys2_Leu3)del;Tyr5Asp]",
            "NP_X.1:p.[Tyr5delinsAsp;(Lys2_Leu3)del]"
        ),
        EquivalenceLevel::Analogous
    );
}

#[test]
fn exact_members_that_differ_beside_an_uncertain_one_differ() {
    assert_eq!(
        genomic_level(
            "NC_TEST.1:g.[(1003_1004)del;1007G>A]",
            "NC_TEST.1:g.[(1003_1004)del;1007G>T]"
        ),
        EquivalenceLevel::Different
    );
}

#[test]
fn a_selenoprotein_predicts_the_protein_the_record_names() {
    let hdp = record_provider(SEC_AT_3, "MKUAYR");
    let mapper = VariantMapper::new(&hdp);
    let SequenceVariant::Coding(c) = parse_hgvs_variant("NM_C.1:c.12dup").unwrap() else {
        panic!("not c.");
    };
    assert_eq!(
        mapper
            .predicted_protein(&c, None)
            .unwrap()
            .map(|p| p.residues)
            .as_deref(),
        Some("MKUALSLTV")
    );
}

#[test]
fn a_change_to_a_non_aug_start_codon_is_an_unknown_change_from_met1() {
    // c.1C>T turns the initiator CTG into TTG.
    let hdp = record_provider(CUG_START, "MKLAYR");
    let mapper = VariantMapper::new(&hdp);
    assert_eq!(
        mapper
            .c_to_p(&coding_on("NM_C.1:c.1C>T"), None)
            .unwrap()
            .to_string(),
        "NP_C.1:p.(Met1?)"
    );
}

#[test]
fn a_non_aug_start_codon_made_atg_is_an_unknown_change_from_met1() {
    // c.1C>A turns the initiator CTG into ATG: both read as Met, yet where
    // translation starts is not predicted.
    let hdp = record_provider(CUG_START, "MKLAYR");
    let mapper = VariantMapper::new(&hdp);
    assert_eq!(
        mapper
            .c_to_p(&coding_on("NM_C.1:c.1C>A"), None)
            .unwrap()
            .to_string(),
        "NP_C.1:p.(Met1?)"
    );
}

/// A CDS that starts at ATG, with the protein of [`CUG_START`].
const AUG_START: &str = "ATGAAACTGGCCTATCGCTAA";

#[test]
fn a_start_codon_change_is_identical_to_the_unknown_change_from_met1_stated() {
    let hdp = record_provider(AUG_START, "MKLAYR");
    assert_eq!(
        level_over(&hdp, "NM_C.1:c.1A>G", "NP_C.1:p.Met1?"),
        EquivalenceLevel::Identity
    );
}

#[test]
fn a_start_codon_change_is_identical_to_the_unknown_change_from_met1_predicted() {
    let hdp = record_provider(AUG_START, "MKLAYR");
    assert_eq!(
        level_over(&hdp, "NM_C.1:c.1A>G", "NP_C.1:p.(Met1?)"),
        EquivalenceLevel::Identity
    );
}

#[test]
fn a_start_codon_change_is_analogous_to_a_substitution_at_met1() {
    // p.Met1Val changes the protein from residue 1, where p.(Met1?) says a
    // change starts.
    let hdp = record_provider(AUG_START, "MKLAYR");
    assert_eq!(
        level_over(&hdp, "NM_C.1:c.1A>G", "NP_C.1:p.Met1Val"),
        EquivalenceLevel::Analogous
    );
}

/// The coding variant `c`.
fn coding_on(c: &str) -> hgvs_weaver::CVariant {
    match parse_hgvs_variant(c).unwrap() {
        SequenceVariant::Coding(c) => c,
        other => panic!("{other} is not c."),
    }
}

#[test]
fn a_change_at_a_selenocysteine_is_named_at_sec() {
    // c.7T>C turns the TGA of codon 3 into CGA, Arg.
    let hdp = record_provider(SEC_AT_3, "MKUAYR");
    let mapper = VariantMapper::new(&hdp);
    assert_eq!(
        mapper
            .c_to_p(&coding_on("NM_C.1:c.7T>C"), None)
            .unwrap()
            .to_string(),
        "NP_C.1:p.(Sec3Arg)"
    );
}

#[test]
fn a_change_at_a_selenocysteine_is_its_prediction_confirmed() {
    let hdp = record_provider(SEC_AT_3, "MKUAYR");
    assert_eq!(
        level_over(&hdp, "NM_C.1:c.7T>C", "NP_C.1:p.Sec3Arg"),
        EquivalenceLevel::Identity
    );
}

#[test]
fn a_selenocysteine_unchanged_is_named_at_sec() {
    let hdp = record_provider(SEC_AT_3, "MKUAYR");
    let mapper = VariantMapper::new(&hdp);
    assert_eq!(
        mapper
            .c_to_p(&coding_on("NM_C.1:c.7_9="), None)
            .unwrap()
            .to_string(),
        "NP_C.1:p.(Sec3=)"
    );
}

#[test]
fn a_non_aug_initiator_unchanged_is_met1_unchanged() {
    let hdp = record_provider(CUG_START, "MKLAYR");
    assert_eq!(
        level_over(&hdp, "NM_C.1:c.1_3=", "NP_C.1:p.Met1="),
        EquivalenceLevel::Identity
    );
}

#[test]
fn a_non_aug_start_has_a_protein_allele() {
    let hdp = record_provider(CUG_START, "MKLAYR");
    let mapper = VariantMapper::new(&hdp);
    assert!(mapper
        .protein_allele(&coding_on("NM_C.1:c.6dup"), None)
        .is_ok());
}

#[test]
fn an_aaa_annotated_start_reads_as_met() {
    // c.6dup reads AAA AAA ACT GGC CTA TCG CTA ACC GTA TAA; the record starts with M.
    let hdp = record_provider("AAAAAACTGGCCTATCGCTAA", "MKLAYR");
    let mapper = VariantMapper::new(&hdp);
    assert_eq!(
        mapper
            .predicted_protein(&coding_on("NM_C.1:c.6dup"), None)
            .unwrap()
            .map(|p| p.residues)
            .as_deref(),
        Some("MKTGLSLTV")
    );
}

#[test]
fn a_stop_made_of_a_selenocysteine_codon_is_a_stop() {
    // c.8G>A turns the TGA of codon 3 into TAA.
    let hdp = record_provider(SEC_AT_3, "MKUAYR");
    let mapper = VariantMapper::new(&hdp);
    assert_eq!(
        mapper
            .c_to_p(&coding_on("NM_C.1:c.8G>A"), None)
            .unwrap()
            .to_string(),
        "NP_C.1:p.(Sec3Ter)"
    );
}

#[test]
fn a_stop_inserted_over_a_selenocysteine_codon_ends_the_protein() {
    // AAA TGA becomes GGG TAA.
    let hdp = record_provider(SEC_AT_3, "MKUAYR");
    let mapper = VariantMapper::new(&hdp);
    assert_eq!(
        mapper
            .predicted_protein(&coding_on("NM_C.1:c.4_9delinsGGGTAA"), None)
            .unwrap()
            .map(|p| p.residues)
            .as_deref(),
        Some("MG")
    );
}

#[test]
fn a_selenocysteine_after_an_in_frame_deletion_stays_one() {
    let hdp = record_provider(SEC_AT_3, "MKUAYR");
    let mapper = VariantMapper::new(&hdp);
    assert_eq!(
        mapper
            .predicted_protein(&coding_on("NM_C.1:c.4_6del"), None)
            .unwrap()
            .map(|p| p.residues)
            .as_deref(),
        Some("MUAYR")
    );
}

#[test]
fn a_description_that_does_not_read_on_its_protein_is_an_error() {
    // Residue 3 of MKLAYR is Leu, not Arg; c.7C>G predicts p.(Leu3Val).
    let hdp = coding_provider();
    let mapper = VariantMapper::new(&hdp);
    let result = VariantEquivalence::new(&mapper, &hdp).equivalent_level(
        &parse_hgvs_variant("NM_X.1:c.7C>G").unwrap(),
        &parse_hgvs_variant("NP_X.1:p.Arg3Gly").unwrap(),
    );
    assert!(result.is_err(), "{result:?}");
}

/// `NC_E.1` holding `NM_X.1`'s sequence from g.1, with two transcripts on it:
/// `NM_E.1`, and `NM_F.1`, the same with its CDS open at the 3' end, on which
/// no frameshift is predicted. Both encode `NP_E.1`.
fn two_transcript_provider() -> Provider {
    let transcript = "GGGGGATGAAACTGGCCTATCGCTAACCGTATAAGTAAGG";
    let len = transcript.len() as i32;
    Provider::new()
        .sequence("NC_E.1", transcript)
        .sequence("NM_E.1", transcript)
        .sequence("NM_F.1", transcript)
        .sequence("NP_E.1", "MKLAYR")
        .transcript(single_exon_transcript(
            "NM_E.1",
            "NC_E.1",
            0,
            Strand::Plus,
            5,
            25,
            len,
        ))
        .transcript(TranscriptData {
            cds_end_open: true,
            ..single_exon_transcript("NM_F.1", "NC_E.1", 0, Strand::Plus, 5, 25, len)
        })
        .protein_for("NM_E.1", "NP_E.1")
        .protein_for("NM_F.1", "NP_E.1")
}

#[test]
fn a_transcript_that_cannot_predict_does_not_hide_one_that_agrees() {
    // g.12del is c.7del on both: p.(Leu3TrpfsTer?) on NM_E.1, nothing on NM_F.1.
    let hdp = two_transcript_provider();
    assert_eq!(
        level_over(&hdp, "NC_E.1:g.12del", "NP_E.1:p.Leu3fs"),
        EquivalenceLevel::Analogous
    );
}

#[test]
fn a_transcript_that_cannot_predict_is_the_error_when_none_agrees() {
    let hdp = two_transcript_provider();
    let mapper = VariantMapper::new(&hdp);
    let SequenceVariant::Genomic(g) = parse_hgvs_variant("NC_E.1:g.12del").unwrap() else {
        panic!("not g.");
    };
    let on = mapper.g_to_c_all(&g, &hdp).unwrap();
    assert_eq!(on.len(), 2, "g.12del lies on both transcripts: {on:?}");
    let result = VariantEquivalence::new(&mapper, &hdp).equivalent_level(
        &SequenceVariant::Genomic(g),
        &parse_hgvs_variant("NP_E.1:p.Ala4fs").unwrap(),
    );
    assert!(result.is_err(), "{result:?}");
}

#[test]
fn a_codon_deleted_before_a_selenocysteine_in_any_spelling_is_a_deletion() {
    // In ATG TCT TGA, c.5_7del (CTT) and c.4_6del (TCT) leave the same ATG TGA.
    let hdp = record_provider("ATGTCTTGAGCCTATCGCTAA", "MSUAYR");
    assert_eq!(
        level_over(&hdp, "NM_C.1:c.5_7del", "NP_C.1:p.Ser2del"),
        EquivalenceLevel::Identity
    );
}

#[test]
fn a_delins_that_rewrites_a_selenocysteine_codon_as_itself_keeps_it() {
    // c.9_12delinsAGCT is c.12C>T: the codons it touches, Sec3 and Ala4, unchanged.
    let hdp = record_provider(SEC_AT_3, "MKUAYR");
    let mapper = VariantMapper::new(&hdp);
    assert_eq!(
        mapper
            .c_to_p(&coding_on("NM_C.1:c.9_12delinsAGCT"), None)
            .unwrap()
            .to_string(),
        "NP_C.1:p.(Sec3_Ala4=)"
    );
}

#[test]
fn a_duplicated_selenocysteine_codon_is_a_selenocysteine() {
    let hdp = record_provider(SEC_AT_3, "MKUAYR");
    let mapper = VariantMapper::new(&hdp);
    assert_eq!(
        mapper
            .c_to_p(&coding_on("NM_C.1:c.7_9dup"), None)
            .unwrap()
            .to_string(),
        "NP_C.1:p.(Sec3dup)"
    );
}

#[test]
fn a_selenocysteine_codon_inserted_after_its_copy_is_the_same_duplication() {
    // c.9_10insTGA is c.7_9dup written at the other end of the copy.
    let hdp = record_provider(SEC_AT_3, "MKUAYR");
    assert_eq!(
        level_over(&hdp, "NM_C.1:c.9_10insTGA", "NP_C.1:p.Sec3dup"),
        EquivalenceLevel::Identity
    );
}
