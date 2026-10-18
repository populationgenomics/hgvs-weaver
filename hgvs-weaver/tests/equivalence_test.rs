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
        EquivalenceLevel::Analogous
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
fn a_coding_change_is_analogous_to_its_prediction_written_unbracketed() {
    // c.7C>G turns CTG into GTG.
    assert_eq!(
        coding_level("NM_X.1:c.7C>G", "NP_X.1:p.Leu3Val"),
        EquivalenceLevel::Analogous
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
