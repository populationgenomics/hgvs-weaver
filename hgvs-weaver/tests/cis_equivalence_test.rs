//! Equivalence of alleles in cis, `c.[a;b]`: two alleles, or an allele and a
//! plain variant, are the same change when the members make the same sequence
//! on the genome and state the same bases unchanged; choices.md, "A cis allele
//! is the sequence its members make".

mod support;

use hgvs_weaver::data::Strand;
use hgvs_weaver::equivalence::{EquivalenceLevel, VariantEquivalence};
use hgvs_weaver::error::HgvsError;
use hgvs_weaver::mapper::VariantMapper;
use hgvs_weaver::parse_hgvs_variant;
use support::{single_exon_transcript, Provider};

use EquivalenceLevel::{Analogous, Different, Identity};

const UTR5: &str = "GGGGG";
/// M K L A Y R *
const CDS: &str = "ATGAAACTGGCCTATCGCTAA";
const UTR3: &str = "CCGTATAAGTAAGG";

/// `NM_X.1` is `UTR5 CDS UTR3` at genomic index 10 of `NC_X.1`, so c.7 is
/// g.22 and c.4_6 is the run `AAA`; `NR_XR.1` reads the same forty bases on
/// the minus strand, so n.23 is g.28 and n.29 is g.22; `NM_S.1` is `NM_X.1`
/// with exon 1 ending at c.5 and the intron `GTAAGTCCAG` after it in `NC_S.1`,
/// so c.5 is g.20, c.5+1 is g.21 and c.6 is g.31; `NC_R.1` is `TTCAGCAGTT`
/// `NC_T.1` ten Ts; `NC_TEST.1` is `ACGT` repeated, so g.1007 is a G.
fn provider() -> Provider {
    let transcript = format!("{UTR5}{CDS}{UTR3}");
    let len = transcript.len() as i32;
    let cds = (UTR5.len() as i32, (UTR5.len() + CDS.len()) as i32 - 1);
    Provider::new()
        .sequence("NC_X.1", &format!("TTTTTTTTTT{transcript}CCCCCCCCCC"))
        .sequence(
            "NC_S.1",
            &format!(
                "TTTTTTTTTT{}GTAAGTCCAG{}CCCCCCCCCC",
                &transcript[..10],
                &transcript[10..]
            ),
        )
        .sequence("NC_R.1", "TTCAGCAGTT")
        .sequence("NC_T.1", "TTTTTTTTTT")
        .sequence("NC_TEST.1", &"ACGT".repeat(300))
        .sequence("NM_X.1", &transcript)
        .sequence("NM_S.1", &transcript)
        .sequence(
            "NR_XR.1",
            &hgvs_weaver::utils::reverse_complement(&transcript),
        )
        .sequence("NP_X.1", "MKLAYR")
        .transcript(single_exon_transcript(
            "NM_X.1",
            "NC_X.1",
            10,
            Strand::Plus,
            cds.0,
            cds.1,
            len,
        ))
        .transcript(support::transcript(
            "NM_S.1",
            "NC_S.1",
            Strand::Plus,
            Some(cds),
            vec![
                support::exon((0, 10), (10, 19), Strand::Plus),
                support::exon((10, len), (30, 30 + len - 10 - 1), Strand::Plus),
            ],
        ))
        .transcript(support::transcript(
            "NR_XR.1",
            "NC_X.1",
            Strand::Minus,
            None,
            vec![support::exon((0, len), (10, 10 + len - 1), Strand::Minus)],
        ))
        .protein_for("NM_X.1", "NP_X.1")
}

fn judge(a: &str, b: &str) -> Result<EquivalenceLevel, HgvsError> {
    let hdp = provider();
    let mapper = VariantMapper::with_refget(&hdp, &hdp);
    let equivalence = VariantEquivalence::new(&mapper, &hdp);
    equivalence.equivalent_level(
        &parse_hgvs_variant(a).unwrap(),
        &parse_hgvs_variant(b).unwrap(),
    )
}

fn level(a: &str, b: &str) -> EquivalenceLevel {
    judge(a, b).unwrap_or_else(|e| panic!("{a} vs {b}: {e}"))
}

// --- The sequence the members make ---

#[test]
fn adjacent_substitutions_are_the_delins_they_make() {
    assert_eq!(
        level("NM_X.1:c.[7C>T;8T>G]", "NM_X.1:c.7_8delinsTG"),
        Analogous
    );
}

#[test]
fn a_delins_is_the_adjacent_substitutions_it_makes() {
    assert_eq!(
        level("NM_X.1:c.7_8delinsTG", "NM_X.1:c.[7C>T;8T>G]"),
        Analogous
    );
}

#[test]
fn adjacent_substitutions_are_not_a_delins_making_another_sequence() {
    assert_eq!(
        level("NM_X.1:c.[7C>T;8T>G]", "NM_X.1:c.7_8delinsTA"),
        Different
    );
}

#[test]
fn an_allele_of_one_delins_is_the_substitutions_it_makes() {
    assert_eq!(
        level("NM_X.1:c.[7_8delinsTG;13T>G]", "NM_X.1:c.[7C>T;8T>G;13T>G]"),
        Analogous
    );
}

#[test]
fn a_delins_over_an_unchanged_base_is_its_changed_members() {
    // c.4_6 is AAA: the delins leaves c.5 as it is and states nothing of it.
    assert_eq!(
        level("NM_X.1:c.[4A>T;6A>T]", "NM_X.1:c.4_6delinsTAT"),
        Analogous
    );
}

#[test]
fn a_change_split_differently_between_members_is_one_change() {
    // On TTCAGCAGTT both make TTCTGTT, though 7A>T and 4A>T differ alone.
    assert_eq!(
        level("NC_R.1:g.[3_5del;7A>T]", "NC_R.1:g.[4A>T;6_8del]"),
        Analogous
    );
}

#[test]
fn a_duplication_in_an_allele_is_the_insertion_it_makes() {
    // c.6 is the last A of the run c.4_6.
    assert_eq!(
        level("NM_X.1:c.[6dup;13T>G]", "NM_X.1:c.[6_7insA;13T>G]"),
        Analogous
    );
}

#[test]
fn a_duplication_in_an_allele_is_one_written_elsewhere_in_its_run() {
    assert_eq!(
        level("NM_X.1:c.[4dup;13T>G]", "NM_X.1:c.[6_7insA;13T>G]"),
        Analogous
    );
}

#[test]
fn a_transcript_allele_is_its_genomic_spelling_unnormalised() {
    // c.4del is g.19del, in the run g.19_21.
    assert_eq!(
        level("NM_X.1:c.[4del;13T>G]", "NC_X.1:g.[21del;28T>G]"),
        Analogous
    );
}

#[test]
fn a_mitochondrial_allele_is_identical_to_the_genomic_one_it_spells() {
    assert_eq!(
        level("NC_X.1:m.[22C>T;28T>G]", "NC_X.1:g.[28T>G;22C>T]"),
        Identity
    );
}

#[test]
fn a_mitochondrial_allele_is_identical_to_its_coding_projection() {
    assert_eq!(
        level("NC_X.1:m.[22C>T;28T>G]", "NM_X.1:c.[7C>T;13T>G]"),
        Identity
    );
}

#[test]
fn a_mitochondrial_allele_is_the_genomic_delins_it_makes() {
    assert_eq!(
        level("NC_X.1:m.[22C>T;23T>G]", "NC_X.1:g.22_23delinsTG"),
        Analogous
    );
}

#[test]
fn an_allele_is_not_one_of_its_members() {
    assert_eq!(level("NM_X.1:c.[7C>T;13T>G]", "NM_X.1:c.7C>T"), Different);
}

#[test]
fn an_allele_is_not_one_with_a_member_more() {
    assert_eq!(
        level("NM_X.1:c.[7C>T;13T>G]", "NM_X.1:c.[7C>T;13T>G;20del]"),
        Different
    );
}

#[test]
fn an_allele_is_not_one_with_a_member_changed() {
    assert_eq!(
        level("NM_X.1:c.[7C>T;13T>G]", "NM_X.1:c.[7C>T;13T>A]"),
        Different
    );
}

#[test]
fn an_allele_of_one_member_is_that_member() {
    assert_eq!(level("NM_X.1:c.[7C>T]", "NM_X.1:c.7C>T"), Identity);
}

#[test]
fn an_allele_of_one_member_is_its_members_genomic_allele() {
    // c.4del is g.19del, in the run g.19_21.
    assert_eq!(level("NM_X.1:c.[4del]", "NC_X.1:g.21del"), Analogous);
}

#[test]
fn an_allele_of_one_member_is_identical_to_its_exact_projection() {
    assert_eq!(level("NC_X.1:g.22C>T", "NM_X.1:c.[7C>T]"), Identity);
}

#[test]
fn an_allele_of_one_member_is_not_another_variant() {
    assert_eq!(level("NM_X.1:c.[7C>T]", "NM_X.1:c.13T>G"), Different);
}

#[test]
fn a_plain_variant_is_the_allele_making_its_sequence() {
    assert_eq!(
        level("NM_X.1:c.7_8delinsTG", "NC_X.1:g.[23T>G;22C>T]"),
        Analogous
    );
}

#[test]
fn a_mitochondrial_allele_of_one_member_is_identical_to_its_projection() {
    assert_eq!(level("NC_X.1:m.[22C>T]", "NM_X.1:c.7C>T"), Identity);
}

#[test]
fn an_allele_of_one_member_against_several_is_judged_by_its_member() {
    assert_eq!(
        level("NM_X.1:c.[7_8delinsTG]", "NM_X.1:c.[7C>T;8T>G]"),
        Analogous
    );
}

// --- Identity ---

#[test]
fn the_same_members_in_another_order_are_one_spelling() {
    assert_eq!(
        level("NM_X.1:c.[13T>G;7C>T]", "NM_X.1:c.[7C>T;13T>G]"),
        Identity
    );
}

#[test]
fn a_transcript_allele_is_identical_to_its_exact_projection() {
    assert_eq!(
        level("NM_X.1:c.[7C>T;13T>G]", "NC_X.1:g.[22C>T;28T>G]"),
        Identity
    );
}

#[test]
fn a_transcript_allele_is_identical_to_its_projection_in_another_order() {
    assert_eq!(
        level("NC_X.1:g.[28T>G;22C>T]", "NM_X.1:c.[7C>T;13T>G]"),
        Identity
    );
}

// --- The minus strand ---

#[test]
fn a_minus_strand_allele_is_identical_to_its_projection() {
    assert_eq!(
        level("NR_XR.1:n.[23A>C;29G>A]", "NC_X.1:g.[22C>T;28T>G]"),
        Identity
    );
}

#[test]
fn adjacent_minus_strand_members_are_the_genomic_delins_they_make() {
    // n.28 and n.29 are g.23 and g.22: A>C and G>A there are T>G and C>T.
    assert_eq!(
        level("NR_XR.1:n.[28A>C;29G>A]", "NC_X.1:g.22_23delinsTG"),
        Analogous
    );
}

#[test]
fn a_minus_strand_allele_is_the_plus_strand_change_it_makes() {
    assert_eq!(
        level("NR_XR.1:n.[28A>C;29G>A]", "NM_X.1:c.7_8delinsTG"),
        Analogous
    );
}

// --- Intronic members ---

#[test]
fn intronic_members_that_meet_are_the_delins_they_make() {
    // c.5+1_5+2 is the GT of the intron, g.21_22.
    assert_eq!(
        level("NM_S.1:c.[5+1G>A;5+2T>C]", "NM_S.1:c.5+1_5+2delinsAC"),
        Analogous
    );
}

#[test]
fn intronic_members_that_meet_are_the_genomic_delins() {
    assert_eq!(
        level("NM_S.1:c.[5+1G>A;5+2T>C]", "NC_S.1:g.21_22delinsAC"),
        Analogous
    );
}

#[test]
fn an_exonic_member_meeting_an_intronic_one_is_the_delins_across_the_border() {
    assert_eq!(
        level("NM_S.1:c.[5A>T;5+1G>A]", "NM_S.1:c.5_5+1delinsTA"),
        Analogous
    );
}

#[test]
fn members_either_side_of_an_intron_are_not_a_delins_over_it() {
    // c.4_6delinsTAT projects over the intron, g.19_31delinsTAT, deleting it.
    assert_eq!(
        level("NM_S.1:c.[4A>T;6A>T]", "NM_S.1:c.4_6delinsTAT"),
        Different
    );
}

// --- Identity members ---

#[test]
fn an_identity_member_states_its_bases_beside_a_change() {
    assert_eq!(level("NM_X.1:c.[7C>T;10=]", "NM_X.1:c.7C>T"), Different);
}

#[test]
fn an_identity_member_adds_no_change() {
    // The change beside it is compared by the sequence it makes.
    assert_eq!(
        level("NM_X.1:c.[7C>T;10=]", "NM_X.1:c.[7delinsT;10=]"),
        Analogous
    );
}

#[test]
fn an_identity_member_inside_a_delins_states_what_the_delins_does_not() {
    assert_eq!(
        level("NM_X.1:c.[4A>T;5=;6A>T]", "NM_X.1:c.4_6delinsTAT"),
        Different
    );
}

#[test]
fn an_identity_member_is_identical_to_its_projection() {
    assert_eq!(
        level("NM_X.1:c.[7C>T;10=]", "NC_X.1:g.[22C>T;25=]"),
        Identity
    );
}

// --- Members that cancel ---

#[test]
fn cancelling_members_are_the_identity_of_the_bases_they_name() {
    assert_eq!(level("NM_X.1:c.[4_5insA;5del]", "NM_X.1:c.4_5="), Analogous);
}

#[test]
fn cancelling_members_do_not_name_the_rest_of_their_run() {
    // c.6 is in the run c.4_6, but no member names it.
    assert_eq!(level("NM_X.1:c.[4_5insA;5del]", "NM_X.1:c.4_6="), Different);
}

#[test]
fn where_cancelling_members_are_written_is_what_they_state() {
    assert_eq!(
        level("NM_X.1:c.[4_5insA;5del]", "NM_X.1:c.[5_6insA;6del]"),
        Different
    );
}

#[test]
fn cancelling_members_beside_a_change_still_state_their_bases() {
    assert_eq!(
        level("NM_X.1:c.[4_5insA;5del;13T>G]", "NM_X.1:c.13T>G"),
        Different
    );
}

#[test]
fn cancelling_members_beside_a_change_are_the_identity_beside_it() {
    assert_eq!(
        level("NM_X.1:c.[4_5insA;5del;13T>G]", "NM_X.1:c.[4_5=;13T>G]"),
        Analogous
    );
}

#[test]
fn whether_members_cancel_is_read_from_the_reference() {
    // On ten Ts the same pair changes g.5.
    assert_eq!(level("NC_T.1:g.[4_5insA;5del]", "NC_T.1:g.5T>A"), Analogous);
}

// --- Members no sequence is built from ---

#[test]
fn a_copy_number_member_compares_by_its_projection() {
    assert_eq!(
        level(
            "NM_X.1:c.[7_8delinsTG;10_12copy3]",
            "NC_X.1:g.[22C>T;23T>G;25_27copy3]"
        ),
        Analogous
    );
}

#[test]
fn unstated_inserted_bases_compare_by_their_genomic_form() {
    // c.15_16 is g.30_31.
    assert_eq!(
        level(
            "NM_X.1:c.[7_8delinsTG;15_16insN[5]]",
            "NC_X.1:g.[22C>T;23T>G;30_31insN[5]]"
        ),
        Analogous
    );
}

#[test]
fn unstated_inserted_bases_in_a_mitochondrial_allele_compare_on_the_genome() {
    assert_eq!(
        level(
            "NC_X.1:m.[22C>T;23T>G;30_31insN[5]]",
            "NC_X.1:g.[22_23delinsTG;30_31insN[5]]"
        ),
        Analogous
    );
}

#[test]
fn a_copy_number_member_is_not_another_count() {
    assert_eq!(
        level("NC_X.1:g.[22C>T;30_35copy3]", "NC_X.1:g.[22C>T;30_35copy4]"),
        Different
    );
}

#[test]
fn the_members_beside_a_copy_number_compare_by_the_sequence_they_make() {
    assert_eq!(
        level(
            "NM_X.1:c.[7C>T;8T>G;10_12copy3]",
            "NM_X.1:c.[7_8delinsTG;10_12copy3]"
        ),
        Analogous
    );
}

#[test]
fn a_member_without_a_sequence_in_another_order_is_one_spelling() {
    assert_eq!(
        level("NC_X.1:g.[22C>T;30_35copy3]", "NC_X.1:g.[30_35copy3;22C>T]"),
        Identity
    );
}

// --- Alleles not read into a sequence ---

#[test]
fn an_allele_with_an_uncertain_member_compares_its_members_as_a_set() {
    // Read as one sequence the two would be the same; with an uncertain
    // member each exact member is compared by its own allele.
    assert_eq!(
        level(
            "NC_X.1:g.[(30_31)del;22C>T;23T>G]",
            "NC_X.1:g.[(30_31)del;22_23delinsTG]"
        ),
        Different
    );
}

#[test]
fn an_allele_with_an_uncertain_member_matches_the_same_members() {
    assert_eq!(
        level(
            "NC_X.1:g.[(30_31)del;22C>T;23T>G]",
            "NC_X.1:g.[23T>G;(30_31)del;22delinsT]"
        ),
        Analogous
    );
}

#[test]
fn protein_alleles_compare_their_members_as_a_set() {
    assert_eq!(
        level(
            "NP_X.1:p.[Lys2Leu;Ala4del]",
            "NP_X.1:p.[Ala4del;Lys2delinsLeu]"
        ),
        Analogous
    );
}

#[test]
fn protein_alleles_in_one_letter_code_are_one_spelling() {
    assert_eq!(
        level("NP_X.1:p.[Lys2Leu;Ala4del]", "NP_X.1:p.[K2L;A4del]"),
        Identity
    );
}

#[test]
fn protein_alleles_with_a_member_changed_differ() {
    assert_eq!(
        level("NP_X.1:p.[Lys2Leu;Ala4del]", "NP_X.1:p.[Lys2Leu;Tyr5del]"),
        Different
    );
}

// --- A nucleotide allele against a protein description ---

#[test]
fn a_coding_allele_is_its_prediction_written_unbracketed() {
    assert_eq!(level("NM_X.1:c.[4A>T;6A>T]", "NP_X.1:p.Lys2Tyr"), Identity);
}

#[test]
fn a_coding_allele_is_not_what_one_member_predicts_alone() {
    // c.6A>T alone is p.(Lys2Asn); with c.4A>T the codon is TAT.
    assert_eq!(level("NM_X.1:c.[4A>T;6A>T]", "NP_X.1:p.Lys2Asn"), Different);
}

#[test]
fn a_coding_allele_is_another_description_of_the_protein_it_leaves() {
    // It predicts p.(Lys2_Leu3delinsAsnArg): MNRAYR, as the wider delins leaves.
    assert_eq!(
        level("NM_X.1:c.[6A>T;8T>G]", "NP_X.1:p.Lys2_Ala4delinsAsnArgAla"),
        Analogous
    );
}

#[test]
fn a_frameshifting_allele_is_the_frameshift_written_without_its_length() {
    // AAA becomes CA, read as CAC: His, which neither member gives alone.
    assert_eq!(
        level("NM_X.1:c.[4A>C;5del]", "NP_X.1:p.Lys2Hisfs"),
        Analogous
    );
}

#[test]
fn a_protein_description_against_a_coding_allele_is_judged_the_same_either_way() {
    assert_eq!(
        level("NP_X.1:p.[Ala4=;Tyr5Asp]", "NM_X.1:c.[13T>G;10_12=]"),
        Identity
    );
}

#[test]
fn a_single_change_is_not_a_protein_allele_stating_more() {
    assert_eq!(
        level("NM_X.1:c.13T>G", "NP_X.1:p.[Ala4=;Tyr5Asp]"),
        Different
    );
}

#[test]
fn a_non_coding_allele_is_judged_by_its_prediction_on_the_coding_transcript() {
    // NR_XR.1's n.30 and n.32 are c.6 and c.4 of NM_X.1, on the other strand.
    assert_eq!(
        level("NR_XR.1:n.[30T>A;32T>A]", "NP_X.1:p.Lys2Tyr"),
        Identity
    );
}

#[test]
fn an_rna_allele_is_judged_by_its_coding_prediction() {
    assert_eq!(level("NM_X.1:r.[4a>u;6a>u]", "NP_X.1:p.Lys2Tyr"), Identity);
}

#[test]
fn a_protein_allele_naming_a_residue_its_protein_lacks_is_an_error() {
    assert!(matches!(
        judge("NM_X.1:c.[4A>T;6A>T]", "NP_X.1:p.[Leu2Ile;Ala4=]"),
        Err(HgvsError::ReferenceMismatch { .. })
    ));
}

#[test]
fn a_genomic_allele_is_judged_by_its_prediction_on_the_transcript() {
    // g.19 and g.21 are c.4 and c.6.
    assert_eq!(
        level("NC_X.1:g.[19A>T;21A>T]", "NP_X.1:p.Lys2Tyr"),
        Identity
    );
}

#[test]
fn an_allele_stating_residues_beside_a_change_is_not_the_change_alone() {
    assert_eq!(
        level("NM_X.1:c.[13T>G;10_12=]", "NP_X.1:p.(Tyr5Asp)"),
        Different
    );
}

#[test]
fn an_allele_stating_residues_beside_a_change_is_its_predicted_protein_allele() {
    assert_eq!(
        level("NM_X.1:c.[13T>G;10_12=]", "NP_X.1:p.[Ala4=;Tyr5Asp]"),
        Identity
    );
}

#[test]
fn an_allele_with_an_uncertain_member_matches_no_protein_description() {
    assert_eq!(
        level("NM_X.1:c.[(4_9)del;13T>G]", "NP_X.1:p.Tyr5Asp"),
        Different
    );
}

// --- Alleles with no one reading ---

#[test]
fn a_contradictory_allele_is_refused() {
    assert!(matches!(
        judge("NM_X.1:c.[7C>T;7C>G]", "NM_X.1:c.7C>T"),
        Err(HgvsError::ValidationError(_))
    ));
}

#[test]
fn an_ambiguous_allele_is_refused_with_its_readings_where_written() {
    match judge("NM_X.1:c.[7_8insA;7_8insG]", "NM_X.1:c.7_8insAG") {
        Err(HgvsError::Ambiguous { readings, .. }) => assert_eq!(
            readings.iter().map(ToString::to_string).collect::<Vec<_>>(),
            ["NM_X.1:c.7_8insAG", "NM_X.1:c.7_8insGA"]
        ),
        other => panic!("{other:?}"),
    }
}

#[test]
fn a_contradictory_allele_is_refused_against_its_exact_projection() {
    assert!(matches!(
        judge("NM_X.1:c.[7C>T;7C>G]", "NC_X.1:g.[22C>T;22C>G]"),
        Err(HgvsError::ValidationError(_))
    ));
}

#[test]
fn an_ambiguous_allele_is_refused_against_its_exact_projection() {
    assert!(matches!(
        judge(
            "NM_X.1:c.[7_8insA;7_8insG]",
            "NC_X.1:g.[22_23insA;22_23insG]"
        ),
        Err(HgvsError::Ambiguous { .. })
    ));
}

#[test]
fn a_contradictory_allele_is_refused_against_itself_in_another_order() {
    assert!(matches!(
        judge("NM_X.1:c.[7C>T;7C>G]", "NM_X.1:c.[7C>G;7C>T]"),
        Err(HgvsError::ValidationError(_))
    ));
}

#[test]
fn an_allele_whose_readings_are_not_listed_is_its_own_spelling() {
    // Each inversion changes bases the other reads (#84).
    assert_eq!(
        level(
            "NC_X.1:g.[20_25inv;22_28inv]",
            "NC_X.1:g.[22_28inv;20_25inv]"
        ),
        Identity
    );
}

#[test]
fn an_allele_whose_readings_are_not_listed_is_not_compared_with_another() {
    assert!(matches!(
        judge("NC_X.1:g.[20_25inv;22_28inv]", "NC_X.1:g.[20_25inv]"),
        Err(HgvsError::UnsupportedOperation(_))
    ));
}

#[test]
fn an_allele_with_no_reading_on_the_genome_is_refused() {
    // Read on the transcript, c.5_6insG inserts at the exon junction, which
    // has no one place on the genome.
    assert!(matches!(
        judge("NM_S.1:c.[5_6insG;6A>T]", "NM_S.1:c.6delinsGT"),
        Err(HgvsError::CisMember { .. })
    ));
}

#[test]
fn exact_members_beside_an_uncertain_one_that_contradict_are_refused() {
    assert!(matches!(
        judge(
            "NC_TEST.1:g.[(1003_1004)del;1007G>A;1007G>C]",
            "NC_TEST.1:g.[(1003_1004)del;1007G>A]"
        ),
        Err(HgvsError::ValidationError(_))
    ));
}

#[test]
fn exact_members_beside_an_uncertain_one_that_contradict_are_refused_against_themselves() {
    assert!(matches!(
        judge(
            "NC_TEST.1:g.[(1003_1004)del;1007G>A;1007G>C]",
            "NC_TEST.1:g.[1007G>C;(1003_1004)del;1007G>A]"
        ),
        Err(HgvsError::ValidationError(_))
    ));
}

#[test]
fn an_exact_member_inside_an_uncertain_one_is_compared_member_by_member() {
    // g.1007 lies inside the bases (1001_1010)del may remove; that overlap is
    // not refused in equivalence.
    assert_eq!(
        level(
            "NC_TEST.1:g.[(1001_1010)del;1007G>A]",
            "NC_TEST.1:g.[1007G>A;(1001_1010)del]"
        ),
        Identity
    );
}

#[test]
fn an_allele_against_a_contradictory_one_is_refused() {
    assert!(matches!(
        judge("NM_X.1:c.[7C>T;13T>G]", "NM_X.1:c.[7C>T;7C>G]"),
        Err(HgvsError::ValidationError(_))
    ));
}
