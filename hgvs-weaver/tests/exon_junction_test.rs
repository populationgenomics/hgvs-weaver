//! Normalisation on a transcript, beside an exon/exon junction whose
//! flanking bases are identical, and at the CDS bounds within one exon.

mod support;

use hgvs_weaver::data::Strand;
use hgvs_weaver::mapper::VariantMapper;
use hgvs_weaver::{parse_hgvs_variant, SequenceVariant};
use support::{exon, single_exon_transcript, transcript, Provider};

/// The transcript, in transcript order: exon 1 is `ATGAA` (c.1 to c.5),
/// exon 2 is `AGGCCTGTAACC` (c.6 to c.*2), so c.4, c.5 and c.6 are all A and
/// the junction falls between c.5 and c.6. The CDS is transcript indices
/// 0..=14 and c.N is index N - 1.
const EXON_1: &str = "ATGAA";
const EXON_2: &str = "AGGCCTGTAACC";
/// The intron in the transcript's sense.
const INTRON: &str = "GTCCCCCCAG";
const FLANK: &str = "CCCCCCCCCCCCCCC";

fn revcomp(s: &str) -> String {
    s.bytes()
        .rev()
        .map(|b| match b {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            other => panic!("not a base: {}", other as char),
        })
        .collect()
}

/// `NM_J.1` and the non-coding `NR_J.1` on the plus strand of `NC_J.1`:
/// exon 1 is g.16 to g.20 (c.5 is g.20), the intron g.21 to g.30, exon 2
/// g.31 to g.42 (c.6 is g.31).
fn plus() -> Provider {
    let genome = format!("{FLANK}{EXON_1}{INTRON}{EXON_2}{FLANK}");
    let exons = || {
        vec![
            exon((0, 5), (15, 19), Strand::Plus),
            exon((5, 17), (30, 41), Strand::Plus),
        ]
    };
    Provider::new()
        .sequence("NC_J.1", &genome)
        .sequence("NM_J.1", &format!("{EXON_1}{EXON_2}"))
        .sequence("NR_J.1", &format!("{EXON_1}{EXON_2}"))
        .transcript(transcript(
            "NM_J.1",
            "NC_J.1",
            Strand::Plus,
            Some((0, 14)),
            exons(),
        ))
        .transcript(transcript("NR_J.1", "NC_J.1", Strand::Plus, None, exons()))
}

/// `NM_JM.1`, the same transcript on the minus strand of `NC_JM.1`: exon 2 is
/// g.16 to g.27 (c.6 is g.27), the intron g.28 to g.37, exon 1 g.38 to g.42
/// (c.5 is g.38, c.4 is g.39).
fn minus() -> Provider {
    let genome = format!(
        "{FLANK}{}{}{}{FLANK}",
        revcomp(EXON_2),
        revcomp(INTRON),
        revcomp(EXON_1)
    );
    Provider::new()
        .sequence("NC_JM.1", &genome)
        .sequence("NM_JM.1", &format!("{EXON_1}{EXON_2}"))
        .transcript(transcript(
            "NM_JM.1",
            "NC_JM.1",
            Strand::Minus,
            Some((0, 14)),
            vec![
                exon((0, 5), (37, 41), Strand::Minus),
                exon((5, 17), (15, 26), Strand::Minus),
            ],
        ))
}

/// `NM_K.1` on `NC_K.1`: exon 1 is `ATGCA` (c.1 to c.5, g.16 to g.20), exon 2
/// `CACACGTAACC` (c.6 to c.*1, g.31 to g.41), so c.4 to c.9 read `CACACA`
/// across the junction.
fn spanning() -> Provider {
    let exon_2 = "CACACGTAACC";
    Provider::new()
        .sequence("NC_K.1", &format!("{FLANK}ATGCA{INTRON}{exon_2}{FLANK}"))
        .sequence("NM_K.1", &format!("ATGCA{exon_2}"))
        .transcript(transcript(
            "NM_K.1",
            "NC_K.1",
            Strand::Plus,
            Some((0, 14)),
            vec![
                exon((0, 5), (15, 19), Strand::Plus),
                exon((5, 16), (30, 40), Strand::Plus),
            ],
        ))
}

/// `NM_C.1`, one exon: `CCAA ATGGCCTAA AACC`, so c.-2, c.-1 and c.1 are A, and
/// so are c.8, c.9, c.*1 and c.*2.
fn one_exon() -> Provider {
    let seq = "CCAAATGGCCTAAAACC";
    Provider::new()
        .sequence("NC_C.1", &format!("{FLANK}{seq}{FLANK}"))
        .sequence("NM_C.1", seq)
        .transcript(single_exon_transcript(
            "NM_C.1",
            "NC_C.1",
            15,
            Strand::Plus,
            4,
            12,
            17,
        ))
}

fn parse(s: &str) -> SequenceVariant {
    parse_hgvs_variant(s).unwrap_or_else(|e| panic!("{s}: {e}"))
}

fn normalized(hdp: &Provider, s: &str) -> String {
    VariantMapper::new(hdp)
        .normalize_variant(parse(s))
        .unwrap_or_else(|e| panic!("normalising {s}: {e}"))
        .to_string()
}

/// The variant projected to the genome and normalised there.
fn genomic(hdp: &Provider, s: &str) -> String {
    let mapper = VariantMapper::new(hdp);
    let g = mapper
        .as_genomic(&parse(s))
        .expect("a nucleotide variant")
        .unwrap_or_else(|e| panic!("projecting {s}: {e}"));
    mapper
        .normalize_variant(SequenceVariant::Genomic(g))
        .unwrap()
        .to_string()
}

#[test]
fn a_repeat_at_its_own_count_across_the_junction_is_not_written_as_an_identity() {
    // From c.4 the run is c.4_6, AA|A: an identity over it would name the
    // intron between c.5 and c.6.
    assert_eq!(normalized(&plus(), "NM_J.1:c.4A[3]"), "NM_J.1:c.4A[3]");
}

#[test]
fn a_deletion_shifts_to_the_last_base_of_its_exon_and_no_further() {
    assert_eq!(normalized(&plus(), "NM_J.1:c.4del"), "NM_J.1:c.5del");
}

#[test]
fn a_deletion_normalised_short_of_the_junction_projects_where_it_was_written() {
    let hdp = plus();
    assert_eq!(genomic(&hdp, "NM_J.1:c.4del"), "NC_J.1:g.20del");
    assert_eq!(
        genomic(&hdp, &normalized(&hdp, "NM_J.1:c.4del")),
        "NC_J.1:g.20del"
    );
}

#[test]
fn a_duplication_shifts_to_the_last_base_of_its_exon_and_no_further() {
    let hdp = plus();
    assert_eq!(normalized(&hdp, "NM_J.1:c.4dup"), "NM_J.1:c.5dup");
    assert_eq!(
        genomic(&hdp, &normalized(&hdp, "NM_J.1:c.4dup")),
        "NC_J.1:g.20dup"
    );
}

#[test]
fn an_insertion_that_slides_to_the_junction_duplicates_the_last_base_of_its_exon() {
    let hdp = plus();
    assert_eq!(normalized(&hdp, "NM_J.1:c.3_4insA"), "NM_J.1:c.5dup");
    assert_eq!(
        genomic(&hdp, &normalized(&hdp, "NM_J.1:c.3_4insA")),
        "NC_J.1:g.20dup"
    );
}

#[test]
fn an_insertion_does_not_slide_onto_the_junction() {
    // Sliding one base would write c.5_6insCA, between two exons: on the
    // genome that is either the end of exon 1 or the start of exon 2.
    let hdp = plus();
    assert_eq!(normalized(&hdp, "NM_J.1:c.4_5insAC"), "NM_J.1:c.4_5insAC");
    assert_eq!(
        genomic(&hdp, &normalized(&hdp, "NM_J.1:c.4_5insAC")),
        "NC_J.1:g.20_21insCA"
    );
}

#[test]
fn an_insertion_is_not_rewritten_as_a_duplication_across_the_junction() {
    // c.4_6 is AAA, so the insertion repeats the bases before it, but they
    // lie in two exons. On the genome the intron's last base, G, precedes
    // c.6, so nothing there repeats the insert.
    let hdp = plus();
    assert_eq!(normalized(&hdp, "NM_J.1:c.6_7insAAA"), "NM_J.1:c.6_7insAAA");
    assert_eq!(
        genomic(&hdp, &normalized(&hdp, "NM_J.1:c.6_7insAAA")),
        "NC_J.1:g.31_32insAAA"
    );
}

#[test]
fn a_deletion_spanning_the_junction_shifts_while_each_end_stays_in_its_exon() {
    // c.3_6del and c.4_7del each remove the whole intron; the bases either
    // side of it they take differ only in c.3 and c.7, which are both G, so on
    // the genome they are one change.
    let hdp = plus();
    assert_eq!(normalized(&hdp, "NM_J.1:c.3_6del"), "NM_J.1:c.4_7del");
    assert_eq!(genomic(&hdp, "NM_J.1:c.3_6del"), "NC_J.1:g.19_32del");
    assert_eq!(
        genomic(&hdp, &normalized(&hdp, "NM_J.1:c.3_6del")),
        "NC_J.1:g.19_32del"
    );
}

#[test]
fn a_deletion_spanning_the_junction_shifts_no_further_than_its_start_can_stay_in_its_exon() {
    // c.4_7del (CACA) matches as far as c.7_10del, which would leave exon 1;
    // c.5_8del is as far as its start can go.
    let hdp = spanning();
    assert_eq!(normalized(&hdp, "NM_K.1:c.4_7del"), "NM_K.1:c.5_8del");
    assert_eq!(genomic(&hdp, "NM_K.1:c.4_7del"), "NC_K.1:g.21_34del");
    assert_eq!(
        genomic(&hdp, &normalized(&hdp, "NM_K.1:c.4_7del")),
        "NC_K.1:g.21_34del"
    );
}

#[test]
fn a_deletion_beside_the_junction_does_not_shift_to_span_it() {
    // c.4_5 is CA and c.6 is C, so the sequence allows c.5_6del, which would
    // take the intron with it on the genome.
    assert_eq!(
        normalized(&spanning(), "NM_K.1:c.4_5del"),
        "NM_K.1:c.4_5del"
    );
}

#[test]
fn a_duplication_spanning_the_junction_shifts_while_each_end_stays_in_its_exon() {
    assert_eq!(normalized(&plus(), "NM_J.1:c.3_6dup"), "NM_J.1:c.4_7dup");
}

#[test]
fn an_insertion_written_on_the_junction_is_left_as_written() {
    // Between c.5 and c.6 is the end of exon 1 or the start of exon 2 on the
    // genome; neither a slide nor a duplication of c.4_5 may pick one.
    assert_eq!(
        normalized(&plus(), "NM_J.1:c.5_6insAA"),
        "NM_J.1:c.5_6insAA"
    );
}

#[test]
fn a_shift_crosses_the_cds_start_within_one_exon() {
    assert_eq!(normalized(&one_exon(), "NM_C.1:c.-2del"), "NM_C.1:c.1del");
}

#[test]
fn a_shift_crosses_the_cds_end_within_one_exon() {
    assert_eq!(normalized(&one_exon(), "NM_C.1:c.8del"), "NM_C.1:c.*2del");
}

#[test]
fn an_insertion_does_not_slide_past_the_last_base_of_the_transcript() {
    // c.*2 is the last base: sliding one base would write c.*2_*3insAC,
    // after a base that does not exist.
    assert_eq!(
        normalized(&plus(), "NM_J.1:c.*1_*2insCA"),
        "NM_J.1:c.*1_*2insCA"
    );
}

#[test]
fn a_non_coding_deletion_stops_at_the_junction() {
    assert_eq!(normalized(&plus(), "NR_J.1:n.4del"), "NR_J.1:n.5del");
}

#[test]
fn a_non_coding_insertion_that_slides_to_the_junction_duplicates_the_last_base_of_its_exon() {
    assert_eq!(normalized(&plus(), "NR_J.1:n.3_4insA"), "NR_J.1:n.5dup");
}

#[test]
fn an_rna_insertion_that_slides_to_the_junction_duplicates_the_last_base_of_its_exon() {
    assert_eq!(normalized(&plus(), "NM_J.1:r.3_4insa"), "NM_J.1:r.5dup");
}

#[test]
fn an_rna_deletion_stops_at_the_junction() {
    assert_eq!(normalized(&plus(), "NM_J.1:r.4del"), "NM_J.1:r.5del");
}

#[test]
fn on_the_minus_strand_a_deletion_stops_at_the_junction() {
    assert_eq!(normalized(&minus(), "NM_JM.1:c.4del"), "NM_JM.1:c.5del");
}

#[test]
fn on_the_minus_strand_a_deletion_normalised_short_of_the_junction_projects_where_it_was_written() {
    let hdp = minus();
    assert_eq!(genomic(&hdp, "NM_JM.1:c.4del"), "NC_JM.1:g.39del");
    assert_eq!(
        genomic(&hdp, &normalized(&hdp, "NM_JM.1:c.4del")),
        "NC_JM.1:g.39del"
    );
}

#[test]
fn on_the_minus_strand_a_duplication_stops_at_the_junction() {
    let hdp = minus();
    assert_eq!(normalized(&hdp, "NM_JM.1:c.4dup"), "NM_JM.1:c.5dup");
    assert_eq!(
        genomic(&hdp, &normalized(&hdp, "NM_JM.1:c.4dup")),
        "NC_JM.1:g.39dup"
    );
}

#[test]
fn on_the_minus_strand_an_insertion_that_slides_to_the_junction_duplicates_the_last_base_of_its_exon(
) {
    let hdp = minus();
    assert_eq!(normalized(&hdp, "NM_JM.1:c.3_4insA"), "NM_JM.1:c.5dup");
    assert_eq!(
        genomic(&hdp, &normalized(&hdp, "NM_JM.1:c.3_4insA")),
        "NC_JM.1:g.39dup"
    );
}

#[test]
fn on_the_minus_strand_an_insertion_does_not_slide_onto_the_junction() {
    let hdp = minus();
    assert_eq!(normalized(&hdp, "NM_JM.1:c.4_5insAC"), "NM_JM.1:c.4_5insAC");
    assert_eq!(
        genomic(&hdp, &normalized(&hdp, "NM_JM.1:c.4_5insAC")),
        "NC_JM.1:g.38_39insGT"
    );
}

#[test]
fn a_delins_that_is_an_insertion_shifts_to_the_last_base_of_its_exon_and_no_further() {
    // G becomes GA: an A inserted after c.3 slides past c.4 and c.5 and
    // duplicates c.5, as c.3_4insA does.
    let hdp = plus();
    assert_eq!(normalized(&hdp, "NM_J.1:c.3delinsGA"), "NM_J.1:c.5dup");
    assert_eq!(genomic(&hdp, "NM_J.1:c.3delinsGA"), "NC_J.1:g.20dup");
}

#[test]
fn a_delins_reduced_to_an_insertion_at_the_end_of_its_exon_keeps_the_exons_last_base() {
    // AA becomes AAC: a C after c.5, which written as c.5_6insC would sit on
    // the junction, at either end of the intron on the genome.
    let hdp = plus();
    assert_eq!(
        normalized(&hdp, "NM_J.1:c.4_5delinsAAC"),
        "NM_J.1:c.5delinsAC"
    );
    assert_eq!(genomic(&hdp, "NM_J.1:c.4_5delinsAAC"), "NC_J.1:g.20_21insC");
    assert_eq!(
        genomic(&hdp, &normalized(&hdp, "NM_J.1:c.4_5delinsAAC")),
        "NC_J.1:g.20_21insC"
    );
}

#[test]
fn a_delins_reduced_to_an_insertion_at_the_start_of_its_exon_keeps_the_exons_first_base() {
    let hdp = plus();
    assert_eq!(normalized(&hdp, "NM_J.1:c.6delinsCA"), "NM_J.1:c.6delinsCA");
    assert_eq!(genomic(&hdp, "NM_J.1:c.6delinsCA"), "NC_J.1:g.30_31insC");
}

#[test]
fn a_delins_spanning_the_junction_is_trimmed_while_each_end_stays_in_its_exon() {
    // AAAG becomes AG: trimming c.4 and c.7 leaves c.5_6del, which still
    // removes the intron; trimming further would leave one of the exons.
    let hdp = plus();
    assert_eq!(normalized(&hdp, "NM_J.1:c.4_7delinsAG"), "NM_J.1:c.5_6del");
    assert_eq!(genomic(&hdp, "NM_J.1:c.4_7delinsAG"), "NC_J.1:g.20_31del");
    assert_eq!(
        genomic(&hdp, &normalized(&hdp, "NM_J.1:c.4_7delinsAG")),
        "NC_J.1:g.20_31del"
    );
}

#[test]
fn a_delins_of_the_two_bases_either_side_of_the_junction_stays_where_it_is() {
    // c.5 and c.6 are both A; c.5_6delinsAGA removes the intron between them,
    // which c.5_6insG would not.
    assert_eq!(
        normalized(&plus(), "NM_J.1:c.5_6delinsAGA"),
        "NM_J.1:c.5_6delinsAGA"
    );
}

#[test]
fn on_the_minus_strand_a_delins_reduced_to_an_insertion_at_the_end_of_its_exon_keeps_the_exons_last_base(
) {
    let hdp = minus();
    assert_eq!(
        normalized(&hdp, "NM_JM.1:c.4_5delinsAAC"),
        "NM_JM.1:c.5delinsAC"
    );
    assert_eq!(
        genomic(&hdp, &normalized(&hdp, "NM_JM.1:c.4_5delinsAAC")),
        "NC_JM.1:g.37_38insG"
    );
}

#[test]
fn on_the_minus_strand_a_delins_that_is_an_insertion_stops_at_the_junction() {
    let hdp = minus();
    assert_eq!(normalized(&hdp, "NM_JM.1:c.3delinsGA"), "NM_JM.1:c.5dup");
    assert_eq!(genomic(&hdp, "NM_JM.1:c.3delinsGA"), "NC_JM.1:g.39dup");
}

#[test]
fn a_delins_with_an_uncertain_position_is_left_as_written() {
    assert_eq!(
        normalized(&plus(), "NM_J.1:c.(3_4)delinsGA"),
        "NM_J.1:c.(3_4)delinsGA"
    );
}
