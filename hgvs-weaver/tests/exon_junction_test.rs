//! Normalisation on a transcript, beside an exon/exon junction whose
//! flanking bases are identical, and at the CDS bounds within one exon.

mod support;

use hgvs_weaver::coords::{GenomicPos, TranscriptPos};
use hgvs_weaver::data::{ExonData, Strand};
use hgvs_weaver::equivalence::{EquivalenceLevel, VariantEquivalence};
use hgvs_weaver::error::HgvsError;
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
/// g.31 to g.42 (c.6 is g.31). `NM_JD.1` is `NM_J.1` with a record G at
/// c.4 over the genome's A. `NR_JDEL.1` lacks the genome's g.19 (cigar
/// `3=1D1=`), so n.3 is g.18 and n.4 g.20; `NR_JINS.1` holds a T the genome
/// lacks after g.19 (cigar `4=1I1=`), so n.4 is g.19, n.5 the T and n.6 g.20.
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
        .sequence("NM_JD.1", &format!("ATGGA{EXON_2}"))
        .transcript(transcript(
            "NM_JD.1",
            "NC_J.1",
            Strand::Plus,
            Some((0, 14)),
            exons(),
        ))
        .sequence("NR_JDEL.1", &format!("ATGA{EXON_2}"))
        .transcript(transcript(
            "NR_JDEL.1",
            "NC_J.1",
            Strand::Plus,
            None,
            vec![
                ExonData::new(
                    TranscriptPos(0),
                    TranscriptPos(4),
                    GenomicPos(15),
                    GenomicPos(19),
                    Strand::Plus,
                    "3=1D1=",
                ),
                exon((4, 16), (30, 41), Strand::Plus),
            ],
        ))
        .sequence("NR_JINS.1", &format!("ATGATA{EXON_2}"))
        .transcript(transcript(
            "NR_JINS.1",
            "NC_J.1",
            Strand::Plus,
            None,
            vec![
                ExonData::new(
                    TranscriptPos(0),
                    TranscriptPos(6),
                    GenomicPos(15),
                    GenomicPos(19),
                    Strand::Plus,
                    "4=1I1=",
                ),
                exon((6, 18), (30, 41), Strand::Plus),
            ],
        ))
}

/// `NM_JM.1`, the same transcript on the minus strand of `NC_JM.1`: exon 2 is
/// g.16 to g.27 (c.6 is g.27), the intron g.28 to g.37, exon 1 g.38 to g.42
/// (c.5 is g.38, c.4 is g.39). `NM_JMD.1` is `NM_JM.1` with a record G at
/// c.4 over the genome's T, an A on the transcript's strand.
fn minus() -> Provider {
    let genome = format!(
        "{FLANK}{}{}{}{FLANK}",
        revcomp(EXON_2),
        revcomp(INTRON),
        revcomp(EXON_1)
    );
    let exons = || {
        vec![
            exon((0, 5), (37, 41), Strand::Minus),
            exon((5, 17), (15, 26), Strand::Minus),
        ]
    };
    Provider::new()
        .sequence("NC_JM.1", &genome)
        .sequence("NM_JM.1", &format!("{EXON_1}{EXON_2}"))
        .transcript(transcript(
            "NM_JM.1",
            "NC_JM.1",
            Strand::Minus,
            Some((0, 14)),
            exons(),
        ))
        .sequence("NM_JMD.1", &format!("ATGGA{EXON_2}"))
        .transcript(transcript(
            "NM_JMD.1",
            "NC_JM.1",
            Strand::Minus,
            Some((0, 14)),
            exons(),
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

/// The variant projected to the genome, as projected.
fn projected(hdp: &Provider, s: &str) -> Result<String, HgvsError> {
    VariantMapper::new(hdp)
        .as_genomic(&parse(s))
        .expect("a nucleotide variant")
        .map(|g| g.to_string())
}

/// The genomic variant `s` projected onto transcript `tx`.
fn on_transcript(hdp: &Provider, s: &str, tx: &str) -> String {
    let mapper = VariantMapper::new(hdp);
    match parse(s) {
        SequenceVariant::Genomic(g) => mapper.g_to_c(&g, tx).map(|v| v.to_string()),
        SequenceVariant::CisPhased(c) => mapper.g_to_c_cis(&c, tx).map(|v| v.to_string()),
        other => panic!("{other} is not genomic"),
    }
    .unwrap_or_else(|e| panic!("projecting {s} onto {tx}: {e}"))
}

/// The `c.` or `n.` variant `s` projected to the genome and back.
fn round_trip(hdp: &Provider, s: &str) -> String {
    let g = projected(hdp, s).unwrap_or_else(|e| panic!("projecting {s}: {e}"));
    let tx = s.split(':').next().expect("an accession");
    on_transcript(hdp, &g, tx)
}

/// Whether `result` is the refusal of an insertion at an exon junction.
fn inserts_at_the_junction<T: std::fmt::Debug>(result: &Result<T, HgvsError>) -> bool {
    matches!(result, Err(HgvsError::UnsupportedOperation(m)) if m.contains("inserts at an exon junction"))
}

/// The message of `result`'s refusal.
fn refusal<T: std::fmt::Debug>(result: Result<T, HgvsError>) -> String {
    match result {
        Err(HgvsError::UnsupportedOperation(m)) => m,
        other => panic!("not a refusal: {other:?}"),
    }
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

#[test]
fn an_identity_across_the_junction_is_the_identity_of_the_genomic_range() {
    // c.4_6 numbers c.4, c.5, the intron and c.6: g.19 to g.31.
    assert_eq!(
        projected(&plus(), "NM_J.1:c.4_6=").unwrap(),
        "NC_J.1:g.19_31="
    );
}

#[test]
fn on_the_minus_strand_an_identity_across_the_junction_is_the_identity_of_the_genomic_range() {
    assert_eq!(
        projected(&minus(), "NM_JM.1:c.4_6=").unwrap(),
        "NC_JM.1:g.27_39="
    );
}

#[test]
fn a_non_coding_identity_across_the_junction_is_the_identity_of_the_genomic_range() {
    assert_eq!(
        projected(&plus(), "NR_J.1:n.4_6=").unwrap(),
        "NC_J.1:g.19_31="
    );
}

#[test]
fn an_identity_stating_its_bases_across_the_junction_is_the_identity_of_the_genomic_range() {
    // AAA are the record's bases at c.4, c.5 and c.6; the genomic range
    // holds the intron between them too, so they are not carried.
    assert_eq!(
        projected(&plus(), "NM_J.1:c.4_6AAA=").unwrap(),
        "NC_J.1:g.19_31="
    );
}

#[test]
fn an_identity_stating_bases_the_record_lacks_across_the_junction_is_a_mismatch() {
    let result = projected(&plus(), "NM_J.1:c.4_6GGG=");
    assert!(
        matches!(result, Err(HgvsError::ReferenceMismatch { .. })),
        "{result:?}"
    );
}

#[test]
fn a_delins_of_the_records_bases_across_the_junction_replaces_the_genomic_range() {
    // c.5 and c.6 are both A, but c.5_6 numbers the intron too: AA replaces
    // all twelve bases, so the intron goes.
    assert_eq!(
        projected(&plus(), "NM_J.1:c.5_6delinsAA").unwrap(),
        "NC_J.1:g.20_31delinsAA"
    );
}

#[test]
fn an_identity_across_the_junction_where_the_record_differs_writes_the_records_base() {
    // NM_JD.1 holds G at c.4 over the genome's A.
    assert_eq!(
        projected(&plus(), "NM_JD.1:c.4_6=").unwrap(),
        "NC_J.1:g.19A>G"
    );
}

#[test]
fn a_duplication_across_the_junction_duplicates_the_genomic_range() {
    assert_eq!(
        projected(&plus(), "NM_J.1:c.4_6dup").unwrap(),
        "NC_J.1:g.19_31dup"
    );
}

#[test]
fn on_the_minus_strand_a_duplication_across_the_junction_duplicates_the_genomic_range() {
    assert_eq!(
        projected(&minus(), "NM_JM.1:c.4_6dup").unwrap(),
        "NC_JM.1:g.27_39dup"
    );
}

#[test]
fn an_inversion_across_the_junction_inverts_the_genomic_range() {
    assert_eq!(
        projected(&plus(), "NM_J.1:c.4_6inv").unwrap(),
        "NC_J.1:g.19_31inv"
    );
}

#[test]
fn on_the_minus_strand_an_inversion_across_the_junction_inverts_the_genomic_range() {
    assert_eq!(
        projected(&minus(), "NM_JM.1:c.5_6inv").unwrap(),
        "NC_JM.1:g.27_38inv"
    );
}

#[test]
fn a_deletion_across_the_junction_deletes_the_genomic_range() {
    assert_eq!(
        projected(&plus(), "NM_J.1:c.4_6del").unwrap(),
        "NC_J.1:g.19_31del"
    );
}

#[test]
fn a_delins_across_the_junction_replaces_the_genomic_range() {
    assert_eq!(
        projected(&plus(), "NM_J.1:c.5_6delinsGC").unwrap(),
        "NC_J.1:g.20_31delinsGC"
    );
}

#[test]
fn on_the_minus_strand_a_delins_across_the_junction_replaces_the_genomic_range() {
    assert_eq!(
        projected(&minus(), "NM_JM.1:c.5_6delinsGT").unwrap(),
        "NC_JM.1:g.27_38delinsAC"
    );
}

#[test]
fn a_delins_of_bases_known_by_number_across_the_junction_replaces_the_genomic_range() {
    assert_eq!(
        projected(&plus(), "NM_J.1:c.4_7delinsN[3]").unwrap(),
        "NC_J.1:g.19_32delinsN[3]"
    );
}

#[test]
fn an_insertion_at_the_junction_has_no_genomic_form() {
    // Between c.5 and c.6 is g.20 or g.31, at either end of the intron.
    let result = projected(&plus(), "NM_J.1:c.5_6insT");
    assert!(inserts_at_the_junction(&result), "{result:?}");
}

#[test]
fn on_the_minus_strand_an_insertion_at_the_junction_has_no_genomic_form() {
    let result = projected(&minus(), "NM_JM.1:c.5_6insT");
    assert!(inserts_at_the_junction(&result), "{result:?}");
}

#[test]
fn a_non_coding_insertion_of_a_length_at_the_junction_has_no_genomic_form() {
    let result = projected(&plus(), "NR_J.1:n.5_6insN[3]");
    assert!(inserts_at_the_junction(&result), "{result:?}");
}

#[test]
fn an_identity_ending_on_the_last_base_of_its_exon_projects() {
    assert_eq!(
        projected(&plus(), "NM_J.1:c.4_5=").unwrap(),
        "NC_J.1:g.19_20="
    );
}

#[test]
fn on_the_minus_strand_an_identity_starting_on_the_first_base_of_its_exon_projects() {
    assert_eq!(
        projected(&minus(), "NM_JM.1:c.6_7=").unwrap(),
        "NC_JM.1:g.26_27="
    );
}

#[test]
fn a_duplication_of_the_last_base_of_its_exon_projects() {
    assert_eq!(
        projected(&plus(), "NM_J.1:c.5dup").unwrap(),
        "NC_J.1:g.20dup"
    );
}

#[test]
fn an_insertion_inside_one_exon_projects() {
    assert_eq!(
        projected(&plus(), "NM_J.1:c.4_5insT").unwrap(),
        "NC_J.1:g.19_20insT"
    );
}

#[test]
fn an_identity_across_the_junction_is_its_genomic_identity() {
    let hdp = plus();
    let mapper = VariantMapper::new(&hdp);
    let level = VariantEquivalence::new(&mapper, &hdp)
        .equivalent_level(&parse("NM_J.1:c.4_6="), &parse("NC_J.1:g.19_31="))
        .unwrap();
    assert_eq!(level, EquivalenceLevel::Identity);
}

#[test]
fn an_identity_across_the_junction_is_not_the_deletion_of_the_intron() {
    let hdp = plus();
    let mapper = VariantMapper::new(&hdp);
    let level = VariantEquivalence::new(&mapper, &hdp)
        .equivalent_level(&parse("NM_J.1:c.4_6="), &parse("NC_J.1:g.21_30del"))
        .unwrap();
    assert_eq!(level, EquivalenceLevel::Different);
}

#[test]
fn the_genomic_deletion_of_the_intron_written_over_the_exons_bases_is_its_projection() {
    let hdp = plus();
    let mapper = VariantMapper::new(&hdp);
    let level = VariantEquivalence::new(&mapper, &hdp)
        .equivalent_level(
            &parse("NC_J.1:g.19_31delinsAAA"),
            &parse("NM_J.1:c.4_6delinsAAA"),
        )
        .unwrap();
    assert_eq!(level, EquivalenceLevel::Identity);
}

#[test]
fn a_genomic_identity_over_the_intron_is_the_identity_of_the_transcript_range() {
    // g.19_31 is c.4, c.5, the intron, and c.6: HGVS numbers c.4_6 over the
    // same bases.
    assert_eq!(
        on_transcript(&plus(), "NC_J.1:g.19_31=", "NM_J.1"),
        "NM_J.1:c.4_6="
    );
}

#[test]
fn a_genomic_identity_stating_the_intron_is_the_identity_of_the_transcript_range() {
    assert_eq!(
        on_transcript(&plus(), "NC_J.1:g.19_31AAGTCCCCCCAGA=", "NM_J.1"),
        "NM_J.1:c.4_6="
    );
}

#[test]
fn a_genomic_delins_of_the_bases_there_over_the_intron_is_the_delins_of_the_transcript_range() {
    assert_eq!(
        on_transcript(&plus(), "NC_J.1:g.19_31delinsAAGTCCCCCCAGA", "NM_J.1"),
        "NM_J.1:c.4_6delinsAAGTCCCCCCAGA"
    );
}

#[test]
fn on_the_minus_strand_a_genomic_identity_over_the_intron_is_the_identity_of_the_transcript_range()
{
    assert_eq!(
        on_transcript(&minus(), "NC_JM.1:g.27_39=", "NM_JM.1"),
        "NM_JM.1:c.4_6="
    );
}

#[test]
fn a_genomic_identity_over_the_intron_as_a_cis_member_is_the_identity_of_the_transcript_range() {
    assert_eq!(
        on_transcript(&plus(), "NC_J.1:g.[19_31=;33G>T]", "NM_J.1"),
        "NM_J.1:c.[4_6=;8G>T]"
    );
}

#[test]
fn the_genomic_deletion_of_the_intron_written_over_the_exons_bases_is_the_delins_of_the_transcript_range(
) {
    assert_eq!(
        on_transcript(&plus(), "NC_J.1:g.19_31delinsAAA", "NM_J.1"),
        "NM_J.1:c.4_6delinsAAA"
    );
}

#[test]
fn on_the_minus_strand_the_genomic_deletion_of_the_intron_written_over_the_exons_bases_is_the_delins_of_the_transcript_range(
) {
    assert_eq!(
        on_transcript(&minus(), "NC_JM.1:g.27_38delinsTT", "NM_JM.1"),
        "NM_JM.1:c.5_6delinsAA"
    );
}

#[test]
fn a_genomic_duplication_over_the_intron_is_the_duplication_of_the_transcript_range() {
    assert_eq!(
        on_transcript(&plus(), "NC_J.1:g.19_31dup", "NM_J.1"),
        "NM_J.1:c.4_6dup"
    );
}

#[test]
fn on_the_minus_strand_a_genomic_inversion_over_the_intron_is_the_inversion_of_the_transcript_range(
) {
    assert_eq!(
        on_transcript(&minus(), "NC_JM.1:g.27_39inv", "NM_JM.1"),
        "NM_JM.1:c.4_6inv"
    );
}

#[test]
fn a_genomic_identity_over_the_intron_writes_the_genomes_base_where_the_record_differs() {
    assert_eq!(
        on_transcript(&plus(), "NC_J.1:g.19_31=", "NM_JD.1"),
        "NM_JD.1:c.4G>A"
    );
}

#[test]
fn on_the_minus_strand_a_genomic_identity_over_the_intron_writes_the_genomes_base_where_the_record_differs(
) {
    assert_eq!(
        on_transcript(&minus(), "NC_JM.1:g.27_39=", "NM_JMD.1"),
        "NM_JMD.1:c.4G>A"
    );
}

#[test]
fn a_genomic_duplication_over_the_intron_where_the_record_differs_is_refused() {
    // The record's G at c.4 is not the genome's A, and the genomic range
    // holds the intron, which the transcript cannot write the change over.
    let hdp = plus();
    let SequenceVariant::Genomic(g) = parse("NC_J.1:g.19_31dup") else {
        unreachable!()
    };
    let result = VariantMapper::new(&hdp).g_to_c(&g, "NM_JD.1");
    assert!(
        matches!(&result, Err(HgvsError::UnsupportedOperation(m)) if m.contains("spans an intron where")),
        "{result:?}"
    );
}

#[test]
fn the_genomic_deletion_written_over_the_intron_alone_is_carried_as_intronic() {
    assert_eq!(
        on_transcript(&plus(), "NC_J.1:g.21_30del", "NM_J.1"),
        "NM_J.1:c.5+1_6-1del"
    );
}

#[test]
fn an_identity_across_the_junction_projected_to_the_genome_and_back_is_itself() {
    assert_eq!(round_trip(&plus(), "NM_J.1:c.4_6="), "NM_J.1:c.4_6=");
}

#[test]
fn a_deletion_across_the_junction_projected_to_the_genome_and_back_is_itself() {
    assert_eq!(round_trip(&plus(), "NM_J.1:c.4_6del"), "NM_J.1:c.4_6del");
}

#[test]
fn a_duplication_across_the_junction_projected_to_the_genome_and_back_is_itself() {
    assert_eq!(round_trip(&plus(), "NM_J.1:c.4_6dup"), "NM_J.1:c.4_6dup");
}

#[test]
fn an_inversion_across_the_junction_projected_to_the_genome_and_back_is_itself() {
    assert_eq!(round_trip(&plus(), "NM_J.1:c.4_6inv"), "NM_J.1:c.4_6inv");
}

#[test]
fn a_delins_across_the_junction_projected_to_the_genome_and_back_is_itself() {
    assert_eq!(
        round_trip(&plus(), "NM_J.1:c.5_6delinsGC"),
        "NM_J.1:c.5_6delinsGC"
    );
}

#[test]
fn on_the_minus_strand_an_identity_across_the_junction_projected_to_the_genome_and_back_is_itself()
{
    assert_eq!(round_trip(&minus(), "NM_JM.1:c.4_6="), "NM_JM.1:c.4_6=");
}

#[test]
fn on_the_minus_strand_a_duplication_across_the_junction_projected_to_the_genome_and_back_is_itself(
) {
    assert_eq!(round_trip(&minus(), "NM_JM.1:c.4_6dup"), "NM_JM.1:c.4_6dup");
}

#[test]
fn on_the_minus_strand_an_inversion_across_the_junction_projected_to_the_genome_and_back_is_itself()
{
    assert_eq!(round_trip(&minus(), "NM_JM.1:c.4_6inv"), "NM_JM.1:c.4_6inv");
}

#[test]
fn a_duplication_across_the_junction_where_the_record_differs_duplicates_the_records_reading() {
    // NM_JD.1 reads g.19_31 as G A, the intron, A, and the duplication makes
    // that twice: the genome's A at g.19 becomes G A, the intron, A, G.
    assert_eq!(
        projected(&plus(), "NM_JD.1:c.4_6dup").unwrap(),
        "NC_J.1:g.19delinsGAGTCCCCCCAGAG"
    );
}

#[test]
fn on_the_minus_strand_a_duplication_across_the_junction_where_the_record_differs_duplicates_the_records_reading(
) {
    assert_eq!(
        projected(&minus(), "NM_JMD.1:c.4_6dup").unwrap(),
        "NC_JM.1:g.39delinsCTCTGGGGGGACTC"
    );
}

#[test]
fn a_statement_across_the_junction_over_a_genome_base_the_record_lacks_deletes_it() {
    assert_eq!(
        projected(&plus(), "NR_JDEL.1:n.3_5=").unwrap(),
        "NC_J.1:g.20del"
    );
}

#[test]
fn a_genomic_statement_over_the_intron_and_a_base_the_record_lacks_inserts_it() {
    assert_eq!(
        on_transcript(&plus(), "NC_J.1:g.18_31=", "NR_JDEL.1"),
        "NR_JDEL.1:c.5_6insA"
    );
}

#[test]
fn a_statement_across_the_junction_over_a_record_base_the_genome_lacks_inserts_it() {
    assert_eq!(
        projected(&plus(), "NR_JINS.1:n.4_7=").unwrap(),
        "NC_J.1:g.19_20insT"
    );
}

#[test]
fn a_range_across_the_junction_ending_on_a_record_base_the_genome_lacks_is_refused() {
    let result = projected(&plus(), "NR_JINS.1:n.5_7=");
    assert!(
        matches!(&result, Err(HgvsError::UnsupportedOperation(m)) if m.contains("ends among bases")),
        "{result:?}"
    );
}

#[test]
fn on_the_minus_strand_a_deletion_across_the_junction_projected_to_the_genome_and_back_is_itself() {
    assert_eq!(round_trip(&minus(), "NM_JM.1:c.4_6del"), "NM_JM.1:c.4_6del");
}

#[test]
fn on_the_minus_strand_a_delins_across_the_junction_projected_to_the_genome_and_back_is_itself() {
    assert_eq!(
        round_trip(&minus(), "NM_JM.1:c.5_6delinsGT"),
        "NM_JM.1:c.5_6delinsGT"
    );
}

#[test]
fn the_refusal_of_an_insertion_at_the_junction_names_its_spelling_at_each_end_of_the_intron() {
    let message = refusal(projected(&plus(), "NM_J.1:c.5_6insT"));
    assert!(message.contains("NM_J.1:c.5_5+1insT"), "{message}");
    assert!(message.contains("NM_J.1:c.6-1_6insT"), "{message}");
}

#[test]
fn on_the_minus_strand_the_refusal_of_an_insertion_at_the_junction_names_its_spelling_at_each_end_of_the_intron(
) {
    let message = refusal(projected(&minus(), "NM_JM.1:c.5_6insN[3]"));
    assert!(message.contains("NM_JM.1:c.5_5+1insN[3]"), "{message}");
    assert!(message.contains("NM_JM.1:c.6-1_6insN[3]"), "{message}");
}

#[test]
fn an_insertion_after_the_earlier_exons_last_base_projects_to_the_start_of_the_intron() {
    assert_eq!(
        projected(&plus(), "NM_J.1:c.5_5+1insT").unwrap(),
        "NC_J.1:g.20_21insT"
    );
}

#[test]
fn an_insertion_before_the_later_exons_first_base_projects_to_the_end_of_the_intron() {
    assert_eq!(
        projected(&plus(), "NM_J.1:c.6-1_6insT").unwrap(),
        "NC_J.1:g.30_31insT"
    );
}

#[test]
fn on_the_minus_strand_an_insertion_after_the_earlier_exons_last_base_projects_to_the_start_of_the_intron(
) {
    // The intron runs down from g.37, after c.5 at g.38.
    assert_eq!(
        projected(&minus(), "NM_JM.1:c.5_5+1insT").unwrap(),
        "NC_JM.1:g.37_38insA"
    );
}

#[test]
fn on_the_minus_strand_an_insertion_before_the_later_exons_first_base_projects_to_the_end_of_the_intron(
) {
    assert_eq!(
        projected(&minus(), "NM_JM.1:c.6-1_6insT").unwrap(),
        "NC_JM.1:g.27_28insA"
    );
}

#[test]
fn an_insertion_across_the_junction_between_bases_that_are_not_neighbours_is_invalid() {
    // c.4 and c.6 have c.5 between them, so nothing is inserted between them.
    let result = projected(&plus(), "NM_J.1:c.4_6insT");
    assert!(
        matches!(&result, Err(HgvsError::ValidationError(m)) if m.contains("not neighbours")),
        "{result:?}"
    );
}
