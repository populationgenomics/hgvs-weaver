//! An intronic offset counts from an exon border: `+k` from the last base of
//! an exon, `-k` from the first. Counted from any other base it names a
//! position the transcript does not have, which every operation refuses
//! rather than reading it as a base `k` along the genome. Issue #124.

mod support;

use hgvs_weaver::coords::{GenomicPos, TranscriptPos};
use hgvs_weaver::data::{ExonData, Strand};
use hgvs_weaver::error::HgvsError;
use hgvs_weaver::mapper::VariantMapper;
use hgvs_weaver::{parse_hgvs_variant, CVariant, NVariant, SequenceVariant};
use support::{exon, transcript, Provider};

/// The record: a 5'UTR of three (c.-3 to c.-1), eight codons `ATG AAA CTG GCC
/// TAT CGC GAG TAA` (c.1 to c.24) and a 3'UTR of three (c.*1 to c.*3). Three
/// exons of ten: c.-3 to c.7, c.8 to c.17, c.18 to c.*3.
const RECORD: &str = "GGCATGAAACTGGCCTATCGCGAGTAACCG";
const FLANK_5: &str = "AACCAACCAA";
const INTRON_1: &str = "GTAAGTTCAG";
const INTRON_2: &str = "GTGAGTCTAG";
const FLANK_3: &str = "TTGGTTGGTT";

fn revcomp(s: &str) -> String {
    s.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            other => other,
        })
        .collect()
}

/// Exon 1 at genome [10, 20), exon 2 at [30, 40), exon 3 at [50, 60).
fn plus_genome() -> String {
    format!(
        "{FLANK_5}{}{INTRON_1}{}{INTRON_2}{}{FLANK_3}",
        &RECORD[..10],
        &RECORD[10..20],
        &RECORD[20..]
    )
}

fn gapped_exon(transcript: (i32, i32), reference: (i32, i32), cigar: &str) -> ExonData {
    ExonData::new(
        TranscriptPos(transcript.0),
        TranscriptPos(transcript.1),
        GenomicPos(reference.0),
        GenomicPos(reference.1),
        Strand::Plus,
        cigar,
    )
}

/// `NM_B.1` on the plus strand of `NC_B.1`; `NM_M.1`, the same record on the
/// minus strand of `NC_M.1`; `NR_B.1`, the record as a non-coding transcript.
/// `NR_GAP.1` lacks the genome's g.36 in exon 2 (cigar `5=1D4=`), so its exon
/// 2 is n.11 to n.19 and exon 3 starts at n.20, where the genome's exon 2
/// alone would put it at n.21.
fn provider() -> Provider {
    let plus_exons = || {
        vec![
            exon((0, 10), (10, 19), Strand::Plus),
            exon((10, 20), (30, 39), Strand::Plus),
            exon((20, 30), (50, 59), Strand::Plus),
        ]
    };
    let gapped_record = format!("{}{}", &RECORD[..15], &RECORD[16..]);
    Provider::new()
        .sequence("NC_B.1", &plus_genome())
        .sequence("NC_M.1", &revcomp(&plus_genome()))
        .sequence("NM_B.1", RECORD)
        .sequence("NM_M.1", RECORD)
        .sequence("NR_B.1", RECORD)
        .sequence("NR_GAP.1", &gapped_record)
        .sequence("NP_B.1", "MKLAYRE")
        .transcript(transcript(
            "NM_B.1",
            "NC_B.1",
            Strand::Plus,
            Some((3, 26)),
            plus_exons(),
        ))
        .transcript(transcript(
            "NM_M.1",
            "NC_M.1",
            Strand::Minus,
            Some((3, 26)),
            vec![
                exon((0, 10), (50, 59), Strand::Minus),
                exon((10, 20), (30, 39), Strand::Minus),
                exon((20, 30), (10, 19), Strand::Minus),
            ],
        ))
        .transcript(transcript(
            "NR_B.1",
            "NC_B.1",
            Strand::Plus,
            None,
            plus_exons(),
        ))
        .transcript(transcript(
            "NR_GAP.1",
            "NC_B.1",
            Strand::Plus,
            None,
            vec![
                gapped_exon((0, 10), (10, 19), "10M"),
                gapped_exon((10, 19), (30, 39), "5=1D4="),
                gapped_exon((19, 29), (50, 59), "10M"),
            ],
        ))
        .protein_for("NM_B.1", "NP_B.1")
}

fn parse(s: &str) -> SequenceVariant {
    parse_hgvs_variant(s).unwrap()
}

fn coding(s: &str) -> CVariant {
    match parse(s) {
        SequenceVariant::Coding(c) => c,
        other => panic!("{other} is not c."),
    }
}

fn non_coding(s: &str) -> NVariant {
    match parse(s) {
        SequenceVariant::NonCoding(n) => n,
        other => panic!("{other} is not n."),
    }
}

fn validate(s: &str) -> bool {
    let provider = provider();
    VariantMapper::new(&provider).validate(&parse(s)).unwrap()
}

/// `result` is the refusal of `position` on `accession` as no exon border.
fn assert_not_a_border<T: std::fmt::Debug>(
    result: Result<T, HgvsError>,
    accession: &str,
    position: &str,
) {
    match result {
        Err(HgvsError::NotAnExonBorder {
            accession: ac,
            position: pos,
            ..
        }) => assert_eq!((ac.as_str(), pos.as_str()), (accession, position)),
        other => panic!("expected {accession}:{position} refused as no exon border, got {other:?}"),
    }
}

#[test]
fn validate_answers_false_for_an_offset_before_a_base_inside_an_exon() {
    // c.9-1 would be the genome base of c.8, the first base of exon 2.
    assert!(!validate("NM_B.1:c.9-1del"));
}

#[test]
fn validate_answers_false_for_an_offset_after_a_base_inside_an_exon() {
    assert!(!validate("NM_B.1:c.6+1del"));
}

#[test]
fn validate_answers_false_for_an_offset_after_the_first_base_of_an_exon() {
    assert!(!validate("NM_B.1:c.8+1del"));
}

#[test]
fn validate_answers_false_for_an_offset_before_the_last_base_of_an_exon() {
    assert!(!validate("NM_B.1:c.7-1del"));
}

#[test]
fn an_offset_from_each_exon_border_validates() {
    for border in ["c.7+1", "c.8-1", "c.17+1", "c.18-1"] {
        assert!(
            validate(&format!("NM_B.1:{border}del")),
            "{border} counts from an exon border"
        );
    }
}

#[test]
fn c_to_g_refuses_an_offset_from_a_base_inside_an_exon() {
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    assert_not_a_border(
        mapper.c_to_g(&coding("NM_B.1:c.9-1del"), None),
        "NM_B.1",
        "c.9-1",
    );
}

#[test]
fn c_to_g_refuses_an_offset_at_the_end_of_a_range() {
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    assert_not_a_border(
        mapper.c_to_g(&coding("NM_B.1:c.7+1_9-1del"), None),
        "NM_B.1",
        "c.9-1",
    );
}

#[test]
fn normalize_refuses_an_offset_from_a_base_inside_an_exon() {
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    assert_not_a_border(
        mapper.normalize_variant(parse("NM_B.1:c.9-1del")),
        "NM_B.1",
        "c.9-1",
    );
}

#[test]
fn c_to_p_refuses_an_offset_from_a_base_inside_an_exon() {
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    assert_not_a_border(
        mapper.c_to_p(&coding("NM_B.1:c.9-1del"), None),
        "NM_B.1",
        "c.9-1",
    );
}

#[test]
fn to_spdi_refuses_an_offset_from_a_base_inside_an_exon() {
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    assert_not_a_border(
        mapper.to_spdi(&parse("NM_B.1:c.9-1del"), false),
        "NM_B.1",
        "c.9-1",
    );
}

#[test]
fn an_offset_upstream_of_the_transcripts_first_base_validates() {
    // c.-3 is the first base of the transcript, so c.-3-1 is upstream of it.
    assert!(validate("NM_B.1:c.-3-1del"));
}

#[test]
fn validate_answers_false_for_an_offset_before_a_five_prime_utr_base_inside_an_exon() {
    assert!(!validate("NM_B.1:c.-2-1del"));
}

#[test]
fn an_offset_downstream_of_the_transcripts_last_base_validates() {
    // c.*3 is the last base of the transcript, so c.*3+1 is downstream of it.
    assert!(validate("NM_B.1:c.*3+1del"));
}

#[test]
fn validate_answers_false_for_an_offset_after_a_three_prime_utr_base_inside_an_exon() {
    assert!(!validate("NM_B.1:c.*2+1del"));
}

#[test]
fn validate_answers_false_for_an_offset_from_inside_an_exon_on_the_minus_strand() {
    assert!(!validate("NM_M.1:c.9-1del"));
}

#[test]
fn an_offset_from_an_exon_border_on_the_minus_strand_validates() {
    assert!(validate("NM_M.1:c.8-1del"));
}

#[test]
fn n_to_g_refuses_an_offset_from_a_non_coding_base_inside_an_exon() {
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    assert_not_a_border(
        mapper.n_to_g(&non_coding("NR_B.1:n.12-1del"), None),
        "NR_B.1",
        "n.12-1",
    );
}

#[test]
fn validate_answers_false_for_an_rna_offset_from_a_base_inside_an_exon() {
    assert!(!validate("NM_B.1:r.9-1del"));
}

#[test]
fn an_offset_after_a_gapped_exons_last_transcript_base_validates() {
    // Exon 2 of NR_GAP.1 carries nine bases, n.11 to n.19.
    assert!(validate("NR_GAP.1:n.19+1del"));
}

#[test]
fn an_offset_before_the_first_base_after_a_gapped_exon_validates() {
    assert!(validate("NR_GAP.1:n.20-1del"));
}

#[test]
fn an_offset_after_the_first_base_past_a_gapped_exon_is_refused() {
    // n.20 is the first base of exon 3, though exon 2's genome span alone
    // would make it the last of exon 2.
    assert!(!validate("NR_GAP.1:n.20+1del"));
}
