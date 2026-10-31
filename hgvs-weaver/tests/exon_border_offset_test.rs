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

fn gapped_exon(
    transcript: (i32, i32),
    reference: (i32, i32),
    strand: Strand,
    cigar: &str,
) -> ExonData {
    ExonData::new(
        TranscriptPos(transcript.0),
        TranscriptPos(transcript.1),
        GenomicPos(reference.0),
        GenomicPos(reference.1),
        strand,
        cigar,
    )
}

/// `NM_B.1` on the plus strand of `NC_B.1`; `NM_M.1`, the same record on the
/// minus strand of `NC_M.1`; `NR_B.1`, the record as a non-coding transcript.
/// `NR_GAP.1` lacks the genome's g.36 in exon 2 (cigar `5=1D4=`), so its exon
/// 2 is n.11 to n.19 and exon 3 starts at n.20, where the genome's exon 2
/// alone would put it at n.21.
///
/// Four transcripts have a gap at an edge of exon 2's alignment. `NM_GI.1`
/// ends it with a base the genome lacks (cigar `10=1I`), so exon 2 is c.8 to
/// c.18 and c.18 has no genome base; `NM_GIM.1` is the same on the minus
/// strand of `NC_M.1`. `NM_GL.1` starts it with one (cigar `1I10=`), so c.8
/// has none. `NM_GD.1` lacks the genome's last base of exon 2, g.40 (cigar
/// `9=1D`), so exon 2 is c.8 to c.16.
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
                gapped_exon((0, 10), (10, 19), Strand::Plus, "10M"),
                gapped_exon((10, 19), (30, 39), Strand::Plus, "5=1D4="),
                gapped_exon((19, 29), (50, 59), Strand::Plus, "10M"),
            ],
        ))
        .sequence("NM_GI.1", &format!("{}A{}", &RECORD[..20], &RECORD[20..]))
        .sequence("NM_GIM.1", &format!("{}A{}", &RECORD[..20], &RECORD[20..]))
        .sequence("NM_GL.1", &format!("{}A{}", &RECORD[..10], &RECORD[10..]))
        .sequence("NM_GD.1", &format!("{}{}", &RECORD[..19], &RECORD[20..]))
        .transcript(transcript(
            "NM_GI.1",
            "NC_B.1",
            Strand::Plus,
            Some((3, 26)),
            vec![
                gapped_exon((0, 10), (10, 19), Strand::Plus, "10M"),
                gapped_exon((10, 21), (30, 39), Strand::Plus, "10=1I"),
                gapped_exon((21, 31), (50, 59), Strand::Plus, "10M"),
            ],
        ))
        .transcript(transcript(
            "NM_GIM.1",
            "NC_M.1",
            Strand::Minus,
            Some((3, 26)),
            vec![
                gapped_exon((0, 10), (50, 59), Strand::Minus, "10M"),
                gapped_exon((10, 21), (30, 39), Strand::Minus, "10=1I"),
                gapped_exon((21, 31), (10, 19), Strand::Minus, "10M"),
            ],
        ))
        .transcript(transcript(
            "NM_GL.1",
            "NC_B.1",
            Strand::Plus,
            Some((3, 26)),
            vec![
                gapped_exon((0, 10), (10, 19), Strand::Plus, "10M"),
                gapped_exon((10, 21), (30, 39), Strand::Plus, "1I10="),
                gapped_exon((21, 31), (50, 59), Strand::Plus, "10M"),
            ],
        ))
        .transcript(transcript(
            "NM_GD.1",
            "NC_B.1",
            Strand::Plus,
            Some((3, 26)),
            vec![
                gapped_exon((0, 10), (10, 19), Strand::Plus, "10M"),
                gapped_exon((10, 19), (30, 39), Strand::Plus, "9=1D"),
                gapped_exon((19, 29), (50, 59), Strand::Plus, "10M"),
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
fn an_offset_after_the_last_base_of_the_first_exon_validates() {
    assert!(validate("NM_B.1:c.7+1del"));
}

#[test]
fn an_offset_before_the_first_base_of_a_middle_exon_validates() {
    assert!(validate("NM_B.1:c.8-1del"));
}

#[test]
fn an_offset_after_the_last_base_of_a_middle_exon_validates() {
    assert!(validate("NM_B.1:c.17+1del"));
}

#[test]
fn an_offset_before_the_first_base_of_the_last_exon_validates() {
    assert!(validate("NM_B.1:c.18-1del"));
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
fn normalize_refuses_an_uncertain_position_offset_from_a_base_inside_an_exon() {
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    assert_not_a_border(
        mapper.normalize_variant(parse("NM_B.1:c.(9-1)del")),
        "NM_B.1",
        "c.9-1",
    );
}

#[test]
fn normalize_refuses_an_uncertain_range_with_an_end_offset_from_inside_an_exon() {
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    assert_not_a_border(
        mapper.normalize_variant(parse("NM_B.1:c.(6+1_8-1)del")),
        "NM_B.1",
        "c.6+1",
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
fn r_to_g_names_an_rna_offset_from_inside_an_exon_as_written() {
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    let r = match parse("NM_B.1:r.9-1del") {
        SequenceVariant::Rna(r) => r,
        other => panic!("{other} is not r."),
    };
    assert_not_a_border(mapper.r_to_g(&r, None), "NM_B.1", "r.9-1");
}

#[test]
fn normalize_names_an_rna_offset_from_inside_an_exon_as_written() {
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    assert_not_a_border(
        mapper.normalize_variant(parse("NM_B.1:r.9-1del")),
        "NM_B.1",
        "r.9-1",
    );
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

fn c_to_g(s: &str) -> String {
    let provider = provider();
    VariantMapper::new(&provider)
        .c_to_g(&coding(s), None)
        .unwrap()
        .to_string()
}

#[test]
fn an_offset_after_an_exon_ending_in_a_base_the_genome_lacks_counts_from_the_exons_genome_end() {
    // g.40 is the exon's last genome base, so c.18+1 is g.41.
    assert_eq!(c_to_g("NM_GI.1:c.18+1del"), "NC_B.1:g.41del");
}

#[test]
fn an_offset_after_an_exon_ending_in_a_base_the_genome_lacks_on_the_minus_strand() {
    // Exon 2 runs down NC_M.1 to g.31, so c.18+1 is g.30.
    assert_eq!(c_to_g("NM_GIM.1:c.18+1del"), "NC_M.1:g.30del");
}

#[test]
fn an_offset_before_an_exon_starting_with_a_base_the_genome_lacks_counts_from_the_exons_genome_start(
) {
    // g.31 is the exon's first genome base, so c.8-1 is g.30.
    assert_eq!(c_to_g("NM_GL.1:c.8-1del"), "NC_B.1:g.30del");
}

#[test]
fn an_offset_after_an_exon_the_genome_runs_past_counts_from_the_exons_genome_end() {
    // The exon's alignment ends at g.40, which the transcript lacks: c.16+1
    // is g.41, the first intron base, not g.40.
    assert_eq!(c_to_g("NM_GD.1:c.16+1del"), "NC_B.1:g.41del");
}

#[test]
fn a_range_into_the_intron_after_an_edge_gap_round_trips_through_the_transcript() {
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    for tx in ["NM_GI.1", "NM_GD.1"] {
        let g = match parse("NC_B.1:g.38_41del") {
            SequenceVariant::Genomic(g) => g,
            other => panic!("{other} is not g."),
        };
        let c = mapper.g_to_c(&g, tx).unwrap();
        assert_eq!(
            mapper.c_to_g(&c, None).unwrap().to_string(),
            "NC_B.1:g.38_41del",
            "{tx}: {c}"
        );
    }
}

/// `result` is the refusal of `position` on `accession` as counting past the
/// end of its intron.
fn assert_past_its_intron<T: std::fmt::Debug>(
    result: Result<T, HgvsError>,
    accession: &str,
    position: &str,
) {
    match result {
        Err(HgvsError::OffsetPastIntron {
            accession: ac,
            position: pos,
            ..
        }) => assert_eq!((ac.as_str(), pos.as_str()), (accession, position)),
        other => {
            panic!("expected {accession}:{position} refused as past its intron, got {other:?}")
        }
    }
}

#[test]
fn validate_answers_false_for_an_offset_past_the_end_of_its_intron() {
    // Intron 2 is g.41 to g.50; c.17+11 would be g.51, c.18 itself.
    assert!(!validate("NM_B.1:c.17+11del"));
}

#[test]
fn an_offset_to_the_last_base_of_its_intron_validates() {
    assert!(validate("NM_B.1:c.17+10del"));
}

#[test]
fn c_to_g_refuses_an_offset_past_the_end_of_its_intron() {
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    assert_past_its_intron(
        mapper.c_to_g(&coding("NM_B.1:c.17+11del"), None),
        "NM_B.1",
        "c.17+11",
    );
}

#[test]
fn c_to_g_refuses_a_range_ending_past_its_intron() {
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    assert_past_its_intron(
        mapper.c_to_g(&coding("NM_B.1:c.17+10_17+11del"), None),
        "NM_B.1",
        "c.17+11",
    );
}

#[test]
fn validate_answers_false_for_a_negative_offset_past_the_start_of_its_intron() {
    // c.18-11 would be g.40, c.17.
    assert!(!validate("NM_B.1:c.18-11del"));
}

#[test]
fn validate_answers_false_for_an_offset_past_its_intron_on_the_minus_strand() {
    assert!(!validate("NM_M.1:c.17+11del"));
}

#[test]
fn an_intron_after_an_edge_gap_runs_from_the_end_of_the_exons_genomic_span() {
    // NM_GD.1's exon 2 spans g.31 to g.40, so its intron is g.41 to g.50.
    assert!(!validate("NM_GD.1:c.16+11del"));
}

#[test]
fn an_offset_downstream_of_the_transcript_is_not_bounded_by_an_intron() {
    // c.*3 is g.60 and the genome runs to g.70.
    assert!(validate("NM_B.1:c.*3+10del"));
}

#[test]
fn validate_answers_false_for_an_offset_from_a_base_past_the_transcripts_end() {
    // c.*3 is the transcript's last base; c.*4 is not one of its bases.
    assert!(!validate("NM_B.1:c.*4+1del"));
}

#[test]
fn c_to_g_refuses_an_offset_from_a_base_past_the_transcripts_end_as_outside_it() {
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    match mapper.c_to_g(&coding("NM_B.1:c.*4+1del"), None) {
        Err(HgvsError::PositionOutsideSequence {
            accession, length, ..
        }) => assert_eq!((accession.as_str(), length), ("NM_B.1", 30)),
        other => panic!("expected NM_B.1's 30 bases to be named, got {other:?}"),
    }
}

#[test]
fn validate_refuses_an_offset_from_a_base_before_the_transcripts_first_as_that_base_is() {
    // c.-3 is the transcript's first base; c.-4del is a ValidationError too.
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    assert!(matches!(
        mapper.validate(&parse("NM_B.1:c.-4-1del")),
        Err(HgvsError::ValidationError(_))
    ));
}
