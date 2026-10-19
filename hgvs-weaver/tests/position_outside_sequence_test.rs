//! A position past the end of the sequence a description is read against is
//! an error from every operation that reads it: a `g.` or `m.` position on
//! the genome, an exonic `c.`, `n.` or `r.` position on the transcript record,
//! a `p.` position on the protein. An insertion needs both its flanking bases.
//! Each test names one operation at one boundary.

mod support;

use hgvs_weaver::data::Strand;
use hgvs_weaver::equivalence::VariantEquivalence;
use hgvs_weaver::error::HgvsError;
use hgvs_weaver::mapper::VariantMapper;
use hgvs_weaver::{parse_hgvs_variant, CVariant, GVariant, SequenceVariant};
use support::{exon, transcript, Provider};

/// A nine-base amplicon.
const AMPLICON: &str = "TTGGGGCTT";

/// A transcript record of 31 bases: a 5'UTR of five, the CDS `ATG AAA CTG GCC
/// TAT CGC TAA` at indices 5 to 25, and a 3'UTR of five, so c.*5 is its last
/// base.
const RECORD: &str = "GGGGGATGAAACTGGCCTATCGCTAACCGTA";

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

/// The genome the record is spliced from: five bases either side of it, so
/// the genome goes on past the record's end.
fn plus_genome() -> String {
    format!("AAAAA{RECORD}AAAAA")
}

/// The plus genome with G where the record's last base, c.*5, is A.
fn differing_genome() -> String {
    let mut genome = plus_genome().into_bytes();
    assert_eq!(genome[35], b'A');
    genome[35] = b'G';
    String::from_utf8(genome).unwrap()
}

/// `NC_000099.1` is the amplicon. `NM_000099.1` (and `NR_000099.1`, the
/// record as a non-coding transcript) is the record on the plus strand of
/// `NC_000098.1`; `NM_000098.1` is the record on the minus strand of
/// `NC_000097.1`, its reverse complement. `NP_000099.1` is the protein
/// `MKLAYR`, served without its stop; `NP_000098.1` is the same protein served
/// with it. `NM_000096.1` is the record again, on `NC_000096.1`, the plus
/// genome with a G under the record's last base.
fn provider() -> Provider {
    Provider::new()
        .sequence("NC_000099.1", AMPLICON)
        .sequence("NC_000098.1", &plus_genome())
        .sequence("NC_000097.1", &revcomp(&plus_genome()))
        .sequence("NM_000099.1", RECORD)
        .sequence("NR_000099.1", RECORD)
        .sequence("NM_000098.1", RECORD)
        .sequence("NC_000096.1", &differing_genome())
        .sequence("NM_000096.1", RECORD)
        .transcript(transcript(
            "NM_000096.1",
            "NC_000096.1",
            Strand::Plus,
            Some((5, 25)),
            vec![exon((0, 31), (5, 35), Strand::Plus)],
        ))
        .sequence("NP_000099.1", "MKLAYR")
        .sequence("NP_000098.1", "MKLAYR*")
        .transcript(transcript(
            "NM_000099.1",
            "NC_000098.1",
            Strand::Plus,
            Some((5, 25)),
            vec![exon((0, 31), (5, 35), Strand::Plus)],
        ))
        .transcript(transcript(
            "NR_000099.1",
            "NC_000098.1",
            Strand::Plus,
            None,
            vec![exon((0, 31), (5, 35), Strand::Plus)],
        ))
        .transcript(transcript(
            "NM_000098.1",
            "NC_000097.1",
            Strand::Minus,
            Some((5, 25)),
            vec![exon((0, 31), (5, 35), Strand::Minus)],
        ))
        .protein_for("NM_000099.1", "NP_000099.1")
        .protein_for("NM_000098.1", "NP_000098.1")
}

fn parse(s: &str) -> SequenceVariant {
    parse_hgvs_variant(s).unwrap_or_else(|e| panic!("{s}: {e}"))
}

fn coding(s: &str) -> CVariant {
    match parse(s) {
        SequenceVariant::Coding(c) => c,
        other => panic!("{other} is not c."),
    }
}

fn genomic(s: &str) -> GVariant {
    match parse(s) {
        SequenceVariant::Genomic(g) => g,
        other => panic!("{other} is not g."),
    }
}

fn normalized(s: &str) -> Result<String, HgvsError> {
    let hdp = provider();
    VariantMapper::new(&hdp)
        .normalize_variant(parse(s))
        .map(|v| v.to_string())
}

/// `result` is the error naming `accession`, its `length`, and the 0-based
/// half-open range `[start, end)` that runs past it.
#[track_caller]
fn assert_outside<T: std::fmt::Debug>(
    result: Result<T, HgvsError>,
    accession: &str,
    (start, end): (usize, usize),
    length: usize,
) {
    match result {
        Err(HgvsError::PositionOutsideSequence {
            accession: ac,
            start: s,
            end: e,
            length: len,
        }) => assert_eq!((ac.as_str(), s, e, len), (accession, start, end, length)),
        other => {
            panic!("expected [{start}, {end}) outside the {length} of {accession}, got {other:?}")
        }
    }
}

// --- g. and m., on the genome ---

#[test]
fn a_deletion_of_the_last_base_normalises() {
    assert_eq!(
        normalized("NC_000099.1:g.9del").unwrap(),
        "NC_000099.1:g.9del"
    );
}

#[test]
fn normalize_refuses_a_deletion_one_past_the_last_base() {
    assert_outside(normalized("NC_000099.1:g.10del"), "NC_000099.1", (9, 10), 9);
}

#[test]
fn to_spdi_refuses_a_range_ending_one_past_the_last_base() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_outside(
        mapper.to_spdi(&parse("NC_000099.1:g.8_10del"), false),
        "NC_000099.1",
        (7, 10),
        9,
    );
}

#[test]
fn an_insertion_between_the_last_two_bases_normalises() {
    assert_eq!(
        normalized("NC_000099.1:g.8_9insA").unwrap(),
        "NC_000099.1:g.8_9insA"
    );
}

#[test]
fn normalize_refuses_an_insertion_after_the_last_base() {
    assert_outside(
        normalized("NC_000099.1:g.9_10insA"),
        "NC_000099.1",
        (8, 10),
        9,
    );
}

#[test]
fn canonical_allele_refuses_an_insertion_after_the_last_base() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_outside(
        mapper.canonical_allele(&parse("NC_000099.1:g.9_10insA")),
        "NC_000099.1",
        (8, 10),
        9,
    );
}

#[test]
fn to_spdi_refuses_a_deletion_wholly_past_the_end() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_outside(
        mapper.to_spdi(&parse("NC_000099.1:g.20del"), false),
        "NC_000099.1",
        (19, 20),
        9,
    );
}

#[test]
fn to_vrs_refuses_a_deletion_wholly_past_the_end() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_outside(
        mapper.to_vrs(&parse("NC_000099.1:g.20del")),
        "NC_000099.1",
        (19, 20),
        9,
    );
}

#[test]
fn a_stated_base_past_the_end_is_refused_for_its_position_not_its_base() {
    assert_outside(
        normalized("NC_000099.1:g.20A>T"),
        "NC_000099.1",
        (19, 20),
        9,
    );
}

#[test]
fn validate_answers_false_for_a_deletion_wholly_past_the_end() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert!(!mapper.validate(&parse("NC_000099.1:g.20del")).unwrap());
}

#[test]
fn equivalence_refuses_a_range_past_the_end_beside_the_part_that_exists() {
    // Read as far as the sequence goes, g.8_12del deletes what g.8_9del does.
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let eq = VariantEquivalence::new(&mapper, &hdp);
    assert_outside(
        eq.equivalent_level(
            &parse("NC_000099.1:g.8_12del"),
            &parse("NC_000099.1:g.8_9del"),
        ),
        "NC_000099.1",
        (7, 12),
        9,
    );
}

#[test]
fn an_uncertain_position_past_the_end_is_not_found_on_the_sequence() {
    // It names no one range; only the bases it states would be checked.
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert!(mapper.validate(&parse("NC_000099.1:g.(20_25)del")).unwrap());
}

#[test]
fn normalize_refuses_a_mitochondrial_position_one_past_the_last_base() {
    assert_outside(normalized("NC_000099.1:m.10del"), "NC_000099.1", (9, 10), 9);
}

#[test]
fn g_to_c_refuses_a_position_past_the_end_of_the_genome() {
    // g.42 would be c.*5+6 had the genome a base there.
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_outside(
        mapper.g_to_c(&genomic("NC_000098.1:g.42del"), "NM_000099.1"),
        "NC_000098.1",
        (41, 42),
        41,
    );
}

#[test]
fn validate_answers_false_for_an_intronic_position_past_the_end_of_the_genome() {
    // c.*5+6 is g.42, one past the genome's last base.
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert!(!mapper.validate(&parse("NM_000099.1:c.*5+6del")).unwrap());
}

#[test]
fn normalize_refuses_an_intronic_position_past_the_end_of_the_genome() {
    assert_outside(
        normalized("NM_000099.1:c.*5+6del"),
        "NC_000098.1",
        (41, 42),
        41,
    );
}

#[test]
fn g_to_c_writes_an_insertion_after_the_records_last_base_as_a_delins_of_it() {
    // On the genome g.36 is G; on the record c.*5 is A, then the record ends.
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_eq!(
        mapper
            .g_to_c(&genomic("NC_000096.1:g.36delinsAT"), "NM_000096.1")
            .unwrap()
            .to_string(),
        "NM_000096.1:c.*5delinsAT"
    );
}

// --- c., n. and r., on the transcript record ---

#[test]
fn a_deletion_of_the_records_last_base_normalises() {
    assert_eq!(
        normalized("NM_000099.1:c.*5del").unwrap(),
        "NM_000099.1:c.*5del"
    );
}

#[test]
fn normalize_refuses_a_coding_position_one_past_the_records_end() {
    assert_outside(
        normalized("NM_000099.1:c.*6del"),
        "NM_000099.1",
        (31, 32),
        31,
    );
}

#[test]
fn normalize_refuses_a_coding_range_ending_one_past_the_records_end() {
    assert_outside(
        normalized("NM_000099.1:c.*4_*6del"),
        "NM_000099.1",
        (29, 32),
        31,
    );
}

#[test]
fn normalize_refuses_an_insertion_after_the_records_last_base() {
    assert_outside(
        normalized("NM_000099.1:c.*5_*6insA"),
        "NM_000099.1",
        (30, 32),
        31,
    );
}

#[test]
fn normalize_refuses_a_coding_deletion_wholly_past_the_records_end() {
    assert_outside(
        normalized("NM_000099.1:c.*20del"),
        "NM_000099.1",
        (45, 46),
        31,
    );
}

#[test]
fn c_to_p_refuses_a_three_prime_utr_position_past_the_records_end() {
    // Inside the record, a 3'UTR deletion predicts p.(=), reading no codon.
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_outside(
        mapper.c_to_p(&coding("NM_000099.1:c.*6del"), None),
        "NM_000099.1",
        (31, 32),
        31,
    );
}

#[test]
fn validate_answers_false_for_a_coding_position_past_the_records_end() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert!(!mapper.validate(&parse("NM_000099.1:c.*6del")).unwrap());
}

#[test]
fn a_deletion_of_the_minus_strand_records_last_base_normalises() {
    assert_eq!(
        normalized("NM_000098.1:c.*5del").unwrap(),
        "NM_000098.1:c.*5del"
    );
}

#[test]
fn normalize_refuses_a_position_one_past_the_minus_strand_records_end() {
    assert_outside(
        normalized("NM_000098.1:c.*6del"),
        "NM_000098.1",
        (31, 32),
        31,
    );
}

#[test]
fn normalize_refuses_an_insertion_after_the_minus_strand_records_last_base() {
    assert_outside(
        normalized("NM_000098.1:c.*5_*6insA"),
        "NM_000098.1",
        (30, 32),
        31,
    );
}

#[test]
fn normalize_refuses_a_non_coding_position_one_past_the_records_end() {
    assert_outside(
        normalized("NR_000099.1:n.32del"),
        "NR_000099.1",
        (31, 32),
        31,
    );
}

#[test]
fn normalize_refuses_an_rna_position_one_past_the_records_end() {
    assert_outside(
        normalized("NM_000099.1:r.*6del"),
        "NM_000099.1",
        (31, 32),
        31,
    );
}

// --- p., on the protein ---

#[test]
fn a_protein_served_without_its_stop_has_the_stop_after_its_last_residue() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert!(mapper.validate(&parse("NP_000099.1:p.Ter7Gln")).unwrap());
}

#[test]
fn normalize_refuses_a_protein_position_one_past_the_stop_it_is_served_without() {
    assert_outside(
        normalized("NP_000099.1:p.Gln8del"),
        "NP_000099.1",
        (7, 8),
        6,
    );
}

#[test]
fn normalize_refuses_a_protein_position_one_past_the_stop_it_is_served_with() {
    assert_outside(
        normalized("NP_000098.1:p.Gln8del"),
        "NP_000098.1",
        (7, 8),
        7,
    );
}
