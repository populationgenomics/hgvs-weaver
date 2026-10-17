//! A repeat projected between a transcript and its genome: the count is the
//! total of the unit in the run read 3' from the repeat's first position, and
//! the target reads its own run, which can hold more copies than the source's.

mod support;

use hgvs_weaver::data::Strand;
use hgvs_weaver::equivalence::VariantEquivalence;
use hgvs_weaver::error::HgvsError;
use hgvs_weaver::mapper::VariantMapper;
use hgvs_weaver::{parse_hgvs_variant, CVariant, GVariant, SequenceVariant};
use support::{exon, single_exon_transcript, transcript, Provider};

/// c.1 to c.23, CDS c.1 to c.21: c.4_6 is `AAA`, c.7_15 is `CAG` three times.
const TX: &str = "ATGAAACAGCAGCAGTTCTAAGG";
const FLANK: &str = "CCCCCGCCCC";

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

/// `NM_REV.1`, [`TX`] on the minus strand of `NC_000099.1` at g.11 to g.33:
/// c.N is g.(34 - N), so c.4_6 (`AAA`) is g.28_30 (`TTT`) and c.7_15 is
/// g.19_27 (`CTGCTGCTG`).
fn minus() -> Provider {
    Provider::new()
        .sequence("NC_000099.1", &format!("{FLANK}{}{FLANK}", revcomp(TX)))
        .sequence("NM_REV.1", TX)
        .transcript(single_exon_transcript(
            "NM_REV.1",
            "NC_000099.1",
            10,
            Strand::Minus,
            0,
            20,
            TX.len() as i32,
        ))
}

/// `NM_FWD.1`, [`TX`] on the plus strand of `NC_000098.1` at g.11 to g.33:
/// c.N is g.(N + 10), so c.4_6 is g.14_16 and c.7_15 is g.17_25.
fn plus() -> Provider {
    Provider::new()
        .sequence("NC_000098.1", &format!("{FLANK}{TX}{FLANK}"))
        .sequence("NM_FWD.1", TX)
        .transcript(single_exon_transcript(
            "NM_FWD.1",
            "NC_000098.1",
            10,
            Strand::Plus,
            0,
            20,
            TX.len() as i32,
        ))
}

/// `NM_EDGE.1` on the plus strand of `NC_000097.1`: exon 1 is `ATGCAGTGT`
/// (c.1 to c.9, g.11 to g.19), whose last four bases are `GTGT`; the intron
/// that follows opens with one more `GT`; exon 2 is `GCCTAAGG` (c.10 to
/// c.*2, g.32 to g.39). The genome reads `GTGTGT` over g.16_21.
fn exon_edge() -> Provider {
    let (exon_1, intron, exon_2) = ("ATGCAGTGT", "GTAAGCCCCCAG", "GCCTAAGG");
    Provider::new()
        .sequence(
            "NC_000097.1",
            &format!("{FLANK}{exon_1}{intron}{exon_2}{FLANK}"),
        )
        .sequence("NM_EDGE.1", &format!("{exon_1}{exon_2}"))
        .transcript(transcript(
            "NM_EDGE.1",
            "NC_000097.1",
            Strand::Plus,
            Some((0, 14)),
            vec![
                exon((0, 9), (10, 18), Strand::Plus),
                exon((9, 17), (31, 38), Strand::Plus),
            ],
        ))
}

/// `NM_DIFF.1`, [`TX`] on the plus strand of `NC_000096.1` at g.11 to g.33,
/// except that the genome holds G where the record holds c.6's A: c.4_6 is
/// `AAA` on the record and g.14_16 `AAG` on the genome.
fn record_differs() -> Provider {
    let mut genome = format!("{FLANK}{TX}{FLANK}").into_bytes();
    genome[15] = b'G';
    Provider::new()
        .sequence("NC_000096.1", &String::from_utf8(genome).unwrap())
        .sequence("NM_DIFF.1", TX)
        .transcript(single_exon_transcript(
            "NM_DIFF.1",
            "NC_000096.1",
            10,
            Strand::Plus,
            0,
            20,
            TX.len() as i32,
        ))
}

/// `NM_JX.1` on the plus strand of `NC_000095.1`: exon 1 is `ATGCAA` (c.1
/// to c.6, g.11 to g.16), exon 2 `AAGCCTAAGG` (c.7 to c.16, g.27 to g.36),
/// so the transcript's run of A is c.5_8, across the junction, and the
/// genome's from g.15 is g.15_16.
fn junction() -> Provider {
    let (exon_1, intron, exon_2) = ("ATGCAA", "GTCCCCCCAG", "AAGCCTAAGG");
    Provider::new()
        .sequence(
            "NC_000095.1",
            &format!("{FLANK}{exon_1}{intron}{exon_2}{FLANK}"),
        )
        .sequence("NM_JX.1", &format!("{exon_1}{exon_2}"))
        .transcript(transcript(
            "NM_JX.1",
            "NC_000095.1",
            Strand::Plus,
            Some((0, 14)),
            vec![
                exon((0, 6), (10, 15), Strand::Plus),
                exon((6, 16), (26, 35), Strand::Plus),
            ],
        ))
}

/// `NM_JXM.1`, `NM_JX.1`'s exons on the minus strand of `NC_000093.1`: exon 2
/// (c.7 to c.16) is g.11 to g.20, the intron g.21 to g.30, exon 1 (c.1 to
/// c.6) g.31 to g.36. c.7_8 (`AA`) is g.19_20 (`TT`).
fn junction_minus() -> Provider {
    let (exon_1, intron, exon_2) = ("ATGCAA", "GTCCCCCCAG", "AAGCCTAAGG");
    let genome = format!(
        "{FLANK}{}{}{}{FLANK}",
        revcomp(exon_2),
        revcomp(intron),
        revcomp(exon_1)
    );
    Provider::new()
        .sequence("NC_000093.1", &genome)
        .sequence("NM_JXM.1", &format!("{exon_1}{exon_2}"))
        .transcript(transcript(
            "NM_JXM.1",
            "NC_000093.1",
            Strand::Minus,
            Some((0, 14)),
            vec![
                exon((0, 6), (30, 35), Strand::Minus),
                exon((6, 16), (10, 19), Strand::Minus),
            ],
        ))
}

/// `NM_INT.1` on the minus strand of `NC_000094.1`: exon 1 `ATGCCC` (c.1 to
/// c.6), an intron that reads `GTCAAAAGCAG` in the transcript's sense, and
/// exon 2 `GCCTAAGG` (c.7 to c.14). c.6+N is g.(30 - N), so c.6+4_7-5
/// (`AAAA`) is g.23_26 (`TTTT`).
fn intronic_minus() -> Provider {
    let pre_mrna = "ATGCCCGTCAAAAGCAGGCCTAAGG";
    Provider::new()
        .sequence(
            "NC_000094.1",
            &format!("{FLANK}{}{FLANK}", revcomp(pre_mrna)),
        )
        .sequence("NM_INT.1", "ATGCCCGCCTAAGG")
        .transcript(transcript(
            "NM_INT.1",
            "NC_000094.1",
            Strand::Minus,
            Some((0, 11)),
            vec![
                exon((0, 6), (29, 34), Strand::Minus),
                exon((6, 14), (10, 17), Strand::Minus),
            ],
        ))
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

fn spdi(hdp: &Provider, s: &str) -> String {
    VariantMapper::new(hdp)
        .canonical_allele(&parse(s))
        .unwrap_or_else(|e| panic!("{s}: {e}"))
        .spdi()
}

fn equivalent(hdp: &Provider, a: &str, b: &str) -> bool {
    let mapper = VariantMapper::new(hdp);
    VariantEquivalence::new(&mapper, hdp)
        .equivalent(&parse(a), &parse(b))
        .unwrap_or_else(|e| panic!("{a} vs {b}: {e}"))
}

fn to_genome(hdp: &Provider, s: &str) -> Result<String, HgvsError> {
    VariantMapper::new(hdp)
        .c_to_g(&coding(s), None)
        .map(|g| g.to_string())
}

#[test]
fn on_the_minus_strand_a_repeat_written_mid_run_inserts_on_the_genome() {
    // From c.5 the run is c.5_6 (two A); three copies insert one A into the
    // genome's TTT at g.28_30.
    assert_eq!(
        spdi(&minus(), "NM_REV.1:c.5A[3]"),
        "NC_000099.1:27:TTT:TTTT"
    );
}

#[test]
fn on_the_minus_strand_a_repeat_written_mid_run_at_its_own_count_is_no_change() {
    let allele = VariantMapper::new(&minus())
        .canonical_allele(&parse("NM_REV.1:c.5A[2]"))
        .unwrap();
    assert_eq!(allele.alternate, allele.reference);
}

#[test]
fn on_the_minus_strand_a_multi_base_repeat_written_mid_run_inserts_one_unit() {
    // From c.10 the run is c.10_15 (two CAG); three copies add one, a CTG in
    // the genome's CTGCTGCTG at g.19_27.
    assert_eq!(
        spdi(&minus(), "NM_REV.1:c.10CAG[3]"),
        "NC_000099.1:18:CTGCTGCTG:CTGCTGCTGCTG"
    );
}

#[test]
fn on_the_plus_strand_a_repeat_written_mid_run_inserts_on_the_genome() {
    assert_eq!(spdi(&plus(), "NM_FWD.1:c.5A[3]"), "NC_000098.1:13:AAA:AAAA");
}

#[test]
fn on_the_plus_strand_a_multi_base_repeat_written_mid_run_inserts_one_unit() {
    assert_eq!(
        spdi(&plus(), "NM_FWD.1:c.10CAG[3]"),
        "NC_000098.1:16:CAGCAGCAG:CAGCAGCAGCAG"
    );
}

#[test]
fn on_the_minus_strand_a_repeat_projects_as_the_genomes_whole_run() {
    // The genome's run from g.28 is three T, the transcript's two A and the
    // A at c.4 the repeat leaves alone.
    assert_eq!(
        to_genome(&minus(), "NM_REV.1:c.5A[3]").unwrap(),
        "NC_000099.1:g.28_30T[4]"
    );
}

#[test]
fn on_the_minus_strand_a_multi_base_repeat_projects_as_the_genomes_whole_run() {
    assert_eq!(
        to_genome(&minus(), "NM_REV.1:c.10CAG[3]").unwrap(),
        "NC_000099.1:g.19_27CTG[4]"
    );
}

#[test]
fn on_the_minus_strand_no_copies_from_mid_run_leaves_the_copies_before_it() {
    // c.5A[0] removes c.5_6 and keeps c.4: the genome's TTT becomes T.
    assert_eq!(spdi(&minus(), "NM_REV.1:c.5A[0]"), "NC_000099.1:27:TTT:T");
}

#[test]
fn an_uncertain_count_grows_at_both_ends() {
    assert_eq!(
        to_genome(&minus(), "NM_REV.1:c.5A[(3_5)]").unwrap(),
        "NC_000099.1:g.28_30T[(4_6)]"
    );
}

#[test]
fn on_the_minus_strand_a_genomic_repeat_written_mid_run_projects_as_the_transcripts_whole_run() {
    // From g.29 the genome's run is g.29_30 (two T); the transcript reads its
    // run from c.4, where it is three A.
    let back = VariantMapper::new(&minus())
        .g_to_c(&genomic("NC_000099.1:g.29T[3]"), "NM_REV.1")
        .unwrap();
    assert_eq!(back.to_string(), "NM_REV.1:c.4_6A[4]");
}

#[test]
fn a_repeat_written_mid_run_is_the_transcripts_duplication() {
    assert!(equivalent(&minus(), "NM_REV.1:c.5A[3]", "NM_REV.1:c.6dup"));
}

#[test]
fn a_repeat_written_mid_run_is_the_genomes_duplication() {
    assert!(equivalent(
        &minus(),
        "NM_REV.1:c.5A[3]",
        "NC_000099.1:g.30dup"
    ));
}

#[test]
fn a_repeat_written_mid_run_is_the_repeat_written_from_the_runs_start() {
    assert!(equivalent(&minus(), "NM_REV.1:c.5A[3]", "NM_REV.1:c.4A[4]"));
}

#[test]
fn a_repeat_at_an_exon_end_inserts_one_unit_into_the_run_the_intron_continues() {
    // On the transcript the run from c.6 is GTGT; on the genome the intron
    // carries it on to GTGTGT over g.16_21.
    assert_eq!(
        spdi(&exon_edge(), "NM_EDGE.1:c.6GT[3]"),
        "NC_000097.1:15:GTGTGT:GTGTGTGT"
    );
}

#[test]
fn a_repeat_at_an_exon_end_counts_the_copies_the_intron_continues() {
    assert_eq!(
        to_genome(&exon_edge(), "NM_EDGE.1:c.6GT[3]").unwrap(),
        "NC_000097.1:g.16_21GT[4]"
    );
}

#[test]
fn a_repeat_over_a_base_the_genome_does_not_share_is_written_as_the_genomes_change() {
    // c.4_6 is AAA on the record and g.14_16 AAG on the genome; four A over
    // them is the G becoming AA.
    assert_eq!(
        to_genome(&record_differs(), "NM_DIFF.1:c.4A[4]").unwrap(),
        "NC_000096.1:g.16delinsAA"
    );
}

#[test]
fn an_uncertain_count_over_a_base_the_genome_does_not_share_is_refused() {
    let err = to_genome(&record_differs(), "NM_DIFF.1:c.4A[(3_5)]").unwrap_err();
    assert!(
        matches!(&err, HgvsError::UnsupportedOperation(m) if m.contains("no single edit")),
        "expected UnsupportedOperation, got {err:?}"
    );
}

#[test]
fn a_repeat_without_a_unit_counts_the_unit_its_range_names() {
    // c.7_9 is one CAG of the three at c.7_15; four copies add one.
    assert_eq!(
        spdi(&minus(), "NM_REV.1:c.7_9[4]"),
        "NC_000099.1:18:CTGCTGCTG:CTGCTGCTGCTG"
    );
}

#[test]
fn a_repeat_whose_run_crosses_an_exon_junction_has_no_genomic_form() {
    // From c.7 the run of TG is c.7_10, across the junction after c.9.
    let err = to_genome(&exon_edge(), "NM_EDGE.1:c.7TG[3]").unwrap_err();
    assert!(
        matches!(&err, HgvsError::UnsupportedOperation(m) if m.contains("exon junction")),
        "expected UnsupportedOperation naming the junction, got {err:?}"
    );
}

#[test]
fn a_genomic_repeat_the_transcript_continues_across_a_junction_has_no_coding_form() {
    // The genome's run from g.15 is g.15_16; the transcript's carries on into
    // exon 2 at c.7_8.
    let err = VariantMapper::new(&junction())
        .g_to_c(&genomic("NC_000095.1:g.15A[3]"), "NM_JX.1")
        .unwrap_err();
    assert!(
        matches!(&err, HgvsError::UnsupportedOperation(m) if m.contains("exon junction")),
        "expected UnsupportedOperation naming the junction, got {err:?}"
    );
}

#[test]
fn on_the_minus_strand_a_repeat_whose_run_crosses_an_exon_junction_has_no_genomic_form() {
    let err = to_genome(&junction_minus(), "NM_JXM.1:c.5A[3]").unwrap_err();
    assert!(
        matches!(&err, HgvsError::UnsupportedOperation(m) if m.contains("exon junction")),
        "expected UnsupportedOperation naming the junction, got {err:?}"
    );
}

#[test]
fn a_repeat_from_an_exons_first_base_is_not_refused_when_the_exon_before_ends_in_its_unit() {
    // c.6 is A too, but the run from c.7 is c.7_8; three copies add one T to
    // the genome's TT at g.19_20.
    assert_eq!(
        spdi(&junction_minus(), "NM_JXM.1:c.7A[3]"),
        "NC_000093.1:18:TT:TTT"
    );
}

#[test]
fn on_the_minus_strand_an_intronic_repeat_written_mid_run_inserts_on_the_genome() {
    // From c.6+5 the run is c.6+5_7-5 (three A); four copies add one T to the
    // genome's TTTT at g.23_26.
    assert_eq!(
        spdi(&intronic_minus(), "NM_INT.1:c.6+5A[4]"),
        "NC_000094.1:22:TTTT:TTTTT"
    );
}

#[test]
fn on_the_minus_strand_an_intronic_repeat_projects_as_the_genomes_whole_run() {
    assert_eq!(
        to_genome(&intronic_minus(), "NM_INT.1:c.6+5A[4]").unwrap(),
        "NC_000094.1:g.23_26T[5]"
    );
}

#[test]
fn on_the_minus_strand_a_genomic_repeat_in_an_intron_projects_as_the_transcripts_whole_run() {
    // From g.25 the genome's run is g.25_26 (two T); the transcript's, from
    // c.6+4, is c.6+4_7-5 (four A).
    let back = VariantMapper::new(&intronic_minus())
        .g_to_c(&genomic("NC_000094.1:g.25T[4]"), "NM_INT.1")
        .unwrap();
    assert_eq!(back.to_string(), "NM_INT.1:c.6+4_7-5A[6]");
}

#[test]
fn a_genomic_repeat_in_an_intron_whose_unit_is_absent_is_refused() {
    // g.21 is the intron's G, not T.
    let err = VariantMapper::new(&intronic_minus())
        .g_to_c(&genomic("NC_000094.1:g.21T[2]"), "NM_INT.1")
        .unwrap_err();
    assert!(
        matches!(&err, HgvsError::ReferenceMismatch { found, .. } if found == "G"),
        "expected ReferenceMismatch finding G, got {err:?}"
    );
}
