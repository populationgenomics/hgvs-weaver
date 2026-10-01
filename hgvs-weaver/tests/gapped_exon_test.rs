//! A transcript base the genome lacks (an `I` in its exon's cigar) has no
//! genomic position of its own. The genome can describe only the whole run
//! of such bases, so an edit touching a run is written over all of it: the
//! genome gets the run, as changed, between its flanking bases. At the edge
//! of an exon's alignment (a soft-clipped end supplied as `I`) nothing places
//! it and it is an error. Issue #43. The mirror holds for a genome base the
//! transcript lacks (a `D`): a genomic edit touching the run is written on
//! the transcript over the whole run. Issue #45.

mod support;

use hgvs_weaver::coords::{GenomicPos, TranscriptPos};
use hgvs_weaver::data::{ExonData, Strand};
use hgvs_weaver::mapper::VariantMapper;
use hgvs_weaver::parse_hgvs_variant;
use hgvs_weaver::SequenceVariant;
use support::{transcript, Provider};

const GENOME_AC: &str = "NC_000099.1";
const GENOME: &str = "ACGTTGCAAGCTAGCTTACGGATCCATGCAAGTCGATCGGCTAGCTAAGGCCTTAAGCTT";

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

fn exon(transcript: (i32, i32), reference: (i32, i32), strand: Strand, cigar: &str) -> ExonData {
    ExonData {
        transcript_start: TranscriptPos(transcript.0),
        transcript_end: TranscriptPos(transcript.1),
        reference_start: GenomicPos(reference.0),
        reference_end: GenomicPos(reference.1),
        alt_strand: strand,
        cigar: cigar.to_string(),
    }
}

/// Four transcripts on a 60-base genome.
///
/// - `TX_INS.1`: one exon on genome 10..=39 with one base, a T, inserted after
///   genome index 24. Its cigar is `15=1I15=`, so n.16 is the inserted base.
///   `TX_CODING.1` is the same with a CDS from n.11, so c.6 is that base.
/// - `TX_MINUS.1`: the same exon read on the minus strand, with the inserted
///   base again at transcript index 15, between genome 24 and 25.
/// - `TX_GAP3.1`: the same exon with three bases, `TCG`, inserted after
///   genome index 24 (cigar `15=3I15=`), so n.16_18 is the run.
/// - `TX_CLIP.1`: one exon whose first five bases are absent from the genome
///   (cigar `5I20=`), aligned to genome 30..=49.
/// - `TX_DEL.1`: the same exon as `TX_INS.1` but lacking genome index 25, the A
///   at g.26 (cigar `15=1D14=`), with a CDS from its first base so c.15 and c.16
///   flank the gap. `TX_DEL3.1` lacks genome 25..=27, `ATG` at g.26_28 (cigar
///   `15=3D12=`). `TX_DELMINUS.1` reads the exon on the minus strand and lacks
///   genome index 24, the C at g.25. `TX_DELEDGE.1` lacks the exon's first three
///   genome bases (cigar `3D27=`).
fn provider() -> Provider {
    let ins = format!("{}T{}", &GENOME[10..25], &GENOME[25..40]);
    let gap3 = format!("{}TCG{}", &GENOME[10..25], &GENOME[25..40]);
    let del = format!("{}{}", &GENOME[10..25], &GENOME[26..40]);
    let del3 = format!("{}{}", &GENOME[10..25], &GENOME[28..40]);
    let del_minus = format!("{}{}", revcomp(&GENOME[25..40]), revcomp(&GENOME[10..24]));
    let del_edge = GENOME[13..40].to_string();
    let minus = format!("{}T{}", revcomp(&GENOME[25..40]), revcomp(&GENOME[10..25]));
    let clip = format!("GGGCC{}", &GENOME[30..50]);
    Provider::new()
        .sequence(GENOME_AC, GENOME)
        .sequence("TX_INS.1", &ins)
        .sequence("TX_CODING.1", &ins)
        .sequence("TX_MINUS.1", &minus)
        .sequence("TX_GAP3.1", &gap3)
        .sequence("TX_CLIP.1", &clip)
        .sequence("TX_DEL.1", &del)
        .sequence("TX_DEL3.1", &del3)
        .sequence("TX_DELMINUS.1", &del_minus)
        .sequence("TX_DELEDGE.1", &del_edge)
        .transcript(transcript(
            "TX_INS.1",
            GENOME_AC,
            Strand::Plus,
            None,
            vec![exon((0, 31), (10, 39), Strand::Plus, "15=1I15=")],
        ))
        .transcript(transcript(
            "TX_CODING.1",
            GENOME_AC,
            Strand::Plus,
            Some((10, 30)),
            vec![exon((0, 31), (10, 39), Strand::Plus, "15=1I15=")],
        ))
        .transcript(transcript(
            "TX_MINUS.1",
            GENOME_AC,
            Strand::Minus,
            None,
            vec![exon((0, 31), (10, 39), Strand::Minus, "15=1I15=")],
        ))
        .transcript(transcript(
            "TX_GAP3.1",
            GENOME_AC,
            Strand::Plus,
            None,
            vec![exon((0, 33), (10, 39), Strand::Plus, "15=3I15=")],
        ))
        .transcript(transcript(
            "TX_CLIP.1",
            GENOME_AC,
            Strand::Plus,
            None,
            vec![exon((0, 25), (30, 49), Strand::Plus, "5I20=")],
        ))
        .transcript(transcript(
            "TX_DEL.1",
            GENOME_AC,
            Strand::Plus,
            Some((0, 28)),
            vec![exon((0, 29), (10, 39), Strand::Plus, "15=1D14=")],
        ))
        .transcript(transcript(
            "TX_DEL3.1",
            GENOME_AC,
            Strand::Plus,
            Some((0, 26)),
            vec![exon((0, 27), (10, 39), Strand::Plus, "15=3D12=")],
        ))
        .transcript(transcript(
            "TX_DELMINUS.1",
            GENOME_AC,
            Strand::Minus,
            Some((0, 28)),
            vec![exon((0, 29), (10, 39), Strand::Minus, "15=1D14=")],
        ))
        .transcript(transcript(
            "TX_DELEDGE.1",
            GENOME_AC,
            Strand::Plus,
            Some((0, 26)),
            vec![exon((0, 27), (10, 39), Strand::Plus, "3D27=")],
        ))
}

fn to_c(mapper: &VariantMapper, s: &str, tx: &str) -> Result<String, String> {
    match parse_hgvs_variant(s).unwrap() {
        SequenceVariant::Genomic(g) => mapper
            .g_to_c(&g, tx)
            .map(|c| c.to_string())
            .map_err(|e| e.to_string()),
        other => panic!("{other} is not g."),
    }
}

fn to_g(mapper: &VariantMapper, s: &str) -> Result<String, String> {
    let projected = match parse_hgvs_variant(s).unwrap() {
        SequenceVariant::NonCoding(n) => mapper.n_to_g(&n, None),
        SequenceVariant::Coding(c) => mapper.c_to_g(&c, None),
        other => panic!("{other} is neither n. nor c."),
    };
    projected.map(|g| g.to_string()).map_err(|e| e.to_string())
}

#[test]
fn a_base_the_genome_lacks_projects_as_an_insertion_between_its_neighbours() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    // The transcript reads C C T A over n.14..17; the genome reads C C A over
    // g.24..26. n.16 is the T the genome lacks.
    assert_eq!(
        to_g(&mapper, "TX_INS.1:n.14C>G").unwrap(),
        "NC_000099.1:g.24C>G"
    );
    assert_eq!(
        to_g(&mapper, "TX_INS.1:n.15C>G").unwrap(),
        "NC_000099.1:g.25C>G"
    );
    assert_eq!(
        to_g(&mapper, "TX_INS.1:n.16T>G").unwrap(),
        "NC_000099.1:g.25_26insG"
    );
    assert_eq!(
        to_g(&mapper, "TX_INS.1:n.17A>G").unwrap(),
        "NC_000099.1:g.26A>G"
    );
    // The genome already lacks the deleted base: nothing changes on it. It
    // lacks an unchanged one too: no change on the record puts it there.
    assert_eq!(
        to_g(&mapper, "TX_INS.1:n.16del").unwrap(),
        "NC_000099.1:g.25="
    );
    assert_eq!(
        to_g(&mapper, "TX_INS.1:n.16=").unwrap(),
        "NC_000099.1:g.25_26insT"
    );
    assert_eq!(
        to_g(&mapper, "TX_INS.1:n.16dup").unwrap(),
        "NC_000099.1:g.25_26insTT"
    );
    // A range across the gap covers the genome bases either side of it.
    assert_eq!(
        to_g(&mapper, "TX_INS.1:n.15_17del").unwrap(),
        "NC_000099.1:g.25_26del"
    );
    assert_eq!(
        to_g(&mapper, "TX_INS.1:n.15_17delinsGGG").unwrap(),
        "NC_000099.1:g.25_26delinsGGG"
    );
    // An insertion next to the gap carries the gap's base with it: the genome
    // gets what the record holds between the flanking bases.
    assert_eq!(
        to_g(&mapper, "TX_INS.1:n.15_16insA").unwrap(),
        "NC_000099.1:g.25_26insAT"
    );
    assert_eq!(
        to_g(&mapper, "TX_INS.1:n.16_17insA").unwrap(),
        "NC_000099.1:g.25_26insTA"
    );
}

#[test]
fn an_edit_touching_a_run_of_bases_the_genome_lacks_is_written_over_the_whole_run() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    // The record reads C [T C G] A over n.15..19; the genome reads C A over
    // g.25..26. No position inside the run exists on the genome, so a change
    // to one of its bases is the run, as changed, inserted where the run is.
    assert_eq!(
        to_g(&mapper, "TX_GAP3.1:n.16T>A").unwrap(),
        "NC_000099.1:g.25_26insACG"
    );
    assert_eq!(
        to_g(&mapper, "TX_GAP3.1:n.17C>A").unwrap(),
        "NC_000099.1:g.25_26insTAG"
    );
    assert_eq!(
        to_g(&mapper, "TX_GAP3.1:n.17del").unwrap(),
        "NC_000099.1:g.25_26insTG"
    );
    assert_eq!(
        to_g(&mapper, "TX_GAP3.1:n.17dup").unwrap(),
        "NC_000099.1:g.25_26insTCCG"
    );
    assert_eq!(
        to_g(&mapper, "TX_GAP3.1:n.16_18inv").unwrap(),
        "NC_000099.1:g.25_26insCGA"
    );
    assert_eq!(
        to_g(&mapper, "TX_GAP3.1:n.17=").unwrap(),
        "NC_000099.1:g.25_26insTCG"
    );
    // Deleting the whole run leaves the genome as it is.
    assert_eq!(
        to_g(&mapper, "TX_GAP3.1:n.16_18del").unwrap(),
        "NC_000099.1:g.25="
    );
    // A range over a flanking base and part of the run: the genome loses the
    // flanking base and gains what remains of the run.
    assert_eq!(
        to_g(&mapper, "TX_GAP3.1:n.15_17del").unwrap(),
        "NC_000099.1:g.25C>G"
    );
    // An insertion into or beside the run carries the run with it.
    assert_eq!(
        to_g(&mapper, "TX_GAP3.1:n.15_16insA").unwrap(),
        "NC_000099.1:g.25_26insATCG"
    );
    assert_eq!(
        to_g(&mapper, "TX_GAP3.1:n.17_18insA").unwrap(),
        "NC_000099.1:g.25_26insTCAG"
    );
    assert_eq!(
        to_g(&mapper, "TX_GAP3.1:n.18_19insA").unwrap(),
        "NC_000099.1:g.25_26insTCGA"
    );
    // The flanking bases themselves do not touch the run.
    assert_eq!(
        to_g(&mapper, "TX_GAP3.1:n.15C>G").unwrap(),
        "NC_000099.1:g.25C>G"
    );
    assert_eq!(
        to_g(&mapper, "TX_GAP3.1:n.19A>G").unwrap(),
        "NC_000099.1:g.26A>G"
    );
    // Applying each projection to the genome gives the record's exon with
    // the change, which is what makes them one allele.
}

#[test]
fn the_same_through_c_numbering_and_on_the_minus_strand() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_eq!(
        to_g(&mapper, "TX_CODING.1:c.6T>G").unwrap(),
        "NC_000099.1:g.25_26insG"
    );
    assert_eq!(
        to_g(&mapper, "TX_CODING.1:c.5C>G").unwrap(),
        "NC_000099.1:g.25C>G"
    );
    // On the minus strand the inserted T sits between genome 24 and 25 as
    // well; its alternate is complemented.
    assert_eq!(
        to_g(&mapper, "TX_MINUS.1:n.16T>G").unwrap(),
        "NC_000099.1:g.25_26insC"
    );
    assert_eq!(
        to_g(&mapper, "TX_MINUS.1:n.15T>G").unwrap(),
        "NC_000099.1:g.26A>C"
    );
    assert_eq!(
        to_g(&mapper, "TX_MINUS.1:n.17G>A").unwrap(),
        "NC_000099.1:g.25C>T"
    );
}

#[test]
fn a_soft_clipped_base_has_no_genomic_position() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    for s in ["TX_CLIP.1:n.1G>C", "TX_CLIP.1:n.3G>C", "TX_CLIP.1:n.5C>G"] {
        let err = to_g(&mapper, s).unwrap_err();
        assert!(err.contains("NC_000099.1 lacks"), "{s}: {err}");
        assert!(err.contains("edge of an exon"), "{s}: {err}");
    }
    // The first aligned base projects as before.
    assert_eq!(
        to_g(&mapper, "TX_CLIP.1:n.6A>G").unwrap(),
        "NC_000099.1:g.31A>G"
    );
    // A range that reaches into the clipped bases has no genomic form either.
    assert!(to_g(&mapper, "TX_CLIP.1:n.4_7del").is_err());
}

#[test]
fn a_genome_base_the_transcript_lacks_projects_as_an_insertion_between_its_neighbours() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    // The genome reads C A T over g.25..27; the record reads C T over c.15..16.
    // g.26 is the A the record lacks.
    assert_eq!(
        to_c(&mapper, "NC_000099.1:g.26A>G", "TX_DEL.1").unwrap(),
        "TX_DEL.1:c.15_16insG"
    );
    // The record already lacks the deleted base; an unchanged one is put there.
    assert_eq!(
        to_c(&mapper, "NC_000099.1:g.26del", "TX_DEL.1").unwrap(),
        "TX_DEL.1:c.15="
    );
    assert_eq!(
        to_c(&mapper, "NC_000099.1:g.26=", "TX_DEL.1").unwrap(),
        "TX_DEL.1:c.15_16insA"
    );
    assert_eq!(
        to_c(&mapper, "NC_000099.1:g.26dup", "TX_DEL.1").unwrap(),
        "TX_DEL.1:c.15_16insAA"
    );
    // A range across the gap covers the record bases either side of it.
    assert_eq!(
        to_c(&mapper, "NC_000099.1:g.25_27del", "TX_DEL.1").unwrap(),
        "TX_DEL.1:c.15_16del"
    );
    // An insertion beside the gap carries the gap's base with it.
    assert_eq!(
        to_c(&mapper, "NC_000099.1:g.25_26insG", "TX_DEL.1").unwrap(),
        "TX_DEL.1:c.15_16insGA"
    );
    assert_eq!(
        to_c(&mapper, "NC_000099.1:g.26_27insG", "TX_DEL.1").unwrap(),
        "TX_DEL.1:c.15_16insAG"
    );
    // The flanking bases do not touch the run.
    assert_eq!(
        to_c(&mapper, "NC_000099.1:g.25C>G", "TX_DEL.1").unwrap(),
        "TX_DEL.1:c.15C>G"
    );
    assert_eq!(
        to_c(&mapper, "NC_000099.1:g.27T>G", "TX_DEL.1").unwrap(),
        "TX_DEL.1:c.16T>G"
    );
}

#[test]
fn a_genomic_edit_touching_a_run_the_transcript_lacks_is_written_over_the_whole_run() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    // The genome reads C [A T G] C over g.25..29; the record reads C C over c.15..16.
    assert_eq!(
        to_c(&mapper, "NC_000099.1:g.27T>A", "TX_DEL3.1").unwrap(),
        "TX_DEL3.1:c.15_16insAAG"
    );
    assert_eq!(
        to_c(&mapper, "NC_000099.1:g.27del", "TX_DEL3.1").unwrap(),
        "TX_DEL3.1:c.15_16insAG"
    );
    assert_eq!(
        to_c(&mapper, "NC_000099.1:g.26_28del", "TX_DEL3.1").unwrap(),
        "TX_DEL3.1:c.15="
    );
    assert_eq!(
        to_c(&mapper, "NC_000099.1:g.26_27insC", "TX_DEL3.1").unwrap(),
        "TX_DEL3.1:c.15_16insACTG"
    );
    assert_eq!(
        to_c(&mapper, "NC_000099.1:g.25_26insC", "TX_DEL3.1").unwrap(),
        "TX_DEL3.1:c.15_16insCATG"
    );
    assert_eq!(
        to_c(&mapper, "NC_000099.1:g.25C>G", "TX_DEL3.1").unwrap(),
        "TX_DEL3.1:c.15C>G"
    );
    assert_eq!(
        to_c(&mapper, "NC_000099.1:g.29C>G", "TX_DEL3.1").unwrap(),
        "TX_DEL3.1:c.16C>G"
    );
}

#[test]
fn a_genome_base_the_transcript_lacks_on_the_minus_strand() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    // The record lacks g.25, a C; c.15 is g.26 and c.16 is g.24, complemented.
    assert_eq!(
        to_c(&mapper, "NC_000099.1:g.25C>G", "TX_DELMINUS.1").unwrap(),
        "TX_DELMINUS.1:c.15_16insC"
    );
    assert_eq!(
        to_c(&mapper, "NC_000099.1:g.25del", "TX_DELMINUS.1").unwrap(),
        "TX_DELMINUS.1:c.15="
    );
    assert_eq!(
        to_c(&mapper, "NC_000099.1:g.26A>G", "TX_DELMINUS.1").unwrap(),
        "TX_DELMINUS.1:c.15T>C"
    );
    assert_eq!(
        to_c(&mapper, "NC_000099.1:g.24C>A", "TX_DELMINUS.1").unwrap(),
        "TX_DELMINUS.1:c.16G>T"
    );
}

#[test]
fn a_genome_base_the_transcript_lacks_at_the_edge_of_an_exon_has_no_transcript_position() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    for s in ["NC_000099.1:g.11A>C", "NC_000099.1:g.13T>C"] {
        let err = to_c(&mapper, s, "TX_DELEDGE.1").unwrap_err();
        assert!(err.contains("TX_DELEDGE.1 lacks"), "{s}: {err}");
        assert!(err.contains("edge of an exon"), "{s}: {err}");
    }
    assert_eq!(
        to_c(&mapper, "NC_000099.1:g.14G>C", "TX_DELEDGE.1").unwrap(),
        "TX_DELEDGE.1:c.1G>C"
    );
}
