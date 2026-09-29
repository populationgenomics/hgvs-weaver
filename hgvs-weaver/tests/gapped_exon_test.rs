//! A transcript base the genome lacks (an `I` in its exon's cigar) has no
//! genomic substitution equivalent. Between two aligned bases it projects as
//! an insertion between them; at the edge of an exon's alignment (a
//! soft-clipped end supplied as `I`) nothing places it and it is an error.
//! Issue #43.

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
/// - `TX_CLIP.1`: one exon whose first five bases are absent from the genome
///   (cigar `5I20=`), aligned to genome 30..=49.
fn provider() -> Provider {
    let ins = format!("{}T{}", &GENOME[10..25], &GENOME[25..40]);
    let minus = format!("{}T{}", revcomp(&GENOME[25..40]), revcomp(&GENOME[10..25]));
    let clip = format!("GGGCC{}", &GENOME[30..50]);
    Provider::new()
        .sequence(GENOME_AC, GENOME)
        .sequence("TX_INS.1", &ins)
        .sequence("TX_CODING.1", &ins)
        .sequence("TX_MINUS.1", &minus)
        .sequence("TX_CLIP.1", &clip)
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
            "TX_CLIP.1",
            GENOME_AC,
            Strand::Plus,
            None,
            vec![exon((0, 25), (30, 49), Strand::Plus, "5I20=")],
        ))
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
    // An insertion next to the gap goes between the genome's flanking pair.
    assert_eq!(
        to_g(&mapper, "TX_INS.1:n.15_16insA").unwrap(),
        "NC_000099.1:g.25_26insA"
    );
    assert_eq!(
        to_g(&mapper, "TX_INS.1:n.16_17insA").unwrap(),
        "NC_000099.1:g.25_26insA"
    );
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
