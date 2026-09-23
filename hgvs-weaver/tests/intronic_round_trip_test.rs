//! Issue #36: an intronic genomic position must be anchored on the right
//! exon boundary on the minus strand, so that `g_to_c(c_to_g(v)) == v`.

mod support;

use hgvs_weaver::data::Strand;
use hgvs_weaver::mapper::VariantMapper;
use hgvs_weaver::parse_hgvs_variant;
use hgvs_weaver::SequenceVariant;
use support::{exon, transcript, Provider};

/// Three 30-base exons; `c.1` is transcript index 10, the stop ends at 69.
const TRANSCRIPT: &str =
    "GCTAGCTAGCATGGCTGGATCCAAGTTCCTGCACGATGAAGTCATCAACCCTCGATACTGGAGCACATAAACGTTGCAAGGTCCATGACC";

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

/// The transcript spliced into a genome at 1000 with two 70-base introns,
/// served on the plus strand, or mirrored onto the minus strand of the
/// reverse complement. Exon `reference_end` is inclusive.
fn provider(strand: Strand) -> Provider {
    let intron = format!("GT{}ACAG", "ACGT".repeat(16)); // 70 bases
    let plus = format!(
        "{}{}{intron}{}{intron}{}{}",
        "N".repeat(1000),
        &TRANSCRIPT[..30],
        &TRANSCRIPT[30..60],
        &TRANSCRIPT[60..],
        "N".repeat(100)
    );
    let len = plus.len() as i32;
    let plus_exons = [
        (0, 30, 1000, 1029),
        (30, 60, 1100, 1129),
        (60, 90, 1200, 1229),
    ];
    let (genome, exons) = match strand {
        Strand::Plus => (
            plus.clone(),
            plus_exons
                .iter()
                .map(|&(t0, t1, g0, g1)| exon((t0, t1), (g0, g1), strand))
                .collect(),
        ),
        Strand::Minus => (
            revcomp(&plus),
            plus_exons
                .iter()
                .map(|&(t0, t1, g0, g1)| exon((t0, t1), (len - 1 - g1, len - 1 - g0), strand))
                .collect(),
        ),
    };
    Provider::new()
        .sequence("NC_TEST.1", &genome)
        .sequence("NM_TEST.1", TRANSCRIPT)
        .transcript(transcript(
            "NM_TEST.1",
            "NC_TEST.1",
            strand,
            Some((10, 69)),
            exons,
        ))
}

fn coding(s: &str) -> hgvs_weaver::structs::CVariant {
    match parse_hgvs_variant(s).unwrap() {
        SequenceVariant::Coding(c) => c,
        other => panic!("{other} is not c."),
    }
}

#[test]
fn intronic_positions_round_trip_on_both_strands() {
    // Intron 1 reads GT ACGT ACGT ..., so c.20+5 is a G on the transcript
    // strand, whichever strand the genome serves it on.
    let cases = [
        (
            "NM_TEST.1:c.20+5_20+7del",
            "NC_TEST.1:g.1035_1037del",
            "NC_TEST.1:g.294_296del",
        ),
        (
            "NM_TEST.1:c.21-5_21-3del",
            "NC_TEST.1:g.1096_1098del",
            "NC_TEST.1:g.233_235del",
        ),
        (
            "NM_TEST.1:c.20+5G>C",
            "NC_TEST.1:g.1035G>C",
            "NC_TEST.1:g.296C>G",
        ),
        (
            "NM_TEST.1:c.50+1_51-1del",
            "NC_TEST.1:g.1131_1200del",
            "NC_TEST.1:g.131_200del",
        ),
    ];
    for (strand, column) in [(Strand::Plus, 1), (Strand::Minus, 2)] {
        let hdp = provider(strand);
        let mapper = VariantMapper::new(&hdp);
        for case in &cases {
            let c = case.0;
            let expected_g = if column == 1 { case.1 } else { case.2 };
            let g = mapper.c_to_g(&coding(c), Some("NC_TEST.1")).unwrap();
            assert_eq!(g.to_string(), expected_g, "{c} on {strand:?}");
            let back = mapper.g_to_c(&g, "NM_TEST.1").unwrap();
            assert_eq!(back.to_string(), c, "{c} on {strand:?} did not round-trip");
        }
    }
}

#[test]
fn an_intron_base_is_counted_from_its_nearer_exon_on_either_strand() {
    // Intron 1 spans c.20+1 .. c.21-1 over 70 bases: base 35 is +35, base 36
    // is -35, and the genomic positions run the other way on the minus strand.
    for (strand, plus35, minus35) in [
        (Strand::Plus, "NC_TEST.1:g.1065", "NC_TEST.1:g.1066"),
        (Strand::Minus, "NC_TEST.1:g.266", "NC_TEST.1:g.265"),
    ] {
        let hdp = provider(strand);
        let mapper = VariantMapper::new(&hdp);
        let g_to_c = |g: &str| {
            let SequenceVariant::Genomic(v) = parse_hgvs_variant(&format!("{g}del")).unwrap()
            else {
                panic!()
            };
            mapper.g_to_c(&v, "NM_TEST.1").unwrap().to_string()
        };
        assert_eq!(g_to_c(plus35), "NM_TEST.1:c.20+35del", "{strand:?}");
        assert_eq!(g_to_c(minus35), "NM_TEST.1:c.21-35del", "{strand:?}");
    }
}
