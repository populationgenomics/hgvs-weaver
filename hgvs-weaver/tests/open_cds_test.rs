//! A CDS that runs off an end of the record has no start or stop codon there. A position
//! numbered from that end has no meaning on the record, in either direction, and a protein
//! consequence that depends on it cannot be predicted. Issue #47.

mod support;

use hgvs_weaver::data::{Strand, TranscriptData};
use hgvs_weaver::mapper::VariantMapper;
use hgvs_weaver::parse_hgvs_variant;
use hgvs_weaver::SequenceVariant;
use support::{exon, transcript, Provider};

const GENOME_AC: &str = "NC_000099.1";

/// One exon of 43 bases at genome 0..=42. The record reads GG, then ATG and eight codons, then
/// TAA, then eleven bases of 3'UTR: the CDS is indices 2..=31, ten codons.
///
/// - `NM_WHOLE.1` states the CDS whole.
/// - `NM_OPEN5.1` says the CDS runs off the 5' end: index 2 is its first coding base, not a start.
/// - `NM_OPEN3.1` says the CDS runs off the 3' end: index 31 is its last coding base, not a stop.
fn provider() -> Provider {
    let record = format!(
        "GG{}{}{}{}",
        "ATG", "GCTGAACAACCACTTTCTGAAGCT", "TAA", "GGCCGGCCGGC"
    );
    assert_eq!(record.len(), 43);
    let genome = format!("{record}ACGTACGTACGTACGTACGT");
    let mut provider = Provider::new()
        .sequence(GENOME_AC, &genome)
        .identifier_type(
            GENOME_AC,
            hgvs_weaver::data::IdentifierType::GenomicAccession,
        );
    for (ac, start_open, end_open) in [
        ("NM_WHOLE.1", false, false),
        ("NM_OPEN5.1", true, false),
        ("NM_OPEN3.1", false, true),
    ] {
        let mut model: TranscriptData = transcript(
            ac,
            GENOME_AC,
            Strand::Plus,
            Some((2, 31)),
            vec![exon((0, 43), (0, 42), Strand::Plus)],
        );
        model.cds_start_open = start_open;
        model.cds_end_open = end_open;
        provider = provider
            .sequence(ac, &record)
            .sequence(&ac.replace("NM_", "NP_"), "MAEQPLSEA")
            .transcript(model)
            .protein_for(ac, &ac.replace("NM_", "NP_"));
    }
    provider
}

fn c_to_g(mapper: &VariantMapper, s: &str) -> Result<String, String> {
    match parse_hgvs_variant(s).unwrap() {
        SequenceVariant::Coding(c) => mapper
            .c_to_g(&c, None)
            .map(|g| g.to_string())
            .map_err(|e| e.to_string()),
        SequenceVariant::NonCoding(n) => mapper
            .n_to_g(&n, None)
            .map(|g| g.to_string())
            .map_err(|e| e.to_string()),
        other => panic!("{other} is not c. or n."),
    }
}

fn g_to_c(mapper: &VariantMapper, s: &str, tx: &str) -> Result<String, String> {
    match parse_hgvs_variant(s).unwrap() {
        SequenceVariant::Genomic(g) => mapper
            .g_to_c(&g, tx)
            .map(|c| c.to_string())
            .map_err(|e| e.to_string()),
        other => panic!("{other} is not g."),
    }
}

fn c_to_p(mapper: &VariantMapper, s: &str) -> Result<String, String> {
    match parse_hgvs_variant(s).unwrap() {
        SequenceVariant::Coding(c) => mapper
            .c_to_p(&c, None)
            .map(|p| p.to_string())
            .map_err(|e| e.to_string()),
        other => panic!("{other} is not c."),
    }
}

#[test]
fn a_whole_cds_numbers_from_both_codons() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_eq!(
        c_to_g(&mapper, "NM_WHOLE.1:c.1A>T").unwrap(),
        "NC_000099.1:g.3A>T"
    );
    assert_eq!(
        c_to_g(&mapper, "NM_WHOLE.1:c.-1G>T").unwrap(),
        "NC_000099.1:g.2G>T"
    );
    assert_eq!(
        c_to_g(&mapper, "NM_WHOLE.1:c.*1G>T").unwrap(),
        "NC_000099.1:g.33G>T"
    );
}

#[test]
fn a_5_open_cds_refuses_every_position_numbered_from_the_start_codon() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    for name in [
        "NM_OPEN5.1:c.1A>T",
        "NM_OPEN5.1:c.-1G>T",
        "NM_OPEN5.1:c.10C>T",
        "NM_OPEN5.1:c.4+2A>T",
    ] {
        let err = c_to_g(&mapper, name).unwrap_err();
        assert!(
            err.contains("open at the 5' end") && err.contains("start codon"),
            "{name}: {err}"
        );
    }
}

#[test]
fn a_5_open_cds_still_numbers_from_the_stop_codon_and_as_n() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_eq!(
        c_to_g(&mapper, "NM_OPEN5.1:c.*1G>T").unwrap(),
        "NC_000099.1:g.33G>T"
    );
    assert_eq!(
        c_to_g(&mapper, "NM_OPEN5.1:n.3A>T").unwrap(),
        "NC_000099.1:g.3A>T"
    );
}

#[test]
fn a_3_open_cds_refuses_every_position_numbered_from_the_stop_codon() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let err = c_to_g(&mapper, "NM_OPEN3.1:c.*1G>T").unwrap_err();
    assert!(
        err.contains("open at the 3' end") && err.contains("stop codon"),
        "{err}"
    );
    assert_eq!(
        c_to_g(&mapper, "NM_OPEN3.1:c.1A>T").unwrap(),
        "NC_000099.1:g.3A>T"
    );
    assert_eq!(
        c_to_g(&mapper, "NM_OPEN3.1:c.-1G>T").unwrap(),
        "NC_000099.1:g.2G>T"
    );
}

#[test]
fn a_genome_position_is_refused_when_its_transcript_position_would_be_numbered_from_an_open_end() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    // genome 33 is the first base after the CDS: c.*1 on a whole CDS, undefined on a 3'-open one
    assert_eq!(
        g_to_c(&mapper, "NC_000099.1:g.33G>T", "NM_WHOLE.1").unwrap(),
        "NM_WHOLE.1:c.*1G>T"
    );
    assert!(g_to_c(&mapper, "NC_000099.1:g.33G>T", "NM_OPEN3.1")
        .unwrap_err()
        .contains("open at the 3' end"));
    assert_eq!(
        g_to_c(&mapper, "NC_000099.1:g.33G>T", "NM_OPEN5.1").unwrap(),
        "NM_OPEN5.1:c.*1G>T"
    );
    // genome 3 is c.1: defined on a 3'-open CDS, undefined on a 5'-open one
    assert_eq!(
        g_to_c(&mapper, "NC_000099.1:g.3A>T", "NM_OPEN3.1").unwrap(),
        "NM_OPEN3.1:c.1A>T"
    );
    assert!(g_to_c(&mapper, "NC_000099.1:g.3A>T", "NM_OPEN5.1")
        .unwrap_err()
        .contains("open at the 5' end"));
}

#[test]
fn a_5_open_cds_predicts_no_protein_consequence() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    for name in ["NM_OPEN5.1:c.4G>A", "NM_OPEN5.1:c.*1G>T"] {
        let err = c_to_p(&mapper, name).unwrap_err();
        assert!(
            err.contains("open at the 5' end") && err.contains("no protein consequence"),
            "{name}: {err}"
        );
    }
}

#[test]
fn a_3_open_cds_predicts_an_in_frame_change_inside_the_coding_bases_it_carries() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    // codon 2, GCT (Ala), becomes GTT (Val); the same on the whole CDS
    assert_eq!(
        c_to_p(&mapper, "NM_OPEN3.1:c.5C>T").unwrap(),
        "NP_OPEN3.1:p.(Ala2Val)"
    );
    assert_eq!(
        c_to_p(&mapper, "NM_WHOLE.1:c.5C>T").unwrap(),
        "NP_WHOLE.1:p.(Ala2Val)"
    );
    // codon 3, GAA (Glu), becomes TAA: a stop the record carries, needing no stop codon of its own
    assert_eq!(
        c_to_p(&mapper, "NM_OPEN3.1:c.7G>T").unwrap(),
        "NP_OPEN3.1:p.(Glu3Ter)"
    );
    // an in-frame deletion of codon 4
    assert_eq!(
        c_to_p(&mapper, "NM_OPEN3.1:c.10_12del").unwrap(),
        "NP_OPEN3.1:p.(Gln4del)"
    );
}

#[test]
fn a_3_open_cds_refuses_a_frameshift_or_a_change_reaching_its_last_codon() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    // a frameshift, an insertion, a change in the last codon the record carries, a deletion reaching it
    for name in [
        "NM_OPEN3.1:c.5del",
        "NM_OPEN3.1:c.5_6insA",
        "NM_OPEN3.1:c.28T>A",
        "NM_OPEN3.1:c.25_30del",
    ] {
        let err = c_to_p(&mapper, name).unwrap_err();
        assert!(err.contains("open at the 3' end"), "{name}: {err}");
    }
    // the whole CDS predicts the same changes
    assert_eq!(
        c_to_p(&mapper, "NM_WHOLE.1:c.28T>A").unwrap(),
        "NP_WHOLE.1:p.(Ter10LysextTer?)"
    );
    assert!(c_to_p(&mapper, "NM_WHOLE.1:c.5del").unwrap().contains("fs"));
    // and the 3'-open CDS predicts codon 9, the last whole codon before the one it stops at
    assert_eq!(
        c_to_p(&mapper, "NM_OPEN3.1:c.25G>T").unwrap(),
        "NP_OPEN3.1:p.(Ala9Ser)"
    );
}

#[test]
fn a_provider_that_says_nothing_about_open_ends_means_both_are_whole() {
    let model: TranscriptData = serde_json::from_str(
        r#"{"ac":"NM_X.1","gene":"X","cds_start_index":2,"cds_end_index":31,"strand":1,"reference_accession":"NC_000099.1","exons":[]}"#,
    )
    .unwrap();
    assert!(!model.cds_start_open && !model.cds_end_open);
}
