//! Protein predictions on real RefSeq transcripts: the biocommons `hgvs`
//! extension, whole-gene and start-codon cases (`tests/data/ext.tsv`), with
//! transcript sequences from `tests/data/real_data.json` and CDS bounds from
//! the GenBank records. Each transcript is served as a single exon on a
//! placeholder reference; c. to p. needs only the transcript.

mod support;

use hgvs_weaver::data::Strand;
use hgvs_weaver::mapper::VariantMapper;
use hgvs_weaver::{parse_hgvs_variant, SequenceVariant};
use support::{single_exon_transcript, Provider};

/// `(accession, 0-based CDS start index, 0-based index of the stop codon's last base)`,
/// from the `CDS` feature of each GenBank record (1-based, inclusive).
const CDS: &[(&str, i32, i32)] = &[
    ("NM_000051.3", 385, 9555),  // ATM, CDS 386..9556
    ("NM_000249.3", 198, 2468),  // MLH1, CDS 199..2469
    ("NM_000425.3", 108, 3881),  // L1CAM, CDS 109..3882
    ("NM_000488.3", 119, 1513),  // SERPINC1, CDS 120..1514
    ("NM_002386.3", 1380, 2333), // MC1R, CDS 1381..2334
    ("NM_007199.2", 102, 1892),  // IRAK3, CDS 103..1893
    ("NM_020451.2", 55, 1827),   // SEPN1, CDS 56..1828, Sec at codons 127 and 462
    ("NM_022051.2", 3156, 4436), // EGLN1, CDS 3157..4437
    ("NM_152263.2", 115, 972),   // TPM3, CDS 116..973
];

/// The real_data.json sequences, each transcript one exon on a placeholder reference.
fn provider() -> Provider {
    let mut provider = Provider::from_json_file("../tests/data/real_data.json");
    for &(ac, cds_start, cds_end) in CDS {
        let len = provider.sequences()[ac].len() as i32;
        provider = provider.transcript(single_exon_transcript(
            ac,
            "NC_PLACEHOLDER.1",
            0,
            Strand::Plus,
            cds_start,
            cds_end,
            len,
        ));
    }
    provider
}

/// `(id, c. variant, biocommons hgvs answer, weaver's answer)`. The two agree
/// unless the last column says otherwise; a comment explains each difference.
const CASES: &[(&str, &str, &str, &str)] = &[
    (
        "EXT01",
        "NM_000051.3:c.9170_9171delGA",
        "NP_000042.3:p.(Ter3057Pheext*4)",
        "",
    ),
    (
        "EXT02",
        "NM_000425.3:c.3772dupT",
        "NP_000416.1:p.(Ter1258Leuext*96)",
        "",
    ),
    (
        "EXT03",
        "NM_000249.3:c.2266_2269dupTGTT",
        "NP_000240.1:p.(Ter757Leuext*34)",
        "",
    ),
    (
        "EXT04",
        "NM_000249.3:c.2269dupT",
        "NP_000240.1:p.(Ter757Leuext*33)",
        "",
    ),
    (
        "EXT05",
        "NM_000488.3:c.1391dupA",
        "NP_000479.1:p.(Ter465Valext*18)",
        "",
    ),
    (
        "EXT06",
        "NM_152263.2:c.855delA",
        "NP_689476.2:p.(Ter286Asnext*73)",
        "",
    ),
    // biocommons gives up on a transcript with two in-frame TGA codons; weaver
    // reads them as selenocysteine because the declared CDS end is the stop.
    (
        "MULTISTOP01",
        "NM_020451.2:c.943G>A",
        "NP_065184.2:p.?",
        "NP_065184.2:p.(Gly315Ser)",
    ),
    (
        "WHOLEGENE01",
        "NM_000249.3:c.-7_*46del",
        "NP_000240.1:p.0?",
        "",
    ),
    (
        "WHOLEGENE02",
        "NM_000249.3:c.-20_*20del",
        "NP_000240.1:p.0?",
        "",
    ),
    // weaver keeps the specific start-codon change; `StartCodonConvention::HgvsQuestion`
    // in transform.rs rewrites it to p.Met1? on request.
    (
        "INITMET01",
        "NM_007199.2:c.1A>G",
        "NP_009130.2:p.Met1?",
        "NP_009130.2:p.(Met1Val)",
    ),
    (
        "INITMET02",
        "NM_022051.2:c.-1_1insGCC",
        "NP_071334.1:p.Met1?",
        "",
    ),
    (
        "INITMET03",
        "NM_002386.3:c.-11_19del",
        "NP_002377.4:p.Met1?",
        "",
    ),
];

#[test]
fn biocommons_real_transcript_cases() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let mut failures = Vec::new();
    for (id, c, biocommons, weaver) in CASES {
        // biocommons writes the new stop of an extension as `*`; weaver as `Ter`.
        let normalised = biocommons.replace("ext*", "extTer");
        let expected = if weaver.is_empty() {
            normalised.as_str()
        } else {
            weaver
        };
        let protein_ac = biocommons.split(':').next().unwrap();
        let actual = match parse_hgvs_variant(c) {
            Ok(SequenceVariant::Coding(v)) => match mapper.c_to_p(&v, Some(protein_ac)) {
                Ok(p) => p.to_string(),
                Err(e) => format!("error: {e}"),
            },
            other => format!("not a c. variant: {other:?}"),
        };
        if actual != *expected {
            failures.push(format!(
                "{id}\t{c}\n  expected {expected}\n  actual   {actual}"
            ));
        }
    }
    assert!(failures.is_empty(), "\n{}", failures.join("\n"));
}
