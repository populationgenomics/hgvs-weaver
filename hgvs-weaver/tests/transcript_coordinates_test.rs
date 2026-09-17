//! Regression tests for c./n. -> genomic coordinate resolution on the paths
//! that used to route through `DataProvider::c_to_g`.
//!
//! Those paths (SPDI intervals and c.-vs-c. equivalence) handed the provider a
//! bare index with no anchor, so `c.*1` was indistinguishable from `c.1` and no
//! adapter could get the minus strand or a non-zero CDS start right. The
//! fixtures here deliberately use a CDS that does not start at index 0 and a
//! minus-strand transcript.

use hgvs_weaver::coords::{GenomicPos, TranscriptPos};
use hgvs_weaver::data::{
    DataProvider, ExonData, IdentifierKind, IdentifierType, Strand, TranscriptData,
    TranscriptSearch,
};
use hgvs_weaver::equivalence::{EquivalenceLevel, VariantEquivalence};
use hgvs_weaver::error::HgvsError;
use hgvs_weaver::mapper::VariantMapper;
use hgvs_weaver::parse_hgvs_variant;
use hgvs_weaver::structs::IntervalSpdi;
use hgvs_weaver::SequenceVariant;

const REF_AC: &str = "NC_TEST.1";

/// One exon: transcript indices 0..=99 <-> genomic indices 1000..=1099.
fn transcript(ac: &str, strand: Strand, cds_start: i32, cds_end: i32) -> TranscriptData {
    TranscriptData {
        ac: ac.to_string(),
        gene: "TEST".to_string(),
        cds_start_index: Some(TranscriptPos(cds_start)),
        cds_end_index: Some(TranscriptPos(cds_end)),
        strand,
        reference_accession: REF_AC.to_string(),
        exons: vec![ExonData {
            transcript_start: TranscriptPos(0),
            transcript_end: TranscriptPos(100),
            reference_start: GenomicPos(1000),
            reference_end: GenomicPos(1099),
            alt_strand: strand,
            cigar: "100M".to_string(),
        }],
    }
}

struct Provider;

impl Provider {
    fn genome() -> String {
        "ACGT".repeat(500)
    }
}

impl DataProvider for Provider {
    fn get_transcript(&self, ac: &str, _ref_ac: Option<&str>) -> Result<TranscriptData, HgvsError> {
        match ac {
            // CDS at transcript indices 10..=39 on the plus strand.
            "NM_PLUS10.1" => Ok(transcript(ac, Strand::Plus, 10, 39)),
            // CDS spanning the whole transcript on the plus strand.
            "NM_PLUS0.1" => Ok(transcript(ac, Strand::Plus, 0, 99)),
            // CDS at transcript indices 10..=39 on the minus strand.
            "NM_MINUS10.1" => Ok(transcript(ac, Strand::Minus, 10, 39)),
            _ => Err(HgvsError::DataProviderError(format!("unknown {}", ac))),
        }
    }

    fn get_seq(
        &self,
        ac: &str,
        start: i32,
        end: i32,
        _kind: IdentifierType,
    ) -> Result<String, HgvsError> {
        let genome = Self::genome();
        let seq: String = if ac == REF_AC {
            genome
        } else {
            let tx = self.get_transcript(ac, None)?;
            let exonic = &genome[1000..1100];
            match tx.strand {
                Strand::Plus => exonic.to_string(),
                Strand::Minus => exonic
                    .chars()
                    .rev()
                    .map(|c| match c {
                        'A' => 'T',
                        'C' => 'G',
                        'G' => 'C',
                        'T' => 'A',
                        other => other,
                    })
                    .collect(),
            }
        };
        let s = start.max(0) as usize;
        let e = if end < 0 {
            seq.len()
        } else {
            (end as usize).min(seq.len())
        };
        Ok(seq[s.min(e)..e].to_string())
    }

    fn get_symbol_accessions(
        &self,
        _symbol: &str,
        _source: IdentifierKind,
        _target: IdentifierKind,
    ) -> Result<Vec<(IdentifierType, String)>, HgvsError> {
        Ok(vec![])
    }

    fn get_identifier_type(&self, id: &str) -> Result<IdentifierType, HgvsError> {
        Ok(if id.starts_with("NC_") {
            IdentifierType::GenomicAccession
        } else {
            IdentifierType::TranscriptAccession
        })
    }
}

impl TranscriptSearch for Provider {
    fn get_transcripts_for_region(
        &self,
        _chrom: &str,
        _start: i32,
        _end: i32,
    ) -> Result<Vec<String>, HgvsError> {
        Ok(vec![])
    }
}

fn coding_interval(hgvs: &str) -> hgvs_weaver::structs::BaseOffsetInterval {
    match parse_hgvs_variant(hgvs).unwrap() {
        SequenceVariant::Coding(v) => v.posedit.pos.unwrap(),
        other => panic!("expected a c. variant, got {:?}", other),
    }
}

#[test]
fn spdi_interval_honours_cds_end_anchor() {
    // c.*1 is the base after the stop codon: transcript index 40, genomic 1040.
    let iv = coding_interval("NM_PLUS10.1:c.*1A>G");
    let got = iv.spdi_interval("NM_PLUS10.1", &Provider).unwrap();
    assert_eq!(got, (1040, 1041, REF_AC.to_string()));
}

#[test]
fn spdi_interval_honours_non_zero_cds_start() {
    // c.1 is transcript index 10, genomic 1010.
    let iv = coding_interval("NM_PLUS10.1:c.1A>G");
    let got = iv.spdi_interval("NM_PLUS10.1", &Provider).unwrap();
    assert_eq!(got, (1010, 1011, REF_AC.to_string()));
}

#[test]
fn spdi_interval_on_minus_strand() {
    // Minus strand: transcript index 0 is genomic 1099, so c.1 (index 10) is 1089.
    let iv = coding_interval("NM_MINUS10.1:c.1A>G");
    let got = iv.spdi_interval("NM_MINUS10.1", &Provider).unwrap();
    assert_eq!(got, (1089, 1090, REF_AC.to_string()));

    // A multi-base interval is reported low-to-high on the genome.
    let iv = coding_interval("NM_MINUS10.1:c.1_3del");
    let got = iv.spdi_interval("NM_MINUS10.1", &Provider).unwrap();
    assert_eq!(got, (1087, 1090, REF_AC.to_string()));
}

#[test]
fn spdi_interval_applies_intronic_offset_by_strand() {
    // On the minus strand a +5 intronic offset moves towards lower genomic indices.
    let iv = coding_interval("NM_MINUS10.1:c.1+5A>G");
    let got = iv.spdi_interval("NM_MINUS10.1", &Provider).unwrap();
    assert_eq!(got, (1084, 1085, REF_AC.to_string()));
}

#[test]
fn coding_variants_on_transcripts_with_different_cds_starts_are_equivalent() {
    // Both name genomic index 1040 (reference base A): c.*1 on the 10..=39 CDS
    // and c.41 on the whole-transcript CDS.
    let v1 = parse_hgvs_variant("NM_PLUS10.1:c.*1A>G").unwrap();
    let v2 = parse_hgvs_variant("NM_PLUS0.1:c.41A>G").unwrap();
    let eq = VariantEquivalence::new(&Provider, &Provider);
    assert_eq!(
        eq.equivalent_level(&v1, &v2).unwrap(),
        EquivalenceLevel::Analogous
    );

    // And a different base is not.
    let v3 = parse_hgvs_variant("NM_PLUS0.1:c.42C>G").unwrap();
    assert_eq!(
        eq.equivalent_level(&v1, &v3).unwrap(),
        EquivalenceLevel::Different
    );
}

#[test]
fn validate_checks_stated_reference_through_transcript_coordinates() {
    let mapper = VariantMapper::new(&Provider);
    let ok = |hgvs: &str| mapper.validate(&parse_hgvs_variant(hgvs).unwrap()).unwrap();

    // Plus strand, CDS at index 10: c.1 is transcript index 10, genome[1010] = G.
    assert!(ok("NM_PLUS10.1:c.1G>A"));
    assert!(!ok("NM_PLUS10.1:c.1A>G"));
    // c.*1 is transcript index 40, genome[1040] = A.
    assert!(ok("NM_PLUS10.1:c.*1A>G"));
    // Minus strand: transcript index 10 is the complement of genome[1089] = C.
    assert!(ok("NM_MINUS10.1:c.1G>A"));
    assert!(!ok("NM_MINUS10.1:c.1C>A"));
    // Intronic positions are accepted unchecked; genomic goes straight to the reference.
    assert!(ok("NM_MINUS10.1:c.1+5T>A"));
    assert!(ok("NC_TEST.1:g.1011G>A"));
    assert!(!ok("NC_TEST.1:g.1011A>G"));
}
