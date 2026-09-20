//! One in-memory [`DataProvider`] for the integration tests, so each test file
//! declares its sequences and transcripts instead of re-implementing the trait.
//!
//! Include it from a test binary with `mod support;` (or, from a directory
//! test such as `properties/`, `#[path = "../support/mod.rs"] mod support;`).
//!
//! The fixture serves named sequences over any `[start, end)` range, clamped to
//! the sequence as the [`DataProvider::get_seq`] contract requires; serves the
//! transcripts it was given; maps transcripts to proteins and back through
//! [`DataProvider::get_symbol_accessions`]; classifies accessions by prefix
//! unless told otherwise; finds transcripts by reference for
//! [`TranscriptSearch`]; and looks sequences up by refget digest for
//! [`Refget`], so `from_vrs` round-trips work.

#![allow(dead_code)]

use hgvs_weaver::coords::{GenomicPos, TranscriptPos};
use hgvs_weaver::data::{
    DataProvider, ExonData, IdentifierKind, IdentifierType, Strand, TranscriptData,
    TranscriptSearch,
};
use hgvs_weaver::error::HgvsError;
use hgvs_weaver::refget::Refget;
use hgvs_weaver::vrs::refget_accession;
use serde::Deserialize;
use std::collections::BTreeMap;
use std::fs;

/// In-memory provider over the sequences and transcripts it is given.
#[derive(Debug, Clone, Default)]
pub struct Provider {
    sequences: BTreeMap<String, String>,
    transcripts: Vec<TranscriptData>,
    /// `(transcript accession, protein accession)` pairs.
    proteins: Vec<(String, String)>,
    identifier_types: BTreeMap<String, IdentifierType>,
    /// When set, every `get_seq` fails: a provider that has no sequences at all.
    failing_sequences: bool,
}

impl Provider {
    pub fn new() -> Self {
        Self::default()
    }

    /// A provider over `{"sequences": {ac: seq}, "transcripts": [...] | {ac: ...}}`.
    pub fn from_json(json: &str) -> Self {
        let data: JsonData = serde_json::from_str(json).expect("provider JSON");
        let mut provider = Self::new();
        for (ac, seq) in data.sequences {
            provider = provider.sequence(&ac, &seq);
        }
        for transcript in data.transcripts.into_vec() {
            provider = provider.transcript(transcript);
        }
        provider
    }

    /// [`from_json`](Self::from_json) on the contents of `path`, relative to
    /// the crate root (`hgvs-weaver/`), as cargo runs the tests.
    pub fn from_json_file(path: &str) -> Self {
        let text = fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
        Self::from_json(&text)
    }

    /// Serves `seq` as the sequence of `ac`.
    pub fn sequence(mut self, ac: &str, seq: &str) -> Self {
        self.sequences.insert(ac.to_string(), seq.to_string());
        self
    }

    /// Serves `transcript` under its own accession.
    pub fn transcript(mut self, transcript: TranscriptData) -> Self {
        self.transcripts.push(transcript);
        self
    }

    /// Maps `tx_ac` to `np_ac` (Transcript -> Protein) and back.
    pub fn protein_for(mut self, tx_ac: &str, np_ac: &str) -> Self {
        self.proteins.push((tx_ac.to_string(), np_ac.to_string()));
        self
    }

    /// Classifies `ac` as `kind`, whatever its prefix says.
    pub fn identifier_type(mut self, ac: &str, kind: IdentifierType) -> Self {
        self.identifier_types.insert(ac.to_string(), kind);
        self
    }

    /// Makes every `get_seq` fail, for tests of what happens without a sequence.
    pub fn failing_sequences(mut self) -> Self {
        self.failing_sequences = true;
        self
    }

    /// The sequences this provider serves, by accession.
    pub fn sequences(&self) -> &BTreeMap<String, String> {
        &self.sequences
    }
}

/// The identifier type an accession's prefix implies.
pub fn identifier_type_from_prefix(ac: &str) -> IdentifierType {
    match ac.get(..3) {
        Some("NC_" | "NT_" | "NW_" | "NG_") => IdentifierType::GenomicAccession,
        Some("NM_" | "NR_" | "XM_" | "XR_") => IdentifierType::TranscriptAccession,
        Some("NP_" | "XP_") => IdentifierType::ProteinAccession,
        _ => IdentifierType::Unknown,
    }
}

/// An exon aligning transcript indices `[transcript.0, transcript.1)` to
/// reference indices `reference.0..=reference.1` with a match-only cigar.
pub fn exon(transcript: (i32, i32), reference: (i32, i32), strand: Strand) -> ExonData {
    ExonData {
        transcript_start: TranscriptPos(transcript.0),
        transcript_end: TranscriptPos(transcript.1),
        reference_start: GenomicPos(reference.0),
        reference_end: GenomicPos(reference.1),
        alt_strand: strand,
        cigar: format!("{}M", transcript.1 - transcript.0),
    }
}

/// A transcript on `reference_ac` with the given exons; `cds` is the 0-based
/// index of the first CDS base and of the stop codon's last base, `None` for
/// a non-coding transcript.
pub fn transcript(
    ac: &str,
    reference_ac: &str,
    strand: Strand,
    cds: Option<(i32, i32)>,
    exons: Vec<ExonData>,
) -> TranscriptData {
    TranscriptData {
        ac: ac.to_string(),
        gene: "TEST".to_string(),
        cds_start_index: cds.map(|(start, _)| TranscriptPos(start)),
        cds_end_index: cds.map(|(_, end)| TranscriptPos(end)),
        strand,
        reference_accession: reference_ac.to_string(),
        exons,
    }
}

/// A `length`-base transcript that is one exon at reference indices
/// `genomic_start..=genomic_start + length - 1`, with the CDS at transcript
/// indices `cds_start..=cds_end`.
pub fn single_exon_transcript(
    ac: &str,
    reference_ac: &str,
    genomic_start: i32,
    strand: Strand,
    cds_start: i32,
    cds_end: i32,
    length: i32,
) -> TranscriptData {
    transcript(
        ac,
        reference_ac,
        strand,
        Some((cds_start, cds_end)),
        vec![exon(
            (0, length),
            (genomic_start, genomic_start + length - 1),
            strand,
        )],
    )
}

impl DataProvider for Provider {
    fn get_transcript(&self, ac: &str, _: Option<&str>) -> Result<TranscriptData, HgvsError> {
        self.transcripts
            .iter()
            .find(|t| t.ac == ac)
            .cloned()
            .ok_or_else(|| HgvsError::DataProviderError(format!("no transcript {ac}")))
    }

    fn get_seq(
        &self,
        ac: &str,
        start: i32,
        end: Option<i32>,
        _: IdentifierType,
    ) -> Result<String, HgvsError> {
        if self.failing_sequences {
            return Err(HgvsError::Other("Missing sequence".into()));
        }
        let seq = self
            .sequences
            .get(ac)
            .ok_or_else(|| HgvsError::DataProviderError(format!("no sequence {ac}")))?;
        let start = (start.max(0) as usize).min(seq.len());
        let end = end.map_or(seq.len(), |e| (e.max(0) as usize).min(seq.len()));
        Ok(seq[start..end.max(start)].to_string())
    }

    fn get_symbol_accessions(
        &self,
        symbol: &str,
        _: IdentifierKind,
        target: IdentifierKind,
    ) -> Result<Vec<(IdentifierType, String)>, HgvsError> {
        Ok(match target {
            IdentifierKind::Protein => self
                .proteins
                .iter()
                .filter(|(tx, _)| tx == symbol)
                .map(|(_, np)| (IdentifierType::ProteinAccession, np.clone()))
                .collect(),
            IdentifierKind::Transcript => self
                .proteins
                .iter()
                .filter(|(_, np)| np == symbol)
                .map(|(tx, _)| (IdentifierType::TranscriptAccession, tx.clone()))
                .collect(),
            IdentifierKind::Genomic => vec![],
        })
    }

    fn get_identifier_type(&self, ac: &str) -> Result<IdentifierType, HgvsError> {
        Ok(self
            .identifier_types
            .get(ac)
            .copied()
            .unwrap_or_else(|| identifier_type_from_prefix(ac)))
    }
}

impl TranscriptSearch for Provider {
    /// The transcripts on `chrom`; all of them if none is.
    fn get_transcripts_for_region(
        &self,
        chrom: &str,
        _: i32,
        _: i32,
    ) -> Result<Vec<String>, HgvsError> {
        let on_chrom: Vec<String> = self
            .transcripts
            .iter()
            .filter(|t| t.reference_accession == chrom)
            .map(|t| t.ac.clone())
            .collect();
        if !on_chrom.is_empty() {
            return Ok(on_chrom);
        }
        Ok(self.transcripts.iter().map(|t| t.ac.clone()).collect())
    }
}

impl Refget for Provider {
    /// `None`: the store computes the digest from the sequence.
    fn refget_accession(&self, _: &str) -> Result<Option<String>, HgvsError> {
        Ok(None)
    }

    fn accession_for_refget(&self, refget: &str) -> Result<Option<String>, HgvsError> {
        Ok(self
            .sequences
            .iter()
            .find(|(_, seq)| refget_accession(seq) == refget)
            .map(|(ac, _)| ac.clone()))
    }
}

#[derive(Deserialize)]
struct JsonData {
    #[serde(default)]
    sequences: BTreeMap<String, String>,
    #[serde(default)]
    transcripts: JsonTranscripts,
}

/// Transcripts as a list, or keyed by accession as `toy_data.json` has them.
#[derive(Deserialize)]
#[serde(untagged)]
enum JsonTranscripts {
    List(Vec<TranscriptData>),
    ByAccession(BTreeMap<String, TranscriptData>),
}

impl JsonTranscripts {
    fn into_vec(self) -> Vec<TranscriptData> {
        match self {
            JsonTranscripts::List(list) => list,
            JsonTranscripts::ByAccession(map) => map.into_values().collect(),
        }
    }
}

impl Default for JsonTranscripts {
    fn default() -> Self {
        JsonTranscripts::List(Vec::new())
    }
}
