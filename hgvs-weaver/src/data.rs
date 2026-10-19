use crate::error::HgvsError;
use crate::structs::{GenomicPos, TranscriptPos};
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use std::convert::TryFrom;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Plus,
    Minus,
}

impl Strand {
    pub fn value(&self) -> i32 {
        match self {
            Strand::Plus => 1,
            Strand::Minus => -1,
        }
    }
}

impl TryFrom<i32> for Strand {
    type Error = HgvsError;

    fn try_from(value: i32) -> Result<Self, Self::Error> {
        match value {
            1 => Ok(Strand::Plus),
            -1 => Ok(Strand::Minus),
            _ => Err(HgvsError::Other(format!("Invalid strand value: {}", value))),
        }
    }
}

impl Serialize for Strand {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        serializer.serialize_i32(self.value())
    }
}

impl<'de> Deserialize<'de> for Strand {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let value = i32::deserialize(deserializer)?;
        Strand::try_from(value).map_err(serde::de::Error::custom)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExonData {
    pub transcript_start: TranscriptPos,
    pub transcript_end: TranscriptPos,
    pub reference_start: GenomicPos,
    pub reference_end: GenomicPos,
    pub alt_strand: Strand,
    pub cigar: String,
}

/// What the codon at a transcript's `cds_start_index` is.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[non_exhaustive]
pub enum CdsStart {
    /// The start codon: translation initiates there, so residue 1 is Met
    /// whatever the codon, and a change to the codon is `p.(Met1?)`.
    #[default]
    Initiator,
    /// A codon inside the protein, from which the coding sequence is
    /// numbered without translation initiating there: a target that is part
    /// of a protein, as MAVE-HGVS numbers it. Residue 1 is the codon's own
    /// residue, and a change to it is read as a change to any other codon is.
    Internal,
    /// The CDS runs off the transcript's 5' end: `cds_start_index` is the
    /// first coding base the record carries, not the A of a start codon. No
    /// position is numbered from the start codon, and no protein consequence
    /// is predicted.
    Open,
}

/// A transcript model as a provider writes it. `cds_start` is two optional
/// keys, `cds_start_open` (default false) and `cds_start_initiates` (default
/// true), and a start that is open and initiates is refused.
#[derive(Serialize, Deserialize)]
struct TranscriptRecord {
    ac: String,
    gene: String,
    cds_start_index: Option<TranscriptPos>,
    cds_end_index: Option<TranscriptPos>,
    #[serde(default)]
    cds_start_open: bool,
    #[serde(
        default,
        deserialize_with = "present_bool",
        skip_serializing_if = "Option::is_none"
    )]
    cds_start_initiates: Option<bool>,
    #[serde(default)]
    cds_end_open: bool,
    strand: Strand,
    reference_accession: String,
    exons: Vec<ExonData>,
}

/// A key that, when present, holds a boolean: `null` is refused, not read as absent.
fn present_bool<'de, D: Deserializer<'de>>(deserializer: D) -> Result<Option<bool>, D::Error> {
    bool::deserialize(deserializer).map(Some)
}

impl TryFrom<TranscriptRecord> for TranscriptData {
    type Error = String;

    fn try_from(r: TranscriptRecord) -> Result<Self, Self::Error> {
        let cds_start = match (r.cds_start_open, r.cds_start_initiates) {
            (true, Some(true)) => {
                return Err(format!(
                    "{}: cds_start_open and cds_start_initiates are both true, but a CDS that runs \
                     off the record's 5' end has no codon there for translation to initiate at",
                    r.ac
                ))
            }
            (true, _) => CdsStart::Open,
            (false, Some(false)) => CdsStart::Internal,
            (false, _) => CdsStart::Initiator,
        };
        Ok(TranscriptData {
            ac: r.ac,
            gene: r.gene,
            cds_start_index: r.cds_start_index,
            cds_end_index: r.cds_end_index,
            cds_start,
            cds_end_open: r.cds_end_open,
            strand: r.strand,
            reference_accession: r.reference_accession,
            exons: r.exons,
        })
    }
}

impl From<TranscriptData> for TranscriptRecord {
    fn from(t: TranscriptData) -> Self {
        let (cds_start_open, cds_start_initiates) = match t.cds_start {
            CdsStart::Initiator => (false, None),
            CdsStart::Internal => (false, Some(false)),
            CdsStart::Open => (true, None),
        };
        TranscriptRecord {
            ac: t.ac,
            gene: t.gene,
            cds_start_index: t.cds_start_index,
            cds_end_index: t.cds_end_index,
            cds_start_open,
            cds_start_initiates,
            cds_end_open: t.cds_end_open,
            strand: t.strand,
            reference_accession: t.reference_accession,
            exons: t.exons,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(try_from = "TranscriptRecord", into = "TranscriptRecord")]
pub struct TranscriptData {
    pub ac: String,
    pub gene: String,
    pub cds_start_index: Option<TranscriptPos>,
    pub cds_end_index: Option<TranscriptPos>,
    /// What the codon at `cds_start_index` is: the start codon, a codon
    /// inside the protein, or a codon past a 5' end the CDS runs off.
    pub cds_start: CdsStart,
    /// The CDS runs off the transcript's 3' end: `cds_end_index` is the last
    /// coding base the record carries, not the last base of a stop codon. No
    /// position is numbered from the stop codon, and no protein consequence is
    /// predicted that needs a base past the last whole codon the record carries.
    pub cds_end_open: bool,
    pub strand: Strand,
    pub reference_accession: String,
    pub exons: Vec<ExonData>,
}

/// Interface for retrieving transcript and sequence data.
pub trait DataProvider {
    fn get_transcript(
        &self,
        transcript_ac: &str,
        reference_ac: Option<&str>,
    ) -> Result<TranscriptData, HgvsError>;
    /// Returns the bases of `ac` in the 0-based half-open range `[start, end)`,
    /// or from `start` to the end of the sequence when `end` is `None`.
    ///
    /// A range that extends past the end must return the bases that exist, not
    /// an error or an empty string: the core pages through sequences in
    /// fixed-size blocks and relies on a short final block to learn where the
    /// sequence ends.
    fn get_seq(
        &self,
        ac: &str,
        start: i32,
        end: Option<i32>,
        kind: IdentifierType,
    ) -> Result<String, HgvsError>;
    fn get_symbol_accessions(
        &self,
        symbol: &str,
        source_kind: IdentifierKind,
        target_kind: IdentifierKind,
    ) -> Result<Vec<(IdentifierType, String)>, HgvsError>;
    fn get_identifier_type(&self, identifier: &str) -> Result<IdentifierType, HgvsError>;
}

/// Interface for discovering transcripts by region.
pub trait TranscriptSearch {
    fn get_transcripts_for_region(
        &self,
        chrom: &str,
        start: i32,
        end: i32,
    ) -> Result<Vec<String>, HgvsError>;
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum IdentifierKind {
    Genomic,
    Transcript,
    Protein,
}

impl IdentifierKind {
    pub fn into_identifier_type(&self) -> IdentifierType {
        match self {
            IdentifierKind::Genomic => IdentifierType::GenomicAccession,
            IdentifierKind::Transcript => IdentifierType::TranscriptAccession,
            IdentifierKind::Protein => IdentifierType::ProteinAccession,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum IdentifierType {
    GenomicAccession,
    TranscriptAccession,
    ProteinAccession,
    GeneSymbol,
    Unknown,
}
