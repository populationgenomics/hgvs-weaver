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

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TranscriptData {
    pub ac: String,
    pub gene: String,
    pub cds_start_index: Option<TranscriptPos>,
    pub cds_end_index: Option<TranscriptPos>,
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
    fn get_seq(
        &self,
        ac: &str,
        start: i32,
        end: i32,
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

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum IdentifierType {
    GenomicAccession,
    TranscriptAccession,
    ProteinAccession,
    GeneSymbol,
    Unknown,
}
