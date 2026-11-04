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
#[non_exhaustive]
pub struct ExonData {
    pub transcript_start: TranscriptPos,
    pub transcript_end: TranscriptPos,
    pub reference_start: GenomicPos,
    pub reference_end: GenomicPos,
    pub alt_strand: Strand,
    pub cigar: String,
}

impl ExonData {
    /// An exon over the transcript indices `[transcript_start, transcript_end)`
    /// and the genomic indices `[reference_start, reference_end]`, aligned by
    /// `cigar` on `alt_strand`.
    pub fn new(
        transcript_start: TranscriptPos,
        transcript_end: TranscriptPos,
        reference_start: GenomicPos,
        reference_end: GenomicPos,
        alt_strand: Strand,
        cigar: impl Into<String>,
    ) -> Self {
        ExonData {
            transcript_start,
            transcript_end,
            reference_start,
            reference_end,
            alt_strand,
            cigar: cigar.into(),
        }
    }
}

/// Where a transcript lies on a genome: the genomic sequence
/// `reference_accession`, the strand the transcript reads, and its exons,
/// each aligned to that sequence.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[non_exhaustive]
pub struct Placement {
    pub reference_accession: String,
    pub strand: Strand,
    pub exons: Vec<ExonData>,
}

impl Placement {
    /// The exons `exons`, on `strand` of `reference_accession`.
    pub fn new(
        reference_accession: impl Into<String>,
        strand: Strand,
        exons: Vec<ExonData>,
    ) -> Self {
        Placement {
            reference_accession: reference_accession.into(),
            strand,
            exons,
        }
    }
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
/// true), and a start that is open and initiates is refused. The placement is
/// three keys, `strand`, `reference_accession` and `exons`, all present or
/// all absent.
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
        deserialize_with = "present",
        skip_serializing_if = "Option::is_none"
    )]
    cds_start_initiates: Option<bool>,
    #[serde(default)]
    cds_end_open: bool,
    #[serde(
        default,
        deserialize_with = "present",
        skip_serializing_if = "Option::is_none"
    )]
    strand: Option<Strand>,
    #[serde(
        default,
        deserialize_with = "present",
        skip_serializing_if = "Option::is_none"
    )]
    reference_accession: Option<String>,
    #[serde(
        default,
        deserialize_with = "present",
        skip_serializing_if = "Option::is_none"
    )]
    exons: Option<Vec<ExonData>>,
}

/// A key that, when present, holds a value: `null` is refused, not read as absent.
fn present<'de, D: Deserializer<'de>, T: Deserialize<'de>>(
    deserializer: D,
) -> Result<Option<T>, D::Error> {
    T::deserialize(deserializer).map(Some)
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
        let placement = match (r.reference_accession, r.strand, r.exons) {
            (Some(reference_accession), Some(strand), Some(exons)) => Some(Placement {
                reference_accession,
                strand,
                exons,
            }),
            (None, None, None) => None,
            (reference_accession, strand, exons) => {
                let missing: Vec<&str> = [
                    ("reference_accession", reference_accession.is_none()),
                    ("strand", strand.is_none()),
                    ("exons", exons.is_none()),
                ]
                .into_iter()
                .filter_map(|(key, absent)| absent.then_some(key))
                .collect();
                return Err(format!(
                    "{}: a placement is reference_accession, strand and exons together, \
                     and {} absent",
                    r.ac,
                    missing.join(" and ")
                ));
            }
        };
        Ok(TranscriptData {
            ac: r.ac,
            gene: r.gene,
            cds_start_index: r.cds_start_index,
            cds_end_index: r.cds_end_index,
            cds_start,
            cds_end_open: r.cds_end_open,
            placement,
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
        let (reference_accession, strand, exons) = match t.placement {
            Some(p) => (Some(p.reference_accession), Some(p.strand), Some(p.exons)),
            None => (None, None, None),
        };
        TranscriptRecord {
            ac: t.ac,
            gene: t.gene,
            cds_start_index: t.cds_start_index,
            cds_end_index: t.cds_end_index,
            cds_start_open,
            cds_start_initiates,
            cds_end_open: t.cds_end_open,
            strand,
            reference_accession,
            exons,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(try_from = "TranscriptRecord", into = "TranscriptRecord")]
#[non_exhaustive]
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
    /// Where the transcript lies on a genome; `None` for a record no
    /// assembly places, such as a GenBank mRNA. An unplaced record is read on
    /// its own sequence alone: an operation that needs the genome is
    /// [`HgvsError::Unplaced`].
    pub placement: Option<Placement>,
}

impl TranscriptData {
    /// A transcript model whose CDS, where it has one, starts at a start codon
    /// and ends at a stop codon the record carries. `cds` is the CDS's first
    /// and last transcript index, inclusive, and `None` for a non-coding
    /// transcript; [`with_cds_start`](Self::with_cds_start) and
    /// [`with_cds_end_open`](Self::with_cds_end_open) set the other cases.
    ///
    /// ```
    /// use hgvs_weaver::coords::{GenomicPos, TranscriptPos};
    /// use hgvs_weaver::data::{CdsStart, ExonData, Strand, TranscriptData};
    ///
    /// // One 40-base exon; the CDS is transcript indices 5 to 34 and runs on
    /// // past the record's 3' end.
    /// let model = TranscriptData::new(
    ///     "NM_000001.1",
    ///     "GENE",
    ///     Some((TranscriptPos(5), TranscriptPos(34))),
    ///     Strand::Plus,
    ///     "NC_000001.11",
    ///     vec![ExonData::new(
    ///         TranscriptPos(0),
    ///         TranscriptPos(40),
    ///         GenomicPos(1000),
    ///         GenomicPos(1039),
    ///         Strand::Plus,
    ///         "40=",
    ///     )],
    /// )
    /// .with_cds_end_open(true);
    /// assert_eq!(model.cds_start, CdsStart::Initiator);
    /// assert!(model.cds_end_open);
    /// ```
    pub fn new(
        ac: impl Into<String>,
        gene: impl Into<String>,
        cds: Option<(TranscriptPos, TranscriptPos)>,
        strand: Strand,
        reference_accession: impl Into<String>,
        exons: Vec<ExonData>,
    ) -> Self {
        TranscriptData {
            ac: ac.into(),
            gene: gene.into(),
            cds_start_index: cds.map(|(start, _)| start),
            cds_end_index: cds.map(|(_, end)| end),
            cds_start: CdsStart::Initiator,
            cds_end_open: false,
            placement: Some(Placement::new(reference_accession, strand, exons)),
        }
    }

    /// A transcript model as [`new`](Self::new) makes one, but placed on no
    /// genome.
    pub fn unplaced(
        ac: impl Into<String>,
        gene: impl Into<String>,
        cds: Option<(TranscriptPos, TranscriptPos)>,
    ) -> Self {
        TranscriptData {
            ac: ac.into(),
            gene: gene.into(),
            cds_start_index: cds.map(|(start, _)| start),
            cds_end_index: cds.map(|(_, end)| end),
            cds_start: CdsStart::Initiator,
            cds_end_open: false,
            placement: None,
        }
    }

    /// Where the transcript lies on a genome, for an operation that needs it:
    /// [`HgvsError::Unplaced`] for a record placed on none.
    pub fn placed(&self) -> Result<&Placement, HgvsError> {
        self.placement.as_ref().ok_or_else(|| HgvsError::Unplaced {
            accession: self.ac.clone(),
        })
    }

    /// The same model with the codon at `cds_start_index` read as `cds_start`.
    pub fn with_cds_start(self, cds_start: CdsStart) -> Self {
        TranscriptData { cds_start, ..self }
    }

    /// The same model with the CDS running off the record's 3' end, or not.
    pub fn with_cds_end_open(self, cds_end_open: bool) -> Self {
        TranscriptData {
            cds_end_open,
            ..self
        }
    }
}

/// Interface for retrieving transcript and sequence data.
pub trait DataProvider {
    /// The model of `transcript_ac`, placed on `reference_ac` where one is
    /// given. A record placed on no genome is returned unplaced whatever
    /// `reference_ac` is, so the caller is told the placement is missing
    /// ([`HgvsError::Unplaced`]) rather than that the lookup failed.
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

/// Interface for what a store holds under an accession: the versions of an
/// unversioned one, and the status and tags of a versioned one. Optional, as
/// [`TranscriptSearch`] is: a provider that cannot list versions is still a
/// [`DataProvider`].
pub trait AccessionCatalogue {
    /// Every versioned accession held under the unversioned `ac`, newest
    /// first: `ac` followed by `.` and its version, or `ac` itself, alone, for
    /// a record with no version of its own, such as a UniProt accession. Empty
    /// when the store holds none.
    fn get_accession_versions(&self, ac: &str) -> Result<Vec<String>, HgvsError>;
    /// The status and tags of the versioned accession `ac`, `None` when the
    /// store holds nothing under it. A protein carries the tags of the
    /// transcripts that encode it.
    fn get_accession_info(&self, ac: &str) -> Result<Option<AccessionInfo>, HgvsError>;
}

/// What a store records about a versioned accession.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct AccessionInfo {
    pub status: AccessionStatus,
    pub tags: Vec<AccessionTag>,
}

/// Whether an accession is the version its publisher's current release names.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[non_exhaustive]
pub enum AccessionStatus {
    /// The version the current release names.
    Current,
    /// Replaced by a later version.
    Superseded,
    /// Withdrawn.
    Suppressed,
}

/// A designation a curated set gives an accession.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[non_exhaustive]
pub enum AccessionTag {
    /// MANE Select: the one transcript MANE chose for a protein-coding gene.
    ManeSelect,
    /// MANE Plus Clinical: a transcript MANE adds for clinical reporting.
    ManePlusClinical,
    /// RefSeq Select: NCBI's representative transcript for a gene.
    RefSeqSelect,
    /// Ensembl canonical: Ensembl's representative transcript for a gene.
    EnsemblCanonical,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[non_exhaustive]
pub enum IdentifierKind {
    Genomic,
    Transcript,
    Protein,
    /// A gene symbol: `get_symbol_accessions(gene, Gene, Protein)` names the
    /// gene's proteins.
    Gene,
}

impl IdentifierKind {
    pub fn into_identifier_type(&self) -> IdentifierType {
        match self {
            IdentifierKind::Genomic => IdentifierType::GenomicAccession,
            IdentifierKind::Transcript => IdentifierType::TranscriptAccession,
            IdentifierKind::Protein => IdentifierType::ProteinAccession,
            IdentifierKind::Gene => IdentifierType::GeneSymbol,
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
