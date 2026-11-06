//! Error types for HGVS operations.
//!
//! The `Display` and `Error` impls are written by hand rather than derived. A
//! derive macro is shorter, but every proc-macro crate is a build script that
//! runs at compile time, and that is the surface supply-chain attacks on the
//! Rust ecosystem have used. `pest_derive` and `serde`'s derive already oblige
//! this crate to carry that toolchain, so dropping `thiserror` does not remove
//! it — it removes one more crate whose compromise would execute here, for the
//! price of a page of formatting code.

use std::fmt;

use crate::coords::SequenceVariant;
use crate::data::IdentifierType;

/// What an accession was read as, for an accession its provider says is the
/// other kind of sequence ([`HgvsError::WrongReferenceKind`]).
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum ReferenceUse {
    /// The accession a description of this coordinate type is written on:
    /// `"c"`, `"p"`.
    Description(String),
    /// A protein an operation was given, named as its caller knows it.
    Protein(String),
    /// A nucleotide sequence an operation was given, named as its caller
    /// knows it.
    Nucleotide(String),
}

impl ReferenceUse {
    /// Whether this use reads residues, so needs a protein.
    pub fn needs_protein(&self) -> bool {
        match self {
            ReferenceUse::Description(ct) => ct == "p",
            ReferenceUse::Protein(_) => true,
            ReferenceUse::Nucleotide(_) => false,
        }
    }
}

impl fmt::Display for ReferenceUse {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ReferenceUse::Description(ct) => {
                // The letter's name, not the letter, takes the article: an m., a c.
                let article = if ct.starts_with(['f', 'h', 'l', 'm', 'n', 'o', 'r', 's', 'x']) {
                    "an"
                } else {
                    "a"
                };
                write!(f, "{article} {ct}. description")
            }
            ReferenceUse::Protein(what) | ReferenceUse::Nucleotide(what) => write!(f, "{what}"),
        }
    }
}

/// A residue, or a run of residues, as a description writes it: `text` is
/// `description[start..end]`, by byte offset, which is also the character
/// offset since a description that parses is ASCII. A run the grammar reads
/// as one sequence, the `DPA` of `insDPA` or the `HIS` of `insHIS`, is one
/// token; a lone residue, the `X` of `Arg111X`, is a run of one.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ResidueToken {
    pub text: String,
    pub start: usize,
    pub end: usize,
}

/// Error types for HGVS operations.
#[derive(Debug)]
#[non_exhaustive]
pub enum HgvsError {
    /// Failure during parsing of an HGVS string or CIGAR string.
    PestError(String),
    /// A protein description that writes residues in both amino acid codes:
    /// `one_letter` and `three_letter` are every residue or run written in
    /// each, in the order they appear in `description`, so a caller can
    /// rewrite one code into the other and parse again. `*`, `Xaa` and `X`,
    /// which stand in either code, are in neither, as is a run of nothing
    /// else; one inside a run is part of that run's token.
    MixedResidueCodes {
        description: String,
        one_letter: Vec<ResidueToken>,
        three_letter: Vec<ResidueToken>,
    },
    /// An `X` in a protein description, in either code, where a stop and any
    /// residue could both stand: the new residue of a substitution at a
    /// residue that is not a stop, or one an insertion or deletion-insertion
    /// names. `ambiguous` is each such `X` as written in `description`, one
    /// token per `X`, in order.
    AmbiguousResidue {
        description: String,
        ambiguous: Vec<ResidueToken>,
    },
    /// Failure during validation of transcript or exon metadata.
    ValidationError(String),
    /// Failure during data retrieval from a `DataProvider`.
    DataProviderError(String),
    /// Attempted an operation that is not yet supported.
    UnsupportedOperation(String),
    /// Error specifically related to CIGAR string mapping.
    CigarError(String),
    /// A base or residue the description states is not what the reference
    /// holds: `expected` is the stated bases, `found` the reference's over
    /// the 0-based half-open `[start, end)` of `accession`, read on the strand
    /// the description is written on.
    ReferenceMismatch {
        accession: String,
        expected: String,
        found: String,
        start: usize,
        end: usize,
    },
    /// A description names a position `accession` does not have: the 0-based
    /// half-open `[start, end)` runs past the end of its `length` bases or
    /// residues. An insertion needs both its flanking bases, so one after the
    /// last base runs past the end too.
    PositionOutsideSequence {
        accession: String,
        start: usize,
        end: usize,
        length: usize,
    },
    /// An intronic offset counted from a base that is not the exon border its
    /// sign needs: `+k` counts from the last base of an exon, `-k` from the
    /// first, and the transcript's first and last bases count outward. The
    /// transcript `accession` has no position `position`, written as given
    /// (`c.9-1`).
    NotAnExonBorder { accession: String, position: String },
    /// An intronic offset that counts past the end of its intron, into the
    /// next exon: the transcript `accession` has no position `position`,
    /// written as given (`c.17+11`), since its intron has `intron_length`
    /// bases. An offset from the transcript's first or last base has no
    /// intron to stay inside.
    OffsetPastIntron {
        accession: String,
        position: String,
        intron_length: usize,
    },
    /// An accession read as one kind of sequence, `read_as`, that its
    /// provider says is the other, `kind`: a `p.` description on a
    /// transcript, a `c.` one on a protein, a genome given as the protein to
    /// predict on.
    WrongReferenceKind {
        accession: String,
        read_as: ReferenceUse,
        kind: IdentifierType,
    },
    /// An operation that needs the genome, asked of the transcript record
    /// `accession`, which its provider places on no genome: a GenBank mRNA
    /// read on its own sequence. A projection, an intronic position and a
    /// genomic form all need a placement.
    Unplaced { accession: String },
    /// A description that makes no statement about the protein: an identity
    /// that covers no whole codon. HGVS has no `p.` form for it; `p.(=)`
    /// would claim the whole protein, and `p.?` an effect.
    NoProteinStatement(String),
    /// Catch-all for other error types.
    Other(String),
    /// A member of an allele in cis failed, so the allele did: `source` is
    /// why, for the member written `member`, the `index`th (from 0) of
    /// `allele`.
    CisMember {
        allele: String,
        index: usize,
        member: String,
        source: Box<HgvsError>,
    },
    /// The description has several readings that make different sequences
    /// or state different bases: `readings` is each, written as a description
    /// with one reading, for a caller to choose from. A transcript allele
    /// with an intronic member is read on the genome, and its readings are
    /// written in `g.`.
    Ambiguous {
        description: String,
        readings: Vec<SequenceVariant>,
    },
}

impl fmt::Display for HgvsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            HgvsError::PestError(msg) => write!(f, "Parse error: {msg}"),
            HgvsError::MixedResidueCodes {
                description,
                one_letter,
                three_letter,
            } => {
                // The one-letter token named is the first holding an X, the
                // residue whose meaning the mix leaves open, else the first.
                let x = one_letter.iter().find(|t| t.text.contains('X'));
                let one = x.or(one_letter.first()).map_or("", |t| t.text.as_str());
                let three = three_letter.first().map_or("", |t| t.text.as_str());
                write!(
                    f,
                    "Mixed residue codes: {description} writes {three} in the three-letter \
                     amino acid code and {one} in the one-letter code; write every residue \
                     in one code"
                )?;
                if x.is_some() {
                    write!(f, " (a stop is Ter or *, an unknown residue Xaa)")?;
                }
                Ok(())
            }
            HgvsError::AmbiguousResidue {
                description,
                ambiguous,
            } => {
                write!(f, "Ambiguous residue: {description} writes X at")?;
                for (i, t) in ambiguous.iter().enumerate() {
                    let sep = if i == 0 { " " } else { ", " };
                    write!(f, "{sep}{}..{}", t.start, t.end)?;
                }
                write!(
                    f,
                    ", which could be a stop or any residue; write * for a stop and Xaa for \
                     any residue"
                )
            }
            HgvsError::ValidationError(msg) => write!(f, "Validation error: {msg}"),
            HgvsError::DataProviderError(msg) => write!(f, "Data provider error: {msg}"),
            HgvsError::UnsupportedOperation(msg) => write!(f, "Unsupported operation: {msg}"),
            HgvsError::CigarError(msg) => write!(f, "CIGAR error: {msg}"),
            HgvsError::ReferenceMismatch {
                accession,
                expected,
                found,
                start,
                end,
            } => write!(
                f,
                "Reference sequence mismatch: expected {expected}, found {found} on {accession} at indices {start}..{end}"
            ),
            HgvsError::PositionOutsideSequence {
                accession,
                start,
                end,
                length,
            } => write!(
                f,
                "Position outside the sequence: {accession} has {length} positions, and the 0-based range [{start}, {end}) runs past its end"
            ),
            HgvsError::NotAnExonBorder {
                accession,
                position,
            } => write!(
                f,
                "Not an exon border: {accession}:{position} counts an intronic offset from a base \
                 that is not an exon border; +k counts from the last base of an exon and -k from \
                 the first"
            ),
            HgvsError::OffsetPastIntron {
                accession,
                position,
                intron_length,
            } => write!(
                f,
                "Offset past its intron: {accession}:{position} counts past the end of an intron of \
                 {intron_length} bases"
            ),
            HgvsError::WrongReferenceKind {
                accession,
                read_as,
                kind,
            } => {
                let named = match kind {
                    IdentifierType::TranscriptAccession => "a transcript",
                    IdentifierType::GenomicAccession => "a genomic",
                    IdentifierType::ProteinAccession => "a protein",
                    IdentifierType::GeneSymbol => "a gene symbol",
                    IdentifierType::Unknown => "an unclassified",
                };
                let needs = if read_as.needs_protein() {
                    "a protein"
                } else {
                    "a nucleotide"
                };
                write!(
                    f,
                    "Wrong reference kind: {accession} is {named} accession, and {read_as} needs {needs} reference"
                )
            }
            HgvsError::Unplaced { accession } => write!(
                f,
                "Unplaced: {accession} has no genome placement, so nothing that needs the genome \
                 can be read on it"
            ),
            HgvsError::NoProteinStatement(msg) => write!(f, "No protein statement: {msg}"),
            HgvsError::Other(msg) => write!(f, "Other error: {msg}"),
            HgvsError::CisMember {
                allele,
                index,
                member,
                source,
            } => write!(f, "Member {} of {allele}, {member}: {source}", index + 1),
            HgvsError::Ambiguous {
                description,
                readings,
            } => {
                write!(f, "Ambiguous: {description} reads as ")?;
                for (i, r) in readings.iter().enumerate() {
                    if i > 0 {
                        write!(f, " | ")?;
                    }
                    write!(f, "{r}")?;
                }
                Ok(())
            }
        }
    }
}

/// A member's error is part of the message, so it is not also a `source`,
/// which a reporter walking the chain would print twice.
impl std::error::Error for HgvsError {}
