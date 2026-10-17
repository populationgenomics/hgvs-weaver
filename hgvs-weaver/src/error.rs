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

/// Error types for HGVS operations.
#[derive(Debug)]
pub enum HgvsError {
    /// Failure during parsing of an HGVS string or CIGAR string.
    PestError(String),
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
    TranscriptMismatch {
        accession: String,
        expected: String,
        found: String,
        start: usize,
        end: usize,
    },
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
}

impl fmt::Display for HgvsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            HgvsError::PestError(msg) => write!(f, "Parse error: {msg}"),
            HgvsError::ValidationError(msg) => write!(f, "Validation error: {msg}"),
            HgvsError::DataProviderError(msg) => write!(f, "Data provider error: {msg}"),
            HgvsError::UnsupportedOperation(msg) => write!(f, "Unsupported operation: {msg}"),
            HgvsError::CigarError(msg) => write!(f, "CIGAR error: {msg}"),
            HgvsError::TranscriptMismatch {
                accession,
                expected,
                found,
                start,
                end,
            } => write!(
                f,
                "Reference sequence mismatch: expected {expected}, found {found} on {accession} at indices {start}..{end}"
            ),
            HgvsError::Other(msg) => write!(f, "Other error: {msg}"),
            HgvsError::CisMember {
                allele,
                index,
                member,
                source,
            } => write!(f, "Member {} of {allele}, {member}: {source}", index + 1),
        }
    }
}

/// A member's error is part of the message, so it is not also a `source`,
/// which a reporter walking the chain would print twice.
impl std::error::Error for HgvsError {}
