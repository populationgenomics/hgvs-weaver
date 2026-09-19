//! The refget seam: a sequence's refget accession (`SQ.` + sha512t24u of its
//! bases) and, the other way, the accession of the sequence behind a refget
//! accession. Separate from [`crate::data::DataProvider`] because a sequence
//! source need know nothing about digests; a [`crate::reference::ReferenceStore`]
//! without a `Refget` computes accessions from the whole sequence and cannot
//! answer the reverse lookup.

use crate::error::HgvsError;

pub trait Refget {
    /// The refget accession of `ac`, if known. `None` lets the store compute it.
    fn refget_accession(&self, ac: &str) -> Result<Option<String>, HgvsError>;
    /// The accession of the sequence whose refget accession is `refget`, if known.
    fn accession_for_refget(&self, refget: &str) -> Result<Option<String>, HgvsError>;
}
