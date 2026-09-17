//! Random-access, cached view of the reference sequences a [`DataProvider`] serves.
//!
//! This is the one place in the crate that fetches sequence. Callers ask the
//! store for a [`Reference`] by accession and kind, then slice it, read single
//! bases, take the whole thing, or walk it while bases match a pattern. Paging
//! and caching sit behind that interface: a provider only ever answers
//! `get_seq(ac, start, end, kind)` for a block at a time, or once for the whole
//! sequence when a caller genuinely needs all of it.

use std::cell::RefCell;
use std::collections::{BTreeMap, HashMap};

use crate::data::{DataProvider, IdentifierType};
use crate::error::HgvsError;

/// Bases fetched per provider call when paging. A provider call may cross into
/// Python, so this is deliberately far larger than any variant.
pub const BLOCK_SIZE: usize = 4096;

#[derive(Default)]
struct Cached {
    /// Fetched blocks by block index. A block shorter than the block size is
    /// the last one and fixes `len`.
    blocks: BTreeMap<usize, String>,
    /// Sequence length, once a short block or a whole fetch has revealed it.
    len: Option<usize>,
    /// The complete sequence, if a caller asked for all of it.
    whole: Option<String>,
    /// The refget accession, once looked up or computed.
    refget: Option<String>,
}

/// Owns the sequence cache for one [`DataProvider`].
pub struct ReferenceStore<'a> {
    hdp: &'a dyn DataProvider,
    block_size: usize,
    cache: RefCell<HashMap<(String, IdentifierType), Cached>>,
}

impl<'a> ReferenceStore<'a> {
    pub fn new(hdp: &'a dyn DataProvider) -> Self {
        Self::with_block_size(hdp, BLOCK_SIZE)
    }

    /// As [`new`](Self::new) with an explicit block size. Mainly for tests.
    pub fn with_block_size(hdp: &'a dyn DataProvider, block_size: usize) -> Self {
        assert!(block_size > 0, "block size must be positive");
        ReferenceStore {
            hdp,
            block_size,
            cache: RefCell::new(HashMap::new()),
        }
    }

    /// The provider behind this store, for callers that need transcript models.
    pub fn provider(&self) -> &'a dyn DataProvider {
        self.hdp
    }

    /// A handle onto one sequence. Cheap; nothing is fetched until it is used.
    pub fn reference<'s>(&'s self, ac: &str, kind: IdentifierType) -> Reference<'s, 'a> {
        Reference {
            store: self,
            ac: ac.to_string(),
            kind,
        }
    }

    fn fetch(
        &self,
        ac: &str,
        kind: IdentifierType,
        start: usize,
        end: Option<usize>,
    ) -> Result<String, HgvsError> {
        let end_i32 = end
            .map(|e| {
                i32::try_from(e).map_err(|_| {
                    HgvsError::ValidationError(format!("Sequence end {} out of range", e))
                })
            })
            .transpose()?;
        let start_i32 = i32::try_from(start).map_err(|_| {
            HgvsError::ValidationError(format!("Sequence start {} out of range", start))
        })?;
        self.hdp.get_seq(ac, start_i32, end_i32, kind)
    }

    /// Ensures block `index` of `ac` is cached (or known to lie past the end).
    fn ensure_block(&self, ac: &str, kind: IdentifierType, index: usize) -> Result<(), HgvsError> {
        let key = (ac.to_string(), kind);
        {
            let cache = self.cache.borrow();
            if let Some(c) = cache.get(&key) {
                if c.whole.is_some() || c.blocks.contains_key(&index) {
                    return Ok(());
                }
                if let Some(len) = c.len {
                    if index * self.block_size >= len {
                        return Ok(());
                    }
                }
            }
        }
        let start = index * self.block_size;
        let block = self.fetch(ac, kind, start, Some(start + self.block_size))?;
        let mut cache = self.cache.borrow_mut();
        let entry = cache.entry(key).or_default();
        if block.len() < self.block_size {
            entry.len = Some(start + block.len());
        }
        entry.blocks.insert(index, block);
        Ok(())
    }
}

/// One reference sequence, addressed by accession and kind.
pub struct Reference<'s, 'a> {
    store: &'s ReferenceStore<'a>,
    ac: String,
    kind: IdentifierType,
}

impl<'s, 'a> Reference<'s, 'a> {
    pub fn accession(&self) -> &str {
        &self.ac
    }

    /// Bases in `[start, end)`. Shorter than requested if the sequence ends first,
    /// and empty if `start` is at or past the end.
    pub fn slice(&self, start: usize, end: usize) -> Result<String, HgvsError> {
        if end <= start {
            return Ok(String::new());
        }
        {
            let cache = self.store.cache.borrow();
            if let Some(whole) = cache.get(&self.key()).and_then(|c| c.whole.as_ref()) {
                return Ok(clamp_slice(whole, start, end));
            }
        }
        let bs = self.store.block_size;
        let first = start / bs;
        let last = (end - 1) / bs;
        let mut out = String::with_capacity(end - start);
        for index in first..=last {
            self.store.ensure_block(&self.ac, self.kind, index)?;
            let cache = self.store.cache.borrow();
            let entry = cache.get(&self.key());
            let block = match entry.and_then(|c| c.blocks.get(&index)) {
                Some(b) => b,
                None => break, // past the known end
            };
            let block_start = index * bs;
            let lo = start.saturating_sub(block_start);
            let hi = (end - block_start).min(block.len());
            if lo >= hi {
                break;
            }
            out.push_str(&block[lo..hi]);
            if block.len() < bs {
                break;
            }
        }
        Ok(out)
    }

    /// The base at `index`, or `None` past the end.
    pub fn base(&self, index: usize) -> Result<Option<u8>, HgvsError> {
        Ok(self.slice(index, index + 1)?.bytes().next())
    }

    /// The refget accession (`SQ.` + sha512t24u) identifying this sequence:
    /// from the provider if it knows it, otherwise computed from the whole
    /// sequence and cached.
    pub fn refget_accession(&self) -> Result<String, HgvsError> {
        {
            let cache = self.store.cache.borrow();
            if let Some(r) = cache.get(&self.key()).and_then(|c| c.refget.clone()) {
                return Ok(r);
            }
        }
        let refget = match self.store.hdp.get_refget_accession(&self.ac)? {
            Some(r) => r,
            None => crate::vrs::refget_accession(&self.whole()?),
        };
        let mut cache = self.store.cache.borrow_mut();
        cache.entry(self.key()).or_default().refget = Some(refget.clone());
        Ok(refget)
    }

    /// The complete sequence.
    pub fn whole(&self) -> Result<String, HgvsError> {
        {
            let cache = self.store.cache.borrow();
            if let Some(whole) = cache.get(&self.key()).and_then(|c| c.whole.as_ref()) {
                return Ok(whole.clone());
            }
        }
        let whole = self.store.fetch(&self.ac, self.kind, 0, None)?;
        let mut cache = self.store.cache.borrow_mut();
        let entry = cache.entry(self.key()).or_default();
        entry.len = Some(whole.len());
        entry.whole = Some(whole.clone());
        Ok(whole)
    }

    /// How many consecutive positions `from, from + 1, …` carry the base
    /// `pattern[k % pattern.len()]` at step `k`. Stops at the first mismatch or
    /// the end of the sequence. Zero for an empty pattern.
    ///
    /// This is the 3' shift rule: an edit over `[start, end)` with reference
    /// bases `pattern` can move right by exactly this many positions.
    pub fn run_right(&self, from: usize, pattern: &[u8]) -> Result<usize, HgvsError> {
        let n = pattern.len();
        if n == 0 {
            return Ok(0);
        }
        let mut k = 0;
        while let Some(b) = self.base(from + k)? {
            if b != pattern[k % n] {
                break;
            }
            k += 1;
        }
        Ok(k)
    }

    /// How many consecutive positions `from - 1, from - 2, …` carry the base
    /// `pattern[(n - 1 - k) mod n]` at step `k`. Stops at the first mismatch or
    /// position zero. Zero for an empty pattern.
    ///
    /// This is the 5' shift rule, the mirror of [`run_right`](Self::run_right).
    pub fn run_left(&self, from: usize, pattern: &[u8]) -> Result<usize, HgvsError> {
        let n = pattern.len();
        if n == 0 {
            return Ok(0);
        }
        let mut k = 0;
        while k < from {
            match self.base(from - 1 - k)? {
                Some(b) if b == pattern[n - 1 - (k % n)] => k += 1,
                _ => break,
            }
        }
        Ok(k)
    }

    fn key(&self) -> (String, IdentifierType) {
        (self.ac.clone(), self.kind)
    }
}

fn clamp_slice(seq: &str, start: usize, end: usize) -> String {
    let end = end.min(seq.len());
    if start >= end {
        String::new()
    } else {
        seq[start..end].to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::{IdentifierKind, TranscriptData};
    use std::cell::Cell;

    /// Serves one fixed sequence and counts provider calls.
    struct CountingProvider {
        seq: &'static str,
        calls: Cell<usize>,
    }

    impl DataProvider for CountingProvider {
        fn get_transcript(&self, _: &str, _: Option<&str>) -> Result<TranscriptData, HgvsError> {
            unreachable!()
        }
        fn get_seq(
            &self,
            _: &str,
            start: i32,
            end: Option<i32>,
            _: IdentifierType,
        ) -> Result<String, HgvsError> {
            self.calls.set(self.calls.get() + 1);
            let s = (start.max(0) as usize).min(self.seq.len());
            let e = end.map_or(self.seq.len(), |e| (e as usize).min(self.seq.len()));
            Ok(self.seq[s..e.max(s)].to_string())
        }
        fn get_symbol_accessions(
            &self,
            _: &str,
            _: IdentifierKind,
            _: IdentifierKind,
        ) -> Result<Vec<(IdentifierType, String)>, HgvsError> {
            Ok(vec![])
        }
        fn get_identifier_type(&self, _: &str) -> Result<IdentifierType, HgvsError> {
            Ok(IdentifierType::GenomicAccession)
        }
    }

    fn provider(seq: &'static str) -> CountingProvider {
        CountingProvider {
            seq,
            calls: Cell::new(0),
        }
    }

    const G: IdentifierType = IdentifierType::GenomicAccession;

    #[test]
    fn slice_spans_blocks_and_clamps_at_end() {
        let hdp = provider("ACGTACGTAC"); // 10 bases
        let store = ReferenceStore::with_block_size(&hdp, 4);
        let r = store.reference("X", G);
        assert_eq!(r.slice(2, 7).unwrap(), "GTACG");
        assert_eq!(r.slice(8, 20).unwrap(), "AC");
        assert_eq!(r.slice(10, 12).unwrap(), "");
        assert_eq!(r.slice(5, 5).unwrap(), "");
        assert_eq!(r.base(9).unwrap(), Some(b'C'));
        assert_eq!(r.base(10).unwrap(), None);
    }

    #[test]
    fn blocks_are_fetched_once() {
        let hdp = provider("ACGTACGTAC");
        let store = ReferenceStore::with_block_size(&hdp, 4);
        let r = store.reference("X", G);
        r.slice(0, 10).unwrap();
        assert_eq!(hdp.calls.get(), 3);
        r.slice(1, 9).unwrap();
        r.base(3).unwrap();
        r.base(11).unwrap(); // past the known end: no fetch
        assert_eq!(hdp.calls.get(), 3);
        // A different kind is a different sequence.
        store
            .reference("X", IdentifierType::TranscriptAccession)
            .base(0)
            .unwrap();
        assert_eq!(hdp.calls.get(), 4);
    }

    #[test]
    fn whole_is_fetched_once_and_serves_slices() {
        let hdp = provider("ACGTACGTAC");
        let store = ReferenceStore::with_block_size(&hdp, 4);
        let r = store.reference("X", G);
        assert_eq!(r.whole().unwrap(), "ACGTACGTAC");
        assert_eq!(r.whole().unwrap(), "ACGTACGTAC");
        assert_eq!(r.slice(3, 6).unwrap(), "TAC");
        assert_eq!(hdp.calls.get(), 1);
    }

    #[test]
    fn run_right_and_left_follow_a_cyclic_pattern() {
        //                 0123456789012
        let hdp = provider("TTCAGCAGCAGTT");
        let store = ReferenceStore::with_block_size(&hdp, 4);
        let r = store.reference("X", G);
        // A CAG deletion at [2,5) can shift right through the repeat to [8,11).
        assert_eq!(r.run_right(5, b"CAG").unwrap(), 6);
        // And a CAG deletion at [8,11) can shift left by the same amount.
        assert_eq!(r.run_left(8, b"CAG").unwrap(), 6);
        // The pattern cycles: from 2 the bases read CAGCAGCAG then T.
        assert_eq!(r.run_right(2, b"CAG").unwrap(), 9);
        // From 3 the bases read A, G, C: an AG insertion can move two places.
        assert_eq!(r.run_right(3, b"AG").unwrap(), 2);
        // Nothing to walk for an empty pattern or at the edges.
        assert_eq!(r.run_right(5, b"").unwrap(), 0);
        assert_eq!(r.run_left(0, b"T").unwrap(), 0);
        assert_eq!(r.run_right(13, b"T").unwrap(), 0);
    }
}
