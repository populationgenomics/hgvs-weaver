use crate::cigar::Projected;
use crate::data::{ExonData, Strand, TranscriptData};
use crate::error::HgvsError;
use crate::structs::{
    Anchor, BaseOffsetInterval, BaseOffsetPosition, GenomicPos, IntronicOffset, TranscriptPos,
};

/// What a transcript index occupies on the genome.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Footprint {
    /// The half-open genomic range: one base, or none for a base the genome
    /// lacks, whose empty range `[at, at)` sits between genome bases `at - 1`
    /// and `at`.
    pub range: (GenomicPos, GenomicPos),
    /// For a base the genome lacks, the whole run of such bases it belongs
    /// to, as half-open transcript indices. No position inside the run exists
    /// on the genome, so the run is the unit the genome can describe.
    pub gap_run: Option<(TranscriptPos, TranscriptPos)>,
}

/// What a genomic position occupies on the transcript.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GenomicFootprint {
    /// The half-open transcript index range: one base, or none for a base the
    /// transcript lacks, whose empty range `[at, at)` sits between transcript
    /// indices `at - 1` and `at`. An intronic position has the exon base its
    /// offset counts from.
    pub range: (TranscriptPos, TranscriptPos),
    /// For a base the transcript lacks, the whole run of such genome bases it
    /// belongs to, as a half-open genomic range, low to high. No index inside
    /// the run exists on the transcript, so the run is the unit the transcript
    /// can describe.
    pub gap_run: Option<(GenomicPos, GenomicPos)>,
    /// Zero in an exon; otherwise the signed intronic offset from `range`'s base.
    pub offset: IntronicOffset,
}

/// Handles coordinate transformations within a single transcript.
pub struct TranscriptMapper {
    /// The transcript model providing exon and CDS information.
    pub transcript: TranscriptData,
    /// Sorted exons (transcript order).
    pub exons: Vec<ExonData>,
    /// CIGAR mappers for exons with non-trivial alignments.
    pub cigar_mappers: Vec<Option<crate::cigar::CigarMapper>>,
}

impl TranscriptMapper {
    /// Creates a new `TranscriptMapper` for the given transcript.
    pub fn new(transcript: TranscriptData) -> Result<Self, HgvsError> {
        let mut exons = transcript.exons.to_vec();
        if transcript.strand == crate::data::Strand::Plus {
            exons.sort_by_key(|e| e.reference_start.0);
        } else {
            exons.sort_by_key(|e| std::cmp::Reverse(e.reference_start.0));
        }
        let mut cigar_mappers = Vec::with_capacity(exons.len());
        for exon in &exons {
            if !exon.cigar.is_empty()
                && exon.cigar != format!("{}M", (exon.reference_end.0 - exon.reference_start.0) + 1)
                && exon.cigar != format!("{}=", (exon.reference_end.0 - exon.reference_start.0) + 1)
            {
                cigar_mappers.push(Some(crate::cigar::CigarMapper::new(&exon.cigar)?));
            } else {
                cigar_mappers.push(None);
            }
        }
        Ok(TranscriptMapper {
            transcript,
            exons,
            cigar_mappers,
        })
    }

    /// Maps a 0-based genomic position to a 0-based transcript position and intronic offset.
    /// An exonic base the transcript lacks has no index and is an error; see
    /// [`Self::g_to_n_footprint`].
    pub fn g_to_n(&self, g_pos: GenomicPos) -> Result<(TranscriptPos, IntronicOffset), HgvsError> {
        let GenomicFootprint {
            range: (lo, hi),
            offset,
            ..
        } = self.g_to_n_footprint(g_pos)?;
        if lo == hi {
            return Err(HgvsError::ValidationError(format!(
                "{}:g.{} is a base {} lacks; it has no transcript index",
                self.transcript.reference_accession,
                g_pos.0 + 1,
                self.transcript.ac,
            )));
        }
        Ok((lo, offset))
    }

    /// What a 0-based genomic position occupies on the transcript: one index
    /// for a base aligned to the transcript, none for a base the transcript
    /// lacks (a `D` in its exon's cigar), together with the run of such bases
    /// it belongs to; an intronic position is its nearest exon base and offset.
    ///
    /// A base the transcript lacks is placed only between two flanking
    /// transcript bases of its own exon. One at either end of the exon's
    /// alignment has no such pair and is an error.
    pub fn g_to_n_footprint(&self, g_pos: GenomicPos) -> Result<GenomicFootprint, HgvsError> {
        let mut n_pos = 0;
        for (i, exon) in self.exons.iter().enumerate() {
            let (e_start, e_end) = (exon.reference_start.0, exon.reference_end.0);
            // e_start and e_end are 0-based inclusive
            if g_pos.0 < e_start || g_pos.0 > e_end {
                n_pos += self.cigar_mappers[i]
                    .as_ref()
                    .map_or(e_end - e_start + 1, |cm| cm.tgt_len());
                continue;
            }
            let minus = exon.alt_strand == Strand::Minus;
            let g_offset = if minus {
                e_end - g_pos.0
            } else {
                g_pos.0 - e_start
            };
            let exonic = |lo: i32, hi: i32, gap_run| GenomicFootprint {
                range: (TranscriptPos(n_pos + lo), TranscriptPos(n_pos + hi)),
                gap_run,
                offset: IntronicOffset(0),
            };
            let Some(cm) = &self.cigar_mappers[i] else {
                return Ok(exonic(g_offset, g_offset + 1, None));
            };
            return match cm.map_ref_to_tgt(g_offset, true)? {
                Projected::Base(t) => Ok(exonic(t, t + 1, None)),
                Projected::Gap { at, run } if at > 0 && at < cm.tgt_len() => {
                    // The run is in alignment offsets; as genome positions it
                    // runs the other way on the minus strand.
                    let run = if minus {
                        (GenomicPos(e_end + 1 - run.1), GenomicPos(e_end + 1 - run.0))
                    } else {
                        (GenomicPos(e_start + run.0), GenomicPos(e_start + run.1))
                    };
                    Ok(exonic(at, at, Some(run)))
                }
                Projected::Gap { .. } => Err(HgvsError::ValidationError(format!(
                    "{}:g.{} is a base {} lacks at the edge of an exon's alignment; \
                     nothing places it on the transcript",
                    self.transcript.reference_accession,
                    g_pos.0 + 1,
                    self.transcript.ac,
                ))),
                Projected::Intronic { .. } => Err(HgvsError::CigarError(format!(
                    "exon cigar {} puts genome position {} in an N op",
                    exon.cigar, g_pos.0
                ))),
            };
        }

        // Handle intronic positions by finding the nearest exon
        let mut best_n = TranscriptPos(0);
        let mut best_dist = i32::MAX;
        let mut best_offset = 0;
        let mut curr_n = 0;
        for (i, exon) in self.exons.iter().enumerate() {
            let (e_start, e_end) = (exon.reference_start, exon.reference_end);
            let d_start = (g_pos.0 - e_start.0).abs();
            let d_end = (g_pos.0 - e_end.0).abs();
            let d = d_start.min(d_end);
            if d < best_dist {
                best_dist = d;
                let minus = exon.alt_strand == crate::data::Strand::Minus;
                let below = g_pos.0 < e_start.0;
                // The offset is signed in transcript direction: negative before
                // the exon's first base, positive after its last. On the minus
                // strand the genome runs the other way, so a position below
                // the exon's genomic start lies after its last transcript base.
                best_offset = match (below, minus) {
                    (true, false) => g_pos.0 - e_start.0,
                    (true, true) => e_start.0 - g_pos.0,
                    (false, false) => g_pos.0 - e_end.0,
                    (false, true) => e_end.0 - g_pos.0,
                };
                let before_first_base = below != minus;
                best_n = if before_first_base {
                    TranscriptPos(curr_n)
                } else {
                    let e_tgt_len = if let Some(cm) = &self.cigar_mappers[i] {
                        cm.tgt_len()
                    } else {
                        (e_end.0 - e_start.0) + 1
                    };
                    TranscriptPos(curr_n + e_tgt_len - 1)
                };
            }
            curr_n += if let Some(cm) = &self.cigar_mappers[i] {
                cm.tgt_len()
            } else {
                (e_end.0 - e_start.0) + 1
            };
        }
        Ok(GenomicFootprint {
            range: (best_n, TranscriptPos(best_n.0 + 1)),
            gap_run: None,
            offset: IntronicOffset(best_offset),
        })
    }

    /// Resolves a c./n. position to a 0-based transcript index.
    ///
    /// The anchor (transcript start, CDS start, CDS end) is applied here, so callers
    /// never do CDS arithmetic themselves. An intronic offset is rejected: an
    /// intronic base has no transcript index.
    pub fn position_to_n(&self, pos: &BaseOffsetPosition) -> Result<TranscriptPos, HgvsError> {
        if pos.offset.is_some_and(|o| o.0 != 0) {
            return Err(HgvsError::UnsupportedOperation(
                "Intronic position has no transcript index".into(),
            ));
        }
        self.c_to_n(pos.base.to_index(), pos.anchor)
    }

    /// Resolves a c./n. position, including any intronic offset, to a 0-based
    /// genomic position on the transcript's reference. Strand is handled here.
    /// A base the genome lacks has no position and is an error; see
    /// [`Self::position_to_g_range`].
    pub fn position_to_g(&self, pos: &BaseOffsetPosition) -> Result<GenomicPos, HgvsError> {
        let n = self.c_to_n(pos.base.to_index(), pos.anchor)?;
        self.n_to_g(n, pos.offset.unwrap_or(IntronicOffset(0)))
    }

    /// Resolves a c./n. position, including any intronic offset, to the
    /// half-open 0-based genomic range it occupies: one base, or none for an
    /// exonic base the genome lacks (see [`Self::n_to_g_range`]).
    pub fn position_to_g_range(
        &self,
        pos: &BaseOffsetPosition,
    ) -> Result<(GenomicPos, GenomicPos), HgvsError> {
        Ok(self.position_footprint(pos)?.range)
    }

    /// What a c./n. position, including any intronic offset, occupies on the
    /// genome. An intronic base is always one genome base.
    pub fn position_footprint(&self, pos: &BaseOffsetPosition) -> Result<Footprint, HgvsError> {
        let n = self.c_to_n(pos.base.to_index(), pos.anchor)?;
        match pos.offset {
            Some(offset) if offset.0 != 0 => {
                let g = self.n_to_g(n, offset)?;
                let end = g.0.checked_add(1).ok_or_else(|| {
                    HgvsError::ValidationError("Genomic end position overflow".into())
                })?;
                Ok(Footprint {
                    range: (g, GenomicPos(end)),
                    gap_run: None,
                })
            }
            _ => self.n_to_g_footprint(n),
        }
    }

    /// Resolves a c./n. interval to a half-open 0-based transcript index range
    /// `[start, end)`. A single position yields a range of length one.
    pub fn interval_to_n(
        &self,
        interval: &BaseOffsetInterval,
    ) -> Result<(TranscriptPos, TranscriptPos), HgvsError> {
        let start = self.position_to_n(&interval.start)?;
        let last = match &interval.end {
            Some(e) => self.position_to_n(e)?,
            None => start,
        };
        if last.0 < start.0 {
            return Err(HgvsError::ValidationError(format!(
                "Transcript range runs backwards: {} is after {}",
                interval.start,
                interval.end.as_ref().unwrap_or(&interval.start)
            )));
        }
        let end = last
            .0
            .checked_add(1)
            .ok_or_else(|| HgvsError::ValidationError("Transcript end position overflow".into()))?;
        Ok((start, TranscriptPos(end)))
    }

    /// Resolves a c./n. interval to a half-open 0-based genomic range `[start, end)`
    /// on the transcript's reference, ordered low-to-high regardless of strand.
    /// The range holds the genome bases the interval's bases align to; it is
    /// empty when every base in the interval is one the genome lacks.
    pub fn interval_to_g(
        &self,
        interval: &BaseOffsetInterval,
    ) -> Result<(GenomicPos, GenomicPos), HgvsError> {
        let (a_lo, a_hi) = self.position_to_g_range(&interval.start)?;
        let (b_lo, b_hi) = match &interval.end {
            Some(e) => self.position_to_g_range(e)?,
            None => (a_lo, a_hi),
        };
        Ok((
            GenomicPos(a_lo.0.min(b_lo.0)),
            GenomicPos(a_hi.0.max(b_hi.0)),
        ))
    }

    /// Maps a 0-based transcript position to a 0-based cDNA position and anchor.
    pub fn n_to_c(
        &self,
        n_pos: TranscriptPos,
    ) -> Result<(TranscriptPos, IntronicOffset, Anchor), HgvsError> {
        if let (Some(cds_start), Some(cds_end)) = (
            self.transcript.cds_start_index,
            self.transcript.cds_end_index,
        ) {
            if n_pos < cds_start {
                Ok((
                    TranscriptPos(n_pos.0 - cds_start.0),
                    IntronicOffset(0),
                    Anchor::CdsStart,
                ))
            } else if n_pos > cds_end {
                Ok((
                    TranscriptPos(n_pos.0 - cds_end.0 - 1),
                    IntronicOffset(0),
                    Anchor::CdsEnd,
                ))
            } else {
                Ok((
                    TranscriptPos(n_pos.0 - cds_start.0),
                    IntronicOffset(0),
                    Anchor::CdsStart,
                ))
            }
        } else {
            Ok((n_pos, IntronicOffset(0), Anchor::TranscriptStart))
        }
    }

    /// Maps a cDNA position and anchor to a 0-based transcript position.
    pub fn c_to_n(&self, c_pos: TranscriptPos, anchor: Anchor) -> Result<TranscriptPos, HgvsError> {
        match anchor {
            Anchor::TranscriptStart => Ok(c_pos),
            Anchor::CdsStart => {
                let cds_start = self
                    .transcript
                    .cds_start_index
                    .ok_or_else(|| HgvsError::ValidationError("Missing CDS start".into()))?;
                Ok(TranscriptPos(cds_start.0 + c_pos.0))
            }
            Anchor::CdsEnd => {
                let cds_end = self
                    .transcript
                    .cds_end_index
                    .ok_or_else(|| HgvsError::ValidationError("Missing CDS end".into()))?;
                Ok(TranscriptPos(cds_end.0 + 1 + c_pos.0))
            }
        }
    }

    /// The half-open 0-based genomic range a transcript index occupies on the
    /// reference: one base for a base aligned to the genome, none for a base
    /// the genome lacks (an `I` in its exon's cigar), whose range `[at, at)`
    /// sits between genome bases `at - 1` and `at`.
    ///
    /// A base the genome lacks is placed only between two flanking genome
    /// bases of its own exon. One at either end of the exon's alignment (a
    /// soft-clipped end, supplied as a leading or trailing `I`) has no such
    /// pair and is an error.
    pub fn n_to_g_range(
        &self,
        n_pos: TranscriptPos,
    ) -> Result<(GenomicPos, GenomicPos), HgvsError> {
        Ok(self.n_to_g_footprint(n_pos)?.range)
    }

    /// What a transcript index occupies on the genome: its range, and for a
    /// base the genome lacks the run it belongs to. See [`Self::n_to_g_range`].
    pub fn n_to_g_footprint(&self, n_pos: TranscriptPos) -> Result<Footprint, HgvsError> {
        let (_, footprint) = self.place(n_pos)?;
        Ok(footprint)
    }

    /// Maps a 0-based transcript position and offset to a 0-based genomic position.
    /// An exonic base the genome lacks has none and is an error.
    pub fn n_to_g(
        &self,
        n_pos: TranscriptPos,
        offset: IntronicOffset,
    ) -> Result<GenomicPos, HgvsError> {
        let (
            exon,
            Footprint {
                range: (lo, hi), ..
            },
        ) = self.place(n_pos)?;
        if lo == hi {
            return Err(HgvsError::ValidationError(format!(
                "{}:n.{} is a base {} lacks; it has no genomic position",
                self.transcript.ac,
                n_pos.0 + 1,
                self.transcript.reference_accession,
            )));
        }
        Ok(GenomicPos(match exon.alt_strand {
            Strand::Plus => lo.0 + offset.0,
            Strand::Minus => lo.0 - offset.0,
        }))
    }

    /// The exon holding a transcript index, and what the index occupies on
    /// the reference.
    fn place(&self, n_pos: TranscriptPos) -> Result<(&ExonData, Footprint), HgvsError> {
        let mut curr_n = 0;
        for (i, exon) in self.exons.iter().enumerate() {
            let (e_start, e_end) = (exon.reference_start.0, exon.reference_end.0);
            let cm = self.cigar_mappers[i].as_ref();
            let e_tgt_len = cm.map_or(e_end - e_start + 1, |cm| cm.tgt_len());
            if n_pos.0 < curr_n || n_pos.0 >= curr_n + e_tgt_len {
                curr_n += e_tgt_len;
                continue;
            }
            let tgt_offset = n_pos.0 - curr_n;
            // Half-open, as offsets along the exon's alignment.
            let mut gap_run = None;
            let (a, b) = match cm {
                None => (tgt_offset, tgt_offset + 1),
                Some(cm) => match cm.map_tgt_to_ref(tgt_offset, true)? {
                    Projected::Base(r) => (r, r + 1),
                    Projected::Gap { at, run } if at > 0 && at < cm.ref_len() => {
                        gap_run =
                            Some((TranscriptPos(curr_n + run.0), TranscriptPos(curr_n + run.1)));
                        (at, at)
                    }
                    Projected::Gap { .. } => {
                        return Err(HgvsError::ValidationError(format!(
                            "{}:n.{} is a base {} lacks at the edge of an exon's alignment; \
                             nothing places it on the genome",
                            self.transcript.ac,
                            n_pos.0 + 1,
                            self.transcript.reference_accession,
                        )))
                    }
                    Projected::Intronic { .. } => {
                        return Err(HgvsError::CigarError(format!(
                            "exon cigar {} puts transcript index {} in an N op",
                            exon.cigar, n_pos.0
                        )))
                    }
                },
            };
            let range = match exon.alt_strand {
                Strand::Plus => (GenomicPos(e_start + a), GenomicPos(e_start + b)),
                Strand::Minus => (GenomicPos(e_end + 1 - b), GenomicPos(e_end + 1 - a)),
            };
            return Ok((exon, Footprint { range, gap_run }));
        }
        Err(HgvsError::ValidationError(
            "Transcript position out of exon bounds".into(),
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::{ExonData, TranscriptData};

    fn create_mock_transcript(strand: crate::data::Strand, exons: Vec<ExonData>) -> TranscriptData {
        TranscriptData {
            ac: "NM_0001.1".to_string(),
            gene: "TEST".to_string(),
            cds_start_index: None,
            cds_end_index: None,
            strand,
            reference_accession: "NC_000001.1".to_string(),
            exons,
        }
    }

    #[test]
    fn test_g_to_n_minus_strand_order() {
        // Exons in genomic ascending order: [1000, 1100], [2000, 2100]
        // For minus strand, transcript order should be [2000, 2100] then [1000, 1100]
        let exons = vec![
            ExonData {
                transcript_start: TranscriptPos(0), // Placeholder
                transcript_end: TranscriptPos(100),
                reference_start: GenomicPos(1000),
                reference_end: GenomicPos(1100),
                alt_strand: crate::data::Strand::Minus,
                cigar: "101M".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(101),
                transcript_end: TranscriptPos(201),
                reference_start: GenomicPos(2000),
                reference_end: GenomicPos(2100),
                alt_strand: crate::data::Strand::Minus,
                cigar: "101M".to_string(),
            },
        ];
        let tx = create_mock_transcript(crate::data::Strand::Minus, exons);
        let mapper = TranscriptMapper::new(tx).unwrap();

        // Genomic 2100 should be n.0
        let (n_pos, offset) = mapper.g_to_n(GenomicPos(2100)).unwrap();
        // CURRENT BEHAVIOR (BUGGY):
        // It visits [1000, 1100] first, n_pos starts at 0.
        // Then it visits [2000, 2100], n_pos = 101.
        // offset_in_exon = 2100 - 2100 = 0 (for alt_strand = -1)
        // Result: n.101.
        // EXPECTED: n.0
        assert_eq!(n_pos.0, 0, "Minus strand genomic 2100 should be n.0");
        assert_eq!(offset.0, 0);
    }

    #[test]
    fn test_g_to_n_intronic_offset() {
        let exons = vec![ExonData {
            transcript_start: TranscriptPos(0),
            transcript_end: TranscriptPos(10),
            reference_start: GenomicPos(1000),
            reference_end: GenomicPos(1010),
            alt_strand: crate::data::Strand::Plus,
            cigar: "11M".to_string(),
        }];
        let tx = create_mock_transcript(crate::data::Strand::Plus, exons);
        let mapper = TranscriptMapper::new(tx).unwrap();

        // Genomic 999 is 1bp upstream of exon start (1000)
        let (n_pos, offset) = mapper.g_to_n(GenomicPos(999)).unwrap();
        assert_eq!(n_pos.0, 0);
        assert_eq!(offset.0, -1);

        // Genomic 1011 is 1bp downstream of exon end (1010)
        let (n_pos, offset) = mapper.g_to_n(GenomicPos(1011)).unwrap();
        assert_eq!(n_pos.0, 10);
        assert_eq!(offset.0, 1);
    }

    #[test]
    fn test_g_to_n_cigar() {
        // Exon: Genomic [1000, 1010] (11bp)
        // CIGAR: 5=1D5= (Ref 11bp -> Tgt 10bp)
        // Ref index: 0 1 2 3 4 5 6 7 8 9 10
        // Tgt index: 0 1 2 3 4 _ 5 6 7 8 9
        let exons = vec![ExonData {
            transcript_start: TranscriptPos(0),
            transcript_end: TranscriptPos(9),
            reference_start: GenomicPos(1000),
            reference_end: GenomicPos(1010),
            alt_strand: crate::data::Strand::Plus,
            cigar: "5=1D5=".to_string(),
        }];
        let tx = create_mock_transcript(crate::data::Strand::Plus, exons);
        let mapper = TranscriptMapper::new(tx).unwrap();

        // g.1000 -> n.0
        assert_eq!(mapper.g_to_n(GenomicPos(1000)).unwrap().0 .0, 0);
        // g.1004 -> n.4
        assert_eq!(mapper.g_to_n(GenomicPos(1004)).unwrap().0 .0, 4);
        // g.1005 is a base the transcript lacks (the 1D): no index, an empty
        // range between n.4 and n.5, and a run of one genome base.
        let err = mapper.g_to_n(GenomicPos(1005)).unwrap_err();
        assert!(
            err.to_string().contains("g.1006 is a base NM_0001.1 lacks"),
            "{err}"
        );
        assert_eq!(
            mapper.g_to_n_footprint(GenomicPos(1005)).unwrap(),
            GenomicFootprint {
                range: (TranscriptPos(5), TranscriptPos(5)),
                gap_run: Some((GenomicPos(1005), GenomicPos(1006))),
                offset: IntronicOffset(0),
            }
        );
        // g.1006 -> n.5
        assert_eq!(mapper.g_to_n(GenomicPos(1006)).unwrap().0 .0, 5);
    }

    /// One exon on genome 10..=39 with a transcript base the genome lacks
    /// after transcript index 14: `15=1I15=`.
    fn inserted_base_exon(strand: crate::data::Strand) -> TranscriptData {
        create_mock_transcript(
            strand,
            vec![ExonData {
                transcript_start: TranscriptPos(0),
                transcript_end: TranscriptPos(31),
                reference_start: GenomicPos(10),
                reference_end: GenomicPos(39),
                alt_strand: strand,
                cigar: "15=1I15=".to_string(),
            }],
        )
    }

    #[test]
    fn a_base_the_genome_lacks_occupies_an_empty_range_between_its_neighbours() {
        let mapper = TranscriptMapper::new(inserted_base_exon(crate::data::Strand::Plus)).unwrap();
        let range = |n: i32| mapper.n_to_g_range(TranscriptPos(n)).unwrap();
        assert_eq!(range(14), (GenomicPos(24), GenomicPos(25)));
        assert_eq!(range(15), (GenomicPos(25), GenomicPos(25)));
        assert_eq!(range(16), (GenomicPos(25), GenomicPos(26)));
        assert_eq!(range(30), (GenomicPos(39), GenomicPos(40)));
        // The base belongs to a run of one; an aligned base to none.
        let run = |n: i32| mapper.n_to_g_footprint(TranscriptPos(n)).unwrap().gap_run;
        assert_eq!(run(14), None);
        assert_eq!(run(15), Some((TranscriptPos(15), TranscriptPos(16))));
        // A single position is not enough for it.
        assert_eq!(
            mapper.n_to_g(TranscriptPos(14), IntronicOffset(0)).unwrap(),
            GenomicPos(24)
        );
        let err = mapper
            .n_to_g(TranscriptPos(15), IntronicOffset(0))
            .unwrap_err();
        assert!(
            err.to_string().contains("n.16 is a base NC_000001.1 lacks"),
            "{err}"
        );
        assert_eq!(
            mapper.n_to_g(TranscriptPos(16), IntronicOffset(0)).unwrap(),
            GenomicPos(25)
        );
    }

    #[test]
    fn a_base_the_genome_lacks_on_the_minus_strand() {
        let mapper = TranscriptMapper::new(inserted_base_exon(crate::data::Strand::Minus)).unwrap();
        let range = |n: i32| mapper.n_to_g_range(TranscriptPos(n)).unwrap();
        // Transcript index 0 is genome 39; index 14 is genome 25; the gap is
        // between genome 24 and 25; index 16 is genome 24.
        assert_eq!(range(0), (GenomicPos(39), GenomicPos(40)));
        assert_eq!(range(14), (GenomicPos(25), GenomicPos(26)));
        assert_eq!(range(15), (GenomicPos(25), GenomicPos(25)));
        assert_eq!(range(16), (GenomicPos(24), GenomicPos(25)));
        assert_eq!(range(30), (GenomicPos(10), GenomicPos(11)));
    }

    #[test]
    fn a_base_the_genome_lacks_at_the_edge_of_an_exon_is_an_error() {
        // A soft-clipped start, supplied as a leading I, and a clipped end.
        for (cigar, clipped, aligned) in [("5I20=", [0, 4], 5), ("20=5I", [20, 24], 19)] {
            let tx = create_mock_transcript(
                crate::data::Strand::Plus,
                vec![ExonData {
                    transcript_start: TranscriptPos(0),
                    transcript_end: TranscriptPos(25),
                    reference_start: GenomicPos(30),
                    reference_end: GenomicPos(49),
                    alt_strand: crate::data::Strand::Plus,
                    cigar: cigar.to_string(),
                }],
            );
            let mapper = TranscriptMapper::new(tx).unwrap();
            for n in clipped {
                let err = mapper.n_to_g_range(TranscriptPos(n)).unwrap_err();
                assert!(
                    err.to_string().contains("edge of an exon"),
                    "{cigar} n index {n}: {err}"
                );
            }
            assert!(mapper.n_to_g_range(TranscriptPos(aligned)).is_ok());
        }
    }
}
