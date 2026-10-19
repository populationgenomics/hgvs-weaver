//! Normalisation of a nucleotide edit against one reference sequence.
//!
//! HGVS wants an edit at its 3'-most equivalent position, with the bases it
//! deletes, duplicates or inverts left unstated, and an insertion that copies
//! the bases just before it written as a duplication. All of that is decided
//! here, on a 0-based half-open range over a [`Reference`]. The coordinate
//! systems (g., c., n.) only convert their positions in and out; see
//! [`VariantMapper::normalize_variant`](crate::mapper::VariantMapper::normalize_variant).

use crate::data::IdentifierType;
use crate::edits::{NaEdit, ResolvedEdit};
use crate::error::HgvsError;
use crate::reference::Reference;
use crate::structs::strip_common_prefix_suffix;

/// An edit placed on a reference by a 0-based half-open range.
///
/// For an insertion the range is empty: `start == end` is the index of the
/// base the inserted bases go in front of.
#[derive(Debug, Clone, PartialEq)]
pub struct PlacedEdit {
    pub start: usize,
    pub end: usize,
    pub edit: NaEdit,
}

impl PlacedEdit {
    /// Places an edit written over the HGVS range `[start, end)`.
    ///
    /// HGVS names an insertion by its two flanking bases; the edit itself is
    /// the empty range in front of the second one. Every other edit covers
    /// the bases it names.
    pub fn from_hgvs_range(start: usize, end: usize, edit: NaEdit) -> Self {
        if matches!(edit, NaEdit::Ins { .. }) && end > start {
            let anchor = end - 1;
            PlacedEdit {
                start: anchor,
                end: anchor,
                edit,
            }
        } else {
            PlacedEdit { start, end, edit }
        }
    }

    pub fn is_insertion(&self) -> bool {
        matches!(self.edit, NaEdit::Ins { .. })
    }

    /// The half-open range HGVS writes: the two flanking bases for an
    /// insertion, the edited bases for everything else.
    pub fn hgvs_range(&self) -> (usize, usize) {
        if self.is_insertion() {
            (self.start.saturating_sub(1), self.start + 1)
        } else {
            (self.start, self.end)
        }
    }
}

/// The exons of a spliced sequence, laid end to end in transcript order from
/// index 0. A sequence that is not spliced, a genome or a protein, is one
/// stretch with no end a shift has to respect. An index past the last exon is
/// in none.
///
/// On a transcript, the bases either side of a junction are neighbours, but on
/// the genome an intron lies between them: a change slid across a junction is
/// a different change on the genome. An exon aligned to the genome with gaps
/// is one stretch here, so a shift within it can still cross a gap.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ExonBounds {
    /// One past the last index of each exon, increasing.
    ends: Vec<usize>,
}

impl ExonBounds {
    /// A sequence that is not spliced.
    pub fn unspliced() -> Self {
        Self::default()
    }

    /// Exons of these lengths, in transcript order.
    pub fn from_lengths(lengths: impl IntoIterator<Item = usize>) -> Self {
        let ends = lengths
            .into_iter()
            .scan(0usize, |at, len| {
                *at += len;
                Some(*at)
            })
            .collect();
        ExonBounds { ends }
    }

    /// Which exon holds `index`, counting from 0; the exon count for an index
    /// past the last.
    pub(crate) fn exon(&self, index: usize) -> usize {
        self.ends.partition_point(|&e| e <= index)
    }

    /// One past the last index of the exon holding `index`.
    pub(crate) fn exon_end(&self, index: usize) -> usize {
        self.ends
            .get(self.exon(index))
            .copied()
            .unwrap_or(usize::MAX)
    }

    /// The first index of the exon holding `index`.
    pub(crate) fn exon_start(&self, index: usize) -> usize {
        match self.exon(index) {
            0 => 0,
            i => self.ends[i - 1],
        }
    }

    /// These exons with each index in `cuts` also ending one, so that no
    /// shift crosses it.
    pub(crate) fn with_cuts(&self, cuts: impl IntoIterator<Item = usize>) -> Self {
        let mut ends: Vec<usize> = self
            .ends
            .iter()
            .copied()
            .chain(cuts)
            .filter(|&c| c > 0)
            .collect();
        ends.sort_unstable();
        ends.dedup();
        ExonBounds { ends }
    }

    /// Whether an exon ends in front of `index`, so that it and the base
    /// before it are in different exons.
    pub(crate) fn is_junction(&self, index: usize) -> bool {
        index > 0 && !self.same_exon(index - 1, index)
    }

    fn same_exon(&self, a: usize, b: usize) -> bool {
        self.exon(a) == self.exon(b)
    }

    /// Whether an insertion in front of `index` has no single place on the
    /// genome, or none HGVS can write: before the first base, on a junction,
    /// or past the last exon.
    fn is_edge(&self, index: usize) -> bool {
        index == 0 || !self.same_exon(index - 1, index)
    }
}

/// Normalises `placed` against `reference`: reduces a delins to the edit it
/// is, shifts it as far 3' as the sequence allows and rewrites an insertion
/// as a duplication when the inserted bases repeat the bases immediately
/// before it.
///
/// A delins is read from the reference's bases over its range, as its allele
/// is, and loses what its inserted bases share with them at either end; what
/// remains is written as the substitution, deletion, insertion, duplication,
/// inversion or delins it is, and normalised as that. One that changes a base
/// at each end cannot slide and stays where it is; a substitution is one.
///
/// An edit that makes the reference, a delins of the bases already there, an
/// inversion of a range that is its own reverse complement or a repeat of the
/// copies already there, is written as the identity of the bases it names,
/// where those lie in one exon.
///
/// No end of the edit leaves the exon it was written in: the shift stops at
/// the last base of the exon, and an insertion is neither slid onto a junction
/// nor rewritten as a duplication of bases in two exons. One written on a
/// junction is left as written.
///
/// The bases an edit states as its reference are checked against the
/// reference and then not written: a deletion, duplication or inversion, the
/// deleted part of a delins and an identity come out bare, whether or not
/// they moved. HGVS writes them bare; a substitution keeps its base, which is
/// part of how it is written.
///
/// # Errors
///
/// [`HgvsError::PositionOutsideSequence`] when the edit's range runs past the
/// end of the reference, an insertion's included when it goes after the last
/// base; [`HgvsError::ReferenceMismatch`] when the edit states bases the
/// reference does not hold, whether or not the edit would move.
pub fn normalize(
    reference: &Reference<'_, '_>,
    placed: PlacedEdit,
    exons: &ExonBounds,
) -> Result<PlacedEdit, HgvsError> {
    let (hgvs_start, hgvs_end) = placed.hgvs_range();
    let PlacedEdit { start, end, edit } = placed;
    edit.check_against(reference, hgvs_start, hgvs_end, crate::data::Strand::Plus)?;
    let edit = without_stated_bases(edit);
    if let Some((start, end, uncertain)) = statement_made(reference, start, end, &edit, exons)? {
        return Ok(PlacedEdit {
            start,
            end,
            edit: NaEdit::RefAlt {
                ref_: None,
                alt: None,
                uncertain,
            },
        });
    }
    match &edit {
        NaEdit::RefAlt {
            alt: Some(alt),
            uncertain,
            ..
        } if end > start && !edit.is_identity() => {
            reduce_delins(reference, start, end, alt, *uncertain, exons)
        }
        NaEdit::Del { .. } | NaEdit::Dup { .. } if end > start => {
            let room = (exons.exon_end(start) - 1 - start).min(exons.exon_end(end - 1) - end);
            let k = shift_3(reference, &edit.resolve(reference, start, end)?)?.min(room);
            Ok(PlacedEdit {
                start: start + k,
                end: end + k,
                edit,
            })
        }
        NaEdit::Ins { .. } if start > 0 && !exons.same_exon(start - 1, start) => {
            Ok(PlacedEdit { start, end, edit })
        }
        NaEdit::Ins {
            alt: Some(seq),
            uncertain,
        } => {
            // The base before the insertion point may reach the last base of
            // its exon, which puts the point on the junction.
            let room = exons.exon_end(start.saturating_sub(1)) - start;
            let k = shift_3(reference, &edit.resolve(reference, start, end)?)?.min(room);
            place_insertion(reference, start, seq, *uncertain, k, exons)
        }
        _ => Ok(PlacedEdit { start, end, edit }),
    }
}

/// [`normalize`] for an edit placed by interbase coordinates, as SPDI and VRS
/// place one. These can put an insertion after the last base, which HGVS
/// cannot name: it is written as [`insertion_in_exon`] writes it.
pub(crate) fn normalize_interbase(
    reference: &Reference<'_, '_>,
    placed: PlacedEdit,
    exons: &ExonBounds,
) -> Result<PlacedEdit, HgvsError> {
    match &placed.edit {
        NaEdit::Ins {
            alt: Some(seq),
            uncertain,
        } if placed.start > 0 && reference.base(placed.start)?.is_none() => insertion_in_exon(
            reference,
            placed.start,
            seq,
            *uncertain,
            exons.exon(placed.start - 1),
            exons,
        ),
        _ => normalize(reference, placed, exons),
    }
}

/// The bases an inversion or a repeat names, with its uncertainty, when what
/// it makes of them is the reference and they lie in one exon: the inverted
/// range, or the whole run a repeat resolves to. Such an edit states those
/// bases unchanged, and is written as their identity. A delins that makes the
/// reference is reduced to one by [`reduce_delins`]; no other edit can make
/// it, and a repeat whose count is a range names no one sequence.
fn statement_made(
    reference: &Reference<'_, '_>,
    start: usize,
    end: usize,
    edit: &NaEdit,
    exons: &ExonBounds,
) -> Result<Option<(usize, usize, bool)>, HgvsError> {
    let uncertain = match edit {
        NaEdit::Inv { uncertain, .. } => *uncertain,
        NaEdit::Repeat {
            min,
            max,
            uncertain,
            ..
        } if min == max => *uncertain,
        _ => return Ok(None),
    };
    if end <= start {
        return Ok(None);
    }
    if matches!(edit, NaEdit::Inv { .. })
        && !may_be_its_own_reverse_complement(reference, start, end)?
    {
        return Ok(None);
    }
    let resolved = edit.resolve(reference, start, end)?;
    let named = (resolved.start, resolved.end);
    Ok((named.1 > named.0
        && resolved.alt.eq_ignore_ascii_case(&resolved.ref_)
        && exons.same_exon(named.0, named.1 - 1))
    .then_some((named.0, named.1, uncertain)))
}

/// Whether `[start, end)` can be its own reverse complement, read from its
/// length and its two end bases, so that a long inversion that is not is
/// refused without reading the bases between.
fn may_be_its_own_reverse_complement(
    reference: &Reference<'_, '_>,
    start: usize,
    end: usize,
) -> Result<bool, HgvsError> {
    if (end - start) % 2 == 1 {
        return Ok(false);
    }
    Ok(match (reference.base(start)?, reference.base(end - 1)?) {
        (Some(first), Some(last)) => {
            crate::utils::complement_dna_char(first.to_ascii_uppercase() as char)
                == last.to_ascii_uppercase() as char
        }
        _ => true,
    })
}

/// `edit` with the reference bases it states, or their number, dropped,
/// where they are not part of how the edit is written.
fn without_stated_bases(edit: NaEdit) -> NaEdit {
    match edit {
        NaEdit::Del { uncertain, .. } => NaEdit::Del {
            ref_: None,
            uncertain,
        },
        NaEdit::Dup { uncertain, .. } => NaEdit::Dup {
            ref_: None,
            uncertain,
        },
        NaEdit::Inv { uncertain, .. } => NaEdit::Inv {
            ref_: None,
            uncertain,
        },
        NaEdit::DelInsLength {
            min,
            max,
            uncertain,
            ..
        } => NaEdit::DelInsLength {
            ref_: None,
            min,
            max,
            uncertain,
        },
        NaEdit::RefAlt { uncertain, .. } if edit.is_identity() => NaEdit::RefAlt {
            ref_: None,
            alt: None,
            uncertain,
        },
        other => other,
    }
}

/// The delins over `[start, end)` that inserts `alt`, reduced and normalised.
///
/// On a transcript, a delins whose range spans a junction removes the intron
/// on the genome; it keeps at least one base in each exon, so that what it
/// reduces to removes the intron too, and is never an inversion, which would
/// turn the intron round instead.
fn reduce_delins(
    reference: &Reference<'_, '_>,
    start: usize,
    end: usize,
    alt: &str,
    uncertain: bool,
    exons: &ExonBounds,
) -> Result<PlacedEdit, HgvsError> {
    let ref_ = reference.slice(start, end)?;
    let spans = !exons.same_exon(start, end - 1);
    let (max_prefix, max_suffix) = if spans {
        (
            exons.exon_end(start) - 1 - start,
            end - 1 - exons.exon_start(end - 1),
        )
    } else {
        (usize::MAX, usize::MAX)
    };
    let (at, r, a) = trim_ends(start, &ref_, alt, max_prefix, max_suffix);
    if r.is_empty() && a.is_empty() {
        return Ok(PlacedEdit {
            start,
            end,
            edit: NaEdit::RefAlt {
                ref_: None,
                alt: None,
                uncertain,
            },
        });
    }
    if r.is_empty() {
        return insertion_in_exon(reference, at, &a, uncertain, exons.exon(start), exons);
    }
    let placed = replacement(reference, at, r, a, !spans, uncertain);
    if matches!(placed.edit, NaEdit::Del { .. }) {
        normalize(reference, placed, exons)
    } else {
        Ok(placed)
    }
}

/// The insertion of `seq` in front of `at`, in the exon numbered `home`, slid
/// 3' within it. Where it then sits on an edge of its exon, where an insertion
/// names no single place on the genome, or after the sequence's last base,
/// which HGVS cannot name, it moves one base back into the exon if the
/// sequence allows, as a slid insertion stops there, and is otherwise written
/// as a delins of the base beside it.
fn insertion_in_exon(
    reference: &Reference<'_, '_>,
    at: usize,
    seq: &str,
    uncertain: bool,
    home: usize,
    exons: &ExonBounds,
) -> Result<PlacedEdit, HgvsError> {
    let resolved = ResolvedEdit {
        start: at,
        end: at,
        ref_: String::new(),
        alt: seq.to_string(),
    };
    let at_home_start = exons.is_edge(at) && exons.exon(at) == home;
    let room = if at_home_start {
        exons.exon_end(at) - at
    } else {
        exons.exon_end(at - 1) - at
    };
    let k = shift_3(reference, &resolved)?.min(room);
    let placed = if k == 0 && at_home_start {
        // At the start of its exon: the bases before it are another exon's,
        // so it duplicates none of them.
        PlacedEdit {
            start: at,
            end: at,
            edit: NaEdit::Ins {
                alt: Some(seq.to_string()),
                uncertain,
            },
        }
    } else {
        place_insertion(reference, at, seq, uncertain, k, exons)?
    };
    let point = placed.start;
    let NaEdit::Ins {
        alt: Some(inserted),
        ..
    } = &placed.edit
    else {
        return Ok(placed);
    };
    if !exons.is_edge(point) && reference.base(point)?.is_some() {
        return Ok(placed);
    }
    let after_home = point > 0 && exons.exon(point - 1) == home;
    if after_home && !exons.is_edge(point - 1) {
        // One base 5', inside the exon, where an insertion slid from there
        // would also stop.
        let back = ResolvedEdit {
            start: point,
            end: point,
            ref_: String::new(),
            alt: inserted.clone(),
        };
        if shift_5(reference, &back)? > 0 {
            let n = inserted.len();
            return Ok(PlacedEdit {
                start: point - 1,
                end: point - 1,
                edit: NaEdit::Ins {
                    alt: Some(format!("{}{}", &inserted[n - 1..], &inserted[..n - 1])),
                    uncertain,
                },
            });
        }
    }
    let (s, alt) = if after_home {
        (
            point - 1,
            format!("{}{inserted}", reference.slice(point - 1, point)?),
        )
    } else {
        (
            point,
            format!("{inserted}{}", reference.slice(point, point + 1)?),
        )
    };
    Ok(PlacedEdit {
        start: s,
        end: s + 1,
        edit: NaEdit::RefAlt {
            ref_: Some(String::new()),
            alt: Some(alt),
            uncertain,
        },
    })
}

/// An insertion at `at` slid `k` positions 3': a duplication where the bases
/// before it repeat the insert and lie in one exon, else an insertion. `k` may
/// take the insertion point onto the end of its exon or of the sequence, where
/// only a duplication is kept; an insertion there would sit on a junction, or
/// past the last base, and stops one base short.
fn place_insertion(
    reference: &Reference<'_, '_>,
    at: usize,
    seq: &str,
    uncertain: bool,
    k: usize,
    exons: &ExonBounds,
) -> Result<PlacedEdit, HgvsError> {
    let n = seq.len();
    let point = at + k;
    // Sliding 3' by k rotates what is inserted: the k bases it passed over now
    // precede it, and the same k bases of the insert follow.
    let r = if n > 0 { k % n } else { 0 };
    let rotated = format!("{}{}", &seq[r..], &seq[..r]);
    if n > 0
        && point >= n
        && exons.same_exon(point - n, point - 1)
        && reference.slice(point - n, point)? == rotated
    {
        return Ok(PlacedEdit {
            start: point - n,
            end: point,
            edit: NaEdit::Dup {
                ref_: None,
                uncertain,
            },
        });
    }
    if k > 0 && (!exons.same_exon(at, point) || reference.base(point)?.is_none()) {
        return place_insertion(reference, at, seq, uncertain, k - 1, exons);
    }
    Ok(PlacedEdit {
        start: point,
        end: point,
        edit: NaEdit::Ins {
            alt: Some(rotated),
            uncertain,
        },
    })
}

/// The part of a resolved edit that can slide: after dropping the bases the
/// reference and alternate share, either a pure insertion (`[at, at)` plus the
/// inserted bases) or a pure deletion (`[at, at + n)` plus the deleted bases).
/// A substitution, an inversion or a delins with bases on both sides cannot
/// slide and yields `None`.
fn movable(resolved: &ResolvedEdit) -> Option<(usize, usize, String)> {
    if resolved.ref_ == resolved.alt {
        return None;
    }
    let (at, r, a) =
        strip_common_prefix_suffix(resolved.start as i32, &resolved.ref_, &resolved.alt);
    let at = at as usize;
    if r.is_empty() && !a.is_empty() {
        Some((at, at, a))
    } else if a.is_empty() && !r.is_empty() {
        Some((at, at + r.len(), r))
    } else {
        None
    }
}

/// The widest range over which `resolved` is ambiguous: every reference
/// position the change could equally well be written at, 5' and 3'. This is
/// what an unambiguous SPDI spells out. An edit that cannot slide is its own
/// range.
pub fn ambiguous_range(
    reference: &Reference<'_, '_>,
    resolved: &ResolvedEdit,
) -> Result<(usize, usize), HgvsError> {
    let Some((at_start, at_end, pattern)) = movable(resolved) else {
        return Ok((resolved.start, resolved.end));
    };
    let k5 = reference.run_left(at_start, pattern.as_bytes())?;
    let k3 = reference.run_right(at_end, pattern.as_bytes())?;
    Ok((
        resolved.start.min(at_start - k5),
        resolved.end.max(at_end + k3),
    ))
}

/// How many positions the edit can move 3' and still describe the same change,
/// measured at its movable insertion or deletion point.
pub fn shift_3(reference: &Reference<'_, '_>, resolved: &ResolvedEdit) -> Result<usize, HgvsError> {
    match movable(resolved) {
        Some((_, end, pattern)) => reference.run_right(end, pattern.as_bytes()),
        None => Ok(0),
    }
}

/// How many positions the edit can move 5' and still describe the same change,
/// measured at its movable insertion or deletion point. For a duplication that
/// point is the end of the duplicated bases, so the count includes them.
pub fn shift_5(reference: &Reference<'_, '_>, resolved: &ResolvedEdit) -> Result<usize, HgvsError> {
    match movable(resolved) {
        Some((start, _, pattern)) => reference.run_left(start, pattern.as_bytes()),
        None => Ok(0),
    }
}

/// `ref_` and `alt` with what they share at both ends removed, prefix first,
/// at most `max_prefix` and `max_suffix` bases, and the index the remainder
/// starts at. Prefix first so that a fully justified insertion or deletion
/// lands at the 3' end of its run, where HGVS writes it.
fn trim_ends(
    start: usize,
    ref_: &str,
    alt: &str,
    max_prefix: usize,
    max_suffix: usize,
) -> (usize, String, String) {
    let prefix = ref_
        .bytes()
        .zip(alt.bytes())
        .take_while(|(x, y)| x == y)
        .count()
        .min(max_prefix);
    let (r, a) = (&ref_[prefix..], &alt[prefix..]);
    let suffix = r
        .bytes()
        .rev()
        .zip(a.bytes().rev())
        .take_while(|(x, y)| x == y)
        .count()
        .min(max_suffix);
    (
        start + prefix,
        r[..r.len() - suffix].to_string(),
        a[..a.len() - suffix].to_string(),
    )
}

/// The HGVS edit, and the range it is written over before normalisation, for
/// "the bases over `[start, end)` (which are `ref_`) become `alt`": identity,
/// substitution, inversion, deletion, insertion or delins.
pub(crate) fn minimal_edit(
    reference: &Reference<'_, '_>,
    start: usize,
    end: usize,
    ref_: &str,
    alt: &str,
) -> Result<(NaEdit, usize, usize), HgvsError> {
    let (at, r, a) = trim_ends(start, ref_, alt, usize::MAX, usize::MAX);
    let uncertain = false;
    Ok(if r.is_empty() && a.is_empty() {
        let edit = NaEdit::RefAlt {
            ref_: None,
            alt: None,
            uncertain,
        };
        if end > start {
            (edit, start, end)
        } else {
            // Nothing over an empty range: say so of the base before it.
            let s = start.saturating_sub(1);
            (edit, s, s + 1)
        }
    } else if r.is_empty() && at == 0 {
        // HGVS has no insertion before the first base; it is written as a
        // delins of that base.
        let first = reference.slice(0, 1)?;
        let edit = NaEdit::RefAlt {
            ref_: Some(String::new()),
            alt: Some(format!("{a}{first}")),
            uncertain,
        };
        (edit, 0, 1)
    } else if r.is_empty() {
        let edit = NaEdit::Ins {
            alt: Some(a),
            uncertain,
        };
        (edit, at - 1, at + 1)
    } else {
        let placed = replacement(reference, at, r, a, true, uncertain);
        (placed.edit, placed.start, placed.end)
    })
}

/// The edit that replaces the bases `r` at `at` by `a`, the two sharing
/// neither their first nor their last base: a deletion, a substitution, an
/// inversion where `inversion` allows one and the sequence is nucleotides, or
/// a delins written without the deleted bases.
fn replacement(
    reference: &Reference<'_, '_>,
    at: usize,
    r: String,
    a: String,
    inversion: bool,
    uncertain: bool,
) -> PlacedEdit {
    let end = at + r.len();
    let edit = if a.is_empty() {
        NaEdit::Del {
            ref_: None,
            uncertain,
        }
    } else if r.len() == 1 && a.len() == 1 {
        NaEdit::RefAlt {
            ref_: Some(r),
            alt: Some(a),
            uncertain,
        }
    } else if inversion
        && reference.kind() != IdentifierType::ProteinAccession
        && r.len() > 1
        && a == crate::utils::reverse_complement(&r)
    {
        NaEdit::Inv {
            ref_: None,
            uncertain,
        }
    } else {
        NaEdit::RefAlt {
            ref_: Some(String::new()),
            alt: Some(a),
            uncertain,
        }
    };
    PlacedEdit {
        start: at,
        end,
        edit,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::{DataProvider, IdentifierKind, IdentifierType, TranscriptData};
    use crate::reference::ReferenceStore;

    struct Fixed(&'static str);

    impl DataProvider for Fixed {
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
            let s = (start.max(0) as usize).min(self.0.len());
            let e = end.map_or(self.0.len(), |e| (e as usize).min(self.0.len()));
            Ok(self.0[s..e.max(s)].to_string())
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

    fn del(ref_: Option<&str>) -> NaEdit {
        NaEdit::Del {
            ref_: ref_.map(str::to_string),
            uncertain: false,
        }
    }
    fn ins(alt: &str) -> NaEdit {
        NaEdit::Ins {
            alt: Some(alt.to_string()),
            uncertain: false,
        }
    }

    fn run(seq: &'static str, placed: PlacedEdit) -> PlacedEdit {
        let hdp = Fixed(seq);
        let store = ReferenceStore::with_block_size(&hdp, 4);
        normalize(
            &store.reference("X", IdentifierType::GenomicAccession),
            placed,
            &ExonBounds::unspliced(),
        )
        .unwrap()
    }

    #[test]
    fn insertion_places_at_the_second_flanking_base() {
        let p = PlacedEdit::from_hgvs_range(4, 6, ins("A"));
        assert_eq!((p.start, p.end), (5, 5));
        assert_eq!(p.hgvs_range(), (4, 6));
        let d = PlacedEdit::from_hgvs_range(4, 6, del(None));
        assert_eq!((d.start, d.end), (4, 6));
    }

    #[test]
    fn deletion_shifts_3_prime_and_states_no_bases() {
        //             0123456789
        let out = run("TTCAGCAGTT", PlacedEdit::from_hgvs_range(2, 5, del(None)));
        assert_eq!(
            out,
            PlacedEdit {
                start: 5,
                end: 8,
                edit: del(None)
            }
        );
    }

    #[test]
    fn insertion_that_repeats_preceding_bases_becomes_duplication() {
        // Insert CAG in front of index 5 (between the two CAGs): shifts to the
        // end of the run and duplicates the last CAG.
        let out = run("TTCAGCAGTT", PlacedEdit::from_hgvs_range(4, 6, ins("CAG")));
        assert_eq!(
            out,
            PlacedEdit {
                start: 5,
                end: 8,
                edit: NaEdit::Dup {
                    ref_: None,
                    uncertain: false
                }
            }
        );
        // A genuine insertion stays one.
        let out = run("TTCAGCAGTT", PlacedEdit::from_hgvs_range(1, 3, ins("GGG")));
        assert_eq!(
            out,
            PlacedEdit {
                start: 2,
                end: 2,
                edit: ins("GGG")
            }
        );
    }

    #[test]
    fn a_slid_insertion_is_rotated() {
        // Inserting GA before index 2 of GGGGGGGGG slides one base (the next
        // base is G) and becomes an insertion of AG before index 3: the same
        // molecule, GGGAGGGGGGG.
        let out = run("GGGGGGGGG", PlacedEdit::from_hgvs_range(1, 3, ins("GA")));
        assert_eq!(
            out,
            PlacedEdit {
                start: 3,
                end: 3,
                edit: ins("AG")
            }
        );
        // Sliding a whole number of periods leaves the insert as written; here
        // it then reads as a duplication of the last CAG copy, at [5, 8).
        let out = run("TTCAGCAGTT", PlacedEdit::from_hgvs_range(1, 3, ins("CAG")));
        assert_eq!((out.start, out.end), (5, 8));
        assert_eq!(
            out.edit,
            NaEdit::Dup {
                ref_: None,
                uncertain: false
            }
        );
    }

    #[test]
    fn substitution_is_left_alone() {
        let sub = NaEdit::RefAlt {
            ref_: Some("C".into()),
            alt: Some("T".into()),
            uncertain: false,
        };
        let out = run("TTCAGCAGTT", PlacedEdit::from_hgvs_range(2, 3, sub.clone()));
        assert_eq!(
            out,
            PlacedEdit {
                start: 2,
                end: 3,
                edit: sub
            }
        );
    }

    #[test]
    fn shift_5_mirrors_shift_3() {
        let hdp = Fixed("TTCAGCAGTT");
        let store = ReferenceStore::with_block_size(&hdp, 4);
        let r = store.reference("X", IdentifierType::GenomicAccession);
        let res = |e: NaEdit, s, t| e.resolve(&r, s, t).unwrap();
        assert_eq!(shift_3(&r, &res(del(None), 2, 5)).unwrap(), 3);
        assert_eq!(shift_5(&r, &res(del(None), 5, 8)).unwrap(), 3);
        assert_eq!(shift_3(&r, &res(ins("CAG"), 5, 5)).unwrap(), 3);
        assert_eq!(shift_5(&r, &res(ins("CAG"), 5, 5)).unwrap(), 3);
        // A duplication slides like an insertion of its bases at its end, so
        // its 5' count includes the duplicated bases themselves.
        let dup = NaEdit::Dup {
            ref_: None,
            uncertain: false,
        };
        assert_eq!(shift_3(&r, &res(dup.clone(), 2, 5)).unwrap(), 3);
        assert_eq!(shift_5(&r, &res(dup.clone(), 5, 8)).unwrap(), 6);
        // The ambiguous range covers the whole CAG run either way.
        assert_eq!(ambiguous_range(&r, &res(dup, 5, 8)).unwrap(), (2, 8));
        assert_eq!(ambiguous_range(&r, &res(del(None), 2, 5)).unwrap(), (2, 8));
        assert_eq!(ambiguous_range(&r, &res(ins("CAG"), 5, 5)).unwrap(), (2, 8));
        // A substitution is not ambiguous.
        let sub = NaEdit::RefAlt {
            ref_: Some("C".into()),
            alt: Some("T".into()),
            uncertain: false,
        };
        assert_eq!(ambiguous_range(&r, &res(sub, 2, 3)).unwrap(), (2, 3));
        // A delins that is really an insertion (delCinsCA at index 2, before A)
        // slides like the insertion of A it is.
        let delins = NaEdit::RefAlt {
            ref_: Some("C".into()),
            alt: Some("CA".into()),
            uncertain: false,
        };
        assert_eq!(shift_3(&r, &res(delins, 2, 3)).unwrap(), 1);
    }
}
