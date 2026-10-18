//! Reading an allele in cis, `c.[a;b]`, into the one sequence its members
//! make, and what it states unchanged.
//!
//! Every member is numbered against the reference, so the allele is the
//! reference with each member's replacement applied, in any order. Members
//! that reach for the same reference base are read every way the description
//! allows: with no reading the allele is contradictory, a
//! [`HgvsError::ValidationError`] naming two members; with readings that make
//! different sequences, or that state different bases, it is
//! [`HgvsError::Ambiguous`], carrying each reading as a description; with one
//! it is a [`CisReading`].
//!
//! A reading is the reference plus the stretches the members replace, each
//! replaced by [`Segment`]s that are either bases the description writes or
//! pointers into the reference: a duplication's copy, an inversion's reverse
//! complement. The sequence itself is never built. Two readings are compared
//! by walking their segments over the stretch from the first base either
//! touches to the last, where two pointers to the same reference bases are
//! equal without reading them.

use std::collections::BTreeMap;

use crate::allele::StatedBases;
use crate::data::{IdentifierType, Strand};
use crate::edits::{NaEdit, ResolvedEdit};
use crate::error::HgvsError;
use crate::mapper::VariantMapper;
use crate::normalize::{self, ExonBounds, PlacedEdit};
use crate::reference::Reference;
use crate::structs::{
    BaseOffsetInterval, CVariant, CisPhasedVariant, GVariant, LinearVariant, MVariant, NVariant,
    PosEdit, SimpleInterval, SimplePosition, TranscriptVariant, Variant,
};
use crate::transcript_mapper::TranscriptMapper;
use crate::utils::reverse_complement;
use crate::SequenceVariant as SV;

mod normalise;

use normalise::Writing;

/// Where the bases of a stretch of an allele's sequence come from. Only a
/// reading builds one, so a range is never backwards.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Segment {
    /// The reference's bases over the 0-based half-open `[start, end)`.
    #[non_exhaustive]
    Reference { start: usize, end: usize },
    /// The reverse complement of the reference's bases over `[start, end)`.
    #[non_exhaustive]
    ReverseComplement { start: usize, end: usize },
    /// Bases the description writes, in upper case.
    #[non_exhaustive]
    Literal(String),
}

impl Segment {
    /// The number of bases the segment stands for.
    pub fn len(&self) -> usize {
        match self {
            Segment::Reference { start, end } | Segment::ReverseComplement { start, end } => {
                end - start
            }
            Segment::Literal(s) => s.len(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The `n` bases of the segment from its `from`th.
    ///
    /// # Errors
    ///
    /// [`HgvsError::ValidationError`] when the segment has fewer than
    /// `from + n` bases, or the reference fewer than the segment names.
    pub fn read(
        &self,
        reference: &Reference<'_, '_>,
        from: usize,
        n: usize,
    ) -> Result<String, HgvsError> {
        if from.checked_add(n).is_none_or(|last| last > self.len()) {
            return Err(HgvsError::ValidationError(format!(
                "a segment of {} bases has none at [{from}, {from} + {n})",
                self.len()
            )));
        }
        let slice = |s: usize, e: usize| {
            let bases = reference.slice(s, e)?;
            if bases.len() != e - s {
                return Err(HgvsError::ValidationError(format!(
                    "{} has no bases over [{s}, {e})",
                    reference.accession()
                )));
            }
            Ok(bases.to_ascii_uppercase())
        };
        match self {
            Segment::Reference { start, .. } => slice(start + from, start + from + n),
            Segment::ReverseComplement { end, .. } => {
                Ok(reverse_complement(&slice(end - from - n, end - from)?))
            }
            Segment::Literal(s) => Ok(s[from..from + n].to_string()),
        }
    }

    fn reverse_complement(&self) -> Segment {
        match self {
            Segment::Reference { start, end } => Segment::ReverseComplement {
                start: *start,
                end: *end,
            },
            Segment::ReverseComplement { start, end } => Segment::Reference {
                start: *start,
                end: *end,
            },
            Segment::Literal(s) => Segment::Literal(reverse_complement(s)),
        }
    }
}

/// A stretch of the reference, `[start, end)`, replaced by the bases of
/// [`with`](Self::with); an insertion replaces the empty stretch at the
/// junction in front of `start`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Replacement {
    start: usize,
    end: usize,
    with: Vec<Segment>,
}

impl Replacement {
    /// The first reference index replaced.
    pub fn start(&self) -> usize {
        self.start
    }

    /// One past the last reference index replaced; `start` for an insertion.
    pub fn end(&self) -> usize {
        self.end
    }

    /// What the stretch becomes.
    pub fn with(&self) -> &[Segment] {
        &self.with
    }
}

/// An allele in cis read into one sequence: the reference with
/// [`replacements`](Self::replacements) applied, the reference bases it
/// [states](Self::stated) unchanged, and the members that
/// [cannot be applied](Self::unapplied) to a sequence.
///
/// The replacements are sorted and do not overlap, and two that touch are
/// one; so the sequence over any stretch is read by walking them.
/// [`segments`](Self::segments) gives that walk, which equivalence compares,
/// normalisation splits into changes and `c_to_p` translates.
#[derive(Debug, Clone)]
pub struct CisReading {
    accession: String,
    kind: IdentifierType,
    replacements: Vec<Replacement>,
    stated: StatedBases,
    unapplied: Vec<SV>,
}

impl CisReading {
    /// The sequence the allele is read on: its own, or the genome a
    /// transcript allele is placed on, when it has an intronic member or is
    /// read with [`read_cis_on_genome`](VariantMapper::read_cis_on_genome).
    pub fn accession(&self) -> &str {
        &self.accession
    }

    /// Whether [`accession`](Self::accession) is a transcript or a genome.
    pub fn kind(&self) -> IdentifierType {
        self.kind
    }

    pub fn replacements(&self) -> &[Replacement] {
        &self.replacements
    }

    /// The bases the allele states unchanged, as sorted 0-based half-open
    /// ranges that neither overlap nor touch: those its identity members name,
    /// and those named by members that together make the reference.
    pub fn stated(&self) -> &[(usize, usize)] {
        self.stated.ranges()
    }

    /// The members no sequence can be built from: a copy number, an uncertain
    /// position, a conversion, inserted bases not stated or written `N`, a
    /// statement. They take part in no replacement; each is held to the
    /// overlap rule on its outer bounds. They are listed as written; identical
    /// ones other than insertions and duplications are one member written
    /// twice, and are listed once.
    pub fn unapplied(&self) -> &[SV] {
        &self.unapplied
    }

    /// The stretch from the first reference base the replacements and
    /// statements touch to the last, half-open; `None` when they touch none.
    pub fn span(&self) -> Option<(usize, usize)> {
        let ends = self
            .replacements
            .iter()
            .map(|r| (r.start, r.end))
            .chain(self.stated.ranges().iter().copied());
        ends.reduce(|a, b| (a.0.min(b.0), a.1.max(b.1)))
    }

    /// The allele's sequence over the reference stretch `[start, end)`, as
    /// segments. A replacement that reaches outside the stretch is taken
    /// whole, so the stretch is widened to hold it.
    pub fn segments(&self, start: usize, end: usize) -> Vec<Segment> {
        let inside: Vec<&Replacement> = self
            .replacements
            .iter()
            .filter(|r| {
                r.start < end && start < r.end
                    || (r.start == r.end && (start..=end).contains(&r.start))
            })
            .collect();
        let lo = inside.iter().map(|r| r.start).fold(start, usize::min);
        let hi = inside.iter().map(|r| r.end).fold(end, usize::max);
        pieces(inside.into_iter(), lo, hi)
    }

    /// The allele's bases over the reference stretch `[start, end)`, widened
    /// as [`segments`](Self::segments) widens it. This builds them: it is for
    /// short stretches.
    ///
    /// # Errors
    ///
    /// [`HgvsError::ValidationError`] when `reference` is not the sequence
    /// the allele is read on.
    pub fn sequence(
        &self,
        reference: &Reference<'_, '_>,
        start: usize,
        end: usize,
    ) -> Result<String, HgvsError> {
        self.check_reference(reference)?;
        materialise(&self.segments(start, end), reference)
    }

    /// Whether the allele makes the reference's sequence.
    ///
    /// # Errors
    ///
    /// [`HgvsError::ValidationError`] when `reference` is not the sequence
    /// the allele is read on.
    pub fn makes_reference(&self, reference: &Reference<'_, '_>) -> Result<bool, HgvsError> {
        self.check_reference(reference)?;
        same_replacements(&self.replacements, &[], reference)
    }

    /// Whether the two readings make the same sequence, compared over the
    /// stretch either touches; bases outside it are the reference's for both.
    ///
    /// # Errors
    ///
    /// [`HgvsError::ValidationError`] when the two are read on different
    /// sequences, or `reference` is not theirs: indices on one say nothing
    /// about the other.
    pub fn same_sequence(
        &self,
        other: &CisReading,
        reference: &Reference<'_, '_>,
    ) -> Result<bool, HgvsError> {
        if (self.accession.as_str(), self.kind) != (other.accession.as_str(), other.kind) {
            return Err(HgvsError::ValidationError(format!(
                "one reading is on {} and the other on {}, so their sequences cannot be compared",
                self.accession, other.accession
            )));
        }
        self.check_reference(reference)?;
        same_replacements(&self.replacements, &other.replacements, reference)
    }

    fn check_reference(&self, reference: &Reference<'_, '_>) -> Result<(), HgvsError> {
        if (reference.accession(), reference.kind()) == (self.accession.as_str(), self.kind) {
            return Ok(());
        }
        Err(HgvsError::ValidationError(format!(
            "the reading is on {}, not {}",
            self.accession,
            reference.accession()
        )))
    }
}

/// The sequence `replacements` make over `[lo, hi)`, which holds them all.
fn pieces<'r>(
    replacements: impl Iterator<Item = &'r Replacement>,
    lo: usize,
    hi: usize,
) -> Vec<Segment> {
    let mut out = Vec::new();
    let mut at = lo;
    for r in replacements {
        if r.start > at {
            out.push(Segment::Reference {
                start: at,
                end: r.start,
            });
        }
        out.extend(r.with.iter().filter(|s| !s.is_empty()).cloned());
        at = at.max(r.end);
    }
    if hi > at {
        out.push(Segment::Reference { start: at, end: hi });
    }
    out
}

fn materialise(segments: &[Segment], reference: &Reference<'_, '_>) -> Result<String, HgvsError> {
    let mut out = String::new();
    for s in segments {
        out.push_str(&s.read(reference, 0, s.len())?);
    }
    Ok(out)
}

/// Bases read and compared at a time where two walks are not on the same
/// reference bases.
const CHUNK: usize = 4096;

/// Whether two sequences of segments spell the same bases. Where both are at
/// the same reference base, the run they share is equal without reading it.
fn same_bases(
    a: &[Segment],
    b: &[Segment],
    reference: &Reference<'_, '_>,
) -> Result<bool, HgvsError> {
    if a.iter().map(Segment::len).sum::<usize>() != b.iter().map(Segment::len).sum::<usize>() {
        return Ok(false);
    }
    let (mut ia, mut oa, mut ib, mut ob) = (0, 0, 0, 0);
    loop {
        while ia < a.len() && oa == a[ia].len() {
            ia += 1;
            oa = 0;
        }
        while ib < b.len() && ob == b[ib].len() {
            ib += 1;
            ob = 0;
        }
        let (Some(sa), Some(sb)) = (a.get(ia), b.get(ib)) else {
            return Ok(ia == a.len() && ib == b.len());
        };
        let n = (sa.len() - oa).min(sb.len() - ob);
        let aligned = matches!(
            (sa, sb),
            (Segment::Reference { start: x, .. }, Segment::Reference { start: y, .. })
                if x + oa == y + ob
        );
        let n = if aligned {
            n
        } else {
            let n = n.min(CHUNK);
            if !sa
                .read(reference, oa, n)?
                .eq_ignore_ascii_case(&sb.read(reference, ob, n)?)
            {
                return Ok(false);
            }
            n
        };
        oa += n;
        ob += n;
    }
}

/// Whether two sorted lists of replacements make the same sequence, compared
/// over the stretch either touches.
fn same_replacements(
    a: &[Replacement],
    b: &[Replacement],
    reference: &Reference<'_, '_>,
) -> Result<bool, HgvsError> {
    let Some((lo, hi)) = a
        .iter()
        .chain(b)
        .map(|r| (r.start, r.end))
        .reduce(|x, y| (x.0.min(y.0), x.1.max(y.1)))
    else {
        return Ok(true);
    };
    same_bases(
        &pieces(a.iter(), lo, hi),
        &pieces(b.iter(), lo, hi),
        reference,
    )
}

/// What a base-by-base member puts at each base of its stretch.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
enum Value {
    /// Nothing: the base is deleted.
    Delete,
    /// The base itself: a positional identity.
    Keep,
    /// These bases, one per base of the stretch.
    Bases(String),
}

/// A member read as the reference stretch it replaces and what it puts there.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
enum Act {
    /// A substitution, deletion, delins that keeps the length, or identity:
    /// each base of `[start, end)` gets one base, or none.
    Write {
        start: usize,
        end: usize,
        value: Value,
    },
    /// A delins that changes the length, replacing `[start, end)` as a whole.
    Block {
        start: usize,
        end: usize,
        alt: String,
    },
    /// Bases inserted at the junction in front of `at`.
    Insert { at: usize, alt: String },
    /// A copy of `[start, end)` inserted at the junction in front of `at`:
    /// `end`, 3' of the original, or `start` for a minus-strand transcript's
    /// duplication read on the genome, where its 3' is the genome's 5'.
    Dup { start: usize, end: usize, at: usize },
    /// `[start, end)` replaced by its reverse complement.
    Inv { start: usize, end: usize },
}

impl Act {
    /// The bases the member writes, replaces or reads.
    fn bases(&self) -> Option<(usize, usize)> {
        match self {
            Act::Write { start, end, .. }
            | Act::Block { start, end, .. }
            | Act::Dup { start, end, .. }
            | Act::Inv { start, end } => Some((*start, *end)),
            Act::Insert { .. } => None,
        }
    }

    /// The junction the member inserts at, if it inserts.
    fn junction(&self) -> Option<usize> {
        match self {
            Act::Insert { at, .. } => Some(*at),
            Act::Dup { at, .. } => Some(*at),
            _ => None,
        }
    }

    /// The positions the description names: a range every base in it, an
    /// insertion its two flanks, a duplication the bases it copies.
    fn named(&self) -> (usize, usize) {
        match self {
            Act::Insert { at, .. } => (at.saturating_sub(1), at + 1),
            other => other.bases().expect("every other act has bases"),
        }
    }

    /// The stretch the member is written over: an insertion's is empty, at
    /// its junction.
    fn written(&self) -> (usize, usize) {
        match self {
            Act::Insert { at, .. } => (*at, *at),
            other => other.bases().expect("every other act has bases"),
        }
    }

    fn is_identity(&self) -> bool {
        matches!(
            self,
            Act::Write {
                value: Value::Keep,
                ..
            }
        )
    }

    fn is_reader(&self) -> bool {
        matches!(self, Act::Dup { .. } | Act::Inv { .. })
    }

    fn changes_length(&self) -> bool {
        matches!(
            self,
            Act::Insert { .. }
                | Act::Dup { .. }
                | Act::Block { .. }
                | Act::Write {
                    value: Value::Delete,
                    ..
                }
        )
    }
}

/// A member's place on the sequence the allele is read on.
#[derive(Debug, Clone)]
enum Placed {
    Applied(Act),
    /// A member no sequence can be built from, held to the overlap rule on
    /// what it reaches for.
    Unapplied(Reach),
}

/// How a duplication or inversion reads bases another member changes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum Mode {
    /// The reference's bases; the other members' changes then land on their
    /// own positions.
    Reference,
    /// The bases as the other members leave them.
    Edited,
}

/// One way of reading the members: which insertions inside a replaced stretch
/// go with it, how each reader reads, and the order of the insertions at each
/// junction (one written twice appears once).
#[derive(Debug, Clone)]
struct Choice {
    dropped: Vec<bool>,
    mode: Vec<Mode>,
    order: BTreeMap<usize, Vec<usize>>,
}

/// The members of an allele placed on one reference, with what is needed to
/// name them in an error.
struct Allele<'x, 's, 'a> {
    reference: &'x Reference<'s, 'a>,
    /// The allele as written, for errors.
    description: String,
    /// Each member as written, by its index in the allele.
    text: Vec<String>,
    /// The applied members: their index in the allele and their act.
    members: Vec<(usize, Act)>,
    /// The members that cannot be applied.
    outer: Vec<Outer>,
}

/// What a member that cannot be applied reaches for.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Reach {
    /// The bases `[lo, hi)`; an unknown bound reaches the end of the sequence
    /// on its side.
    Stretch(usize, usize),
    /// The junction in front of an index, where unknown bases are inserted.
    Junction(usize),
    /// A duplication at an uncertain position, or a copy number: the bases
    /// `[lo, hi)`, and the junctions at either end, where a copy may go.
    Copy(usize, usize),
}

impl Reach {
    fn bases(self) -> Option<(usize, usize)> {
        match self {
            Reach::Stretch(lo, hi) | Reach::Copy(lo, hi) => Some((lo, hi)),
            Reach::Junction(_) => None,
        }
    }

    fn junctions(self) -> Vec<usize> {
        match self {
            Reach::Stretch(..) => Vec::new(),
            Reach::Junction(at) => vec![at],
            Reach::Copy(lo, hi) => vec![lo, hi],
        }
    }

    /// The reach of `edit` written over `[lo, hi)` at a position, or with
    /// bases, no sequence can be built from.
    fn of(edit: &NaEdit, lo: usize, hi: usize) -> Reach {
        match edit {
            NaEdit::Ins { .. } | NaEdit::InsLength { .. } if hi == lo + 2 => {
                Reach::Junction(lo + 1)
            }
            NaEdit::Dup { .. } | NaEdit::NACopy { .. } => Reach::Copy(lo, hi),
            _ => Reach::Stretch(lo, hi),
        }
    }
}

/// A member that cannot be applied: its index in the allele and its reach.
#[derive(Debug, Clone, Copy)]
struct Outer {
    index: usize,
    reach: Reach,
}

/// A choice described by what its members do: how many read edited bases,
/// how many insertions are dropped, which, and the order at each junction.
type ChoiceKey<'m> = (
    usize,
    usize,
    Vec<(&'m Act, Mode, bool)>,
    Vec<(usize, Vec<&'m Act>)>,
);

/// The most readings listed, and the most insertions ordered at one junction.
const MAX_READINGS: usize = 1024;
const MAX_AT_JUNCTION: usize = 6;
/// The most members that change the length or invert, and the most changes
/// among them, searched for sets that cancel.
const MAX_CANCEL_SEARCH: usize = 16;
/// The most choices a description leaves open besides the order at a
/// junction: an insertion kept or dropped, a copy of the reference or of the
/// edited bases.
const MAX_CHOICES: usize = 10;
/// The most sets of members tried for cancelling, among members that meet.
const MAX_CANDIDATE_SETS: usize = 1 << 16;
/// The most sets of cancelling members that share members with each other.
const MAX_CLASHING_SETS: usize = 16;

fn overlap(a: (usize, usize), b: (usize, usize)) -> Option<(usize, usize)> {
    let lo = a.0.max(b.0);
    let hi = a.1.min(b.1);
    (lo < hi).then_some((lo, hi))
}

fn strictly_inside(j: usize, (s, e): (usize, usize)) -> bool {
    s < j && j < e
}

impl<'x, 's, 'a> Allele<'x, 's, 'a> {
    fn slice(&self, s: usize, e: usize) -> Result<String, HgvsError> {
        Segment::Reference { start: s, end: e }.read(self.reference, 0, e - s)
    }

    fn contradiction(&self, a: usize, b: usize, why: &str) -> HgvsError {
        HgvsError::ValidationError(format!(
            "{} has no reading: {} and {} {why}",
            self.description, self.text[a], self.text[b]
        ))
    }

    /// Whether the inversion over `[s, e)` changes any base of `[a, b)`.
    fn inverted_differs(
        &self,
        (s, e): (usize, usize),
        (a, b): (usize, usize),
    ) -> Result<bool, HgvsError> {
        let inverted = reverse_complement(&self.slice(s + e - b, s + e - a)?);
        Ok(!inverted.eq_ignore_ascii_case(&self.slice(a, b)?))
    }

    /// Whether two base-by-base values agree over `[a, b)`.
    fn writes_agree(
        &self,
        (v1, s1): (&Value, usize),
        (v2, s2): (&Value, usize),
        (a, b): (usize, usize),
    ) -> Result<bool, HgvsError> {
        Ok(match (v1, v2) {
            (Value::Delete, Value::Delete) | (Value::Keep, Value::Keep) => true,
            (Value::Delete, _) | (_, Value::Delete) => false,
            (Value::Keep, Value::Bases(x)) => {
                x[a - s2..b - s2].eq_ignore_ascii_case(&self.slice(a, b)?)
            }
            (Value::Bases(x), Value::Keep) => {
                x[a - s1..b - s1].eq_ignore_ascii_case(&self.slice(a, b)?)
            }
            (Value::Bases(x), Value::Bases(y)) => {
                x[a - s1..b - s1].eq_ignore_ascii_case(&y[a - s2..b - s2])
            }
        })
    }

    /// Refuses the allele if two members leave it no reading, and then if a
    /// member that cannot be applied leaves readings that cannot be listed.
    fn check_contradictions(&self) -> Result<(), HgvsError> {
        if let Some(e) = self.contradiction_in()? {
            return Err(e);
        }
        self.outer_unlisted()
    }

    /// Refuses a member that inserts bases it does not state, or a copy at
    /// an uncertain position, where another member acts at its junction: the
    /// readings exist, and cannot be written down.
    fn outer_unlisted(&self) -> Result<(), HgvsError> {
        for o in &self.outer {
            for j in o.reach.junctions() {
                // Copies read their bases, which an identity leaves as they are.
                let copy = matches!(o.reach, Reach::Copy(..));
                let applied = self.members.iter().find(|(_, act)| {
                    !(copy && act.is_identity())
                        && (act.bases().is_some_and(|r| strictly_inside(j, r))
                            || act.junction() == Some(j))
                });
                let unapplied = self.outer.iter().find(|other| {
                    other.index != o.index
                        && (other.reach.junctions().contains(&j)
                            || other.reach.bases().is_some_and(|r| strictly_inside(j, r)))
                });
                let beside = applied
                    .map(|(k, _)| *k)
                    .or(unapplied.map(|other| other.index));
                if let Some(k) = beside {
                    return Err(HgvsError::UnsupportedOperation(format!(
                        "{}: {} has no sequence to read, and {} acts where it inserts, so its \
                         readings cannot be listed",
                        self.description, self.text[o.index], self.text[k]
                    )));
                }
            }
        }
        Ok(())
    }

    /// The contradiction two members make, if any; `Err` only when the
    /// reference cannot be read.
    fn contradiction_in(&self) -> Result<Option<HgvsError>, HgvsError> {
        for o in &self.outer {
            if let Some(e) = self.outer_contradiction(o)? {
                return Ok(Some(e));
            }
        }
        for (x, &(i, ref a)) in self.members.iter().enumerate() {
            for &(j, ref b) in &self.members[x + 1..] {
                if let Some(e) = self.check_pair((i, a), (j, b))? {
                    return Ok(Some(e));
                }
                if let Some(e) = self.check_pair((j, b), (i, a))? {
                    return Ok(Some(e));
                }
            }
        }
        Ok(None)
    }

    /// A member that cannot be applied, held to its outer bounds against the
    /// others: another member reaching for one of its bases, or inserting
    /// among them, leaves no reading.
    fn outer_contradiction(&self, o: &Outer) -> Result<Option<HgvsError>, HgvsError> {
        if let Reach::Junction(at) = o.reach {
            // An identity, or a delins that changes the length, leaves no
            // junction inside it to insert at, whatever is inserted.
            for &(j, ref act) in &self.members {
                let closed = matches!(act, Act::Block { .. }) || act.is_identity();
                if closed && act.bases().is_some_and(|r| strictly_inside(at, r)) {
                    return Ok(Some(self.contradiction(
                        o.index,
                        j,
                        "conflict: the first inserts where the second leaves no junction",
                    )));
                }
            }
        }
        let Some((lo, hi)) = o.reach.bases() else {
            return Ok(None);
        };
        for &(j, ref act) in &self.members {
            // A copy reads its bases, which an identity leaves as they are.
            if matches!(o.reach, Reach::Copy(..)) && act.is_identity() {
                continue;
            }
            let inside = act.bases().and_then(|r| overlap(r, (lo, hi))).is_some()
                || act.junction().is_some_and(|k| strictly_inside(k, (lo, hi)));
            if inside {
                return Ok(Some(self.contradiction(
                    o.index,
                    j,
                    "reach for one base, and the first has no sequence to read it by",
                )));
            }
        }
        for other in &self.outer {
            if other.index == o.index {
                continue;
            }
            let clash = other
                .reach
                .bases()
                .and_then(|r| overlap(r, (lo, hi)))
                .is_some()
                || other
                    .reach
                    .junctions()
                    .iter()
                    .any(|&k| strictly_inside(k, (lo, hi)));
            if clash {
                return Ok(Some(self.contradiction(
                    o.index,
                    other.index,
                    "reach for one base, and neither has a sequence to read it by",
                )));
            }
        }
        Ok(None)
    }

    /// The contradictions in which `a` plays its part against `b`; each pair
    /// is checked both ways round.
    fn check_pair(
        &self,
        (i, a): (usize, &Act),
        (j, b): (usize, &Act),
    ) -> Result<Option<HgvsError>, HgvsError> {
        match (a, b) {
            (
                Act::Write {
                    start: s1,
                    end: e1,
                    value: v1,
                },
                Act::Write {
                    start: s2,
                    end: e2,
                    value: v2,
                },
            ) if i < j => {
                if let Some(r) = overlap((*s1, *e1), (*s2, *e2)) {
                    if !self.writes_agree((v1, *s1), (v2, *s2), r)? {
                        return Ok(Some(self.contradiction(
                            i,
                            j,
                            "put different bases at one position",
                        )));
                    }
                }
            }
            (Act::Block { start, end, .. }, other) => {
                if other
                    .bases()
                    .and_then(|r| overlap(r, (*start, *end)))
                    .is_some()
                {
                    return Ok(Some(self.contradiction(
                        i,
                        j,
                        "share a base, and the first changes the length over it, so no base is left for the second",
                    )));
                }
                if other
                    .junction()
                    .is_some_and(|k| strictly_inside(k, (*start, *end)))
                {
                    return Ok(Some(self.contradiction(
                        i,
                        j,
                        "conflict: the second inserts inside a stretch the first replaces with another length, where no junction is left",
                    )));
                }
            }
            (
                Act::Write {
                    start,
                    end,
                    value: Value::Keep,
                },
                other,
            ) => {
                if other
                    .junction()
                    .is_some_and(|k| strictly_inside(k, (*start, *end)))
                {
                    return Ok(Some(self.contradiction(
                        i,
                        j,
                        "conflict: the first states its bases unchanged and side by side, and the second inserts between them",
                    )));
                }
                if let Act::Inv { start: s, end: e } = other {
                    if let Some(r) = overlap((*start, *end), (*s, *e)) {
                        if self.inverted_differs((*s, *e), r)? {
                            return Ok(Some(self.contradiction(
                                i,
                                j,
                                "conflict: the first states a base unchanged that the second inverts",
                            )));
                        }
                    }
                }
            }
            _ => {}
        }
        Ok(None)
    }

    /// Whether another member changes, removes or inserts among the bases the
    /// reader `x` reads.
    fn interacts(&self, x: usize) -> Result<bool, HgvsError> {
        let (_, reader) = &self.members[x];
        let range = reader.bases().expect("a reader has bases");
        let mut changed = false;
        for (y, (_, other)) in self.members.iter().enumerate() {
            if y == x {
                continue;
            }
            if other.junction().is_some_and(|k| strictly_inside(k, range)) {
                changed = true;
            }
            match other {
                Act::Write { start, end, value } => {
                    if let Some(r) = overlap(range, (*start, *end)) {
                        changed |= !self.writes_agree((value, *start), (&Value::Keep, r.0), r)?;
                    }
                }
                Act::Inv { start, end } => {
                    let Some(r) = overlap(range, (*start, *end)) else {
                        continue;
                    };
                    if !self.inverted_differs((*start, *end), r)? {
                        continue;
                    }
                    if matches!(reader, Act::Inv { .. }) {
                        return Err(HgvsError::UnsupportedOperation(format!(
                            "{}: {} and {} invert overlapping stretches each changes, \
                             whose readings are not listed",
                            self.description,
                            self.text[self.members[x].0],
                            self.text[self.members[y].0]
                        )));
                    }
                    changed = true;
                }
                _ => {}
            }
        }
        Ok(changed)
    }

    /// The reading of the members under `choice`: over `window` alone when
    /// one is given, leaving out the member `exclude`, as a duplication or
    /// inversion reading its stretch sees it.
    fn assemble(
        &self,
        choice: &Choice,
        window: Option<(usize, usize)>,
        exclude: Option<usize>,
        stack: &mut Vec<usize>,
    ) -> Result<Vec<Replacement>, HgvsError> {
        Ok(merge_touching(
            self.assemble_runs(choice, window, exclude, stack)?,
        ))
    }

    /// As [`assemble`](Self::assemble), with replacements that touch left
    /// apart: no replacement then reaches past the bases of the members that
    /// make it.
    fn assemble_runs(
        &self,
        choice: &Choice,
        window: Option<(usize, usize)>,
        exclude: Option<usize>,
        stack: &mut Vec<usize>,
    ) -> Result<Vec<Replacement>, HgvsError> {
        let live = |k: usize| Some(k) != exclude && !choice.dropped[k];
        let in_window = |r: (usize, usize)| {
            window.is_none_or(|w| overlap(r, w).is_some() || r.0 == r.1 && strictly_inside(r.0, w))
        };
        let straddles = |r: (usize, usize)| {
            window.is_some_and(|(lo, hi)| r.0 < lo && lo < r.1 || r.0 < hi && hi < r.1)
        };

        let mut whole: Vec<Replacement> = Vec::new();
        let mut edited_inv: Vec<(usize, usize)> = Vec::new();
        for (k, (_, act)) in self.members.iter().enumerate() {
            if !live(k) {
                continue;
            }
            let replaced = match act {
                Act::Block { start, end, alt } => {
                    Some((*start, *end, vec![Segment::Literal(alt.clone())]))
                }
                Act::Inv { start, end } if choice.mode[k] == Mode::Edited => {
                    let stretch = self.stretch(choice, k, stack)?;
                    edited_inv.push((*start, *end));
                    Some((
                        *start,
                        *end,
                        stretch
                            .iter()
                            .rev()
                            .map(Segment::reverse_complement)
                            .collect(),
                    ))
                }
                _ => None,
            };
            if let Some((start, end, with)) = replaced {
                if !in_window((start, end)) {
                    continue;
                }
                if straddles((start, end)) {
                    return Err(HgvsError::UnsupportedOperation(format!(
                        "{}: a stretch replaced as a whole reaches across the edge of a \
                         stretch a duplication or inversion reads",
                        self.description
                    )));
                }
                whole.push(Replacement { start, end, with });
            }
        }
        let consumed = |p: usize| edited_inv.iter().any(|&(s, e)| s <= p && p < e);

        let mut junctions: Vec<(usize, Vec<Segment>)> = Vec::new();
        for (&j, slots) in &choice.order {
            if window.is_some_and(|w| !strictly_inside(j, w))
                || edited_inv.iter().any(|&r| strictly_inside(j, r))
            {
                continue;
            }
            let mut with = Vec::new();
            for &k in slots {
                if !live(k) {
                    continue;
                }
                match &self.members[k].1 {
                    Act::Insert { alt, .. } => with.push(Segment::Literal(alt.clone())),
                    Act::Dup { start, end, .. } => match choice.mode[k] {
                        Mode::Reference => with.push(Segment::Reference {
                            start: *start,
                            end: *end,
                        }),
                        Mode::Edited => with.extend(self.stretch(choice, k, stack)?),
                    },
                    _ => unreachable!("only insertions and duplications sit at a junction"),
                }
            }
            if !with.is_empty() {
                junctions.push((j, with));
            }
        }

        // Base by base: the writes, and inversions read on the reference with
        // the writes landing on their positions; split where a junction is.
        let mut cuts: Vec<usize> = Vec::new();
        let mut writes: Vec<(usize, usize, &Value)> = Vec::new();
        let mut inverted: Vec<(usize, usize)> = Vec::new();
        for (k, (_, act)) in self.members.iter().enumerate() {
            if !live(k) {
                continue;
            }
            match act {
                Act::Write { start, end, value } if *value != Value::Keep => {
                    writes.push((*start, *end, value))
                }
                Act::Inv { start, end } if choice.mode[k] == Mode::Reference => {
                    inverted.push((*start, *end))
                }
                _ => continue,
            }
            let (s, e) = act.bases().expect("writes and inversions have bases");
            cuts.extend([s, e]);
        }
        cuts.extend(junctions.iter().map(|(j, _)| *j));
        cuts.extend(edited_inv.iter().flat_map(|&(s, e)| [s, e]));
        if let Some((lo, hi)) = window {
            cuts.extend([lo, hi]);
        }
        cuts.sort_unstable();
        cuts.dedup();
        let mut runs: Vec<Replacement> = Vec::new();
        for pair in cuts.windows(2) {
            let (a, b) = (pair[0], pair[1]);
            if consumed(a) || window.is_some_and(|(lo, hi)| a < lo || b > hi) {
                continue;
            }
            let with = if let Some((s, _, value)) = writes.iter().find(|w| w.0 <= a && b <= w.1) {
                match value {
                    Value::Delete => Vec::new(),
                    Value::Bases(x) => vec![Segment::Literal(x[a - s..b - s].to_string())],
                    Value::Keep => unreachable!("identities make no run"),
                }
            } else if let Some((s, e)) = inverted.iter().find(|r| r.0 <= a && b <= r.1) {
                vec![Segment::ReverseComplement {
                    start: s + e - b,
                    end: s + e - a,
                }]
            } else {
                continue;
            };
            runs.push(Replacement {
                start: a,
                end: b,
                with,
            });
        }

        let mut all: Vec<Replacement> = whole;
        all.extend(runs);
        all.extend(junctions.into_iter().map(|(j, with)| Replacement {
            start: j,
            end: j,
            with,
        }));
        // An insertion at a junction goes before a stretch starting there.
        all.sort_by_key(|r| (r.start, r.end));
        Ok(all)
    }

    /// Whether `replacements`, a reading of every member, change a base an
    /// identity states unchanged or insert between two of its bases: a
    /// duplication or inversion reading the edited bases can.
    fn breaks_an_identity(&self, replacements: &[Replacement]) -> Result<bool, HgvsError> {
        for (_, act) in &self.members {
            let Act::Write {
                start: s,
                end: e,
                value: Value::Keep,
            } = act
            else {
                continue;
            };
            for r in replacements {
                if r.start == r.end {
                    if strictly_inside(r.start, (*s, *e)) {
                        return Ok(true);
                    }
                    continue;
                }
                let Some((a, b)) = overlap((r.start, r.end), (*s, *e)) else {
                    continue;
                };
                let made: usize = r.with.iter().map(Segment::len).sum();
                if made != r.end - r.start {
                    return Ok(true);
                }
                let bases = materialise(&r.with, self.reference)?;
                let kept = &bases[a - r.start..b - r.start];
                if !kept.eq_ignore_ascii_case(&self.slice(a, b)?) {
                    return Ok(true);
                }
            }
        }
        Ok(false)
    }

    /// The bases the reader `k` reads, as the other members leave them.
    fn stretch(
        &self,
        choice: &Choice,
        k: usize,
        stack: &mut Vec<usize>,
    ) -> Result<Vec<Segment>, HgvsError> {
        if stack.contains(&k) {
            return Err(HgvsError::UnsupportedOperation(format!(
                "{}: {} reads bases that a member reading its own bases changes, whose readings are not enumerated",
                self.description, self.text[self.members[k].0]
            )));
        }
        let (s, e) = self.members[k].1.bases().expect("a reader has bases");
        stack.push(k);
        let replacements = self.assemble(choice, Some((s, e)), Some(k), stack)?;
        stack.pop();
        Ok(pieces(replacements.iter(), s, e))
    }

    /// The bases of each insertion at a junction, as `choice` reads it.
    fn item_bases(&self, choice: &Choice, k: usize) -> Result<String, HgvsError> {
        match &self.members[k].1 {
            Act::Insert { alt, .. } => Ok(alt.clone()),
            Act::Dup { start, end, .. } => match choice.mode[k] {
                Mode::Reference => self.slice(*start, *end),
                Mode::Edited => {
                    materialise(&self.stretch(choice, k, &mut Vec::new())?, self.reference)
                }
            },
            _ => unreachable!("only insertions and duplications sit at a junction"),
        }
    }

    /// Every reading the members have, distinct by the sequence it makes.
    fn readings(&self) -> Result<Vec<(Choice, Vec<Replacement>)>, HgvsError> {
        let n = self.members.len();
        let mut droppable = Vec::new();
        let mut readers = Vec::new();
        for (k, (_, act)) in self.members.iter().enumerate() {
            if let Some(j) = act.junction() {
                let inside_replaced = self.members.iter().any(|(_, o)| {
                    matches!(o, Act::Write { start, end, value }
                        if *value != Value::Keep && strictly_inside(j, (*start, *end)))
                });
                if inside_replaced {
                    droppable.push(k);
                }
            }
            if act.is_reader() && self.interacts(k)? {
                readers.push(k);
            }
        }
        let dims = droppable.len() + readers.len();
        if dims > MAX_CHOICES {
            return Err(self.too_many());
        }
        let mut written: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for (k, (_, act)) in self.members.iter().enumerate() {
            if let Some(j) = act.junction() {
                written.entry(j).or_default().push(k);
            }
        }
        let mut out: Vec<(Choice, Vec<Replacement>)> = Vec::new();
        for bits in 0u32..(1 << dims) {
            let mut base = Choice {
                dropped: vec![false; n],
                mode: vec![Mode::Reference; n],
                order: written.clone(),
            };
            for (d, &k) in droppable.iter().enumerate() {
                base.dropped[k] = bits & (1 << d) != 0;
            }
            for (d, &k) in readers.iter().enumerate() {
                if bits & (1 << (droppable.len() + d)) != 0 {
                    base.mode[k] = Mode::Edited;
                }
            }
            let mut per_junction: Vec<(usize, Vec<Vec<usize>>)> = Vec::new();
            for (&j, slots) in &written {
                let present: Vec<usize> = slots
                    .iter()
                    .copied()
                    .filter(|&k| !base.dropped[k])
                    .collect();
                if present.len() < 2 {
                    continue;
                }
                if present.len() > MAX_AT_JUNCTION {
                    return Err(HgvsError::UnsupportedOperation(format!(
                        "{} inserts {} times at one junction, whose orders are not listed",
                        self.description,
                        present.len()
                    )));
                }
                per_junction.push((j, self.orders(&base, &present, &written)?));
            }
            let combos = per_junction
                .iter()
                .try_fold(1usize, |acc, (_, a)| acc.checked_mul(a.len()))
                .filter(|&c| out.len() + c <= MAX_READINGS)
                .ok_or_else(|| self.too_many())?;
            for c in 0..combos {
                let mut choice = base.clone();
                let mut rest = c;
                for (j, options) in &per_junction {
                    choice
                        .order
                        .insert(*j, options[rest % options.len()].clone());
                    rest /= options.len();
                }
                let replacements = self.assemble(&choice, None, None, &mut Vec::new())?;
                if self.breaks_an_identity(&replacements)? {
                    continue;
                }
                // Choices that make one sequence are one reading, kept as the
                // choice whose key is least, so that which one is written does
                // not depend on the order of the members.
                let mut seen = None;
                for (i, (_, r)) in out.iter().enumerate() {
                    if same_replacements(r, &replacements, self.reference)? {
                        seen = Some(i);
                        break;
                    }
                }
                match seen {
                    None => out.push((choice, replacements)),
                    Some(i) if self.choice_key(&choice) < self.choice_key(&out[i].0) => {
                        out[i] = (choice, replacements);
                    }
                    Some(_) => {}
                }
            }
        }
        if out.is_empty() {
            return Err(HgvsError::ValidationError(format!(
                "{} has no reading: each way of reading it changes a base an identity states \
                 unchanged",
                self.description
            )));
        }
        Ok(out)
    }

    /// What a choice does, described by the members' acts rather than their
    /// places in the allele: fewest readers of edited bases first, then
    /// fewest insertions dropped, then the acts themselves.
    fn choice_key(&self, choice: &Choice) -> ChoiceKey<'_> {
        let edited = choice.mode.iter().filter(|&&m| m == Mode::Edited).count();
        let dropped = choice.dropped.iter().filter(|&&d| d).count();
        let mut acts: Vec<(&Act, Mode, bool)> = self
            .members
            .iter()
            .enumerate()
            .filter(|&(k, _)| choice.mode[k] == Mode::Edited || choice.dropped[k])
            .map(|(k, (_, act))| (act, choice.mode[k], choice.dropped[k]))
            .collect();
        acts.sort_unstable();
        let orders = choice
            .order
            .iter()
            .map(|(&j, slots)| (j, slots.iter().map(|&k| &self.members[k].1).collect()))
            .collect();
        (edited, dropped, acts, orders)
    }

    fn too_many(&self) -> HgvsError {
        HgvsError::UnsupportedOperation(format!(
            "{} has too many readings to list",
            self.description
        ))
    }

    /// The orders the insertions `present` at one junction can go in: every
    /// order of all of them, and, where some insert the same bases and so may
    /// be one written twice, every order of one of each. Orders that make one
    /// sequence are told apart later, by the readings they give.
    fn orders(
        &self,
        choice: &Choice,
        present: &[usize],
        written: &BTreeMap<usize, Vec<usize>>,
    ) -> Result<Vec<Vec<usize>>, HgvsError> {
        let mut bases = Vec::with_capacity(present.len());
        for &k in present {
            if let Act::Dup { start, end, .. } = self.members[k].1 {
                // What an edited copy holds depends on the order at a
                // junction inside it, which is being chosen too.
                let inner_order = written
                    .iter()
                    .any(|(&j, slots)| strictly_inside(j, (start, end)) && slots.len() > 1);
                if choice.mode[k] == Mode::Edited && inner_order {
                    return Err(HgvsError::UnsupportedOperation(format!(
                        "{}: {} copies a junction where several members insert, and inserts \
                         beside others, so its readings are not listed",
                        self.description, self.text[self.members[k].0]
                    )));
                }
            }
            bases.push(self.item_bases(choice, k)?);
        }
        // Insertions of the same bases, in the order written.
        let mut classes: Vec<Vec<usize>> = Vec::new();
        for (i, &k) in present.iter().enumerate() {
            match (0..i).find(|&p| bases[p] == bases[i]) {
                Some(p) => {
                    let class = classes
                        .iter_mut()
                        .find(|c| c[0] == present[p])
                        .expect("an earlier insertion has its class");
                    class.push(k);
                }
                None => classes.push(vec![k]),
            }
        }
        // Each class may be any number of its insertions, one written as
        // many times as it appears or several.
        let mut kept: Vec<Vec<usize>> = vec![Vec::new()];
        for class in &classes {
            kept = kept
                .iter()
                .flat_map(|k| {
                    (1..=class.len()).map(move |n| {
                        let mut k = k.clone();
                        k.extend(&class[..n]);
                        k
                    })
                })
                .collect();
        }
        let orders = kept.iter().flat_map(|k| permutations(k));
        // Where no copy reads edited bases, what each inserts is fixed, and
        // orders spelling one sequence at this junction are one.
        let fixed = present.iter().all(|&k| {
            !matches!(self.members[k].1, Act::Dup { .. }) || choice.mode[k] == Mode::Reference
        });
        if !fixed {
            return Ok(orders.collect());
        }
        let mut seen: Vec<String> = Vec::new();
        let mut out = Vec::new();
        for order in orders {
            let text: String = order
                .iter()
                .map(|k| bases[present.iter().position(|p| p == k).expect("present")].as_str())
                .collect();
            if !seen.contains(&text) {
                seen.push(text);
                out.push(order);
            }
        }
        Ok(out)
    }

    /// The members' extents: the stretch each is written over, widened to
    /// every place an insertion or deletion could equally be written.
    fn extents(&self) -> Result<Vec<(usize, usize)>, HgvsError> {
        self.members
            .iter()
            .map(|(_, act)| {
                let (s, e) = act.written();
                let resolved = match act {
                    Act::Insert { at, alt } => ResolvedEdit {
                        start: *at,
                        end: *at,
                        ref_: String::new(),
                        alt: alt.clone(),
                    },
                    Act::Dup { start, end, .. } => {
                        // The copy can go anywhere along the run of its bases.
                        let copied = self.slice(*start, *end)?;
                        let k5 = self.reference.run_left(*end, copied.as_bytes())?;
                        let k3 = self.reference.run_right(*end, copied.as_bytes())?;
                        return Ok(((*end - k5).min(s), (*end + k3).max(e)));
                    }
                    Act::Write {
                        start,
                        end,
                        value: Value::Delete,
                    } => ResolvedEdit {
                        start: *start,
                        end: *end,
                        ref_: self.slice(*start, *end)?,
                        alt: String::new(),
                    },
                    Act::Block { start, end, alt } => ResolvedEdit {
                        start: *start,
                        end: *end,
                        ref_: self.slice(*start, *end)?,
                        alt: alt.clone(),
                    },
                    _ => return Ok((s, e)),
                };
                let (lo, hi) = normalize::ambiguous_range(self.reference, &resolved)?;
                Ok((s.min(lo), e.max(hi)))
            })
            .collect()
    }

    /// The sub-allele of the members at `slots`.
    fn subset(&self, slots: &[usize]) -> Allele<'x, 's, 'a> {
        Allele {
            reference: self.reference,
            description: self.description.clone(),
            text: self.text.clone(),
            members: slots.iter().map(|&k| self.members[k].clone()).collect(),
            outer: Vec::new(),
        }
    }

    /// Whether the members at `slots` together have one reading, and it
    /// makes the reference.
    fn cancels(&self, slots: &[usize]) -> Result<bool, HgvsError> {
        let sub = self.subset(slots);
        if sub.contradiction_in()?.is_some() {
            return Ok(false);
        }
        let readings = sub.readings()?;
        match readings.as_slice() {
            [(_, only)] => same_replacements(only, &[], self.reference),
            _ => Ok(false),
        }
    }

    /// The members that meet, as groups of slots, among `slots`: extents
    /// that overlap or touch, directly or through another member.
    fn groups(&self, slots: &[usize], extents: &[(usize, usize)]) -> Vec<Vec<usize>> {
        let mut sorted: Vec<usize> = slots.to_vec();
        sorted.sort_by_key(|&k| extents[k]);
        let mut groups: Vec<(usize, Vec<usize>)> = Vec::new();
        for k in sorted {
            match groups.last_mut() {
                Some((reach, g)) if extents[k].0 <= *reach => {
                    *reach = (*reach).max(extents[k].1);
                    g.push(k);
                }
                _ => groups.push((extents[k].1, vec![k])),
            }
        }
        groups.into_iter().map(|(_, g)| g).collect()
    }

    /// The sets of members to try for cancelling in one group, smallest
    /// first: each base-by-base change alone, and each set of members that
    /// change the length or invert whose length changes sum to nothing,
    /// with any of the base-by-base changes inside its span. A base-by-base
    /// change outside that span changes a base nothing else touches, so it
    /// cannot be part of a set that makes the reference.
    fn candidate_sets(&self, group: &[usize]) -> Result<Vec<Vec<usize>>, HgvsError> {
        let (movers, writes): (Vec<usize>, Vec<usize>) = group.iter().partition(|&&k| {
            let act = &self.members[k].1;
            act.changes_length() || matches!(act, Act::Inv { .. })
        });
        if movers.len() > MAX_CANCEL_SEARCH {
            return Err(HgvsError::UnsupportedOperation(format!(
                "{}: {} members that change the length or invert meet, too many to search \
                 for members that cancel",
                self.description,
                movers.len()
            )));
        }
        let mut sets: Vec<Vec<usize>> = writes.iter().map(|&w| vec![w]).collect();
        let too_many = || {
            HgvsError::UnsupportedOperation(format!(
                "{}: the members that meet can be combined in too many ways to search for \
                 members that cancel",
                self.description
            ))
        };
        for size in 1..=movers.len() {
            for set in combinations(&movers, size) {
                let net: i64 = set.iter().map(|&k| self.length_change(k)).sum();
                if net != 0 {
                    continue;
                }
                let span = self.span_of(&set);
                let inside: Vec<usize> = writes
                    .iter()
                    .copied()
                    .filter(|&w| {
                        let bases = self.members[w].1.bases().expect("a write has bases");
                        overlap(bases, span).is_some()
                    })
                    .collect();
                if inside.len() > MAX_CANCEL_SEARCH {
                    return Err(HgvsError::UnsupportedOperation(format!(
                        "{}: {} changes lie among members that might cancel, too many to search",
                        self.description,
                        inside.len()
                    )));
                }
                let extras = 1usize
                    .checked_shl(inside.len() as u32)
                    .ok_or_else(too_many)?;
                if sets.len().saturating_add(extras) > MAX_CANDIDATE_SETS {
                    return Err(too_many());
                }
                for extra in 0..=inside.len() {
                    for w in combinations(&inside, extra) {
                        let mut candidate = set.clone();
                        candidate.extend(w);
                        candidate.sort_unstable();
                        sets.push(candidate);
                    }
                }
            }
        }
        sets.sort_by_key(Vec::len);
        sets.dedup();
        Ok(sets)
    }

    fn length_change(&self, k: usize) -> i64 {
        let len = |s: &usize, e: &usize| (e - s) as i64;
        match &self.members[k].1 {
            Act::Insert { alt, .. } => alt.len() as i64,
            Act::Dup { start, end, .. } => len(start, end),
            Act::Block { start, end, alt } => alt.len() as i64 - len(start, end),
            Act::Write {
                start,
                end,
                value: Value::Delete,
            } => -len(start, end),
            Act::Write { .. } | Act::Inv { .. } => 0,
        }
    }

    /// From the first position the members `set` name to the last.
    fn span_of(&self, set: &[usize]) -> (usize, usize) {
        set.iter()
            .map(|&k| self.members[k].1.named())
            .reduce(|a, b| (a.0.min(b.0), a.1.max(b.1)))
            .expect("a set has members")
    }

    /// The sets of members that together make the reference, no smaller set
    /// of them doing so, whose span holds no base another member names and
    /// lies in one of `exons`: members in two exons do not meet, and on the
    /// genome an intron lies between them.
    fn cancelling_sets(
        &self,
        extents: &[(usize, usize)],
        exons: &ExonBounds,
    ) -> Result<Vec<Vec<usize>>, HgvsError> {
        let changes: Vec<usize> = (0..self.members.len())
            .filter(|&k| !self.members[k].1.is_identity())
            .collect();
        let mut found: Vec<Vec<usize>> = Vec::new();
        for group in self.groups(&changes, extents) {
            let mut making: Vec<Vec<usize>> = Vec::new();
            for set in self.candidate_sets(&group)? {
                let span = self.span_of(&set);
                if exons.exon(span.0) != exons.exon(span.1 - 1)
                    || making.iter().any(|m| m.iter().all(|k| set.contains(k)))
                {
                    continue;
                }
                if self.cancels(&set)? {
                    making.push(set);
                }
            }
            found.extend(making);
        }
        let named = |k: usize| self.members[k].1.named();
        Ok(found
            .into_iter()
            .filter(|set| {
                let span = self.span_of(set);
                let clear_of_members = (0..self.members.len())
                    .filter(|k| !set.contains(k))
                    .all(|k| overlap(named(k), span).is_none());
                let clear_of_outer = self.outer.iter().all(|o| {
                    o.reach.bases().is_none_or(|r| overlap(r, span).is_none())
                        && o.reach
                            .junctions()
                            .iter()
                            .all(|&j| !strictly_inside(j, span))
                });
                clear_of_members && clear_of_outer
            })
            .collect())
    }
}

fn permutations(items: &[usize]) -> Vec<Vec<usize>> {
    if items.len() <= 1 {
        return vec![items.to_vec()];
    }
    let mut out = Vec::new();
    for i in 0..items.len() {
        let mut rest = items.to_vec();
        let first = rest.remove(i);
        for mut p in permutations(&rest) {
            p.insert(0, first);
            out.push(p);
        }
    }
    out
}

fn combinations(items: &[usize], size: usize) -> Vec<Vec<usize>> {
    if size == 0 {
        return vec![Vec::new()];
    }
    let mut out = Vec::new();
    for i in 0..items.len() {
        for mut rest in combinations(&items[i + 1..], size - 1) {
            rest.insert(0, items[i]);
            out.push(rest);
        }
    }
    out
}

/// `sorted` with every two replacements that touch made one.
fn merge_touching(sorted: Vec<Replacement>) -> Vec<Replacement> {
    let mut out: Vec<Replacement> = Vec::new();
    for r in sorted {
        match out.last_mut() {
            Some(last) if last.end == r.start => {
                last.end = r.end;
                last.with.extend(r.with);
            }
            _ => out.push(r),
        }
    }
    out
}

/// Sorted, disjoint ranges that neither overlap nor touch.
fn merge_ranges(mut ranges: Vec<(usize, usize)>) -> Vec<(usize, usize)> {
    ranges.sort_unstable();
    let mut out: Vec<(usize, usize)> = Vec::new();
    for r in ranges {
        match out.last_mut() {
            Some(last) if r.0 <= last.1 => last.1 = last.1.max(r.1),
            _ => out.push(r),
        }
    }
    out
}

/// The families of cancelling sets that can be statements together: every
/// set that shares no member with another is in each, and of the sets that
/// do share, each family takes a choice that shares none and to which no
/// other set could be added. One family when no two share a member.
fn families(sets: &[Vec<usize>]) -> Result<Vec<Vec<usize>>, HgvsError> {
    let n = sets.len();
    let clash = |a: usize, b: usize| sets[a].iter().any(|k| sets[b].contains(k));
    let mut component: Vec<usize> = (0..n).collect();
    fn root(component: &mut [usize], mut i: usize) -> usize {
        while component[i] != i {
            component[i] = component[component[i]];
            i = component[i];
        }
        i
    }
    for a in 0..n {
        for b in a + 1..n {
            if clash(a, b) {
                let (ra, rb) = (root(&mut component, a), root(&mut component, b));
                component[ra] = rb;
            }
        }
    }
    let mut components: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for i in 0..n {
        let r = root(&mut component, i);
        components.entry(r).or_default().push(i);
    }
    let mut families: Vec<Vec<usize>> = vec![Vec::new()];
    for members in components.values() {
        if members.len() > MAX_CLASHING_SETS {
            return Err(HgvsError::UnsupportedOperation(format!(
                "{} sets of members that cancel share members, too many to list the readings of",
                members.len()
            )));
        }
        let mut choices: Vec<Vec<usize>> = Vec::new();
        for mask in 1u32..(1 << members.len()) {
            let chosen: Vec<usize> = (0..members.len())
                .filter(|i| mask & (1 << i) != 0)
                .map(|i| members[i])
                .collect();
            let independent = chosen
                .iter()
                .enumerate()
                .all(|(x, &a)| chosen[x + 1..].iter().all(|&b| !clash(a, b)));
            let maximal = members
                .iter()
                .filter(|i| !chosen.contains(i))
                .all(|&i| chosen.iter().any(|&a| clash(a, i)));
            if independent && maximal {
                choices.push(chosen);
            }
        }
        if families.len().saturating_mul(choices.len()) > MAX_READINGS {
            return Err(HgvsError::UnsupportedOperation(
                "the sets of members that cancel can be chosen in too many ways to list".into(),
            ));
        }
        families = families
            .iter()
            .flat_map(|f| {
                choices.iter().map(move |c| {
                    let mut f = f.clone();
                    f.extend(c);
                    f
                })
            })
            .collect();
    }
    Ok(families)
}

/// The system an allele is read in, which decides how a reading is written.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum System {
    Genomic,
    Mitochondrial,
    Coding,
    NonCoding,
    Rna,
}

/// The sequence an allele is read on, and how to write positions on it.
struct Frame {
    accession: String,
    gene: Option<String>,
    kind: IdentifierType,
    system: System,
    exons: ExonBounds,
    mapper: Option<TranscriptMapper>,
    /// Each member placed, with its index in the allele as written, the
    /// member as written in the frame's system, and, for an applied member,
    /// its edit over the HGVS range it is written over.
    placed: Vec<(usize, Placed, SV, Option<PlacedEdit>)>,
}

/// The statements an allele's one reading makes: the bases it states
/// unchanged, the members set aside as statements, and the one reading of
/// the members left beside them.
struct Statement {
    stated: Vec<(usize, usize)>,
    /// The members, as slots of the allele, set aside: sets that together
    /// make the reference.
    chosen: Vec<usize>,
    /// The members left, as slots of the allele, and their one reading.
    rest: Vec<usize>,
    choice: Choice,
}

/// The one reading of an allele: the replacements its members make, and its
/// statements.
struct Settled {
    replacements: Vec<Replacement>,
    statement: Statement,
}

/// Which sequence a transcript allele is read on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum On {
    /// The transcript, unless a member is intronic.
    Written,
    /// The genome the transcript is placed on.
    Genome,
}

impl VariantMapper<'_> {
    /// Reads an allele in cis of nucleotide members into the one sequence
    /// they make, and the bases it states unchanged.
    ///
    /// The members are read on the sequence they are written on; a
    /// transcript allele with an intronic member is read on the genome, each
    /// member projected as [`c_to_g`](Self::c_to_g) projects it, and the
    /// readings of an ambiguous one are then written in `g.`.
    ///
    /// # Errors
    ///
    /// [`HgvsError::ReferenceMismatch`] when a member states a base the
    /// reference does not hold; [`HgvsError::ValidationError`] when two
    /// members leave the allele no reading, naming both;
    /// [`HgvsError::Ambiguous`] when it has several, listing each;
    /// [`HgvsError::UnsupportedOperation`] for an allele of protein members,
    /// or one whose readings weaver does not list: too many, two overlapping
    /// inversions or a duplication and an inversion that change each other's
    /// bases, or unstated inserted bases beside another member.
    pub fn read_cis(&self, cis: &CisPhasedVariant) -> Result<CisReading, HgvsError> {
        self.read_cis_on(cis, On::Written, true)
    }

    /// As [`read_cis`](Self::read_cis), with an ambiguity's readings left
    /// unwritten: for checking that a written reading has one reading.
    fn read_cis_unwritten(&self, cis: &CisPhasedVariant) -> Result<CisReading, HgvsError> {
        self.read_cis_on(cis, On::Written, false)
    }

    /// As [`read_cis`](Self::read_cis), always on the genome: a transcript
    /// allele's members are each projected as [`c_to_g`](Self::c_to_g)
    /// projects one, a minus-strand duplication's copy going 5' of the
    /// original on the genome, where it goes on the transcript.
    pub fn read_cis_on_genome(&self, cis: &CisPhasedVariant) -> Result<CisReading, HgvsError> {
        self.read_cis_on(cis, On::Genome, true)
    }

    fn read_cis_on(
        &self,
        cis: &CisPhasedVariant,
        on: On,
        write: bool,
    ) -> Result<CisReading, HgvsError> {
        let frame = self.checked_frame(cis, on)?;
        let reference = self.refs.reference(&frame.accession, frame.kind);
        let (allele, unapplied) = allele_on(cis, &frame, &reference);
        let settled = self.settle(cis, &frame, &allele, write)?;
        Ok(CisReading {
            stated: StatedBases::of_ranges(&frame.accession, settled.statement.stated),
            accession: frame.accession,
            kind: frame.kind,
            replacements: settled.replacements,
            unapplied,
        })
    }

    /// Where an allele of nucleotide members is read, once each member is
    /// found to state only bases its reference holds.
    fn checked_frame(&self, cis: &CisPhasedVariant, on: On) -> Result<Frame, HgvsError> {
        if cis.coordinate_type() == "p" {
            return Err(HgvsError::UnsupportedOperation(format!(
                "{cis} is an allele of protein members, which is not read into a sequence"
            )));
        }
        // A member stating no bases has none to check, and one whose bound is
        // unknown names no position to check them at.
        for m in cis
            .members
            .iter()
            .filter(|m| !has_ranged_bound(m) || states_bases(m))
        {
            self.check_stated(m)?;
        }
        self.frame(cis, on)
    }

    /// The one reading of `allele` and the statements it makes, or the error
    /// that it has none or several: `write` decides whether an ambiguity's
    /// readings are written out.
    fn settle(
        &self,
        cis: &CisPhasedVariant,
        frame: &Frame,
        allele: &Allele,
        write: bool,
    ) -> Result<Settled, HgvsError> {
        allele.check_contradictions()?;
        let mut readings = allele.readings()?;
        if readings.len() > 1 {
            return Err(self.ambiguity(cis, write, || {
                readings
                    .iter()
                    .map(|(choice, _)| self.describe(frame, allele, choice, &[], Writing::Reading))
                    .collect()
            }));
        }
        let (_, replacements) = readings
            .pop()
            .expect("members that are not contradictory have a reading");
        let extents = allele.extents()?;
        let sets = allele.cancelling_sets(&extents, &frame.exons)?;
        let mut statements: Vec<Statement> = Vec::new();
        for family in families(&sets)? {
            let chosen: Vec<usize> = family
                .iter()
                .flat_map(|&s| sets[s].iter().copied())
                .collect();
            let s = self.statement(frame, allele, chosen)?;
            if !statements.iter().any(|o| o.stated == s.stated) {
                statements.push(s);
            }
        }
        if statements.len() > 1 {
            return Err(self.ambiguity(cis, write, || {
                statements
                    .iter()
                    .map(|s| self.describe_statement(frame, allele, s, Writing::Reading))
                    .collect()
            }));
        }
        let statement = statements.pop().expect("the cancelling sets have a family");
        Ok(Settled {
            replacements,
            statement,
        })
    }

    /// The statements made when the sets of members `chosen` are set aside:
    /// those bases, and what the members left state once normalised, their
    /// identities and each group of them that meets and makes the reference.
    fn statement(
        &self,
        frame: &Frame,
        allele: &Allele,
        chosen: Vec<usize>,
    ) -> Result<Statement, HgvsError> {
        let rest: Vec<usize> = (0..allele.members.len())
            .filter(|k| !chosen.contains(k))
            .collect();
        let sub = allele.subset(&rest);
        let mut readings = sub.readings()?;
        let (Some((choice, _)), true) = (readings.pop(), readings.is_empty()) else {
            return Err(HgvsError::UnsupportedOperation(format!(
                "{}: the members left beside a statement have several readings",
                allele.description
            )));
        };
        let set_aside: Vec<(usize, usize)> = chosen
            .iter()
            .map(|&k| allele.members[k].1.named())
            .collect();
        let stated = self
            .normal_form(frame, &sub, &choice, &set_aside, Writing::Reading)?
            .stated;
        Ok(Statement {
            stated,
            chosen,
            rest,
            choice,
        })
    }

    /// The description of the reading in which `statement`'s sets are set
    /// aside.
    fn describe_statement(
        &self,
        frame: &Frame,
        allele: &Allele,
        statement: &Statement,
        writing: Writing,
    ) -> Result<SV, HgvsError> {
        let set_aside: Vec<(usize, usize)> = statement
            .chosen
            .iter()
            .map(|&k| allele.members[k].1.named())
            .collect();
        self.describe(
            frame,
            &allele.subset(&statement.rest),
            &statement.choice,
            &set_aside,
            writing,
        )
    }

    /// The error for an allele with several readings: listing each, written
    /// by `readings`, when `write` asks for them.
    fn ambiguity(
        &self,
        cis: &CisPhasedVariant,
        write: bool,
        readings: impl FnOnce() -> Result<Vec<SV>, HgvsError>,
    ) -> HgvsError {
        if !write {
            return HgvsError::Ambiguous {
                description: cis.to_string(),
                readings: Vec::new(),
            };
        }
        match readings() {
            Ok(mut readings) => {
                readings.sort_by_key(ToString::to_string);
                HgvsError::Ambiguous {
                    description: cis.to_string(),
                    readings,
                }
            }
            Err(e) => HgvsError::UnsupportedOperation(format!(
                "{cis} has several readings, and one could not be written: {e}"
            )),
        }
    }
    /// A placed edit as a variant of the frame's system.
    fn member_variant(&self, frame: &Frame, placed: &PlacedEdit) -> Result<SV, HgvsError> {
        let (s, e) = crate::mapper::hgvs_positions(placed);
        fn posedit<I>(pos: I, edit: &NaEdit) -> PosEdit<I, NaEdit> {
            PosEdit {
                pos: Some(pos),
                edit: edit.clone(),
                uncertain: false,
                predicted: false,
            }
        }
        fn transcript<V: TranscriptVariant>(
            am: &TranscriptMapper,
            s: usize,
            e: Option<usize>,
        ) -> Result<BaseOffsetInterval, HgvsError> {
            Ok(BaseOffsetInterval {
                start: V::position_from_index(am, s as i32)?,
                end: e
                    .map(|e| V::position_from_index(am, e as i32))
                    .transpose()?,
                uncertain: false,
            })
        }
        let linear = || SimpleInterval {
            start: crate::mapper::simple_position(s),
            end: e.map(crate::mapper::simple_position),
            uncertain: false,
        };
        let am = || {
            frame.mapper.as_ref().ok_or_else(|| {
                HgvsError::ValidationError(format!("{} has no transcript model", frame.accession))
            })
        };
        let edit = &placed.edit;
        let (ac, gene) = (frame.accession.clone(), frame.gene.clone());
        Ok(match frame.system {
            System::Genomic => SV::Genomic(GVariant::from_parts(ac, gene, posedit(linear(), edit))),
            System::Mitochondrial => {
                SV::Mitochondrial(MVariant::from_parts(ac, gene, posedit(linear(), edit)))
            }
            System::Coding => SV::Coding(CVariant::from_parts(
                ac,
                gene,
                posedit(transcript::<CVariant>(am()?, s, e)?, edit),
            )),
            System::NonCoding => SV::NonCoding(NVariant::from_parts(
                ac,
                gene,
                posedit(transcript::<NVariant>(am()?, s, e)?, edit),
            )),
            System::Rna => {
                let am = am()?;
                let coding = am.transcript.cds_start_index.is_some()
                    && am.transcript.cds_end_index.is_some();
                let spelled = if coding {
                    SV::Coding(CVariant::from_parts(
                        ac,
                        gene,
                        posedit(transcript::<CVariant>(am, s, e)?, edit),
                    ))
                } else {
                    SV::NonCoding(NVariant::from_parts(
                        ac,
                        gene,
                        posedit(transcript::<NVariant>(am, s, e)?, edit),
                    ))
                };
                SV::Rna(self.tx_to_r(&spelled)?)
            }
        })
    }

    /// A placed edit as a variant of the frame's system, an `r.` one in the
    /// `c.` or `n.` it is spelled from.
    fn transcript_spelling(&self, frame: &Frame, placed: &PlacedEdit) -> Result<SV, HgvsError> {
        match self.member_variant(frame, placed)? {
            SV::Rna(r) => self.r_as_transcript(&r),
            other => Ok(other),
        }
    }

    /// Where the allele is read, and each member placed there.
    fn frame(&self, cis: &CisPhasedVariant, on: On) -> Result<Frame, HgvsError> {
        let system = match cis.coordinate_type() {
            "g" => System::Genomic,
            "m" => System::Mitochondrial,
            "c" => System::Coding,
            "n" => System::NonCoding,
            "r" => System::Rna,
            other => {
                return Err(HgvsError::UnsupportedOperation(format!(
                    "{cis} is an allele of {other}. members, which is not read into a sequence"
                )))
            }
        };
        if matches!(system, System::Genomic | System::Mitochondrial) {
            return self.linear_frame(cis, system);
        }
        // r. members are read as the c. or n. members they are spelled from.
        let members: Vec<SV> = cis
            .members
            .iter()
            .map(|m| match m {
                SV::Rna(r) => self.r_as_transcript(r),
                other => Ok(other.clone()),
            })
            .collect::<Result<_, _>>()?;
        let posedits: Vec<&PosEdit<BaseOffsetInterval, NaEdit>> = members
            .iter()
            .map(|m| match m {
                SV::Coding(c) => Ok(&c.posedit),
                SV::NonCoding(n) => Ok(&n.posedit),
                other => Err(HgvsError::UnsupportedOperation(format!(
                    "{other} is not a transcript variant"
                ))),
            })
            .collect::<Result<_, _>>()?;
        let intronic = posedits.iter().any(|p| {
            p.pos
                .as_ref()
                .is_some_and(crate::mapper::has_intronic_offset)
        });
        if intronic || on == On::Genome {
            return self.genomic_frame(cis, &members);
        }
        let transcript = self.provider().get_transcript(&cis.ac, None)?;
        let am = TranscriptMapper::new(transcript)?;
        let reference = self
            .refs
            .reference(&cis.ac, IdentifierType::TranscriptAccession);
        let mut placed = Vec::new();
        for (index, posedit) in posedits.iter().enumerate() {
            let mut edit = None;
            let place = match &posedit.pos {
                None => Placed::Unapplied(Reach::Stretch(0, usize::MAX)),
                Some(pos) => {
                    let (n_start, n_end) = am.interval_to_n(pos)?;
                    if n_start.0 < 0 {
                        return Err(HgvsError::ValidationError(format!(
                            "{} names a position before the first base of {}",
                            cis.members[index], cis.ac
                        )));
                    }
                    let (start, end) = (n_start.0 as usize, n_end.0 as usize);
                    let uncertain = !pos.is_exact() || posedit.uncertain;
                    edit = Some(PlacedEdit::from_hgvs_range(
                        start,
                        end,
                        posedit.edit.clone(),
                    ));
                    place(
                        &posedit.edit,
                        start,
                        end,
                        uncertain,
                        &reference,
                        &cis.members[index],
                    )?
                }
            };
            let edit = edit.filter(|_| matches!(place, Placed::Applied(_)));
            placed.push((index, place, cis.members[index].clone(), edit));
        }
        Ok(Frame {
            accession: cis.ac.clone(),
            gene: cis.gene.clone(),
            kind: IdentifierType::TranscriptAccession,
            system,
            exons: am.exon_bounds()?,
            mapper: Some(am),
            placed,
        })
    }

    fn linear_frame(&self, cis: &CisPhasedVariant, system: System) -> Result<Frame, HgvsError> {
        let reference = self
            .refs
            .reference(&cis.ac, IdentifierType::GenomicAccession);
        let mut placed = Vec::new();
        for (index, m) in cis.members.iter().enumerate() {
            let posedit = match m {
                SV::Genomic(g) => &g.posedit,
                SV::Mitochondrial(v) => &v.posedit,
                other => {
                    return Err(HgvsError::UnsupportedOperation(format!(
                        "{other} is not a g. or m. variant"
                    )))
                }
            };
            let (place, edit) = place_linear(posedit, &reference, m)?;
            placed.push((index, place, m.clone(), edit));
        }
        Ok(Frame {
            accession: cis.ac.clone(),
            gene: cis.gene.clone(),
            kind: IdentifierType::GenomicAccession,
            system,
            exons: ExonBounds::unspliced(),
            mapper: None,
            placed,
        })
    }

    /// A transcript allele read on the genome: each member, spelled `members`
    /// in `c.` or `n.`, projected on its own. On the minus strand a
    /// duplication's copy goes 3' of the original on the transcript, which is
    /// 5' of it on the genome, so it is placed at the junction before the
    /// bases it copies there.
    fn genomic_frame(&self, cis: &CisPhasedVariant, members: &[SV]) -> Result<Frame, HgvsError> {
        let transcript = self.provider().get_transcript(&cis.ac, None)?;
        let minus = transcript.strand == Strand::Minus;
        let genome = transcript.reference_accession.clone();
        let reference = self
            .refs
            .reference(&genome, IdentifierType::GenomicAccession);
        let mut placed = Vec::new();
        for (index, m) in members.iter().enumerate() {
            let projected = match m {
                SV::Coding(c) => self.c_to_g(c, None),
                SV::NonCoding(n) => self.n_to_g(n, None),
                other => Err(HgvsError::UnsupportedOperation(format!(
                    "{other} is not a transcript variant"
                ))),
            }
            .map_err(|e| HgvsError::CisMember {
                allele: cis.to_string(),
                index,
                member: cis.members[index].to_string(),
                source: Box::new(e),
            })?;
            if let (NaEdit::Ins { .. } | NaEdit::InsLength { .. }, Some(pos)) =
                (&projected.posedit.edit, &projected.posedit.pos)
            {
                let (s, e) = crate::mapper::simple_interval_range(pos)?;
                if e != s + 2 {
                    return Err(HgvsError::UnsupportedOperation(format!(
                        "{} inserts on an exon junction, which has no one place on the genome",
                        cis.members[index]
                    )));
                }
            }
            let (place, edit) = place_linear(&projected.posedit, &reference, &cis.members[index])?;
            let place = match place {
                Placed::Applied(Act::Dup { start, end, .. }) if minus => {
                    Placed::Applied(Act::Dup {
                        start,
                        end,
                        at: start,
                    })
                }
                other => other,
            };
            placed.push((index, place, SV::Genomic(projected), edit));
        }
        Ok(Frame {
            accession: genome,
            gene: cis.gene.clone(),
            kind: IdentifierType::GenomicAccession,
            system: System::Genomic,
            exons: ExonBounds::unspliced(),
            mapper: None,
            placed,
        })
    }
}

/// The members of `cis` placed in `frame` on `reference`, and those no
/// sequence can be built from as written. Identical members other than
/// insertions and duplications are one edit written twice, and are kept once.
fn allele_on<'x, 's, 'a>(
    cis: &CisPhasedVariant,
    frame: &Frame,
    reference: &'x Reference<'s, 'a>,
) -> (Allele<'x, 's, 'a>, Vec<SV>) {
    let text: Vec<String> = cis.members.iter().map(ToString::to_string).collect();
    let mut members: Vec<(usize, Act)> = Vec::new();
    let mut outer: Vec<Outer> = Vec::new();
    let mut unapplied: Vec<SV> = Vec::new();
    for (index, placed, _, _) in &frame.placed {
        let member = &cis.members[*index];
        match placed {
            Placed::Applied(act) => {
                let twice = !matches!(act, Act::Insert { .. } | Act::Dup { .. })
                    && members.iter().any(|(_, a)| a == act);
                if !twice {
                    members.push((*index, act.clone()));
                }
            }
            Placed::Unapplied(reach) => {
                let inserts = member_edit(member).is_some_and(|e| {
                    matches!(
                        e,
                        NaEdit::Ins { .. } | NaEdit::InsLength { .. } | NaEdit::Dup { .. }
                    )
                });
                let twice = !inserts && unapplied.contains(member);
                if !twice {
                    outer.push(Outer {
                        index: *index,
                        reach: *reach,
                    });
                    unapplied.push(member.clone());
                }
            }
        }
    }
    // In order of position, so that nothing read from the allele depends on
    // the order its members are written in; each keeps its index for errors.
    members.sort_by(|(_, a), (_, b)| {
        (a.written(), a.junction(), a).cmp(&(b.written(), b.junction(), b))
    });
    outer.sort_by_key(|o| (o.reach.bases(), o.reach.junctions()));
    let allele = Allele {
        reference,
        description: cis.to_string(),
        text,
        members,
        outer,
    };
    (allele, unapplied)
}

/// A `g.` or `m.` member placed on its genome, and, where it is applied, its
/// edit over the HGVS range it is written over.
fn place_linear(
    posedit: &PosEdit<SimpleInterval, NaEdit>,
    reference: &Reference<'_, '_>,
    member: &SV,
) -> Result<(Placed, Option<PlacedEdit>), HgvsError> {
    let Some(pos) = &posedit.pos else {
        return Ok((Placed::Unapplied(Reach::Stretch(0, usize::MAX)), None));
    };
    let uncertain = !pos.is_exact() || posedit.uncertain;
    if uncertain {
        let index =
            |b: crate::structs::GenomicBound| b.known().map(|p| p.to_index().0.max(0) as usize);
        let last = pos.end.as_ref().unwrap_or(&pos.start);
        let lo = index(pos.start.base).unwrap_or(0);
        let hi = match last.end {
            Some(b) => index(b).map_or(usize::MAX, |i| i + 1),
            None => index(last.base).map_or(usize::MAX, |i| i + 1),
        };
        return Ok((Placed::Unapplied(Reach::of(&posedit.edit, lo, hi)), None));
    }
    let (start, end) = crate::mapper::simple_interval_range(pos)?;
    let placed = place(&posedit.edit, start, end, false, reference, member)?;
    let edit = matches!(placed, Placed::Applied(_))
        .then(|| PlacedEdit::from_hgvs_range(start, end, posedit.edit.clone()));
    Ok((placed, edit))
}

/// Whether a bound of a `g.` or `m.` member is unknown or known only to lie
/// in a range, so that it names no bases to check.
fn has_ranged_bound(member: &SV) -> bool {
    let ranged = |p: &SimplePosition| p.end.is_some() || p.base.known().is_none();
    let linear = |pos: &Option<SimpleInterval>| {
        pos.as_ref()
            .is_some_and(|p| ranged(&p.start) || p.end.as_ref().is_some_and(ranged))
    };
    match member {
        SV::Genomic(g) => linear(&g.posedit.pos),
        SV::Mitochondrial(m) => linear(&m.posedit.pos),
        _ => false,
    }
}

/// A nucleotide member's edit.
fn member_edit(member: &SV) -> Option<&NaEdit> {
    match member {
        SV::Genomic(g) => Some(&g.posedit.edit),
        SV::Mitochondrial(m) => Some(&m.posedit.edit),
        SV::Coding(c) => Some(&c.posedit.edit),
        SV::NonCoding(n) => Some(&n.posedit.edit),
        SV::Rna(r) => Some(&r.posedit.edit),
        SV::Protein(_) | SV::CisPhased(_) => None,
    }
}

/// Whether a nucleotide member's edit states bases the reference must hold.
fn states_bases(member: &SV) -> bool {
    member_edit(member).is_some_and(NaEdit::states_bases)
}

/// Whether `s` is bases the reading can apply: no `N`, which names a base
/// without saying which.
fn is_bases(s: &str) -> bool {
    !s.is_empty() && s.bytes().all(|b| b"ACGTacgt".contains(&b))
}

/// A member written over the HGVS range `[start, end)` placed as what it
/// replaces there, or as a member that cannot be applied.
fn place(
    edit: &NaEdit,
    start: usize,
    end: usize,
    uncertain: bool,
    reference: &Reference<'_, '_>,
    member: &SV,
) -> Result<Placed, HgvsError> {
    let insertion = matches!(edit, NaEdit::Ins { .. } | NaEdit::InsLength { .. });
    if insertion && end != start + 2 {
        return Err(HgvsError::ValidationError(format!(
            "{member} inserts between two bases that are not side by side on {}",
            reference.accession()
        )));
    }
    let unapplied = Placed::Unapplied(Reach::of(edit, start, end));
    if uncertain {
        return Ok(unapplied);
    }
    if end == 0 || reference.base(end - 1)?.is_none() {
        return Err(HgvsError::ValidationError(format!(
            "{member} names bases past the end of {}",
            reference.accession()
        )));
    }
    let write = |value| Placed::Applied(Act::Write { start, end, value });
    let replace = |start: usize, end: usize, alt: String| {
        Placed::Applied(if alt.is_empty() {
            Act::Write {
                start,
                end,
                value: Value::Delete,
            }
        } else if alt.len() == end - start {
            Act::Write {
                start,
                end,
                value: Value::Bases(alt),
            }
        } else {
            Act::Block { start, end, alt }
        })
    };
    Ok(match edit {
        NaEdit::RefAlt {
            uncertain: true, ..
        }
        | NaEdit::Del {
            uncertain: true, ..
        }
        | NaEdit::Ins {
            uncertain: true, ..
        }
        | NaEdit::Dup {
            uncertain: true, ..
        }
        | NaEdit::Inv {
            uncertain: true, ..
        }
        | NaEdit::Repeat {
            uncertain: true, ..
        } => unapplied,
        NaEdit::RefAlt { alt: None, .. } => write(Value::Keep),
        _ if edit.is_identity() => write(Value::Keep),
        NaEdit::RefAlt { alt: Some(a), .. } if a.is_empty() || is_bases(a) => {
            replace(start, end, a.to_ascii_uppercase())
        }
        NaEdit::Del { .. } => write(Value::Delete),
        NaEdit::Ins { alt: Some(a), .. } if is_bases(a) => Placed::Applied(Act::Insert {
            at: start + 1,
            alt: a.to_ascii_uppercase(),
        }),
        NaEdit::Dup { .. } => Placed::Applied(Act::Dup {
            start,
            end,
            at: end,
        }),
        NaEdit::Inv { .. } => Placed::Applied(Act::Inv { start, end }),
        NaEdit::Repeat { min, max, .. } if min == max => {
            let resolved = edit.resolve(reference, start, end)?;
            replace(
                resolved.start,
                resolved.end,
                resolved.alt.to_ascii_uppercase(),
            )
        }
        _ => unapplied,
    })
}
