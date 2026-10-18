//! Normalising an allele in cis: the members that meet written as one change.
//!
//! The rule is choices.md's "Members of a cis allele that meet are written as
//! one change". Statements are set aside first, by reading the allele; the
//! changes left are then merged, round by round, wherever their extents meet,
//! until no two do. An extent is the bases a change is written over together
//! with the region it could equally be written anywhere in, cut where an exon
//! ends or a statement stands; on `c.` two changes in one codon meet too.

use super::{
    allele_on, materialise, merge_ranges, overlap, pieces, same_replacements, strictly_inside, Act,
    Allele, Choice, Frame, Mode, On, Placed, Reach, Replacement, Segment, Settled, System,
};
use crate::data::IdentifierType;
use crate::edits::NaEdit;
use crate::error::HgvsError;
use crate::mapper::VariantMapper;
use crate::normalize::{self, ExonBounds, PlacedEdit};
use crate::reference::Reference;
use crate::structs::{BaseOffsetInterval, CisPhasedVariant, Variant};
use crate::SequenceVariant as SV;

/// What a reading is normalised for.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Writing {
    /// The description `normalize` returns, which refuses to write members
    /// that share bases across an exon junction as one change.
    Normalised,
    /// The statements a reading makes, or a reading an ambiguity lists, which
    /// write such members as the change across the junction they make.
    Reading,
}

/// A 0-based half-open stretch of the sequence the allele is read on.
type Span = (usize, usize);

/// Where a transcript position falls, in transcript order.
type OrderKey = (u8, i32, i32);

/// One reading normalised: the bases it states unchanged, sorted, and the
/// changes it makes.
pub(super) struct NormalForm {
    pub(super) stated: Vec<Span>,
    changes: Vec<PlacedEdit>,
}

/// What a change does not reach across: the bases stated unchanged, the
/// bases and junctions members no sequence is built from reach for, and the
/// exons, as stretches no shift leaves.
struct Barriers {
    /// What they are built from: the statements, and the windows of units
    /// written across a barrier.
    stated: Vec<(usize, usize)>,
    fixed: Vec<(usize, usize)>,
    ranges: Vec<(usize, usize)>,
    junctions: Vec<usize>,
    /// The exons with a cut at each barrier's edges.
    bounds: ExonBounds,
}

impl Barriers {
    /// The barriers among `units`: the statements, the reach of the members
    /// no sequence is built from, and every unit written across one of
    /// those, which is left as written and so stops a shift as they do.
    fn among(
        exons: &ExonBounds,
        stated: &[(usize, usize)],
        frame: &Frame,
        units: &[Unit],
    ) -> Barriers {
        let base = Barriers::new(exons, stated, frame, &[]);
        let mut fixed = fixed_among(units, &base);
        // A member written across an exon junction meets nothing, and no
        // other shifts onto it.
        let spans = |(a, b): Span| a < b && exons.exon(a) != exons.exon(b - 1);
        fixed.extend(units.iter().filter_map(|u| match &u.form {
            Form::Spanning(p) => Some((p.start, p.end)),
            Form::Member(_) if spans(u.window) => Some(u.window),
            _ => None,
        }));
        if fixed.is_empty() {
            return base;
        }
        Barriers::new(exons, stated, frame, &fixed)
    }

    fn new(
        exons: &ExonBounds,
        stated: &[(usize, usize)],
        frame: &Frame,
        fixed: &[(usize, usize)],
    ) -> Barriers {
        let mut ranges = stated.to_vec();
        ranges.extend_from_slice(fixed);
        let mut junctions = Vec::new();
        for (_, placed, _, _) in &frame.placed {
            if let Placed::Unapplied(reach) = placed {
                ranges.extend(reach.bases());
                junctions.extend(reach.junctions().into_iter().filter(|&j| j != usize::MAX));
            }
        }
        // Not merged: two statements that touch leave a junction between
        // them where something may be inserted.
        ranges.sort_unstable();
        ranges.dedup();
        let cuts: Vec<usize> = ranges
            .iter()
            .flat_map(|&(s, e)| [s, e])
            .chain(junctions.iter().copied())
            .filter(|&c| c != usize::MAX)
            .collect();
        Barriers {
            stated: stated.to_vec(),
            fixed: fixed.to_vec(),
            bounds: exons.with_cuts(cuts),
            ranges,
            junctions,
        }
    }

    /// These barriers as a unit over `own` meets them: a unit written
    /// across a barrier is not one to itself, nor to a unit it is part of.
    fn beside(&self, exons: &ExonBounds, frame: &Frame, own: (usize, usize)) -> Barriers {
        let fixed: Vec<(usize, usize)> = self
            .fixed
            .iter()
            .copied()
            .filter(|&f| overlap(f, own).is_none())
            .collect();
        Barriers::new(exons, &self.stated, frame, &fixed)
    }

    /// Whether the base at `i` is a barrier's.
    fn holds(&self, i: usize) -> bool {
        self.ranges.iter().any(|&(s, e)| s <= i && i < e)
    }

    /// Whether a barrier lies in `[lo, hi)`: one of its bases, or a junction
    /// inside it.
    fn across(&self, span: (usize, usize)) -> bool {
        self.ranges.iter().any(|&r| overlap(r, span).is_some())
            || self.junctions.iter().any(|&j| strictly_inside(j, span))
    }

    /// Whether `p` changes a barrier's base or inserts between two of them,
    /// or inserts at, or reaches across, a junction where a member no
    /// sequence is built from may insert: the overlap rule refuses both.
    fn intrudes(&self, p: &PlacedEdit) -> bool {
        let inserts_at = match p.edit {
            NaEdit::Ins { .. } => Some(p.start),
            NaEdit::Dup { .. } => Some(p.end),
            _ => None,
        };
        let on_junction = self.junctions.iter().any(|&j| {
            inserts_at == Some(j) || (!p.is_insertion() && strictly_inside(j, (p.start, p.end)))
        });
        on_junction
            || self.ranges.iter().any(|&r| {
                if p.is_insertion() {
                    strictly_inside(p.start, r)
                } else {
                    overlap((p.start, p.end), r).is_some()
                }
            })
    }
}

/// Where the units written across one of the `base` barriers lie: a member
/// written over one, or a change that could not be moved off one. Each is left
/// as written, and stops a shift as the barrier does.
fn fixed_among(units: &[Unit], base: &Barriers) -> Vec<(usize, usize)> {
    units
        .iter()
        .filter_map(|u| match &u.form {
            Form::Change(p) | Form::Part(p) => base.intrudes(p).then_some((p.start, p.end)),
            _ => (u.window.0 < u.window.1 && base.across(u.window)).then_some(u.window),
        })
        .collect()
}

/// A change being normalised: one member, part of one, or members merged.
struct Unit {
    /// The slots, in the allele, of the members it holds.
    slots: Vec<usize>,
    /// The replacements that make it.
    reps: Vec<Replacement>,
    /// The stretch its members are written over, holding `reps`.
    window: (usize, usize),
    /// The positions its members name, which it states if it makes the
    /// reference.
    named: Vec<(usize, usize)>,
    form: Form,
}

enum Form {
    /// A member as written, normalised alone if it meets no other.
    Member(usize),
    /// The part of a member between statements, whose change is not yet read.
    Piece,
    /// The part of a member between statements that changes the bases
    /// there, normalised. Its extent holds the bases it is written over.
    Part(PlacedEdit),
    /// What the replacements make over the window, normalised.
    Change(PlacedEdit),
    /// A member written across an exon junction, normalised alone: it meets
    /// nothing, and others stop where it ends up.
    Spanning(PlacedEdit),
}

/// Where a unit reaches.
struct Extent {
    /// The bases it is written over and could equally be written over, in
    /// the stretch it is in.
    extent: (usize, usize),
    /// An insertion on an exon junction, which names no one place and meets
    /// nothing.
    pinned: bool,
    /// The first and last exon of its window.
    exons: (usize, usize),
    /// The first and last codon its extent shares, on `c.` inside the CDS.
    codons: Option<(usize, usize)>,
    /// The junctions its members insert at. Two units that insert at one
    /// junction share the replacement there, and so are one change.
    inserts_at: Vec<usize>,
    /// The bases its members are written over: two units written over one
    /// base are one change too.
    window: (usize, usize),
}

/// What members merged make.
enum Merged {
    /// The reference: the positions the members name, stated unchanged.
    Statement(Vec<(usize, usize)>),
    Change(Unit),
}

/// An allele being normalised in its frame.
struct Context<'f, 'x, 's, 'a> {
    frame: &'f Frame,
    allele: &'f Allele<'x, 's, 'a>,
    /// How the allele's members are read.
    choice: &'f Choice,
    writing: Writing,
    reference: &'x Reference<'s, 'a>,
    /// Each member's extent before cutting: written over or justified over.
    extents: Vec<(usize, usize)>,
    /// The CDS, `[start, end)`, where codons group members.
    cds: Option<(usize, usize)>,
}

impl Context<'_, '_, '_, '_> {
    fn slice(&self, (s, e): (usize, usize)) -> Result<String, HgvsError> {
        Segment::Reference { start: s, end: e }.read(self.reference, 0, e - s)
    }

    /// The units the changing members start as: each member, or each part
    /// of a member that agrees with an identity, between the statements;
    /// and the bases named by each member cut there that makes the reference.
    fn units(
        &self,
        reps: Vec<Replacement>,
        stated: &[(usize, usize)],
        barriers: &Barriers,
    ) -> Result<(Vec<Unit>, Vec<Span>), HgvsError> {
        let (members, choice) = (&self.allele.members, self.choice);
        // A member that makes the reference on its own changes nothing
        // another could read together with it.
        let mut inert = Vec::with_capacity(members.len());
        for (j, (_, act)) in members.iter().enumerate() {
            inert.push(
                matches!(act, Act::Write { .. } | Act::Inv { .. })
                    && !act.is_identity()
                    && self.allele.cancels(&[j])?,
            );
        }
        let mut restated = Vec::new();
        for k in 0..members.len() {
            for j in 0..members.len() {
                if j != k && self.agrees_within(k, j)? {
                    restated.push((k, j));
                }
            }
        }
        let exons = &self.frame.exons;
        // A member is cut where an identity it agrees with stands, unless
        // another change inserts inside it or shares a base with it: the
        // identity is then part of how the two read together.
        let cuttable = |k: usize, b: (usize, usize)| {
            // Across an exon junction it names the intron on the genome,
            // which no part of it would.
            exons.exon(b.0) == exons.exon(b.1 - 1)
        } && {
            matches!(members[k].1, Act::Write { .. })
                || matches!(members[k].1, Act::Inv { .. }) && choice.mode[k] == Mode::Reference
        } && members.iter().enumerate().all(|(j, (_, other))| {
            j == k
                || choice.dropped[j]
                || inert[j]
                || restated.contains(&(k, j))
                || !(other.junction().is_some_and(|at| strictly_inside(at, b))
                    || !other.is_identity() && other.bases().and_then(|o| overlap(o, b)).is_some())
        });
        let cut_members: Vec<(usize, (usize, usize))> = members
            .iter()
            .enumerate()
            .filter_map(|(k, (_, act))| {
                let b = act.bases()?;
                let crossed = !act.is_identity()
                    && stated.iter().any(|&s| overlap(s, b).is_some())
                    && cuttable(k, b);
                crossed.then_some((k, b))
            })
            .collect();
        let cut_ranges: Vec<(usize, usize)> = cut_members.iter().map(|&(_, b)| b).collect();
        let mut cut = Vec::new();
        for r in reps {
            if cut_ranges
                .iter()
                .any(|&(s, e)| s <= r.start && r.end <= e && r.start < r.end)
            {
                cut.extend(cut_at(r, stated)?);
            } else {
                cut.push(r);
            }
        }
        // What lies inside a statement, of a member cut there, agrees with it.
        cut.retain(|r| {
            r.start == r.end
                || !(stated.iter().any(|&(s, e)| s <= r.start && r.end <= e)
                    && cut_ranges.iter().any(|&(s, e)| s <= r.start && r.end <= e))
        });
        let mut units = Vec::new();
        for (k, (_, act)) in members.iter().enumerate() {
            // A member that changes nothing inside what is stated adds
            // nothing to it.
            let stated_whole = act.bases().is_some_and(|b| outside(b, stated).is_empty());
            if act.is_identity() || choice.dropped[k] || inert[k] && stated_whole {
                continue;
            }
            if act.is_reader() && choice.mode[k] == Mode::Edited {
                // What it reads is what the other members leave, so it is
                // written as what it makes: a copy inserted at its junction,
                // or the stretch it turns round.
                let window = match act {
                    Act::Dup { at, .. } => (*at, *at),
                    other => other.written(),
                };
                units.push(Unit {
                    slots: vec![k],
                    reps: Vec::new(),
                    window,
                    named: vec![act.named()],
                    form: Form::Piece,
                });
                continue;
            }
            let crossed = cut_members.iter().find(|&&(j, _)| j == k).map(|&(_, b)| b);
            if let Some(b) = crossed {
                for piece in outside(b, stated) {
                    units.push(Unit {
                        slots: vec![k],
                        reps: Vec::new(),
                        window: piece,
                        named: vec![piece],
                        form: Form::Piece,
                    });
                }
            } else {
                units.push(Unit {
                    slots: vec![k],
                    reps: Vec::new(),
                    window: act.written(),
                    named: vec![act.named()],
                    form: Form::Member(k),
                });
            }
        }
        for r in cut {
            let owner = if r.start == r.end {
                units.iter().position(|u| {
                    !matches!(u.form, Form::Change(_))
                        && u.slots
                            .iter()
                            .any(|&k| members[k].1.junction() == Some(r.start))
                })
            } else {
                let holds =
                    |u: &Unit| u.window.0 <= r.start && r.end <= u.window.1 && self.replaces(u);
                // A part of a member cut at a statement takes what it makes
                // before a member that only restates it does.
                units
                    .iter()
                    .position(|u| matches!(u.form, Form::Piece) && holds(u))
                    .or_else(|| units.iter().position(holds))
            };
            let inside_statement =
                r.start < r.end && stated.iter().any(|&(s, e)| s <= r.start && r.end <= e);
            let Some(i) = owner else {
                if inside_statement {
                    // Made by a member that changes nothing there.
                    continue;
                }
                return Err(HgvsError::Other(format!(
                    "{}: no member makes the replacement over [{}, {})",
                    self.allele.description, r.start, r.end
                )));
            };
            units[i].reps.push(r);
        }
        // Each piece with the bases it makes, or `None` where it makes the
        // reference's.
        let mut made: Vec<Option<String>> = Vec::with_capacity(units.len());
        for u in &units {
            made.push(match u.form {
                Form::Piece => {
                    let alt = materialise(
                        &pieces(u.reps.iter(), u.window.0, u.window.1),
                        self.reference,
                    )?;
                    (!alt.eq_ignore_ascii_case(&self.slice(u.window)?)).then_some(alt)
                }
                _ => None,
            });
        }
        // A member cut at a statement whose every part makes the reference
        // makes the reference, and states what it names; a part of one that
        // also changes something, and makes the reference itself, is dropped.
        let is_cut = |k: usize| cut_members.iter().any(|&(j, _)| j == k);
        let mut changing = vec![false; members.len()];
        for (u, alt) in units.iter().zip(&made) {
            if alt.is_some() {
                changing[u.slots[0]] = true;
            }
        }
        let mut statements = Vec::new();
        let mut out = Vec::with_capacity(units.len());
        for (mut u, alt) in units.into_iter().zip(made) {
            if !matches!(u.form, Form::Piece) {
                out.push(u);
                continue;
            }
            let k = u.slots[0];
            match alt {
                Some(alt) => {
                    u.form = Form::Part(self.change_over(u.window, alt, barriers)?);
                    out.push(u);
                }
                // A reader of edited bases that makes the reference there
                // may still take in another member's change: it is grouped.
                None if !is_cut(k) => {
                    let unchanged = self.slice(u.window)?;
                    u.form = Form::Part(self.change_over(u.window, unchanged, barriers)?);
                    out.push(u);
                }
                None if !changing[k] => statements.push(members[k].1.named()),
                None => {}
            }
        }
        statements.sort_unstable();
        statements.dedup();
        Ok((out, statements))
    }

    /// Whether member `j` writes, base by base, within the bases member `k`
    /// writes, and what it writes there is what `k` does: one edit `k` states
    /// again.
    fn agrees_within(&self, k: usize, j: usize) -> Result<bool, HgvsError> {
        let members = &self.allele.members;
        match (&members[k].1, &members[j].1) {
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
            ) if s1 <= s2 && e2 <= e1 => self.allele.writes_agree((v1, *s1), (v2, *s2), (*s2, *e2)),
            _ => Ok(false),
        }
    }

    /// The member at `slot` as written, over the HGVS range it is written
    /// over.
    fn written_edit(&self, slot: usize) -> Result<PlacedEdit, HgvsError> {
        let index = self.allele.members[slot].0;
        self.frame
            .placed
            .iter()
            .find(|(i, ..)| *i == index)
            .and_then(|(.., edit)| edit.clone())
            .ok_or_else(|| {
                HgvsError::Other(format!(
                    "{}: member {} has no placed edit",
                    self.allele.description,
                    index + 1
                ))
            })
    }

    /// Whether the unit replaces bases, rather than inserting at a junction.
    fn replaces(&self, u: &Unit) -> bool {
        match u.form {
            Form::Member(k) => matches!(
                self.allele.members[k].1,
                Act::Write { .. } | Act::Block { .. } | Act::Inv { .. }
            ),
            Form::Piece | Form::Part(_) | Form::Change(_) | Form::Spanning(_) => true,
        }
    }

    fn extent(&self, u: &Unit, barriers: &Barriers) -> Result<Extent, HgvsError> {
        let (written, reach) = match &u.form {
            Form::Member(k) => (self.allele.members[*k].1.written(), self.extents[*k]),
            Form::Change(p) | Form::Spanning(p) => self.justified(p)?,
            Form::Part(p) => {
                let (w, j) = self.justified(p)?;
                let o = u.window;
                ((w.0.min(o.0), w.1.max(o.1)), (j.0.min(o.0), j.1.max(o.1)))
            }
            Form::Piece => unreachable!("a piece's change is read when it is made"),
        };
        let (extent, pinned) = self.clip(written, reach, barriers);
        let exons = &self.frame.exons;
        let w = u.window;
        let exons = if w.0 < w.1 {
            (exons.exon(w.0), exons.exon(w.1 - 1))
        } else {
            let e = exons.exon(w.0.saturating_sub(1));
            (e, e)
        };
        let inserts_at = u
            .slots
            .iter()
            .filter(|&&k| !self.choice.dropped[k])
            .filter_map(|&k| self.allele.members[k].1.junction())
            .collect();
        Ok(Extent {
            extent,
            pinned,
            exons,
            codons: self.codons(extent),
            inserts_at,
            window: u.window,
        })
    }

    /// The bases `p` is written over, and those together with every place
    /// it could equally be written.
    fn justified(&self, p: &PlacedEdit) -> Result<(Span, Span), HgvsError> {
        let w = if p.is_insertion() {
            (p.start, p.start)
        } else {
            (p.start, p.end)
        };
        let (s, e) = p.hgvs_range();
        let resolved = p.edit.resolve(self.reference, s, e)?;
        let j = normalize::ambiguous_range(self.reference, &resolved)?;
        Ok((w, (w.0.min(j.0), w.1.max(j.1))))
    }

    /// `reach`, around `written`, cut to the stretch between barriers and
    /// exon junctions that `written` is in; and whether `written` is an
    /// insertion on an exon junction.
    fn clip(
        &self,
        written: (usize, usize),
        reach: (usize, usize),
        barriers: &Barriers,
    ) -> ((usize, usize), bool) {
        let bounds = &barriers.bounds;
        let (lo, hi) = if written.0 < written.1 {
            (bounds.exon_start(written.0), bounds.exon_end(written.1 - 1))
        } else {
            let j = written.0;
            if self.frame.exons.is_junction(j) {
                return ((j, j), true);
            }
            if j == 0 {
                (0, bounds.exon_end(0))
            } else if !bounds.is_junction(j) {
                (bounds.exon_start(j - 1), bounds.exon_end(j))
            } else {
                // Beside a barrier: on the side away from it.
                let lo = if barriers.holds(j - 1) {
                    j
                } else {
                    bounds.exon_start(j - 1)
                };
                let hi = if barriers.holds(j) {
                    j
                } else {
                    bounds.exon_end(j)
                };
                (lo, hi)
            }
        };
        (
            (
                reach.0.max(lo).min(written.0),
                reach.1.min(hi).max(written.1),
            ),
            false,
        )
    }

    /// The first and last codon of the CDS that `[lo, hi)` shares bases
    /// with; for an empty extent, the codon whose bases it lies between.
    fn codons(&self, (lo, hi): (usize, usize)) -> Option<(usize, usize)> {
        let (s, e) = self.cds?;
        if lo < hi {
            let (a, b) = (lo.max(s), hi.min(e));
            (a < b).then(|| ((a - s) / 3, (b - 1 - s) / 3))
        } else {
            (s < lo && lo < e && (lo - s) % 3 != 0).then(|| ((lo - s) / 3, (lo - s) / 3))
        }
    }

    /// Whether two units meet: in one exon, with no barrier between, their
    /// extents overlapping or touching, or sharing a codon.
    fn meet(&self, a: &Extent, b: &Extent, barriers: &Barriers) -> bool {
        // A unit written across an exon junction lies in two exons, and
        // meets nothing: a change merged with it would say something else on
        // the genome.
        let inserts_inside =
            |x: &Extent, y: &Extent| x.inserts_at.iter().any(|&j| strictly_inside(j, y.window));
        if a.inserts_at.iter().any(|j| b.inserts_at.contains(j))
            || overlap(a.window, b.window).is_some()
            || inserts_inside(a, b)
            || inserts_inside(b, a)
        {
            return true;
        }
        let one_exon = a.exons.0 == a.exons.1 && a.exons == b.exons;
        if a.pinned || b.pinned || !one_exon {
            return false;
        }
        let hull = (a.extent.0.min(b.extent.0), a.extent.1.max(b.extent.1));
        if barriers.across(hull) {
            return false;
        }
        let touch = a.extent.0 <= b.extent.1 && b.extent.0 <= a.extent.1;
        touch || matches!((a.codons, b.codons), (Some(x), Some(y)) if x.0 <= y.1 && y.0 <= x.1)
    }

    fn merge(&self, parts: Vec<Unit>, barriers: &Barriers) -> Result<Merged, HgvsError> {
        let mut slots = Vec::new();
        let mut reps = Vec::new();
        let mut named = Vec::new();
        let mut window = (usize::MAX, 0);
        for u in parts {
            slots.extend(u.slots);
            reps.extend(u.reps);
            named.extend(u.named);
            window = (window.0.min(u.window.0), window.1.max(u.window.1));
        }
        let exons = &self.frame.exons;
        let across = window.0 < window.1 && exons.exon(window.0) != exons.exon(window.1 - 1);
        if across && self.writing == Writing::Normalised {
            // Members that share bases or a junction are one change, and a
            // change across an exon junction projects over the intron, which
            // none of them did.
            return Err(HgvsError::UnsupportedOperation(format!(
                "{}: members that share bases across an exon junction are not merged, \
                 since one change there would delete the intron",
                self.allele.description
            )));
        }
        reps.sort_by_key(|r| (r.start, r.end));
        let alt = materialise(&pieces(reps.iter(), window.0, window.1), self.reference)?;
        if alt.eq_ignore_ascii_case(&self.slice(window)?) {
            return Ok(Merged::Statement(named));
        }
        let change = self.change_over(
            window,
            alt,
            &barriers.beside(&self.frame.exons, self.frame, window),
        )?;
        Ok(Merged::Change(Unit {
            slots,
            reps,
            window,
            named,
            form: Form::Change(change),
        }))
    }

    /// The change that puts `alt` in place of `window`, normalised.
    fn change_over(
        &self,
        (lo, hi): (usize, usize),
        alt: String,
        barriers: &Barriers,
    ) -> Result<PlacedEdit, HgvsError> {
        let exons = &self.frame.exons;
        if lo == hi || exons.exon(lo) == exons.exon(hi - 1) {
            // Trimmed to what it alters before it is shifted, so that a
            // statement's edge does not hold it as an exon junction would.
            let ref_ = self.slice((lo, hi))?;
            let (edit, s, e) = normalize::minimal_edit(self.reference, lo, hi, &ref_, &alt)?;
            let settled = self.settle_change(PlacedEdit::from_hgvs_range(s, e, edit), barriers)?;
            // An insertion on a junction names no one place on the genome,
            // and one trimmed onto a barrier it cannot leave is refused by
            // the overlap rule; reduced as a delins, either keeps the base
            // beside it.
            let on_junction = settled.is_insertion()
                && (exons.is_junction(settled.start)
                    || self.reference.base(settled.start)?.is_none());
            if !on_junction && !barriers.intrudes(&settled) {
                return Ok(settled);
            }
        }
        // Across an exon junction: reduced keeping a base in each exon.
        let uncertain = false;
        let placed = if lo == hi && lo == 0 {
            // HGVS has no insertion before the first base.
            let first = self.slice((0, 1))?;
            PlacedEdit {
                start: 0,
                end: 1,
                edit: NaEdit::RefAlt {
                    ref_: Some(String::new()),
                    alt: Some(format!("{alt}{first}")),
                    uncertain,
                },
            }
        } else if lo == hi && self.reference.base(lo)?.is_none() {
            // Nor one after the last.
            let last = self.slice((lo - 1, lo))?;
            PlacedEdit {
                start: lo - 1,
                end: lo,
                edit: NaEdit::RefAlt {
                    ref_: Some(String::new()),
                    alt: Some(format!("{last}{alt}")),
                    uncertain,
                },
            }
        } else if lo == hi {
            PlacedEdit {
                start: lo,
                end: lo,
                edit: NaEdit::Ins {
                    alt: Some(alt),
                    uncertain,
                },
            }
        } else if alt.is_empty() {
            PlacedEdit {
                start: lo,
                end: hi,
                edit: NaEdit::Del {
                    ref_: None,
                    uncertain,
                },
            }
        } else {
            PlacedEdit {
                start: lo,
                end: hi,
                edit: NaEdit::RefAlt {
                    ref_: Some(String::new()),
                    alt: Some(alt),
                    uncertain,
                },
            }
        };
        self.settle_change(placed, barriers)
    }

    /// `placed` normalised no further than a barrier's edge, or as written
    /// where it is written across a barrier, which it then cannot leave.
    fn settle_change(
        &self,
        placed: PlacedEdit,
        barriers: &Barriers,
    ) -> Result<PlacedEdit, HgvsError> {
        let normalised = normalize::normalize(self.reference, placed.clone(), &barriers.bounds)?;
        // An identity may overlap a statement: it states the same bases.
        let identity = matches!(normalised.edit, NaEdit::RefAlt { alt: None, .. });
        Ok(if !identity && barriers.intrudes(&normalised) {
            placed
        } else {
            normalised
        })
    }
}

/// `r` cut where a statement starts or ends inside it. Only a replacement
/// that keeps the length can agree with a statement it overlaps, so only
/// such a replacement is cut.
fn cut_at(r: Replacement, stated: &[(usize, usize)]) -> Result<Vec<Replacement>, HgvsError> {
    let mut cuts: Vec<usize> = stated
        .iter()
        .flat_map(|&(s, e)| [s, e])
        .filter(|&c| r.start < c && c < r.end)
        .collect();
    if cuts.is_empty() {
        return Ok(vec![r]);
    }
    if r.with.iter().map(Segment::len).sum::<usize>() != r.end - r.start {
        return Err(HgvsError::Other(format!(
            "a replacement over [{}, {}) that changes the length crosses a statement",
            r.start, r.end
        )));
    }
    cuts.sort_unstable();
    cuts.dedup();
    let mut out = Vec::new();
    let mut rest = r;
    for c in cuts {
        let (head, tail) = split_segments(&rest.with, c - rest.start);
        out.push(Replacement {
            start: rest.start,
            end: c,
            with: head,
        });
        rest = Replacement {
            start: c,
            end: rest.end,
            with: tail,
        };
    }
    out.push(rest);
    Ok(out)
}

/// `segments` split after their first `n` bases.
fn split_segments(segments: &[Segment], mut n: usize) -> (Vec<Segment>, Vec<Segment>) {
    let mut head = Vec::new();
    let mut tail = Vec::new();
    for s in segments {
        if n >= s.len() {
            n -= s.len();
            head.push(s.clone());
        } else if n == 0 {
            tail.push(s.clone());
        } else {
            let (a, b) = match s {
                Segment::Reference { start, end } => (
                    Segment::Reference {
                        start: *start,
                        end: start + n,
                    },
                    Segment::Reference {
                        start: start + n,
                        end: *end,
                    },
                ),
                // The first n bases of a reverse complement are those of the
                // last n reference bases it reads.
                Segment::ReverseComplement { start, end } => (
                    Segment::ReverseComplement {
                        start: end - n,
                        end: *end,
                    },
                    Segment::ReverseComplement {
                        start: *start,
                        end: end - n,
                    },
                ),
                Segment::Literal(x) => (
                    Segment::Literal(x[..n].to_string()),
                    Segment::Literal(x[n..].to_string()),
                ),
            };
            head.push(a);
            tail.push(b);
            n = 0;
        }
    }
    (head, tail)
}

/// The stretches of `(lo, hi)` no statement covers.
fn outside((lo, hi): (usize, usize), stated: &[(usize, usize)]) -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    let mut at = lo;
    for &(s, e) in stated {
        if e <= at || s >= hi {
            continue;
        }
        if s > at {
            out.push((at, s));
        }
        at = at.max(e);
    }
    if at < hi {
        out.push((at, hi));
    }
    out
}

/// The groups of `0..n` joined by `meet`, directly or through another, each
/// in increasing order, in order of their first.
fn components(n: usize, meet: impl Fn(usize, usize) -> bool) -> Vec<Vec<usize>> {
    let mut root: Vec<usize> = (0..n).collect();
    fn find(root: &mut [usize], mut i: usize) -> usize {
        while root[i] != i {
            root[i] = root[root[i]];
            i = root[i];
        }
        i
    }
    for a in 0..n {
        for b in a + 1..n {
            if meet(a, b) {
                let (ra, rb) = (find(&mut root, a), find(&mut root, b));
                if ra != rb {
                    root[ra.max(rb)] = ra.min(rb);
                }
            }
        }
    }
    let mut groups: Vec<Vec<usize>> = Vec::new();
    let mut of_root: Vec<Option<usize>> = vec![None; n];
    for i in 0..n {
        let r = find(&mut root, i);
        match of_root[r] {
            Some(g) => groups[g].push(i),
            None => {
                of_root[r] = Some(groups.len());
                groups.push(vec![i]);
            }
        }
    }
    groups
}

/// `ranges` sorted, with those that overlap or touch written as one where
/// the one lies in one exon. A range across an exon junction projects over
/// the intron, so two written apart there stay apart, and one written across
/// a junction is kept as written.
fn merge_statements(mut ranges: Vec<(usize, usize)>, exons: &ExonBounds) -> Vec<(usize, usize)> {
    ranges.sort_unstable();
    ranges.dedup();
    let mut out: Vec<(usize, usize)> = Vec::with_capacity(ranges.len());
    for r in ranges {
        match out.last_mut() {
            Some(last)
                if r.0 <= last.1 && exons.exon(last.0) == exons.exon(last.1.max(r.1) - 1) =>
            {
                last.1 = last.1.max(r.1);
            }
            _ => out.push(r),
        }
    }
    out
}

/// Where a member no sequence is built from sits, for ordering.
fn reach_key(reach: &Reach) -> (usize, usize) {
    match *reach {
        Reach::Stretch(lo, hi) | Reach::Copy(lo, hi) => (lo, hi),
        Reach::Junction(at) => (at.saturating_sub(1), at + 1),
    }
}

/// Whether the edit inserts, so that writing it twice is not one edit.
fn inserts(member: &SV) -> bool {
    super::member_edit(member).is_some_and(|e| {
        matches!(
            e,
            NaEdit::Ins { .. } | NaEdit::InsLength { .. } | NaEdit::Dup { .. }
        )
    })
}

/// Where a transcript member starts and ends, in transcript order.
fn transcript_key(member: &SV) -> Option<(OrderKey, OrderKey)> {
    let pos: &BaseOffsetInterval = match member {
        SV::Coding(c) => c.posedit.pos.as_ref()?,
        SV::NonCoding(n) => n.posedit.pos.as_ref()?,
        SV::Rna(r) => r.posedit.pos.as_ref()?,
        _ => return None,
    };
    let start = pos.start.order_key();
    Some((start, pos.end.as_ref().map_or(start, |e| e.order_key())))
}

impl VariantMapper<'_> {
    /// An allele in cis of nucleotide members, normalised: the members that
    /// meet written as the one change they make, statements as identities,
    /// each member that meets no other normalised as a single variant is,
    /// in order of position; the plain variant when one member is left.
    ///
    /// # Errors
    ///
    /// As [`read_cis`](Self::read_cis): a contradictory allele is a
    /// [`HgvsError::ValidationError`], an ambiguous one
    /// [`HgvsError::Ambiguous`] listing its readings normalised, and one whose
    /// readings are not listed [`HgvsError::UnsupportedOperation`]. So is a
    /// transcript allele with an intronic member that meets an exon or
    /// another member on the genome.
    pub(crate) fn normalize_cis(&self, cis: &CisPhasedVariant) -> Result<SV, HgvsError> {
        let frame = self.checked_frame(cis, On::Written)?;
        if matches!(cis.coordinate_type(), "c" | "n" | "r")
            && frame.kind == IdentifierType::GenomicAccession
        {
            return self.normalize_cis_with_introns(cis, &frame);
        }
        let reference = self.refs.reference(&frame.accession, frame.kind);
        let (allele, _) = allele_on(cis, &frame, &reference);
        let settled = self.settle(cis, &frame, &allele, true)?;
        let normalised =
            self.describe_statement(&frame, &allele, &settled.statement, Writing::Normalised)?;
        self.read_as(cis, &settled, &reference, normalised)
    }

    /// `normalised`, once it reads as the allele `cis`, whose reading on
    /// `reference` is `settled`, does: one reading, the same sequence, the
    /// same bases stated. Where it would not, the allele is refused rather
    /// than answered wrongly.
    fn read_as(
        &self,
        cis: &CisPhasedVariant,
        settled: &Settled,
        reference: &Reference<'_, '_>,
        normalised: SV,
    ) -> Result<SV, HgvsError> {
        let as_allele = match &normalised {
            SV::CisPhased(c) => c.clone(),
            one => CisPhasedVariant::new(cis.ac.clone(), cis.gene.clone(), vec![one.clone()])?,
        };
        let unsupported = |why: String| {
            HgvsError::UnsupportedOperation(format!(
                "{cis} is not yet normalised: {as_allele}, the form it would be written in, {why}"
            ))
        };
        let after = self
            .read_cis(&as_allele)
            .map_err(|e| unsupported(format!("does not read: {e}")))?;
        if !same_replacements(&settled.replacements, after.replacements(), reference)? {
            return Err(unsupported("makes another sequence".into()));
        }
        if merge_ranges(settled.statement.stated.clone()) != after.stated() {
            return Err(unsupported("states other bases".into()));
        }
        if matches!(cis.coordinate_type(), "c" | "n" | "r") {
            if let Ok(before) = self.read_cis_on_genome(cis) {
                let after = self
                    .read_cis_on_genome(&as_allele)
                    .map_err(|e| unsupported(format!("does not read on the genome: {e}")))?;
                let genome = self.refs.reference(before.accession(), before.kind());
                if !before.same_sequence(&after, &genome)? {
                    return Err(unsupported("makes another sequence on the genome".into()));
                }
            }
        }
        Ok(normalised)
    }

    /// A transcript allele with intronic members, read on the genome in
    /// `frame`. An intronic member and the members it meets on the genome are
    /// written as the one change they make there, in transcript coordinates
    /// with intronic offsets and not shifted, as a single intronic variant is
    /// not; an intronic member that meets nothing is carried as given; the
    /// other exonic members are normalised as an allele on the transcript.
    /// What that writes can reach an intronic member none of its members
    /// reached, so it repeats until nothing changes.
    fn normalize_cis_with_introns(
        &self,
        cis: &CisPhasedVariant,
        frame: &Frame,
    ) -> Result<SV, HgvsError> {
        let reference = self.refs.reference(&frame.accession, frame.kind);
        let (allele, _) = allele_on(cis, frame, &reference);
        let settled = self.settle(cis, frame, &allele, true)?;
        let mut current = cis.clone();
        for _ in 0..=cis.members.len() {
            let next = self.introns_once(&current)?;
            if next == current {
                let normalised = match <[SV; 1]>::try_from(next.members) {
                    Ok([one]) => one,
                    Err(members) => SV::CisPhased(CisPhasedVariant::new(
                        cis.ac.clone(),
                        cis.gene.clone(),
                        members,
                    )?),
                };
                return self.read_as(cis, &settled, &reference, normalised);
            }
            current = next;
        }
        Err(HgvsError::UnsupportedOperation(format!(
            "{cis} is not yet normalised: writing its intronic members with what they meet \
             does not settle"
        )))
    }

    /// One round of [`normalize_cis_with_introns`](Self::normalize_cis_with_introns).
    fn introns_once(&self, cis: &CisPhasedVariant) -> Result<CisPhasedVariant, HgvsError> {
        let frame = self.checked_frame(cis, On::Written)?;
        if frame.kind != IdentifierType::GenomicAccession {
            // No intronic member is left: the allele is a transcript one.
            return Ok(match self.normalize_cis(cis)? {
                SV::CisPhased(c) => c,
                one => CisPhasedVariant::new(cis.ac.clone(), cis.gene.clone(), vec![one])?,
            });
        }
        let genome = self.refs.reference(&frame.accession, frame.kind);
        let (allele, _) = allele_on(cis, &frame, &genome);
        let extents = allele.extents()?;
        let mut intronic = Vec::with_capacity(cis.members.len());
        for m in &cis.members {
            let spelled = match m {
                SV::Rna(r) => self.r_as_transcript(r)?,
                other => other.clone(),
            };
            let pos = match &spelled {
                SV::Coding(c) => c.posedit.pos.as_ref(),
                SV::NonCoding(n) => n.posedit.pos.as_ref(),
                _ => None,
            };
            intronic.push(pos.is_some_and(crate::mapper::has_intronic_offset));
        }
        // The members an intronic member meets on the genome, directly or
        // through another, as groups of slots.
        let n = allele.members.len();
        let meets = |a: Span, b: Span| a.0 <= b.1 && b.0 <= a.1;
        let groups = components(n, |a, b| {
            let (ia, ib) = (allele.members[a].0, allele.members[b].0);
            (intronic[ia] || intronic[ib]) && meets(extents[a], extents[b])
        });
        let mut grouped = vec![false; cis.members.len()];
        let mut members: Vec<SV> = Vec::new();
        for group in groups.iter().filter(|g| g.len() > 1) {
            let indices: Vec<usize> = group.iter().map(|&k| allele.members[k].0).collect();
            for &i in &indices {
                grouped[i] = true;
            }
            let window = group
                .iter()
                .map(|&k| allele.members[k].1.written())
                .reduce(|a, b| (a.0.min(b.0), a.1.max(b.1)))
                .expect("a group has members");
            members.push(self.merged_on_genome(cis, &frame, &genome, &indices, window)?);
        }
        let mut exonic = Vec::new();
        for (i, m) in cis.members.iter().enumerate() {
            if grouped[i] {
                continue;
            }
            if intronic[i] {
                if inserts(m) || !members.contains(m) {
                    members.push(m.clone());
                }
            } else {
                exonic.push(m.clone());
            }
        }
        if !exonic.is_empty() {
            let sub = CisPhasedVariant::new(cis.ac.clone(), cis.gene.clone(), exonic)?;
            match self.normalize_cis(&sub)? {
                SV::CisPhased(c) => members.extend(c.members),
                one => members.push(one),
            }
        }
        members.sort_by_key(transcript_key);
        CisPhasedVariant::new(cis.ac.clone(), cis.gene.clone(), members)
    }

    /// The members of `cis` at `indices`, which meet on the genome, written
    /// as the one change they make over `window` there, trimmed to what it
    /// alters and not shifted, in the allele's own system.
    fn merged_on_genome(
        &self,
        cis: &CisPhasedVariant,
        frame: &Frame,
        genome: &Reference<'_, '_>,
        indices: &[usize],
        (lo, hi): Span,
    ) -> Result<SV, HgvsError> {
        let group = CisPhasedVariant::new(
            cis.ac.clone(),
            cis.gene.clone(),
            indices.iter().map(|&i| cis.members[i].clone()).collect(),
        )?;
        let reading = self.read_cis_on_genome(&group)?;
        let alt = reading.sequence(genome, lo, hi)?;
        let ref_ = Segment::Reference { start: lo, end: hi }.read(genome, 0, hi - lo)?;
        let (edit, s, e) = normalize::minimal_edit(genome, lo, hi, &ref_, &alt)?;
        let placed = PlacedEdit::from_hgvs_range(s, e, edit);
        let (first, last) = crate::mapper::hgvs_positions(&placed);
        let g = crate::structs::GVariant {
            ac: frame.accession.clone(),
            gene: cis.gene.clone(),
            posedit: crate::structs::PosEdit {
                pos: Some(crate::structs::SimpleInterval {
                    start: crate::mapper::simple_position(first),
                    end: last.map(crate::mapper::simple_position),
                    uncertain: false,
                }),
                edit: placed.edit,
                uncertain: false,
                predicted: false,
            },
        };
        let c = self.g_to_c(&g, &cis.ac)?;
        let coding = {
            let t = self.provider().get_transcript(&cis.ac, None)?;
            t.cds_start_index.is_some() && t.cds_end_index.is_some()
        };
        let as_n = || {
            if coding {
                return Err(HgvsError::UnsupportedOperation(format!(
                    "{cis}: an n. allele on a coding transcript is not written back from the genome"
                )));
            }
            Ok(SV::NonCoding(crate::structs::NVariant {
                ac: c.ac.clone(),
                gene: c.gene.clone(),
                posedit: c.posedit.clone(),
            }))
        };
        Ok(match cis.coordinate_type() {
            "c" => SV::Coding(c.clone()),
            "n" => as_n()?,
            _ if coding => SV::Rna(self.tx_to_r(&SV::Coding(c.clone()))?),
            _ => SV::Rna(self.tx_to_r(&SV::NonCoding(crate::structs::NVariant {
                ac: c.ac.clone(),
                gene: c.gene.clone(),
                posedit: c.posedit.clone(),
            }))?),
        })
    }

    /// The description of one reading, normalised: the members of `allele`
    /// under `choice`, with the bases `set_aside` and the identities among
    /// them stated unchanged.
    pub(super) fn describe(
        &self,
        frame: &Frame,
        allele: &Allele,
        choice: &Choice,
        set_aside: &[(usize, usize)],
        writing: Writing,
    ) -> Result<SV, HgvsError> {
        let normalised = self
            .normal_form(frame, allele, choice, set_aside, writing)
            .and_then(|form| self.written(frame, form.stated, form.changes));
        match normalised {
            Err(_) if writing == Writing::Reading => {
                self.written_plainly(frame, allele, choice, set_aside)
            }
            other => other,
        }
    }

    /// One reading written as its replacements make it, cut where a
    /// statement stands and not shifted: for a reading an ambiguity lists
    /// that the rule cannot write, so that the ambiguity still lists it.
    fn written_plainly(
        &self,
        frame: &Frame,
        allele: &Allele,
        choice: &Choice,
        set_aside: &[(usize, usize)],
    ) -> Result<SV, HgvsError> {
        let reference = allele.reference;
        let mut stated: Vec<Span> = set_aside.to_vec();
        stated.extend(
            allele
                .members
                .iter()
                .filter(|(_, a)| a.is_identity())
                .map(|(_, a)| a.named()),
        );
        let stated = merge_statements(stated, &frame.exons);
        let mut changes = Vec::new();
        for r in allele.assemble(choice, None, None, &mut Vec::new())? {
            let pieces_of = if r.start < r.end {
                cut_at(r, &stated)?
            } else {
                vec![r]
            };
            for piece in pieces_of {
                let window = (piece.start, piece.end);
                let alt = materialise(
                    &pieces(std::iter::once(&piece), window.0, window.1),
                    reference,
                )?;
                let ref_ = Segment::Reference {
                    start: window.0,
                    end: window.1,
                }
                .read(reference, 0, window.1 - window.0)?;
                if alt.eq_ignore_ascii_case(&ref_) {
                    continue;
                }
                let (edit, s, e) =
                    normalize::minimal_edit(reference, window.0, window.1, &ref_, &alt)?;
                let mut placed = PlacedEdit::from_hgvs_range(s, e, edit);
                if let NaEdit::Ins {
                    alt: Some(inserted),
                    uncertain,
                } = &placed.edit
                {
                    if placed.start > 0 && reference.base(placed.start)?.is_none() {
                        // No insertion after the last base: a delins of it.
                        let at = placed.start - 1;
                        let last = Segment::Reference {
                            start: at,
                            end: at + 1,
                        }
                        .read(reference, 0, 1)?;
                        placed = PlacedEdit {
                            start: at,
                            end: at + 1,
                            edit: NaEdit::RefAlt {
                                ref_: Some(String::new()),
                                alt: Some(format!("{last}{inserted}")),
                                uncertain: *uncertain,
                            },
                        };
                    }
                }
                changes.push(placed);
            }
        }
        self.written(frame, stated, changes)
    }

    /// One reading normalised, as the bases it states unchanged, sorted and
    /// neither overlapping nor touching, and the changes it makes.
    pub(super) fn normal_form(
        &self,
        frame: &Frame,
        allele: &Allele,
        choice: &Choice,
        set_aside: &[(usize, usize)],
        writing: Writing,
    ) -> Result<NormalForm, HgvsError> {
        let reference = allele.reference;
        let cds = match (frame.system, frame.mapper.as_ref()) {
            (System::Coding | System::Rna, Some(am)) => {
                match (am.transcript.cds_start_index, am.transcript.cds_end_index) {
                    (Some(s), Some(e)) if s.0 >= 0 && e.0 >= s.0 => {
                        Some((s.0 as usize, e.0 as usize + 1))
                    }
                    (Some(s), Some(e)) => {
                        return Err(HgvsError::ValidationError(format!(
                            "{} has a CDS from transcript index {} to {}, which names no codons",
                            am.transcript.ac, s.0, e.0
                        )))
                    }
                    _ => None,
                }
            }
            _ => None,
        };
        let ctx = Context {
            frame,
            allele,
            choice,
            writing,
            reference,
            extents: allele.extents()?,
            cds,
        };
        let mut stated: Vec<(usize, usize)> = set_aside.to_vec();
        stated.extend(
            allele
                .members
                .iter()
                .filter(|(_, a)| a.is_identity())
                .map(|(_, a)| a.named()),
        );
        let reps = allele.assemble_runs(choice, None, None, &mut Vec::new())?;
        let barriers = Barriers::new(&frame.exons, &stated, frame, &[]);
        // A statement across an exon junction projects over the intron, so
        // on the genome it states nothing a member agreeing with it does: it
        // stops a shift, and cuts no member.
        let in_one_exon: Vec<(usize, usize)> = stated
            .iter()
            .copied()
            .filter(|&(s, e)| frame.exons.exon(s) == frame.exons.exon(e - 1))
            .collect();
        let (mut units, statements) = ctx.units(reps, &merge_ranges(in_one_exon), &barriers)?;
        if !statements.is_empty() {
            stated.extend(statements);
            // A member that changes nothing inside what is now stated adds
            // nothing to it.
            let now_stated = merge_ranges(stated.clone());
            let mut kept = Vec::with_capacity(units.len());
            for u in units {
                let inside = outside(u.window, &now_stated).is_empty() && u.window.0 < u.window.1;
                let drop = match u.form {
                    Form::Member(k) => inside && allele.cancels(&[k])?,
                    _ => false,
                };
                if !drop {
                    kept.push(u);
                }
            }
            units = kept;
        }
        let first = Barriers::among(&frame.exons, &stated, frame, &units);
        for u in &mut units {
            if let Form::Member(k) = u.form {
                let w = u.window;
                if w.0 < w.1 && frame.exons.exon(w.0) != frame.exons.exon(w.1 - 1) {
                    u.form = Form::Spanning(self.alone(&ctx, k, &first)?);
                }
            }
        }
        loop {
            let barriers = Barriers::among(&frame.exons, &stated, frame, &units);
            let extents = units
                .iter()
                .map(|u| ctx.extent(u, &barriers))
                .collect::<Result<Vec<_>, _>>()?;
            let groups = components(units.len(), |a, b| {
                ctx.meet(&extents[a], &extents[b], &barriers)
            });
            if groups.iter().all(|g| g.len() == 1) {
                break;
            }
            let mut taken: Vec<Option<Unit>> = units.into_iter().map(Some).collect();
            let mut next = Vec::with_capacity(groups.len());
            for group in groups {
                let parts: Vec<Unit> = group
                    .iter()
                    .map(|&i| taken[i].take().expect("groups partition the units"))
                    .collect();
                if parts.len() == 1 {
                    next.extend(parts);
                    continue;
                }
                match ctx.merge(parts, &barriers)? {
                    Merged::Statement(named) => stated.extend(named),
                    Merged::Change(u) => next.push(u),
                }
            }
            units = next;
        }
        // A lone member settles against what the other units finally are:
        // the statements, and each unit left written across one.
        let barriers = Barriers::among(&frame.exons, &stated, frame, &units);
        let fixed = fixed_among(&units, &Barriers::new(&frame.exons, &stated, frame, &[]));
        let mut changes = Vec::with_capacity(units.len());
        for u in &units {
            let p = match &u.form {
                // Written across a statement: left as written.
                Form::Member(k) if fixed.contains(&u.window) => ctx.written_edit(*k)?,
                Form::Member(k) => self.alone(&ctx, *k, &barriers)?,
                Form::Change(p) | Form::Part(p) | Form::Spanning(p) => p.clone(),
                Form::Piece => unreachable!("a piece's change is read when it is made"),
            };
            if matches!(p.edit, NaEdit::RefAlt { alt: None, .. }) {
                stated.push((p.start, p.end));
            } else {
                changes.push(p);
            }
        }
        Ok(NormalForm {
            stated: merge_statements(stated, &frame.exons),
            changes,
        })
    }

    /// A member that meets no other, normalised as a single variant is,
    /// stopping at a barrier.
    fn alone(
        &self,
        ctx: &Context,
        slot: usize,
        barriers: &Barriers,
    ) -> Result<PlacedEdit, HgvsError> {
        let written = ctx.written_edit(slot)?;
        let own = ctx.allele.members[slot].1.written();
        let normalised = ctx.settle_change(
            written.clone(),
            &barriers.beside(&ctx.frame.exons, ctx.frame, own),
        )?;
        let statement = matches!(normalised.edit, NaEdit::RefAlt { alt: None, .. });
        if statement && matches!(written.edit, NaEdit::Repeat { .. }) {
            if let Some(am) = ctx.frame.mapper.as_ref() {
                let (s, e) = crate::mapper::hgvs_positions(&normalised);
                if crate::mapper::genome_may_count_more(am, s, e)?
                    && self.genome_range(ctx.frame, &written)?
                        != self.genome_range(ctx.frame, &normalised)?
                {
                    // On the genome the repeat states bases the identity
                    // over its transcript run does not.
                    return Ok(written);
                }
            }
        }
        Ok(normalised)
    }

    /// The genomic range a transcript change written at `placed` covers.
    fn genome_range(
        &self,
        frame: &Frame,
        placed: &PlacedEdit,
    ) -> Result<(usize, usize), HgvsError> {
        let g = match self.transcript_spelling(frame, placed)? {
            SV::Coding(c) => self.c_to_g(&c, None)?,
            SV::NonCoding(n) => self.n_to_g(&n, None)?,
            other => {
                return Err(HgvsError::Other(format!(
                    "{other} is not a transcript variant"
                )))
            }
        };
        let pos = g
            .posedit
            .pos
            .as_ref()
            .ok_or_else(|| HgvsError::ValidationError("Missing position".into()))?;
        crate::mapper::simple_interval_range(pos)
    }

    /// The description of a reading: `stated` as identities, cut where
    /// something is inserted between their bases, `changes`,
    /// and the members no sequence is built from as written, in order of
    /// position. The plain variant when there is one member.
    fn written(
        &self,
        frame: &Frame,
        stated: Vec<(usize, usize)>,
        changes: Vec<PlacedEdit>,
    ) -> Result<SV, HgvsError> {
        let inserted: Vec<usize> = changes
            .iter()
            .filter_map(|p| match p.edit {
                NaEdit::Ins { .. } => Some(p.start),
                NaEdit::Dup { .. } => Some(p.end),
                _ => None,
            })
            .chain(frame.placed.iter().filter_map(|(_, p, _, _)| match p {
                Placed::Unapplied(Reach::Junction(at)) => Some(*at),
                _ => None,
            }))
            .collect();
        let mut parts: Vec<PlacedEdit> = Vec::with_capacity(stated.len() + changes.len());
        for (s, e) in stated {
            let mut at = s;
            for cut in (s + 1..e).filter(|c| inserted.contains(c)) {
                parts.push(identity(at, cut));
                at = cut;
            }
            parts.push(identity(at, e));
        }
        parts.extend(changes);
        let mut members: Vec<((usize, usize), SV)> = parts
            .iter()
            .map(|p| Ok((p.hgvs_range(), self.member_variant(frame, p)?)))
            .collect::<Result<_, HgvsError>>()?;
        let mut unapplied: Vec<&SV> = Vec::new();
        for (_, p, written, _) in &frame.placed {
            let Placed::Unapplied(reach) = p else {
                continue;
            };
            // One written twice is one member, unless it inserts.
            if !inserts(written) && unapplied.contains(&written) {
                continue;
            }
            unapplied.push(written);
            members.push((reach_key(reach), written.clone()));
        }
        members.sort_by_key(|(k, _)| *k);
        let members: Vec<SV> = members.into_iter().map(|(_, m)| m).collect();
        let description = match <[SV; 1]>::try_from(members) {
            Ok([one]) => one,
            Err(members) => SV::CisPhased(CisPhasedVariant::new(
                frame.accession.clone(),
                frame.gene.clone(),
                members,
            )?),
        };
        if let SV::CisPhased(cis) = &description {
            if let Err(e) = self.read_cis_unwritten(cis) {
                let why = match e {
                    HgvsError::Ambiguous { .. } => "it has several readings itself".to_string(),
                    other => other.to_string(),
                };
                return Err(HgvsError::UnsupportedOperation(format!(
                    "{cis}, written from one reading, has no one reading of its own: {why}"
                )));
            }
        }
        Ok(description)
    }
}

fn identity(start: usize, end: usize) -> PlacedEdit {
    PlacedEdit {
        start,
        end,
        edit: NaEdit::RefAlt {
            ref_: None,
            alt: None,
            uncertain: false,
        },
    }
}
