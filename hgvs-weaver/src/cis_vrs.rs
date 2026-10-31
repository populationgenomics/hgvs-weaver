//! Reading a VRS `CisPhasedBlock` back as an allele in cis, `g.[a;b]`.
//!
//! A producer that justifies each member Allele alone can give two members
//! overlapping locations where their edits do not overlap: the deletion of
//! `g.[3_5del;7A>T]` on `TTCAGCAGTT` is justified over g.3_8, which holds
//! g.7. Each member is therefore read as every placement of its minimal edit
//! within its location, and the placements of members whose locations meet
//! are tried together. A combination is a reading when the overlap rule
//! accepts it, every member makes its own change (no two placed as one
//! edit, which the rule would read once), and it states only what the
//! block's reference Alleles state. Where some readings render as the block's
//! own members, as a block weaver wrote does, only those count. One sequence
//! is the answer, normalised as one allele, its members' own 3'-most
//! spellings preferred; several are an [`HgvsError::Ambiguous`] listing every
//! combination of each set's readings.

use crate::data::IdentifierType;
use crate::error::HgvsError;
use crate::mapper::{g_variant, interval, VariantMapper};
use crate::normalize;
use crate::reading::CisReading;
use crate::structs::{CisPhasedVariant, Variant};
use crate::vrs::{VrsAllele, VrsBound, VrsCisPhasedBlock, VrsState, VrsVariation};
use crate::SequenceVariant as SV;

/// The most placements of one set of meeting members tried together.
const MAX_PLACEMENTS: usize = 16384;

/// A member read back: the variant its Allele normalises to alone, every
/// placement of its minimal edit within its location (the variant alone
/// when it has only one), and the location, for an exact one.
struct Member {
    alone: SV,
    placements: Vec<SV>,
    location: Option<(usize, usize)>,
    /// The location and the bases the Allele puts there, for an exact
    /// location and a state with its bases: what identifies it.
    key: Option<Key>,
}

/// An Allele as its exact location and the bases it puts there.
type Key = (usize, usize, String);

/// A member a block renders, and the bounds it is justified within.
type Rendered = (SV, (usize, usize));

impl VariantMapper<'_> {
    /// `ac:g.[a;b]` (or `p.`) from a VRS `CisPhasedBlock`, on the sequence
    /// `accession` or the first member's refget accession names; every member
    /// must lie on it. A nucleotide block is read with its members placed
    /// together and normalised as one allele; a protein block member by member.
    pub(crate) fn cis_phased_from_vrs(
        &self,
        json: &str,
        accession: Option<&str>,
    ) -> Result<SV, HgvsError> {
        let block = VrsCisPhasedBlock::from_json(json)?;
        if let (Some(reference), Some(first)) = (&block.sequence_reference, block.members.first()) {
            let stated = &first.location.sequence_reference.refget_accession;
            if reference.refget_accession != *stated {
                return Err(HgvsError::ValidationError(format!(
                    "The CisPhasedBlock is on {} but its members on {stated}",
                    reference.refget_accession
                )));
            }
        }
        let mut ac = accession.map(str::to_string);
        let mut members = Vec::with_capacity(block.members.len());
        for allele in &block.members {
            let alone = self.from_vrs(&allele.to_json(), ac.as_deref())?;
            let ac = ac.get_or_insert_with(|| alone.ac().to_string());
            members.push(self.member(ac, allele, alone)?);
        }
        let ac = ac.ok_or_else(|| {
            HgvsError::ValidationError("A CisPhasedBlock has at least one member".into())
        })?;
        if members.iter().any(|m| matches!(m.alone, SV::Protein(_))) {
            let alone = members.into_iter().map(|m| m.alone).collect();
            return Ok(SV::CisPhased(CisPhasedVariant::new(ac, None, alone)?));
        }
        self.placed_together(&ac, members)
    }

    /// `allele`, which reads back alone as `alone`, with its placements,
    /// 3'-most first; `alone` when it has none to list.
    fn member(&self, ac: &str, allele: &VrsAllele, alone: SV) -> Result<Member, HgvsError> {
        let location = match (allele.location.start, allele.location.end) {
            (VrsBound::Exact(s), VrsBound::Exact(e)) => Some((s, e)),
            _ => None,
        };
        let key = match &alone {
            SV::Genomic(_) => self.key(ac, allele)?,
            _ => None,
        };
        let mut placements = match &key {
            Some((start, end, alt)) => self.placements(ac, *start, *end, alt)?,
            None => Vec::new(),
        };
        placements.reverse();
        if placements.is_empty() {
            placements.push(alone.clone());
        }
        Ok(Member {
            alone,
            placements,
            location,
            key,
        })
    }

    /// `allele` on `ac` as its exact location and the bases it puts there;
    /// `None` for a Range bound or a state without bases.
    fn key(&self, ac: &str, allele: &VrsAllele) -> Result<Option<Key>, HgvsError> {
        let (VrsBound::Exact(start), VrsBound::Exact(end)) =
            (allele.location.start, allele.location.end)
        else {
            return Ok(None);
        };
        let alt = match &allele.state {
            VrsState::Literal { sequence, .. } => sequence.to_ascii_uppercase(),
            VrsState::ReferenceLength {
                length,
                repeat_subunit_length,
                ..
            } => {
                let reference = self.refs.reference(ac, IdentifierType::GenomicAccession);
                let unit = reference.slice(start, start + repeat_subunit_length)?;
                unit.to_ascii_uppercase()
                    .chars()
                    .cycle()
                    .take(*length)
                    .collect()
            }
            _ => return Ok(None),
        };
        Ok(Some((start, end, alt)))
    }

    /// Every spelling of "the bases over `[start, end)` of `ac` become `alt`"
    /// as a minimal edit, 5'-most first; none when they are the reference.
    fn placements(
        &self,
        ac: &str,
        start: usize,
        end: usize,
        alt: &str,
    ) -> Result<Vec<SV>, HgvsError> {
        let reference = self.refs.reference(ac, IdentifierType::GenomicAccession);
        let ref_ = reference.slice(start, end)?;
        if ref_.eq_ignore_ascii_case(alt) {
            return Ok(Vec::new());
        }
        let (r, a) = (ref_.as_bytes(), alt.as_bytes());
        let room = r.len().min(a.len());
        let same = |(x, y): &(&u8, &u8)| x.eq_ignore_ascii_case(y);
        let prefix = r.iter().zip(a).take_while(same).count().min(room);
        let suffix = r.iter().rev().zip(a.iter().rev()).take_while(same).count();
        let trimmed = (prefix + suffix).min(room);
        (trimmed - suffix.min(trimmed)..=prefix.min(trimmed))
            .map(|left| {
                let right = trimmed - left;
                let (edit, s, e) = normalize::minimal_edit(
                    &reference,
                    start + left,
                    end - right,
                    &ref_[left..r.len() - right],
                    &alt[left..a.len() - right],
                )?;
                Ok(SV::Genomic(g_variant(ac, interval(s, e)?, edit)))
            })
            .collect()
    }

    /// The allele `members` make on `ac` with each set of members whose
    /// locations meet placed together, normalised as one allele; an
    /// [`HgvsError::Ambiguous`] listing every combination of the sets'
    /// readings where a set has several.
    fn placed_together(&self, ac: &str, members: Vec<Member>) -> Result<SV, HgvsError> {
        let reference = self.refs.reference(ac, IdentifierType::GenomicAccession);
        let mut sets: Vec<Vec<Vec<SV>>> = Vec::new();
        for set in meeting(&members) {
            let set: Vec<&Member> = set.into_iter().map(|i| &members[i]).collect();
            sets.push(match set.as_slice() {
                [one] => vec![vec![one.alone.clone()]],
                _ => self.readings(ac, &set, &reference)?,
            });
        }
        let combinations = sets
            .iter()
            .try_fold(1usize, |n, s| n.checked_mul(s.len()))
            .filter(|&n| n <= MAX_PLACEMENTS)
            .ok_or_else(|| too_many(ac))?;
        let mut readings = Vec::with_capacity(combinations);
        for n in 0..combinations {
            let mut rest = n;
            let mut all = Vec::new();
            for set in &sets {
                all.extend(set[rest % set.len()].iter().cloned());
                rest /= set.len();
            }
            readings.push(self.normalize_variant(SV::CisPhased(CisPhasedVariant::new(
                ac.to_string(),
                None,
                all,
            )?))?);
        }
        if readings.len() == 1 {
            return Ok(readings.remove(0));
        }
        let alone = members.into_iter().map(|m| m.alone).collect();
        Err(HgvsError::Ambiguous {
            description: CisPhasedVariant::new(ac.to_string(), None, alone)?.to_string(),
            readings,
        })
    }

    /// The readings of `set`, members whose locations meet: one placement
    /// for each sequence they make, the first found in the order the
    /// members' placements are listed, and each reading of a placement the
    /// overlap rule finds ambiguous, such as two insertions at one point.
    /// Where some render as the members of `set` do, only those.
    ///
    /// # Errors
    ///
    /// When no combination is a reading, what reading the first that the
    /// overlap rule refuses raises, else a [`HgvsError::ValidationError`];
    /// [`HgvsError::UnsupportedOperation`] when there are more than
    /// [`MAX_PLACEMENTS`] to try.
    fn readings(
        &self,
        ac: &str,
        set: &[&Member],
        reference: &crate::reference::Reference<'_, '_>,
    ) -> Result<Vec<Vec<SV>>, HgvsError> {
        let count = set
            .iter()
            .try_fold(1usize, |n, m| n.checked_mul(m.placements.len()))
            .filter(|&n| n <= MAX_PLACEMENTS)
            .ok_or_else(|| too_many(ac))?;
        let allele = |members: Vec<SV>| CisPhasedVariant::new(ac.to_string(), None, members);
        let identities: Vec<SV> = set
            .iter()
            .filter(|m| {
                m.key
                    .as_ref()
                    .is_some_and(|k| reference_allele(reference, k))
            })
            .map(|m| m.alone.clone())
            .collect();
        let stated = match identities.is_empty() {
            true => Vec::new(),
            false => self.read_cis(&allele(identities)?)?.stated().to_vec(),
        };
        let mut keys: Vec<&Key> = set.iter().filter_map(|m| m.key.as_ref()).collect();
        keys.sort_unstable();
        let mut refused: Option<HgvsError> = None;
        // Each sequence made, with a placement making it and whether that
        // placement renders as the block's members.
        let mut found: Vec<(Vec<SV>, CisReading, bool)> = Vec::new();
        for n in 0..count {
            let mut rest = n;
            let placed: Vec<SV> = set
                .iter()
                .map(|m| {
                    let p = m.placements[rest % m.placements.len()].clone();
                    rest /= m.placements.len();
                    p
                })
                .collect();
            // Two members placed as one edit would be read once; two
            // insertions at one point are not, and read as one copy or two.
            if (1..placed.len()).any(|i| placed[..i].contains(&placed[i]) && !inserts(&placed[i])) {
                continue;
            }
            let candidates = match self.read_cis(&allele(placed.clone())?) {
                Ok(r) => vec![(placed, r)],
                Err(HgvsError::Ambiguous { readings, .. }) => {
                    let mut each = Vec::with_capacity(readings.len());
                    for r in readings {
                        let members = match r {
                            SV::CisPhased(c) => c.members,
                            one => vec![one],
                        };
                        let reading = self.read_cis(&allele(members.clone())?)?;
                        each.push((members, reading));
                    }
                    each
                }
                Err(e) => {
                    refused.get_or_insert(e);
                    continue;
                }
            };
            for (placed, reading) in candidates {
                if reading.stated() != stated.as_slice() {
                    continue;
                }
                // Rendered as written: an insertion that is a duplication
                // renders as one.
                let placed = self.with_duplications(ac, placed, reference)?;
                let renders = self.renders_as(ac, &placed, &keys)?;
                let mut seen = false;
                for (p, other, r) in found.iter_mut() {
                    if other.same_sequence(&reading, reference)? {
                        if renders && !*r {
                            (*p, *r) = (placed.clone(), true);
                        }
                        seen = true;
                        break;
                    }
                }
                if !seen {
                    found.push((placed, reading, renders));
                }
            }
        }
        if found.iter().any(|(.., r)| *r) {
            found.retain(|(.., r)| *r);
        }
        if found.is_empty() {
            return Err(refused.unwrap_or_else(|| {
                HgvsError::ValidationError(format!(
                    "no placement of the CisPhasedBlock's members on {ac} applies each of them"
                ))
            }));
        }
        Ok(found.into_iter().map(|(placed, ..)| placed).collect())
    }

    /// `members` with each insertion of a copy of the bases just before it
    /// written as their duplication, as HGVS writes it, where the allele
    /// still reads as one and makes what it made: a duplication reads the
    /// bases it copies, which a neighbour may change.
    fn with_duplications(
        &self,
        ac: &str,
        mut members: Vec<SV>,
        reference: &crate::reference::Reference<'_, '_>,
    ) -> Result<Vec<SV>, HgvsError> {
        let allele = |m: &[SV]| CisPhasedVariant::new(ac.to_string(), None, m.to_vec());
        let before = self.read_cis(&allele(&members)?)?;
        for i in 0..members.len() {
            let Some(dup) = as_duplication(reference, &members[i])? else {
                continue;
            };
            let mut respelled = members.clone();
            respelled[i] = dup;
            if let Ok(after) = self.read_cis(&allele(&respelled)?) {
                if after.stated() == before.stated() && after.same_sequence(&before, reference)? {
                    members = respelled;
                }
            }
        }
        Ok(members)
    }

    /// Whether the allele `placed` on `ac` renders as the Alleles `keys`,
    /// sorted.
    fn renders_as(&self, ac: &str, placed: &[SV], keys: &[&Key]) -> Result<bool, HgvsError> {
        let cis = SV::CisPhased(CisPhasedVariant::new(
            ac.to_string(),
            None,
            placed.to_vec(),
        )?);
        let Ok(rendered) = self.to_vrs_variation(&cis) else {
            return Ok(false);
        };
        let alleles = match rendered {
            VrsVariation::Allele(a) => vec![a],
            VrsVariation::CisPhasedBlock(b) => b.members,
            _ => return Ok(false),
        };
        let mut rendered = Vec::with_capacity(alleles.len());
        for a in &alleles {
            match self.key(ac, a)? {
                Some(k) => rendered.push(k),
                None => return Ok(false),
            }
        }
        rendered.sort_unstable();
        Ok(rendered.iter().eq(keys.iter().copied()))
    }
}

impl VariantMapper<'_> {
    /// `cis`, a `g.` allele normalised on its genome, as the members its
    /// block renders, and the bounds each is justified within. Normalisation
    /// keeps a change apart from the bases a positional identity names, and
    /// a change's region of ambiguity stops there too: the identity says
    /// those bases are as they were. So each change is justified within the
    /// stretch between the identities either side of it, and changes whose
    /// Alleles then share a base, or an insertion point, are one change,
    /// until none do. Without the bound,
    /// two changes either side of an identity in one run are each justified
    /// over the whole run, and the block no longer says which comes first:
    /// two alleles, one identifier.
    pub(crate) fn rendered_members(
        &self,
        cis: &CisPhasedVariant,
    ) -> Result<Vec<Rendered>, HgvsError> {
        let whole = (0, usize::MAX);
        let reference = self
            .refs
            .reference(&cis.ac, IdentifierType::GenomicAccession);
        let mut stated = Vec::new();
        for m in &cis.members {
            if let SV::Genomic(g) = m {
                if g.posedit.edit.is_identity() && m.has_exact_position() {
                    if let Some(pos) = g.posedit.pos.as_ref() {
                        stated.push(crate::mapper::simple_interval_range(pos)?);
                    }
                }
            }
        }
        let mut members: Vec<Rendered> = Vec::with_capacity(cis.members.len());
        for m in &cis.members {
            let bounds = match self.change_core(m)? {
                Some(core) if !stated.is_empty() => between(&stated, core),
                _ => whole,
            };
            members.push((m.clone(), bounds));
        }
        if stated.is_empty() {
            return Ok(members);
        }
        loop {
            let mut spans = Vec::with_capacity(members.len());
            for (m, bounds) in &members {
                spans.push(match self.change_core(m)? {
                    Some(_) => {
                        let a = self.canonical_allele_within(m, *bounds)?;
                        Some((a.start, a.end))
                    }
                    None => None,
                });
            }
            // The changes whose Alleles meet, directly or through another,
            // made one change together: one pair at a time, a merged change
            // written at its own minimal spelling could move past a third
            // member it meets.
            let meets = |i: usize, j: usize| match (spans[i], spans[j]) {
                (Some(a), Some(b)) => overlap(a, b),
                _ => false,
            };
            let Some(first) =
                (0..members.len()).find(|&i| (0..members.len()).any(|j| j != i && meets(i, j)))
            else {
                return Ok(members);
            };
            let mut group = vec![first];
            let mut k = 0;
            while k < group.len() {
                for j in 0..members.len() {
                    if !group.contains(&j) && meets(group[k], j) {
                        group.push(j);
                    }
                }
                k += 1;
            }
            group.sort_unstable();
            let spans_of: Vec<(usize, usize)> = group.iter().filter_map(|&g| spans[g]).collect();
            let lo = spans_of.iter().map(|s| s.0).min().unwrap_or(0);
            let hi = spans_of.iter().map(|s| s.1).max().unwrap_or(0);
            let together = group.iter().map(|&g| members[g].0.clone()).collect();
            let reading = self.read_cis(&CisPhasedVariant::new(cis.ac.clone(), None, together)?)?;
            let alt = reading.sequence(&reference, lo, hi)?;
            let ref_ = reference.slice(lo, hi)?;
            let (edit, s, e) = normalize::minimal_edit(&reference, lo, hi, &ref_, &alt)?;
            let one = SV::Genomic(g_variant(&cis.ac, interval(s, e)?, edit));
            let bounds = members[group[0]].1;
            for &g in group.iter().rev() {
                members.remove(g);
            }
            members.insert(group[0], (one, bounds));
        }
    }

    /// The bases `member`, a change at an exact position on the genome,
    /// alters, trimmed of what it writes back at either end: an insertion
    /// as the empty range at its junction. `None` for an identity or a
    /// member without bases to apply.
    fn change_core(&self, member: &SV) -> Result<Option<(usize, usize)>, HgvsError> {
        let SV::Genomic(g) = member else {
            return Ok(None);
        };
        if g.posedit.edit.is_identity() || !member.has_exact_position() {
            return Ok(None);
        }
        let Some(pos) = g.posedit.pos.as_ref() else {
            return Ok(None);
        };
        let reference = self.refs.reference(&g.ac, IdentifierType::GenomicAccession);
        let (s, e) = crate::mapper::simple_interval_range(pos)?;
        let Ok(resolved) = g.posedit.edit.resolve(&reference, s, e) else {
            return Ok(None);
        };
        let actual = reference.slice(resolved.start, resolved.end)?;
        let (at, r, _) = crate::structs::strip_common_prefix_suffix(
            resolved.start as i32,
            &actual,
            &resolved.alt,
        );
        Ok(Some((at as usize, at as usize + r.len())))
    }
}

/// Whether two Allele locations share a base, or one inserts at a point
/// the other holds, borders or also inserts at: two insertions there would
/// leave their order unsaid.
fn overlap(a: (usize, usize), b: (usize, usize)) -> bool {
    let inside = |p: usize, (s, e): (usize, usize)| s <= p && p <= e;
    (a.0 < b.1 && b.0 < a.1) || (a.0 == a.1 && inside(a.0, b)) || (b.0 == b.1 && inside(b.0, a))
}

/// The stretch between the `stated` ranges either side of `core`, which
/// meets none of them inside.
fn between(stated: &[(usize, usize)], (s, e): (usize, usize)) -> (usize, usize) {
    let lo = stated
        .iter()
        .map(|r| r.1)
        .filter(|&end| end <= s)
        .max()
        .unwrap_or(0);
    let hi = stated
        .iter()
        .map(|r| r.0)
        .filter(|&start| start >= e)
        .min()
        .unwrap_or(usize::MAX);
    (lo, hi)
}

/// Whether the Allele `key` puts back the bases its location holds.
fn reference_allele(reference: &crate::reference::Reference<'_, '_>, key: &Key) -> bool {
    reference
        .slice(key.0, key.1)
        .is_ok_and(|bases| bases.eq_ignore_ascii_case(&key.2))
}

/// `member`, an insertion of a copy of the bases just before it, as the
/// duplication of those bases; `None` for any other member.
fn as_duplication(
    reference: &crate::reference::Reference<'_, '_>,
    member: &SV,
) -> Result<Option<SV>, HgvsError> {
    let SV::Genomic(g) = member else {
        return Ok(None);
    };
    let (crate::edits::NaEdit::Ins { alt: Some(ins), .. }, Some(pos)) =
        (&g.posedit.edit, g.posedit.pos.as_ref())
    else {
        return Ok(None);
    };
    // The insertion goes after its first flanking base, g.N of g.N_N+1.
    let at = usize::try_from(pos.start.index()?.0)
        .map_err(|_| HgvsError::Other(format!("{g} starts before its sequence")))?
        + 1;
    let Some(from) = at.checked_sub(ins.len()) else {
        return Ok(None);
    };
    if !reference.slice(from, at)?.eq_ignore_ascii_case(ins) {
        return Ok(None);
    }
    let dup = crate::edits::NaEdit::Dup {
        ref_: None,
        uncertain: false,
    };
    Ok(Some(SV::Genomic(g_variant(
        &g.ac,
        interval(from, at)?,
        dup,
    ))))
}

/// The refusal of a block whose meeting members have too many placements.
fn too_many(ac: &str) -> HgvsError {
    HgvsError::UnsupportedOperation(format!(
        "the CisPhasedBlock's members on {ac} meet with more than {MAX_PLACEMENTS} placements \
         between them, which are not tried"
    ))
}

/// The indices of `members` in sets whose exact locations overlap or touch,
/// directly or through another, each sorted, in order of their first member.
fn meeting(members: &[Member]) -> Vec<Vec<usize>> {
    let mut sets: Vec<Vec<usize>> = Vec::new();
    for (i, m) in members.iter().enumerate() {
        let meets = |j: &usize| match (m.location, members[*j].location) {
            (Some((s1, e1)), Some((s2, e2))) => s1 <= e2 && s2 <= e1,
            _ => false,
        };
        let (mut joined, apart): (Vec<Vec<usize>>, Vec<Vec<usize>>) =
            sets.into_iter().partition(|set| set.iter().any(meets));
        let mut set: Vec<usize> = joined.drain(..).flatten().collect();
        set.push(i);
        set.sort_unstable();
        sets = apart;
        sets.push(set);
    }
    sets.sort_unstable_by_key(|set| set[0]);
    sets
}

/// Whether `member` inserts bases at a point: an insertion or a duplication.
fn inserts(member: &SV) -> bool {
    matches!(
        member,
        SV::Genomic(g) if matches!(
            g.posedit.edit,
            crate::edits::NaEdit::Ins { .. } | crate::edits::NaEdit::Dup { .. }
        )
    )
}
