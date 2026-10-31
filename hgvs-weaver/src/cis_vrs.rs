//! Reading a VRS `CisPhasedBlock` back as an allele in cis, `g.[a;b]`.
//!
//! A producer that justifies each member Allele alone can give two members
//! overlapping locations where their edits do not overlap: the deletion of
//! `g.[3_5del;7A>T]` on `TTCAGCAGTT` is justified over g.3_8, which holds
//! g.7. Each member is therefore read as every placement of its minimal edit
//! within its location, the placements of members whose locations meet are
//! tried together, and those the overlap rule accepts are the block's
//! readings: one sequence is the answer, normalised as one allele; several
//! are an [`HgvsError::Ambiguous`] listing them.

use crate::data::IdentifierType;
use crate::error::HgvsError;
use crate::mapper::{g_variant, interval, VariantMapper};
use crate::normalize;
use crate::reading::CisReading;
use crate::structs::{CisPhasedVariant, Variant};
use crate::vrs::{VrsAllele, VrsBound, VrsCisPhasedBlock, VrsState};
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
}

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

    /// `allele`, which reads back alone as `alone`, with its placements.
    fn member(&self, ac: &str, allele: &VrsAllele, alone: SV) -> Result<Member, HgvsError> {
        let location = match (allele.location.start, allele.location.end) {
            (VrsBound::Exact(s), VrsBound::Exact(e)) => Some((s, e)),
            _ => None,
        };
        let placements = match (&allele.state, location, &alone) {
            (VrsState::Literal { sequence, .. }, Some((start, end)), SV::Genomic(_)) => {
                self.placements(ac, start, end, sequence)?
            }
            (
                VrsState::ReferenceLength {
                    length,
                    repeat_subunit_length,
                    ..
                },
                Some((start, end)),
                SV::Genomic(_),
            ) => {
                // `alone` read, so the subunit fits the location.
                let reference = self.refs.reference(ac, IdentifierType::GenomicAccession);
                let unit = reference.slice(start, start + repeat_subunit_length)?;
                let alt: String = unit.chars().cycle().take(*length).collect();
                self.placements(ac, start, end, &alt)?
            }
            _ => Vec::new(),
        };
        let placements = match placements.is_empty() {
            true => vec![alone.clone()],
            false => placements,
        };
        Ok(Member {
            alone,
            placements,
            location,
        })
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
    /// locations meet placed together, normalised as one allele.
    fn placed_together(&self, ac: &str, members: Vec<Member>) -> Result<SV, HgvsError> {
        let reference = self.refs.reference(ac, IdentifierType::GenomicAccession);
        let mut chosen: Vec<Vec<SV>> = Vec::new();
        let mut ambiguous: Option<(usize, Vec<Vec<SV>>)> = None;
        for set in meeting(&members) {
            let set: Vec<&Member> = set.into_iter().map(|i| &members[i]).collect();
            if let [one] = set.as_slice() {
                chosen.push(vec![one.alone.clone()]);
                continue;
            }
            let mut readings = self.readings(ac, &set, &reference)?;
            let first = readings.remove(0);
            if !readings.is_empty() && ambiguous.is_none() {
                ambiguous = Some((chosen.len(), readings));
            }
            chosen.push(first);
        }
        let normalised = |chosen: &[Vec<SV>]| {
            let all = chosen.iter().flatten().cloned().collect();
            self.normalize_variant(SV::CisPhased(CisPhasedVariant::new(
                ac.to_string(),
                None,
                all,
            )?))
        };
        let Some((at, others)) = ambiguous else {
            return normalised(&chosen);
        };
        let mut readings = vec![normalised(&chosen)?];
        for other in others {
            let mut alternative = chosen.clone();
            alternative[at] = other;
            readings.push(normalised(&alternative)?);
        }
        let alone = members.into_iter().map(|m| m.alone).collect();
        Err(HgvsError::Ambiguous {
            description: CisPhasedVariant::new(ac.to_string(), None, alone)?.to_string(),
            readings,
        })
    }

    /// The placements of `set`, members whose locations meet, that the
    /// overlap rule accepts, one for each sequence and statement they make,
    /// at least one.
    ///
    /// # Errors
    ///
    /// When none is accepted, what reading the members as they read back
    /// alone raises; [`HgvsError::UnsupportedOperation`] when there are more
    /// than [`MAX_PLACEMENTS`] to try.
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
            .ok_or_else(|| {
                HgvsError::UnsupportedOperation(format!(
                    "the CisPhasedBlock's members on {ac} meet with more than {MAX_PLACEMENTS} \
                     placements between them, which are not tried"
                ))
            })?;
        let allele = |members: Vec<SV>| CisPhasedVariant::new(ac.to_string(), None, members);
        let mut found: Vec<(Vec<SV>, CisReading)> = Vec::new();
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
            let Ok(reading) = self.read_cis(&allele(placed.clone())?) else {
                continue;
            };
            let mut seen = false;
            for (_, other) in &found {
                if other.stated() == reading.stated() && other.same_sequence(&reading, reference)? {
                    seen = true;
                    break;
                }
            }
            if !seen {
                found.push((placed, reading));
            }
        }
        if found.is_empty() {
            self.read_cis(&allele(set.iter().map(|m| m.alone.clone()).collect())?)?;
        }
        Ok(found.into_iter().map(|(placed, _)| placed).collect())
    }
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
