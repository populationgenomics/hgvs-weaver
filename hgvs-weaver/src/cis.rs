//! Projecting an allele in cis, `c.[a;b]`, between coordinate systems.
//!
//! Every member is numbered against the reference, so an allele projects
//! member by member, by the same rules as a single variant, and the result is
//! an allele in cis of the projections: one member stays an allele of one
//! member. If any member fails, the projection fails with
//! [`HgvsError::CisMember`], naming the member. Across strands the members
//! are written in reverse, so an allele written in transcript order is in
//! genomic order on the genome and back.
//!
//! One case reads differently member by member. A duplication's copy goes 3'
//! of the original on the strand it is written on, which on the other strand
//! is 5'. Where another member inserts at the junction the duplication's
//! projected form would insert at, the two would be insertions at one
//! junction, whose order the target cannot say; so that duplication is written
//! as the insertion it is, at the junction it uses where written. That
//! insertion can meet another duplication the same way, so the rewriting
//! repeats until no duplication meets an insertion.

use crate::data::{IdentifierType, Strand, TranscriptSearch};
use crate::edits::NaEdit;
use crate::error::HgvsError;
use crate::mapper::VariantMapper;
use crate::structs::{
    CisPhasedVariant, GVariant, GenomicPos, LinearVariant, PosEdit, SimpleInterval, SimplePosition,
    Variant,
};
use crate::SequenceVariant as SV;

/// The first and last 0-based index an exact g. interval names; `None` for
/// an uncertain or unknown one.
fn exact_range(pos: Option<&SimpleInterval>) -> Option<(i32, i32)> {
    let pos = pos?;
    let exact = |p: &SimplePosition| (p.end.is_none() && !p.uncertain).then(|| p.index().ok())?;
    if pos.uncertain {
        return None;
    }
    let first = exact(&pos.start)?;
    let last = match &pos.end {
        Some(e) => exact(e)?,
        None => first,
    };
    Some((first.0, last.0))
}

fn is_insertion(edit: &NaEdit) -> bool {
    matches!(edit, NaEdit::Ins { .. } | NaEdit::InsLength { .. })
}

/// Whether a member is a duplication, whether an insertion, and whether its
/// position is exact, as written.
fn shape(member: &SV) -> (bool, bool, bool) {
    let tx_exact =
        |pos: Option<&crate::structs::BaseOffsetInterval>| pos.is_some_and(|p| p.is_exact());
    let (edit, exact) = match member {
        SV::Genomic(g) => (
            &g.posedit.edit,
            exact_range(g.posedit.pos.as_ref()).is_some(),
        ),
        SV::Coding(c) => (&c.posedit.edit, tx_exact(c.posedit.pos.as_ref())),
        SV::NonCoding(n) => (&n.posedit.edit, tx_exact(n.posedit.pos.as_ref())),
        SV::Rna(r) => (&r.posedit.edit, tx_exact(r.posedit.pos.as_ref())),
        _ => return (false, false, true),
    };
    (is_dup(edit), is_insertion(edit), exact)
}

/// Refuses, across strands, an allele with a duplication and an insertion
/// where either is at an uncertain position: whether the copy meets the
/// insertion, and so how the allele reads on the other strand, is unknown.
fn refuse_uncertain_copies(cis: &CisPhasedVariant) -> Result<(), HgvsError> {
    let shapes: Vec<_> = cis.members.iter().map(shape).collect();
    if !shapes.iter().any(|s| s.0) || !shapes.iter().any(|s| s.1) {
        return Ok(());
    }
    match shapes
        .iter()
        .position(|&(dup, ins, exact)| (dup || ins) && !exact)
    {
        Some(i) => Err(member_error(
            cis,
            i,
            HgvsError::UnsupportedOperation(
                "its position is uncertain, so whether a duplication's copy meets an \
                 insertion on the other strand cannot be told"
                    .into(),
            ),
        )),
        None => Ok(()),
    }
}

fn is_dup(edit: &NaEdit) -> bool {
    matches!(edit, NaEdit::Dup { .. })
}

/// Whether some member other than the `skip`th is an insertion, of stated
/// bases or of a length, between the 0-based indices `left` and `left + 1`.
fn inserts_at(members: &[GVariant], skip: usize, left: i32) -> bool {
    members.iter().enumerate().any(|(j, m)| {
        j != skip
            && is_insertion(&m.posedit.edit)
            && exact_range(m.posedit.pos.as_ref()) == Some((left, left + 1))
    })
}

/// The g. insertion of `bases` between the 0-based indices `left` and
/// `left + 1`, with `like`'s accession, gene and flags.
fn insertion(like: &GVariant, left: i32, bases: String, uncertain: bool) -> GVariant {
    let position = |i: i32| SimplePosition {
        base: GenomicPos(i).to_hgvs().into(),
        end: None,
        uncertain: false,
    };
    GVariant::from_parts(
        like.ac.clone(),
        like.gene.clone(),
        PosEdit {
            pos: Some(SimpleInterval {
                start: position(left),
                end: Some(position(left + 1)),
                uncertain: false,
            }),
            edit: NaEdit::Ins {
                alt: Some(bases),
                uncertain,
            },
            uncertain: like.posedit.uncertain,
            predicted: like.posedit.predicted,
        },
    )
}

fn edit_uncertain(edit: &NaEdit) -> bool {
    matches!(
        edit,
        NaEdit::Dup {
            uncertain: true,
            ..
        }
    )
}

impl VariantMapper<'_> {
    /// Projects an allele in cis of `c.` members to the genome, member by
    /// member as [`c_to_g`](Self::c_to_g) projects one.
    pub fn c_to_g_cis(
        &self,
        cis: &CisPhasedVariant,
        reference_ac: Option<&str>,
    ) -> Result<CisPhasedVariant, HgvsError> {
        expect_system(cis, "c", "c_to_g")?;
        self.cis_tx_to_g(cis, reference_ac, |m| match m {
            SV::Coding(c) => self.c_to_g(c, reference_ac),
            other => Err(wrong_member(other, "c")),
        })
    }

    /// Projects an allele in cis of `n.` members to the genome, member by
    /// member as [`n_to_g`](Self::n_to_g) projects one.
    pub fn n_to_g_cis(
        &self,
        cis: &CisPhasedVariant,
        reference_ac: Option<&str>,
    ) -> Result<CisPhasedVariant, HgvsError> {
        expect_system(cis, "n", "n_to_g")?;
        self.cis_tx_to_g(cis, reference_ac, |m| match m {
            SV::NonCoding(n) => self.n_to_g(n, reference_ac),
            other => Err(wrong_member(other, "n")),
        })
    }

    /// Projects an allele in cis of `r.` members to the genome, member by
    /// member as [`r_to_g`](Self::r_to_g) projects one: a member spanning a
    /// splice junction has no genomic form, so neither has the allele.
    pub fn r_to_g_cis(
        &self,
        cis: &CisPhasedVariant,
        reference_ac: Option<&str>,
    ) -> Result<CisPhasedVariant, HgvsError> {
        expect_system(cis, "r", "r_to_g")?;
        self.cis_tx_to_g(cis, reference_ac, |m| match m {
            SV::Rna(r) => self.r_to_g(r, reference_ac),
            other => Err(wrong_member(other, "r")),
        })
    }

    /// Projects an allele in cis of `g.` members to `transcript_ac`, member
    /// by member as [`g_to_c`](Self::g_to_c) projects one.
    pub fn g_to_c_cis(
        &self,
        cis: &CisPhasedVariant,
        transcript_ac: &str,
    ) -> Result<CisPhasedVariant, HgvsError> {
        expect_system(cis, "g", "g_to_c")?;
        let strand = self
            .provider()
            .get_transcript(transcript_ac, Some(&cis.ac))?
            .strand;
        let mut genomic: Vec<GVariant> = project_members(cis, |m| match m {
            SV::Genomic(g) => Ok(g.clone()),
            other => Err(wrong_member(other, "g")),
        })?;
        let mut projected = project_members(cis, |m| match m {
            SV::Genomic(g) => self.g_to_c(g, transcript_ac),
            other => Err(wrong_member(other, "g")),
        })?;
        if strand == Strand::Minus {
            refuse_uncertain_copies(cis)?;
            let candidates: Vec<bool> = genomic
                .iter()
                .zip(&projected)
                .map(|(g, c)| is_dup(&g.posedit.edit) && is_dup(&c.posedit.edit))
                .collect();
            for i in self.copies_as_insertions(cis, &mut genomic, &candidates, Written::OnGenome)? {
                projected[i] = self
                    .g_to_c(&genomic[i], transcript_ac)
                    .map_err(|e| member_error(cis, i, e))?;
            }
            projected.reverse();
        }
        allele(cis, projected.into_iter().map(SV::Coding).collect())
    }

    /// Projects an allele in cis of `g.` members to every transcript the
    /// searcher finds over the span of its members, as
    /// [`g_to_c_all`](Self::g_to_c_all) does for one variant: the alleles of
    /// the transcripts every member projects to, an error if there are none.
    pub fn g_to_c_all_cis(
        &self,
        cis: &CisPhasedVariant,
        searcher: &dyn TranscriptSearch,
    ) -> Result<Vec<CisPhasedVariant>, HgvsError> {
        expect_system(cis, "g", "g_to_c_all")?;
        let spans = project_members(cis, |m| {
            let SV::Genomic(g) = m else {
                return Err(wrong_member(m, "g"));
            };
            let pos = g
                .posedit
                .pos
                .as_ref()
                .ok_or_else(|| HgvsError::ValidationError("Missing position".into()))?;
            let start = pos.start.index()?.0;
            let end = match &pos.end {
                Some(e) => e.index()?.0 + 1,
                None => start + 1,
            };
            Ok((start, end))
        })?;
        let (Some(start), Some(end)) = (
            spans.iter().map(|s| s.0).min(),
            spans.iter().map(|s| s.1).max(),
        ) else {
            return Err(HgvsError::ValidationError(format!(
                "{cis} has no members to find transcripts for"
            )));
        };
        let transcripts = searcher.get_transcripts_for_region(&cis.ac, start, end)?;
        if transcripts.is_empty() {
            return Err(HgvsError::ValidationError(format!(
                "No transcripts found for region {}:{}-{}",
                cis.ac, start, end
            )));
        }
        let mut results = Vec::new();
        let mut errors = Vec::new();
        for tx_ac in transcripts {
            match self.g_to_c_cis(cis, &tx_ac) {
                Ok(c) => results.push(c),
                Err(e) => errors.push(format!("{tx_ac}: {e}")),
            }
        }
        if results.is_empty() && !errors.is_empty() {
            return Err(HgvsError::ValidationError(format!(
                "All mapping attempts failed: {}",
                errors.join("; ")
            )));
        }
        Ok(results)
    }

    /// Respells an allele in cis of `c.` members in `r.`, member by member.
    pub fn c_to_r_cis(&self, cis: &CisPhasedVariant) -> Result<CisPhasedVariant, HgvsError> {
        expect_system(cis, "c", "c_to_r")?;
        let members = project_members(cis, |m| match m {
            SV::Coding(c) => self.c_to_r(c).map(SV::Rna),
            other => Err(wrong_member(other, "c")),
        })?;
        allele(cis, members)
    }

    /// Respells an allele in cis of `n.` members in `r.`, member by member.
    pub fn n_to_r_cis(&self, cis: &CisPhasedVariant) -> Result<CisPhasedVariant, HgvsError> {
        expect_system(cis, "n", "n_to_r")?;
        let members = project_members(cis, |m| match m {
            SV::NonCoding(n) => self.n_to_r(n).map(SV::Rna),
            other => Err(wrong_member(other, "n")),
        })?;
        allele(cis, members)
    }

    /// Respells an allele in cis of `r.` members in `c.`, member by member.
    pub fn r_to_c_cis(&self, cis: &CisPhasedVariant) -> Result<CisPhasedVariant, HgvsError> {
        expect_system(cis, "r", "r_to_c")?;
        let members = project_members(cis, |m| match m {
            SV::Rna(r) => self.r_to_c(r).map(SV::Coding),
            other => Err(wrong_member(other, "r")),
        })?;
        allele(cis, members)
    }

    /// Respells an allele in cis of `r.` members on a non-coding transcript
    /// in `n.`, member by member.
    pub fn r_to_n_cis(&self, cis: &CisPhasedVariant) -> Result<CisPhasedVariant, HgvsError> {
        expect_system(cis, "r", "r_to_n")?;
        let members = project_members(cis, |m| match m {
            SV::Rna(r) => self.r_to_n(r).map(SV::NonCoding),
            other => Err(wrong_member(other, "r")),
        })?;
        allele(cis, members)
    }

    /// A transcript allele on the genome: each member by `project`, a
    /// minus-strand duplication beside an insertion at its genomic junction
    /// written as the insertion, and the order reversed on the minus strand.
    fn cis_tx_to_g(
        &self,
        cis: &CisPhasedVariant,
        reference_ac: Option<&str>,
        project: impl Fn(&SV) -> Result<GVariant, HgvsError>,
    ) -> Result<CisPhasedVariant, HgvsError> {
        let strand = self
            .provider()
            .get_transcript(&cis.ac, reference_ac)?
            .strand;
        let mut projected = project_members(cis, project)?;
        if strand == Strand::Minus {
            refuse_uncertain_copies(cis)?;
            let candidates: Vec<bool> = cis
                .members
                .iter()
                .zip(&projected)
                .map(|(m, g)| {
                    let source_dup = match m {
                        SV::Coding(c) => is_dup(&c.posedit.edit),
                        SV::NonCoding(n) => is_dup(&n.posedit.edit),
                        SV::Rna(r) => is_dup(&r.posedit.edit),
                        _ => false,
                    };
                    source_dup && is_dup(&g.posedit.edit)
                })
                .collect();
            self.copies_as_insertions(cis, &mut projected, &candidates, Written::OnTranscript)?;
            projected.reverse();
        }
        allele(cis, projected.into_iter().map(SV::Genomic).collect())
    }

    /// Rewrites, in `genome`, the members' genomic forms, each candidate
    /// duplication whose copy the other strand would place where another
    /// member inserts as the insertion of its bases where it was `written`,
    /// round after round until none is left. Returns the members rewritten.
    fn copies_as_insertions(
        &self,
        cis: &CisPhasedVariant,
        genome: &mut [GVariant],
        candidates: &[bool],
        written: Written,
    ) -> Result<Vec<usize>, HgvsError> {
        let mut rewritten = Vec::new();
        loop {
            let mut changed = false;
            for i in 0..genome.len() {
                if !candidates[i] || rewritten.contains(&i) {
                    continue;
                }
                let g = &genome[i];
                let Some((first, last)) = exact_range(g.posedit.pos.as_ref()) else {
                    continue;
                };
                // The junction left of the copy: on the genome's strand it is
                // after `last`, on the minus-strand transcript's before `first`.
                let (keep, other) = match written {
                    Written::OnGenome => (last, first - 1),
                    Written::OnTranscript => (first - 1, last),
                };
                if !inserts_at(genome, i, other) {
                    continue;
                }
                if keep < 0 {
                    return Err(member_error(
                        cis,
                        i,
                        HgvsError::ValidationError(format!(
                            "{g} copies the first base of {}, so its copy has no junction \
                             before it",
                            g.ac
                        )),
                    ));
                }
                let bases = self
                    .genome_bases(&g.ac, first, last)
                    .map_err(|e| member_error(cis, i, e))?;
                genome[i] = insertion(g, keep, bases, edit_uncertain(&g.posedit.edit));
                rewritten.push(i);
                changed = true;
            }
            if !changed {
                return Ok(rewritten);
            }
        }
    }

    /// The bases of `ac` at the 0-based indices `first..=last`.
    fn genome_bases(&self, ac: &str, first: i32, last: i32) -> Result<String, HgvsError> {
        let bases =
            self.provider()
                .get_seq(ac, first, Some(last + 1), IdentifierType::GenomicAccession)?;
        if bases.len() != (last + 1 - first) as usize {
            return Err(HgvsError::ValidationError(format!(
                "{ac} has no bases over [{first}, {})",
                last + 1
            )));
        }
        Ok(bases)
    }
}

/// Which strand a duplication was written on, which decides the side of the
/// original its copy goes.
#[derive(Clone, Copy)]
enum Written {
    OnGenome,
    OnTranscript,
}

/// Refuses an allele whose members are not in the `system` that `method` projects from.
pub(crate) fn expect_system(
    cis: &CisPhasedVariant,
    system: &str,
    method: &str,
) -> Result<(), HgvsError> {
    if cis.coordinate_type() == system {
        return Ok(());
    }
    Err(HgvsError::UnsupportedOperation(format!(
        "{method} projects an allele of {system}. members, not {cis}"
    )))
}

fn wrong_member(member: &SV, system: &str) -> HgvsError {
    HgvsError::UnsupportedOperation(format!("{member} is not a {system}. variant"))
}

pub(crate) fn member_error(cis: &CisPhasedVariant, index: usize, source: HgvsError) -> HgvsError {
    HgvsError::CisMember {
        allele: cis.to_string(),
        index,
        member: cis.members[index].to_string(),
        source: Box::new(source),
    }
}

/// `project` over every member, in order; the first failure names its member.
fn project_members<'c, T>(
    cis: &'c CisPhasedVariant,
    project: impl Fn(&'c SV) -> Result<T, HgvsError>,
) -> Result<Vec<T>, HgvsError> {
    cis.members
        .iter()
        .enumerate()
        .map(|(i, m)| project(m).map_err(|e| member_error(cis, i, e)))
        .collect()
}

/// The allele in cis of `members`, on their accession, with `like`'s gene.
fn allele(like: &CisPhasedVariant, members: Vec<SV>) -> Result<CisPhasedVariant, HgvsError> {
    let ac = members.first().map(|m| m.ac().to_string()).ok_or_else(|| {
        HgvsError::ValidationError("A cis allele needs at least one member".into())
    })?;
    CisPhasedVariant::new(ac, like.gene.clone(), members)
}
