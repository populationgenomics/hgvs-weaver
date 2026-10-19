//! The protein an allele in cis of coding members predicts, `c.[a;b]`.
//!
//! The members are applied to the transcript together and the protein is
//! read from the sequence they make, as a single variant's is from its own,
//! so `c.[4A>T;6A>T]` predicts what `c.4_6delinsTAT` predicts, not what each
//! member predicts alone. The change is then written as HGVS writes changes
//! in cis, those apart as separate members of a protein allele,
//! `p.[(Lys2Gln;Tyr5Asp)]`. The bases the allele states unchanged, by an
//! identity or by members that cancel, add the residues whose whole codon
//! they cover and the changes leave as they were, as members beside them,
//! `p.[(Ala4=;Tyr5Asp)]`.

use crate::data::IdentifierType;
use crate::edits::{AaEdit, ProteinStatement, ResolvedEdit};
use crate::error::HgvsError;
use crate::mapper::{
    cds_bounds, covered_cds_index, has_intronic_offset, place_on_coding, predicted_statement,
    whole_codons, CdsResidues, Placement, PredictedProtein, VariantMapper,
};
use crate::protein::{CodingChange, Reading};
use crate::reading::CisReading;
use crate::structs::{CVariant, CisPhasedVariant, PVariant};
use crate::SequenceVariant as SV;

/// What an allele of more than one coding member predicts.
enum Prediction {
    /// A statement about the whole protein, which nothing stands beside:
    /// `p.?`, `p.(Met1?)`, `p.0?`, `p.(=)`.
    Statement(PVariant),
    /// The change the members make together, read as one; the changes
    /// HGVS writes as separate members, each described; and the residues
    /// stated unchanged beside them.
    Change {
        reading: Reading,
        changes: Vec<PVariant>,
        beside: Vec<PVariant>,
    },
    /// The allele makes the reference: the residues it states unchanged.
    Stated(Vec<PVariant>),
}

/// The CDS an allele is predicted on.
struct Cds<'c> {
    cis: &'c CisPhasedVariant,
    protein_ac: String,
    ref_seq: String,
    /// The 0-based transcript index of the first CDS base.
    start: usize,
    /// Bases from the first CDS base through the stop codon.
    len: usize,
    end_open: bool,
    /// Whether translation initiates at the first codon.
    initiates: bool,
}

impl Cds<'_> {
    fn statement(&self, pos: Option<crate::structs::AaInterval>, s: ProteinStatement) -> PVariant {
        predicted_statement(self.cis.gene.clone(), &self.protein_ac, pos, s)
    }

    fn residues(&self) -> CdsResidues<'_> {
        CdsResidues::new(
            &self.cis.ac,
            &self.ref_seq,
            self.start,
            self.len,
            self.initiates,
        )
    }

    /// The prediction that the residues of the codons `[from, to)` are unchanged.
    fn unchanged(&self, (from, to): (usize, usize)) -> Result<PVariant, HgvsError> {
        self.residues()
            .unchanged(self.protein_ac.clone(), self.cis.gene.clone(), from, to)
    }

    /// The change `edits`, sorted and apart, make, relative to the first CDS base.
    fn change(&self, edits: &[ResolvedEdit]) -> CodingChange {
        let coding = &self.ref_seq[self.start..];
        CodingChange {
            coding: coding.to_string(),
            cds_len: self.len,
            end_open: self.end_open,
            edit: combined(coding, edits),
            protein_ac: self.protein_ac.clone(),
            initiates: self.initiates,
        }
    }
}

/// The one edit over `seq` from the first of `edits` to the last, the bases
/// between them as they are.
fn combined(seq: &str, edits: &[ResolvedEdit]) -> ResolvedEdit {
    let (start, end) = (edits[0].start, edits[edits.len() - 1].end);
    let mut alt = String::new();
    let mut at = start;
    for e in edits {
        alt.push_str(&seq[at..e.start]);
        alt.push_str(&e.alt);
        at = e.end;
    }
    ResolvedEdit {
        start,
        end,
        ref_: seq[start..end].to_string(),
        alt,
    }
}

impl VariantMapper<'_> {
    /// The protein consequence of an allele in cis of `c.` members. The
    /// members are applied to the transcript together and the protein is
    /// read from the sequence they make, by the rules
    /// [`c_to_p`](Self::c_to_p) reads one change by: `c.[4A>T;6A>T]` is
    /// `p.(Lys2Tyr)`, as `c.4_6delinsTAT` is. Changes apart, an untouched
    /// codon between them, are the members of a protein allele,
    /// `c.[4A>C;13T>G]` giving `p.[(Lys2Gln;Tyr5Asp)]`; residues whose whole
    /// codon the allele states unchanged, and the changes leave as they were,
    /// are members beside them, `c.[13T>G;10_12=]` giving
    /// `p.[(Ala4=;Tyr5Asp)]`. One prediction alone is written as a single
    /// variant. An allele of one member predicts what its member does.
    ///
    /// # Errors
    ///
    /// [`HgvsError::NoProteinStatement`] for an allele that makes the
    /// reference and states no whole codon; [`HgvsError::CisMember`] naming
    /// a member no sequence can be built from, an uncertain position among
    /// them; [`HgvsError::UnsupportedOperation`] for two members whose
    /// 3'-most spellings overlap, or a member past the stop that shifts the
    /// frame of the read-through; and what [`read_cis`](Self::read_cis)
    /// raises for a contradictory or ambiguous allele.
    pub fn c_to_p_cis(
        &self,
        cis: &CisPhasedVariant,
        protein_ac: Option<&str>,
    ) -> Result<SV, HgvsError> {
        crate::cis::expect_system(cis, "c", "c_to_p")?;
        if let [SV::Coding(one)] = cis.members.as_slice() {
            return self.c_to_p(one, protein_ac).map(SV::Protein);
        }
        let mut members = match self.cis_prediction(cis, protein_ac)? {
            Prediction::Statement(p) => return Ok(SV::Protein(p)),
            Prediction::Stated(stated) => stated,
            Prediction::Change {
                changes, beside, ..
            } => {
                let mut members = beside;
                for mut change in changes {
                    change.gene = cis.gene.clone();
                    change.posedit.predicted = true;
                    members.push(change);
                }
                members
            }
        };
        members.sort_by_key(|p| {
            p.posedit
                .pos
                .as_ref()
                .map(|pos| pos.start.base.to_index().0)
        });
        if members.len() == 1 {
            return Ok(SV::Protein(members.remove(0)));
        }
        let ac = members[0].ac.clone();
        Ok(SV::CisPhased(CisPhasedVariant::new(
            ac,
            cis.gene.clone(),
            members.into_iter().map(SV::Protein).collect(),
        )?))
    }

    /// The protein an allele in cis of `c.` members leaves, as
    /// [`predicted_protein`](Self::predicted_protein) reads one variant's:
    /// no residues when it deletes the CDS, those its changes leave whatever
    /// it states beside them, and `None` when its prediction is a statement
    /// (`p.?`, `p.(Met1?)`) or it changes no base.
    pub fn predicted_protein_cis(
        &self,
        cis: &CisPhasedVariant,
        protein_ac: Option<&str>,
    ) -> Result<Option<PredictedProtein>, HgvsError> {
        crate::cis::expect_system(cis, "c", "predicted_protein")?;
        if let [SV::Coding(one)] = cis.members.as_slice() {
            return self.predicted_protein(one, protein_ac);
        }
        Ok(match self.cis_prediction(cis, protein_ac)? {
            Prediction::Statement(p)
                if matches!(
                    p.posedit.edit,
                    AaEdit::Statement(ProteinStatement::NoProtein { .. })
                ) =>
            {
                Some(PredictedProtein {
                    residues: String::new(),
                    stops: true,
                })
            }
            Prediction::Change { reading, .. } => Some(PredictedProtein {
                residues: reading.proteins()?.1,
                stops: reading.alternate_stops(),
            }),
            _ => None,
        })
    }

    fn cis_prediction(
        &self,
        cis: &CisPhasedVariant,
        protein_ac: Option<&str>,
    ) -> Result<Prediction, HgvsError> {
        let protein_ac = self.protein_accession(&cis.ac, protein_ac)?;
        let transcript = self.provider().get_transcript(&cis.ac, None)?;
        if transcript.cds_start == crate::data::CdsStart::Open {
            return Err(HgvsError::ValidationError(format!(
                "{}: its CDS is open at the 5' end, so the record carries no start codon and no protein \
                 consequence can be predicted",
                cis.ac
            )));
        }
        let whole = self.read_cis(cis)?;
        refuse_unapplied(cis, &whole)?;
        let ref_seq = self
            .refs
            .reference(&cis.ac, IdentifierType::TranscriptAccession)
            .whole()?;
        let (start, end) = cds_bounds(&transcript, ref_seq.len())?;
        let cds = Cds {
            cis,
            protein_ac,
            ref_seq,
            start,
            len: end + 1 - start,
            end_open: transcript.cds_end_open,
            initiates: transcript.cds_start == crate::data::CdsStart::Initiator,
        };
        let coding: Vec<&CVariant> = cis
            .members
            .iter()
            .filter_map(|m| match m {
                SV::Coding(c) => Some(c),
                _ => None,
            })
            .collect();
        let intronic = |c: &&CVariant| c.posedit.pos.as_ref().is_some_and(has_intronic_offset);
        // A change in an intron is read on no codon, so the protein it leaves
        // is not known, as for one alone.
        if coding
            .iter()
            .any(|c| intronic(c) && !c.posedit.edit.is_identity())
        {
            return Ok(Prediction::Statement(
                cds.statement(None, ProteinStatement::Unknown),
            ));
        }
        let reading = match whole.kind() {
            IdentifierType::TranscriptAccession => Some(whole),
            _ => self.read_exonic(cis, &coding, intronic)?,
        };
        let mut stated = intronic_statements(&coding, intronic, cds.len)?;
        if let Some(r) = &reading {
            stated.extend(
                r.stated()
                    .iter()
                    .map(|&(s, e)| (s as i64 - cds.start as i64, e as i64 - cds.start as i64)),
            );
        }
        let codons = stated_codons(stated, (cds.len / 3) as i64);
        let edits = match &reading {
            Some(r) => self.edits(r)?,
            None => Vec::new(),
        };
        if edits.is_empty() {
            return stated_only(&cds, &codons);
        }
        let (reading, edits) = match self.read_changes(&cds, edits)? {
            Changes::Statement(p) => return Ok(Prediction::Statement(p)),
            Changes::Change { reading, edits } => (reading, edits),
        };
        let Some(described) = described_apart(&cds, &edits)? else {
            return Ok(Prediction::Statement(
                cds.statement(None, ProteinStatement::Unknown),
            ));
        };
        let (silent, changes): (Vec<PVariant>, Vec<PVariant>) = described
            .into_iter()
            .filter(|p| p.posedit.pos.is_some())
            .partition(is_silent);
        if changes.is_empty() {
            // A silent change names the residues of the codons it touches
            // unchanged, as a statement does: the two are one run.
            let mut all = codons;
            all.extend(silent.iter().filter_map(named_residues));
            if all.is_empty() {
                let only = reading.describe()?;
                return Ok(Prediction::Change {
                    reading,
                    changes: vec![only],
                    beside: Vec::new(),
                });
            }
            return Ok(match stated_only(&cds, &merged(all))? {
                Prediction::Stated(beside) => Prediction::Change {
                    reading,
                    changes: Vec::new(),
                    beside,
                },
                other => other,
            });
        }
        // Beside a change that alters residues, a silent one says nothing the
        // change does not: the change predicts every residue it leaves.
        let beside = residues_beside(&changes, &codons)
            .into_iter()
            .map(|c| cds.unchanged(c))
            .collect::<Result<_, _>>()?;
        Ok(Prediction::Change {
            reading,
            changes,
            beside,
        })
    }

    /// The reading on the transcript of the members not in an intron, which
    /// a transcript allele with an intronic identity is read without; `None`
    /// when every member is intronic.
    fn read_exonic(
        &self,
        cis: &CisPhasedVariant,
        coding: &[&CVariant],
        intronic: impl Fn(&&CVariant) -> bool,
    ) -> Result<Option<CisReading>, HgvsError> {
        let exonic: Vec<SV> = coding
            .iter()
            .filter(|c| !intronic(c))
            .map(|&c| SV::Coding(c.clone()))
            .collect();
        if exonic.is_empty() {
            return Ok(None);
        }
        let sub = CisPhasedVariant::new(cis.ac.clone(), cis.gene.clone(), exonic)?;
        self.read_cis(&sub).map(Some)
    }

    /// The replacements of `reading` that change bases, as edits on the
    /// sequence it is read on; members that cancel make none.
    fn edits(&self, reading: &CisReading) -> Result<Vec<ResolvedEdit>, HgvsError> {
        let reference = self.refs.reference(reading.accession(), reading.kind());
        let mut edits = Vec::new();
        for r in reading.replacements() {
            let ref_ = reference.slice(r.start(), r.end())?.to_ascii_uppercase();
            let alt = reading.sequence(&reference, r.start(), r.end())?;
            if ref_ != alt {
                edits.push(ResolvedEdit {
                    start: r.start(),
                    end: r.end(),
                    ref_,
                    alt,
                });
            }
        }
        Ok(edits)
    }

    /// What the edits an allele makes on its transcript, sorted and apart,
    /// predict together: each placed on its 3'-most spelling as a single
    /// change is, then the ones in the CDS read as one change. A change to
    /// an initiator is `p.(Met1?)`, and one wholly in the 5'UTR beside
    /// changes that leave the first codon makes the protein `p.?`. An edit
    /// past the stop codon enters the reading only where translation of the
    /// edited sequence reaches it.
    fn read_changes(&self, cds: &Cds, edits: Vec<ResolvedEdit>) -> Result<Changes, HgvsError> {
        let statement = |pos, s| Ok(Changes::Statement(cds.statement(pos, s)));
        if edits
            .iter()
            .any(|e| e.start < cds.start && e.end >= cds.start + cds.len && e.alt.is_empty())
        {
            return statement(None, ProteinStatement::NoProtein { uncertain: true });
        }
        let mut upstream = false;
        let mut placed: Vec<ResolvedEdit> = Vec::new();
        for edit in edits {
            match place_on_coding(&cds.ref_seq, edit, cds.start, cds.initiates) {
                Placement::AcrossStart => {
                    return statement(Some(cds.residues().first()?), ProteinStatement::Unknown)
                }
                Placement::FivePrimeUtr => upstream = true,
                Placement::Cds(e) => placed.push(ResolvedEdit {
                    start: e.start - cds.start,
                    end: e.end - cds.start,
                    ..e
                }),
            }
        }
        if placed.windows(2).any(|w| w[0].end > w[1].start) {
            return Err(HgvsError::UnsupportedOperation(format!(
                "{}: two members, each written at its 3'-most place, overlap, so the codons they \
                 change together are not read",
                cds.cis
            )));
        }
        if placed.is_empty() {
            let s = match upstream {
                true => ProteinStatement::Unknown,
                false => ProteinStatement::Unchanged,
            };
            return statement(None, s);
        }
        if cds.end_open && !upstream {
            Reading::new(cds.change(&placed)).check_within_record(&cds.cis.ac)?;
        }
        let (reading, read) = read_through(cds, &placed)?;
        if reading.runs_off() {
            return statement(None, ProteinStatement::Unknown);
        }
        if reading.changes_start() {
            let mut var_p = reading.describe()?;
            var_p.gene = cds.cis.gene.clone();
            var_p.posedit.predicted = true;
            return Ok(Changes::Statement(var_p));
        }
        if upstream {
            return statement(None, ProteinStatement::Unknown);
        }
        Ok(Changes::Change {
            reading,
            edits: placed[..read].to_vec(),
        })
    }
}

/// What the edits an allele makes predict together.
enum Changes {
    /// A statement about the whole protein: `p.?`, `p.(Met1?)`, `p.0?`, `p.(=)`.
    Statement(PVariant),
    /// The change the edits make, read together, and the edits, relative to
    /// the first CDS base, that the reading holds.
    Change {
        reading: Reading,
        edits: Vec<ResolvedEdit>,
    },
}

/// The reading of the edits `placed`, relative to the first CDS base: those
/// that start before the end of the stop codon, and those past it only if
/// translation of the sequence the first make reaches them.
///
/// # Errors
///
/// [`HgvsError::UnsupportedOperation`] when an edit translation reaches
/// past the stop codon turns the edits before it from in frame to out of
/// frame, or back: the
/// read-through is then neither the frameshift nor the in-frame change those
/// edits make, and is not predicted. The one exception is a frame those
/// edits shift at the stop itself and the later edit restores: that is an
/// extension, read in frame.
fn read_through(cds: &Cds, placed: &[ResolvedEdit]) -> Result<(Reading, usize), HgvsError> {
    let in_cds = placed.iter().take_while(|e| e.start < cds.len).count();
    if in_cds == 0 {
        return Ok((Reading::new(cds.change(placed)), placed.len()));
    }
    let early = Reading::new(cds.change(&placed[..in_cds]));
    let Some(late) = placed.get(in_cds) else {
        return Ok((early, in_cds));
    };
    if early
        .translation_end()
        .is_some_and(|end| end <= late.start as i64)
    {
        return Ok((early, in_cds));
    }
    let all = Reading::new(cds.change(placed));
    if all.in_frame() != early.in_frame() && (all.changes_before_stop() || !all.in_frame()) {
        return Err(HgvsError::UnsupportedOperation(format!(
            "{}: a member past the stop codon changes the frame the read-through is translated in, \
             so the protein it makes is not predicted",
            cds.cis
        )));
    }
    Ok((all, placed.len()))
}

/// Refuses an allele with a member no sequence can be built from: an
/// uncertain position, a copy number, a conversion, inserted bases not
/// stated. Its codons cannot be read, so neither can the protein.
fn refuse_unapplied(cis: &CisPhasedVariant, reading: &CisReading) -> Result<(), HgvsError> {
    let Some(member) = reading.unapplied().first() else {
        return Ok(());
    };
    let index = cis
        .members
        .iter()
        .position(|m| m == member)
        .ok_or_else(|| HgvsError::Other(format!("{member} is not a member of {cis}")))?;
    let why = match member.has_exact_position() {
        false => "it has an uncertain position, so it names no codons to translate",
        true => "it makes no bases to apply, so the codons the allele leaves cannot be read",
    };
    Err(crate::cis::member_error(
        cis,
        index,
        HgvsError::UnsupportedOperation(why.into()),
    ))
}

/// The coding bases each intronic identity covers, 0-based half-open from
/// the first base of the CDS.
fn intronic_statements(
    coding: &[&CVariant],
    intronic: impl Fn(&&CVariant) -> bool,
    cds_len: usize,
) -> Result<Vec<(i64, i64)>, HgvsError> {
    coding
        .iter()
        .filter(|c| intronic(c))
        .filter_map(|c| c.posedit.pos.as_ref())
        .map(|pos| {
            let first = covered_cds_index(&pos.start, cds_len as i64, true)?;
            let last = covered_cds_index(
                pos.end.as_ref().unwrap_or(&pos.start),
                cds_len as i64,
                false,
            )?;
            Ok((first, last + 1))
        })
        .collect()
}

/// The ranges, sorted, with those that overlap or touch made one.
fn merged<T: Ord + Copy>(mut ranges: Vec<(T, T)>) -> Vec<(T, T)> {
    ranges.sort_unstable();
    let mut merged: Vec<(T, T)> = Vec::with_capacity(ranges.len());
    for (from, to) in ranges {
        match merged.last_mut() {
            Some(last) if from <= last.1 => last.1 = last.1.max(to),
            _ => merged.push((from, to)),
        }
    }
    merged
}

/// The whole codons the stated coding bases cover: the stretches of bases
/// made one first, so that statements which together cover a codon state it.
fn stated_codons(stated: Vec<(i64, i64)>, codons: i64) -> Vec<(usize, usize)> {
    merged(stated)
        .into_iter()
        .filter_map(|(first, end)| whole_codons(first, end, codons))
        .collect()
}

/// What an allele that makes the reference predicts: the residues of the
/// whole codons it states, or `p.(=)` when that is every codon of a CDS that
/// ends in a stop.
fn stated_only(cds: &Cds, codons: &[(usize, usize)]) -> Result<Prediction, HgvsError> {
    if codons.is_empty() {
        return Err(HgvsError::NoProteinStatement(format!(
            "{} states no whole codon, so it states no residue of {}",
            cds.cis, cds.protein_ac
        )));
    }
    if codons == [(0, cds.len / 3)] && !cds.end_open {
        return Ok(Prediction::Statement(
            cds.statement(None, ProteinStatement::Unchanged),
        ));
    }
    codons
        .iter()
        .map(|&c| cds.unchanged(c))
        .collect::<Result<_, _>>()
        .map(Prediction::Stated)
}

/// The residues a described change names, as a 0-based half-open range:
/// those its position names, and every residue after its first when it
/// shifts the frame or ends the protein, since those are read anew or not at
/// all. `None` for a description that names no residue, `p.(=)`.
fn named_residues(change: &PVariant) -> Option<(usize, usize)> {
    let pos = change.posedit.pos.as_ref()?;
    let index = |a: &crate::structs::AAPosition| a.base.to_index().0 as usize;
    let first = index(&pos.start);
    let ends_protein = match &change.posedit.edit {
        AaEdit::Fs { .. } | AaEdit::Ext { .. } => true,
        AaEdit::Subst { alt, .. } | AaEdit::DelIns { alt, .. } | AaEdit::Ins { alt, .. } => {
            alt.contains("Ter")
        }
        AaEdit::Statement(_) => true,
        _ => false,
    };
    Some(match ends_protein {
        true => (first, usize::MAX),
        false => (first, pos.end.as_ref().map_or(first, index) + 1),
    })
}

/// The parts of the stated codon ranges `codons` that none of the predicted
/// `changes` names.
fn residues_beside(changes: &[PVariant], codons: &[(usize, usize)]) -> Vec<(usize, usize)> {
    let mut left: Vec<(usize, usize)> = codons.to_vec();
    for (first, after) in changes.iter().filter_map(named_residues) {
        left = left
            .into_iter()
            .flat_map(|(from, to)| [(from, to.min(first)), (from.max(after), to)])
            .filter(|(from, to)| from < to)
            .collect();
    }
    left
}

/// Whether a described change is silent: it names residues and leaves them
/// as they are.
fn is_silent(change: &PVariant) -> bool {
    matches!(change.posedit.edit, AaEdit::Identity { .. }) && change.posedit.pos.is_some()
}

/// The edits `edits`, sorted and apart and relative to the first CDS base,
/// grouped into the changes HGVS writes as separate members: a group ends
/// where the frame is as the reference's and at least one codon between it
/// and the next edit is untouched.
fn groups(edits: &[ResolvedEdit]) -> Vec<Vec<ResolvedEdit>> {
    let mut groups: Vec<Vec<ResolvedEdit>> = Vec::new();
    let mut net: i64 = 0;
    for e in edits {
        let last_codon = groups
            .last()
            .and_then(|g| g.last())
            .map(|l| (l.end.max(l.start + 1) - 1) / 3);
        let apart = last_codon.is_some_and(|c| e.start / 3 > c + 1);
        match groups.last_mut() {
            Some(g) if !(apart && net % 3 == 0) => g.push(e.clone()),
            _ => groups.push(vec![e.clone()]),
        }
        net += e.alt.len() as i64 - (e.end - e.start) as i64;
    }
    groups
}

/// The changes `edits` make, each group of them read on its own and
/// described, as the members HGVS writes for changes in cis: a group that
/// shifts the frame or ends the protein is read with every edit after it,
/// since what follows is read anew, and groups whose descriptions meet are
/// read as one. `None` when a group leaves no codon where its change
/// starts, so that nothing can be said from there on.
fn described_apart(cds: &Cds, edits: &[ResolvedEdit]) -> Result<Option<Vec<PVariant>>, HgvsError> {
    let describe = |group: &[ResolvedEdit]| {
        let reading = Reading::new(cds.change(group));
        match reading.runs_off() {
            true => Ok(None),
            false => reading.describe().map(Some),
        }
    };
    let mut pending: std::collections::VecDeque<Vec<ResolvedEdit>> = groups(edits).into();
    let mut described: Vec<(Vec<ResolvedEdit>, PVariant)> = Vec::new();
    while let Some(mut group) = pending.pop_front() {
        let Some(mut change) = describe(&group)? else {
            return Ok(None);
        };
        if named_residues(&change).is_some_and(|(_, after)| after == usize::MAX) {
            group.extend(pending.drain(..).flatten());
            let Some(whole) = describe(&group)? else {
                return Ok(None);
            };
            change = whole;
        }
        let meets = |previous: &PVariant| match (named_residues(previous), named_residues(&change))
        {
            (Some((_, after)), Some((first, _))) => first <= after,
            _ => false,
        };
        if let Some((previous, _)) = described.pop_if(|(_, p)| meets(p)) {
            let mut joined = previous;
            joined.extend(group);
            pending.push_front(joined);
            continue;
        }
        described.push((group, change));
    }
    if cds.end_open {
        for (group, _) in &described {
            Reading::new(cds.change(group)).check_within_record(&cds.cis.ac)?;
        }
    }
    Ok(Some(described.into_iter().map(|(_, p)| p).collect()))
}
