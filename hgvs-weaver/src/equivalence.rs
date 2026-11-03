use crate::allele::{CanonicalAllele, StatedBases};
use crate::data::{IdentifierKind, TranscriptSearch};
use crate::edits::{AaEdit, NaEdit};
use crate::error::HgvsError;
use crate::mapper::VariantMapper;
use crate::protein_reading::{first_difference, same_protein, Outcome, ProteinRead, Stated};
use crate::reading::CisReading;
use crate::structs::{
    AaInterval, CVariant, CisPhasedVariant, LinearVariant, PVariant, PosEdit, SequenceVariant,
    Variant,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EquivalenceLevel {
    /// Identical notation after basic normalization.
    Identity,
    /// Biologically identical but different notation (e.g., ins vs dup).
    Analogous,
    /// Definitively different edits/outcomes.
    Different,
    /// Missing data or unsupported variant type for comparison.
    Unknown,
}

impl EquivalenceLevel {
    pub fn is_equivalent(&self) -> bool {
        matches!(self, Self::Identity | Self::Analogous)
    }
}

#[non_exhaustive]
pub struct VariantEquivalence<'a> {
    /// The mapper whose provider, cache and refget lookup the comparison uses.
    pub mapper: &'a VariantMapper<'a>,
    pub searcher: &'a dyn TranscriptSearch,
}

/// A coding variant, or an allele in cis of coding members, whose prediction
/// is compared with a protein description.
enum Coding {
    One(CVariant),
    Cis(CisPhasedVariant),
}

fn mito_as_genomic(var: &SequenceVariant) -> std::borrow::Cow<'_, SequenceVariant> {
    match var {
        SequenceVariant::Mitochondrial(m) => {
            std::borrow::Cow::Owned(SequenceVariant::Genomic(m.to_genomic()))
        }
        SequenceVariant::CisPhased(cis) if cis.coordinate_type() == "m" => {
            std::borrow::Cow::Owned(SequenceVariant::CisPhased(CisPhasedVariant {
                ac: cis.ac.clone(),
                gene: cis.gene.clone(),
                members: cis
                    .members
                    .iter()
                    .map(|m| mito_as_genomic(m).into_owned())
                    .collect(),
            }))
        }
        other => std::borrow::Cow::Borrowed(other),
    }
}

impl<'a> VariantEquivalence<'a> {
    pub fn new(mapper: &'a VariantMapper<'a>, searcher: &'a dyn TranscriptSearch) -> Self {
        VariantEquivalence { mapper, searcher }
    }

    pub fn equivalent(
        &self,
        var1: &SequenceVariant,
        var2: &SequenceVariant,
    ) -> Result<bool, HgvsError> {
        Ok(self.equivalent_level(var1, var2)?.is_equivalent())
    }

    pub fn equivalent_level(
        &self,
        var1: &SequenceVariant,
        var2: &SequenceVariant,
    ) -> Result<EquivalenceLevel, HgvsError> {
        // A mitochondrial variant is a genomic variant on the mitochondrial
        // reference; compare it as one.
        let var1 = mito_as_genomic(var1);
        let var2 = mito_as_genomic(var2);
        // An r. allele of one member is that member, which is read in c. or
        // n. letters below as an allele's members are not.
        let var1 = rna_alone(var1);
        let var2 = rna_alone(var2);
        // A gene symbol is read on each of the gene's accessions of its kind,
        // before an r. description is converted on its transcript.
        let read1 = self.held(self.readings(&var1)?)?;
        let read2 = self.held(self.readings(&var2)?)?;
        // An r. change across a splice junction has no genomic form, so it is
        // compared only with its own spelling or by the protein it leaves. A
        // reading refused so takes no part, and is the error only when no
        // other reading makes the two the same.
        let same_spelling = spelling(&var1) == spelling(&var2);
        let mut refusals = Vec::new();
        let mut comparable = |read: &[(SequenceVariant, SequenceVariant)],
                              other: &SequenceVariant|
         -> Result<Vec<SequenceVariant>, HgvsError> {
            let mut kept = Vec::with_capacity(read.len());
            for (written, converted) in read {
                if same_spelling || other.coordinate_type() == "p" {
                    kept.push(converted.clone());
                    continue;
                }
                match self.rna_on_the_genome(written) {
                    Ok(()) => kept.push(converted.clone()),
                    Err(e) if is_a_refusal(&e) => refusals.push((spelling(written), e)),
                    Err(e) => return Err(e),
                }
            }
            Ok(kept)
        };
        let vars1 = comparable(&read1, &var2)?;
        let vars2 = comparable(&read2, &var1)?;

        for v1 in &vars1 {
            for v2 in &vars2 {
                let lvl = self.equivalent_level_single(v1, v2)?;
                if lvl.is_equivalent() {
                    return Ok(lvl);
                }
            }
        }
        // The refusal names one reading whichever order the two come in.
        if let Some((_, e)) = refusals.into_iter().min_by(|a, b| a.0.cmp(&b.0)) {
            return Err(e);
        }
        Ok(EquivalenceLevel::Different)
    }

    /// Refuses an r. change, alone or a member of an allele, that spans a
    /// splice junction.
    fn rna_on_the_genome(&self, var: &SequenceVariant) -> Result<(), HgvsError> {
        let on_the_genome = |r: &crate::structs::RVariant| {
            if r.posedit.edit.is_whole_transcript() {
                Ok(())
            } else {
                self.mapper.r_on_the_genome(r).map(|_| ())
            }
        };
        match var {
            SequenceVariant::Rna(r) => on_the_genome(r),
            SequenceVariant::CisPhased(cis) => {
                for (i, member) in cis.members.iter().enumerate() {
                    if let SequenceVariant::Rna(r) = member {
                        on_the_genome(r).map_err(|e| crate::cis::member_error(cis, i, e))?;
                    }
                }
                Ok(())
            }
            _ => Ok(()),
        }
    }

    /// An r. variant as the c. or n. variant it is spelled from; any other
    /// variant as it is.
    fn rna_as_transcript(&self, var: &SequenceVariant) -> Result<SequenceVariant, HgvsError> {
        match var {
            SequenceVariant::Rna(r) if !r.posedit.edit.is_whole_transcript() => {
                self.mapper.r_as_transcript(r)
            }
            var => Ok(var.clone()),
        }
    }

    /// `var` on each of a gene's accessions of its kind where it names a
    /// gene symbol, else `var` alone.
    fn on_each_gene_accession(
        &self,
        var: &SequenceVariant,
    ) -> Result<Vec<SequenceVariant>, HgvsError> {
        if self.mapper.refs.identifier_type(var.ac())? == crate::data::IdentifierType::GeneSymbol {
            self.mapper.on_each_of_the_genes_accessions(var)
        } else {
            Ok(vec![var.clone()])
        }
    }

    /// The readings of a description, one per accession a gene symbol expands
    /// to, each as written and as compared: an r. reading in its c. or n.
    /// spelling.
    fn readings(
        &self,
        var: &SequenceVariant,
    ) -> Result<Vec<(SequenceVariant, SequenceVariant)>, HgvsError> {
        self.on_each_gene_accession(var)?
            .into_iter()
            .map(|written| {
                let converted = self.rna_as_transcript(&written)?;
                Ok((written, converted))
            })
            .collect()
    }

    /// The readings whose references hold what they state ([`check_stated`]).
    /// Each accession of a gene numbers its own sequence, so a reading the
    /// reference contradicts is set aside; when every reading is contradicted,
    /// the first contradiction is the error.
    ///
    /// [`check_stated`]: Self::check_stated
    fn held(
        &self,
        readings: Vec<(SequenceVariant, SequenceVariant)>,
    ) -> Result<Vec<(SequenceVariant, SequenceVariant)>, HgvsError> {
        let mut held = Vec::with_capacity(readings.len());
        let mut contradiction = None;
        for reading in readings {
            match self.check_stated(&reading.1) {
                Ok(()) => held.push(reading),
                Err(
                    e @ (HgvsError::ReferenceMismatch { .. }
                    | HgvsError::PositionOutsideSequence { .. }
                    | HgvsError::NotAnExonBorder { .. }
                    | HgvsError::OffsetPastIntron { .. }),
                ) => {
                    contradiction.get_or_insert(e);
                }
                Err(e) => return Err(e),
            }
        }
        match contradiction {
            Some(e) if held.is_empty() => Err(e),
            _ => Ok(held),
        }
    }

    /// Every position `var` names and every base and residue it states, checked
    /// against its reference as [`validate`](VariantMapper::validate) checks
    /// them, so a description the reference contradicts is refused even
    /// against its own spelling. An allele is read whole, an uncertain member
    /// on its outer bounds, and a contradictory or ambiguous one is refused;
    /// one whose readings weaver does not list has had its members checked.
    fn check_stated(&self, var: &SequenceVariant) -> Result<(), HgvsError> {
        match self.mapper.check_stated(var) {
            Err(HgvsError::UnsupportedOperation(_))
                if matches!(var, SequenceVariant::CisPhased(_)) =>
            {
                Ok(())
            }
            checked => checked,
        }
    }

    /// Two variants are `Identity` when they are the same text after spelling
    /// normalisation, `Analogous` when they name the same change (the same
    /// [`Comparand`] on the genome, or the same protein left behind) and
    /// `Different` otherwise.
    fn equivalent_level_single(
        &self,
        var1: &SequenceVariant,
        var2: &SequenceVariant,
    ) -> Result<EquivalenceLevel, HgvsError> {
        if spelling(var1) == spelling(var2) {
            return Ok(EquivalenceLevel::Identity);
        }
        if let Some(level) = self.cis_phased_level(var1, var2)? {
            return Ok(level);
        }
        let same = match (var1, var2) {
            (SequenceVariant::Protein(p1), SequenceVariant::Protein(p2)) => {
                // Read each on the other's version too; Analogous only if both agree.
                let level = self.protein_level(p1, || self.protein_outcome_on(p1, &p1.ac), p2)?;
                if level == EquivalenceLevel::Analogous && p1.ac != p2.ac {
                    return self.protein_level(p2, || self.protein_outcome_on(p2, &p2.ac), p1);
                }
                return Ok(level);
            }
            (protein @ SequenceVariant::Protein(_), nucleotide)
            | (nucleotide, protein @ SequenceVariant::Protein(_)) => {
                return self.nucleotide_vs_protein(nucleotide, protein);
            }
            _ => {
                // A transcript variant that projects to exactly the genomic text.
                if self.exact_projection(var1, var2)? || self.exact_projection(var2, var1)? {
                    return Ok(EquivalenceLevel::Identity);
                }
                match (self.comparand(var1)?, self.comparand(var2)?) {
                    (Some(a), Some(b)) => a == b,
                    _ => false,
                }
            }
        };
        Ok(if same {
            EquivalenceLevel::Analogous
        } else {
            EquivalenceLevel::Different
        })
    }

    /// The level when either side is an allele in cis, `None` when neither
    /// is. A cis allele of one member is that member. Otherwise the two are
    /// compared by the [`AlleleComparand`] each makes on the genome, a plain
    /// variant read as an allele of one member: `Identity` when a transcript
    /// allele projects to the genomic one's text, `Analogous` when the
    /// comparands are equal, `Different` otherwise. An allele that cannot be
    /// read into a sequence, of protein members or with a member that is
    /// uncertain or a statement, is compared by the set of its members'
    /// [`MemberKey`]s, and against a plain variant is `Different`. An allele
    /// of nucleotide members against a protein description is compared by
    /// its prediction, as a single variant is.
    ///
    /// # Errors
    ///
    /// The error reading either allele on the genome gives: a contradictory
    /// allele's [`HgvsError::ValidationError`], an ambiguous one's
    /// [`HgvsError::Ambiguous`], a member that does not project.
    fn cis_phased_level(
        &self,
        var1: &SequenceVariant,
        var2: &SequenceVariant,
    ) -> Result<Option<EquivalenceLevel>, HgvsError> {
        let (cis, other) = match (var1, var2) {
            (SequenceVariant::CisPhased(cis), other) | (other, SequenceVariant::CisPhased(cis)) => {
                (cis, other)
            }
            _ => return Ok(None),
        };
        if let [member] = cis.members.as_slice() {
            return self
                .equivalent_level_single(&mito_as_genomic(member), other)
                .map(Some);
        }
        if let SequenceVariant::CisPhased(one) = other {
            if let [member] = one.members.as_slice() {
                return self
                    .equivalent_level_single(
                        &SequenceVariant::CisPhased(cis.clone()),
                        &mito_as_genomic(member),
                    )
                    .map(Some);
            }
        }
        let allele = SequenceVariant::CisPhased(cis.clone());
        match (cis.coordinate_type() == "p", other.coordinate_type() == "p") {
            (false, true) => return self.nucleotide_vs_protein(&allele, other).map(Some),
            (true, false) => return self.nucleotide_vs_protein(other, &allele).map(Some),
            (true, true) => return self.protein_allele_level(&allele, other).map(Some),
            _ => {}
        }
        let as_allele = |var: &SequenceVariant| -> Result<Option<CisPhasedVariant>, HgvsError> {
            Ok(match var {
                SequenceVariant::CisPhased(c) => Some(c.clone()),
                SequenceVariant::Protein(_) => None,
                one => Some(CisPhasedVariant::new(
                    one.ac().to_string(),
                    one.gene().map(str::to_string),
                    vec![one.clone()],
                )?),
            })
        };
        let (Some(a), Some(b)) = (as_allele(var1)?, as_allele(var2)?) else {
            return Ok(Some(EquivalenceLevel::Different));
        };
        if !(reads_as_a_sequence(&a) && reads_as_a_sequence(&b)) {
            let both = matches!(
                (var1, var2),
                (SequenceVariant::CisPhased(_), SequenceVariant::CisPhased(_))
            );
            let same = both
                && matches!(
                    (self.member_keys(&a)?, self.member_keys(&b)?),
                    (Some(x), Some(y)) if x == y
                );
            return Ok(Some(if same {
                EquivalenceLevel::Analogous
            } else {
                EquivalenceLevel::Different
            }));
        }
        if self.exact_cis_projection(&a, &b)? || self.exact_cis_projection(&b, &a)? {
            return Ok(Some(EquivalenceLevel::Identity));
        }
        let same = self
            .allele_comparand(&a)?
            .same(&self.allele_comparand(&b)?, self)?;
        Ok(Some(if same {
            EquivalenceLevel::Analogous
        } else {
            EquivalenceLevel::Different
        }))
    }

    /// A protein allele against another protein allele or a single p.
    /// description, on one protein accession (any version, read on the
    /// allele's): `Analogous` when the two can leave the same protein and
    /// state the same residues unchanged, however the change is split between
    /// members. Where a member names no stretch of residues (an uncertain
    /// range, `p.Met1?`, `p.?`), two alleles are compared by the set of their
    /// members, and an allele with a single description is `Different`.
    fn protein_allele_level(
        &self,
        allele: &SequenceVariant,
        other: &SequenceVariant,
    ) -> Result<EquivalenceLevel, HgvsError> {
        if base_accession(allele.ac()) != base_accession(other.ac()) {
            return Ok(EquivalenceLevel::Different);
        }
        let differ = allele.ac() != other.ac();
        let same_on = |ac: &str| -> Result<bool, HgvsError> {
            let (Some(a), Some(b)) = (
                read_across_versions(self.mapper.protein_comparand(allele, ac), differ)?,
                read_across_versions(self.mapper.protein_comparand(other, ac), differ)?,
            ) else {
                return Ok(false);
            };
            Ok(match (a, b) {
                (Some(a), Some(b)) => a.same(b),
                _ => match (allele, other) {
                    (SequenceVariant::CisPhased(a), SequenceVariant::CisPhased(b)) => matches!(
                        (self.member_keys(a)?, self.member_keys(b)?),
                        (Some(x), Some(y)) if x == y
                    ),
                    _ => false,
                },
            })
        };
        // On two versions, the same change on each.
        let same = same_on(allele.ac())? && (!differ || same_on(other.ac())?);
        Ok(if same {
            EquivalenceLevel::Analogous
        } else {
            EquivalenceLevel::Different
        })
    }

    /// Whether `tx`, a transcript allele, projected member by member onto
    /// `g`'s genome, is `g`'s text, in any order of the members.
    fn exact_cis_projection(
        &self,
        g: &CisPhasedVariant,
        tx: &CisPhasedVariant,
    ) -> Result<bool, HgvsError> {
        if g.coordinate_type() != "g" {
            return Ok(false);
        }
        let projected = match tx.coordinate_type() {
            "c" => self.mapper.c_to_g_cis(tx, Some(&g.ac)),
            "n" => self.mapper.n_to_g_cis(tx, Some(&g.ac)),
            "r" => self.mapper.r_to_g_cis(tx, Some(&g.ac)),
            _ => return Ok(false),
        };
        Ok(projected.is_ok_and(|p| {
            spelling(&SequenceVariant::CisPhased(p))
                == spelling(&SequenceVariant::CisPhased(g.clone()))
        }))
    }

    /// What an allele of nucleotide members is compared by, read on the
    /// genome.
    ///
    /// # Errors
    ///
    /// As [`read_cis_on_genome`](VariantMapper::read_cis_on_genome).
    fn allele_comparand(&self, cis: &CisPhasedVariant) -> Result<AlleleComparand, HgvsError> {
        let reading = self.mapper.read_cis_on_genome(cis)?;
        let mut unapplied = reading
            .unapplied()
            .iter()
            .map(|m| self.unapplied_key(m))
            .collect::<Result<Vec<_>, _>>()?;
        unapplied.sort_unstable();
        unapplied.dedup();
        Ok(AlleleComparand { reading, unapplied })
    }

    /// The value a member no sequence is built from is compared by: a copy
    /// number by its projection to the genome, as alone; anything else by
    /// the [`spelling`] of its genomic form.
    fn unapplied_key(&self, member: &SequenceVariant) -> Result<MemberKey, HgvsError> {
        match self.mapper.as_genomic(member) {
            None => Ok(MemberKey::Spelling(spelling(member))),
            Some(g) => {
                let g = SequenceVariant::Genomic(g?);
                match self.comparand(&g)? {
                    Some(Comparand::Projection(ac, posedit)) => {
                        Ok(MemberKey::Projection(ac, posedit))
                    }
                    _ => Ok(MemberKey::Spelling(spelling(&g))),
                }
            }
        }
    }

    /// The [`MemberKey`]s of a cis allele's members, as a sorted set; `None`
    /// when an exact member has no canonical allele.
    fn member_keys(&self, cis: &CisPhasedVariant) -> Result<Option<Vec<MemberKey>>, HgvsError> {
        let mut keys = Vec::with_capacity(cis.members.len());
        for m in &cis.members {
            if !m.has_exact_position() {
                let m = self.rna_as_transcript(&mito_as_genomic(m))?;
                keys.push(MemberKey::Spelling(spelling(&m)));
                continue;
            }
            match self.nucleotide_allele(m)? {
                Some(a) => keys.push(MemberKey::Allele(a.spdi())),
                None => return Ok(None),
            }
        }
        keys.sort_unstable();
        keys.dedup();
        Ok(Some(keys))
    }

    /// Whether `tx` is a c. or n. variant whose projection onto `g`'s
    /// reference is `g` to the letter.
    fn exact_projection(
        &self,
        g: &SequenceVariant,
        tx: &SequenceVariant,
    ) -> Result<bool, HgvsError> {
        let SequenceVariant::Genomic(g) = g else {
            return Ok(false);
        };
        let projected = match tx {
            SequenceVariant::Coding(c) => self.mapper.tx_to_g(c, Some(&g.ac)),
            SequenceVariant::NonCoding(n) => self.mapper.tx_to_g(n, Some(&g.ac)),
            _ => return Ok(false),
        };
        Ok(projected.is_ok_and(|p| p.to_string() == g.to_string()))
    }

    /// What a single nucleotide variant is compared by, read on the genome:
    /// the change it makes, the bases it states when it makes the reference,
    /// or, for a copy number, its projection. `None` for a variant with none
    /// of these (a conversion, an insertion of bases known only by number).
    fn comparand(&self, var: &SequenceVariant) -> Result<Option<Comparand>, HgvsError> {
        let g = match self.mapper.as_genomic(var) {
            None | Some(Err(HgvsError::UnsupportedOperation(_))) => return Ok(None),
            Some(g) => g?,
        };
        if matches!(g.posedit.edit, NaEdit::NACopy { .. }) {
            let mut posedit = g.posedit;
            // A parenthesised single base, `(5)`, is exact: written bare, it
            // compares as the base it names.
            if let Some(pos) = posedit.pos.as_mut().filter(|p| p.is_exact()) {
                pos.start.uncertain = false;
                if let Some(end) = pos.end.as_mut() {
                    end.uncertain = false;
                }
            }
            return Ok(Some(Comparand::Projection(g.ac, posedit.to_string())));
        }
        Ok(self
            .nucleotide_allele(&SequenceVariant::Genomic(g))?
            .map(|a| {
                if a.is_reference() {
                    Comparand::Statement(StatedBases::of_range(&a.accession, a.start, a.end))
                } else {
                    Comparand::Change(a)
                }
            }))
    }

    /// The canonical allele of a variant, or `None` for one that has none (a
    /// conversion, a copy number, an allele in cis).
    fn nucleotide_allele(
        &self,
        var: &SequenceVariant,
    ) -> Result<Option<CanonicalAllele>, HgvsError> {
        match self.mapper.canonical_allele(var) {
            Ok(a) => Ok(Some(a)),
            Err(HgvsError::UnsupportedOperation(_)) => Ok(None),
            Err(e) => Err(e),
        }
    }

    /// How the protein description `p1`, which leaves the protein `left`, and
    /// `p2` compare: `Identity` as one spelling, `Analogous` on two versions of
    /// one accession as one spelling or when they can leave the same protein,
    /// read on `p1`'s accession. A `p2` that does not read on that version,
    /// having been checked on its own, names a change that version does not
    /// have, and is `Different`. `left` is asked for only when the spellings
    /// differ, since it needs the protein sequence.
    fn protein_level(
        &self,
        p1: &PVariant,
        left: impl FnOnce() -> Result<Option<Outcome>, HgvsError>,
        p2: &PVariant,
    ) -> Result<EquivalenceLevel, HgvsError> {
        if protein_spelling(p1) == protein_spelling(p2) {
            return Ok(EquivalenceLevel::Identity);
        }
        let same = base_accession(&p1.ac) == base_accession(&p2.ac)
            && (posedit_spelling(&p1.posedit) == posedit_spelling(&p2.posedit)
                || match read_across_versions(self.protein_outcome_on(p2, &p1.ac), p1.ac != p2.ac)?
                {
                    Some(outcome) => same_protein(left()?, outcome),
                    None => false,
                });
        Ok(if same {
            EquivalenceLevel::Analogous
        } else {
            EquivalenceLevel::Different
        })
    }

    /// A nucleotide variant against a protein description: through every
    /// coding transcript the variant lies on that encodes the description's
    /// protein (see `encodes`), its prediction compared with the description
    /// as two protein descriptions are, the prediction leaving the protein the
    /// edited transcript encodes. A transcript of another protein takes no
    /// part. The level through the best transcript; a transcript on which the
    /// variant does not read, whose protein cannot be told, or whose
    /// prediction fails, is the error only when no transcript makes the two
    /// equivalent. A nucleotide allele predicts
    /// one description or a protein allele, `p.[(Ala4=;Tyr5Asp)]`; a protein
    /// allele on either side is compared as two protein alleles are.
    fn nucleotide_vs_protein(
        &self,
        nucleotide: &SequenceVariant,
        described: &SequenceVariant,
    ) -> Result<EquivalenceLevel, HgvsError> {
        // An uncertain position names no codons, so it predicts no protein.
        if !nucleotide.has_exact_position() {
            return Ok(EquivalenceLevel::Different);
        }
        let ac = described.ac();
        let coding = self.coding_forms(nucleotide)?;
        let mut best = EquivalenceLevel::Different;
        let mut failed = None;
        for c in &coding {
            let transcript = match c {
                Coding::One(one) => &one.ac,
                Coding::Cis(cis) => &cis.ac,
            };
            // The variant is checked on its transcript before the transcript
            // is set aside; a transcript of another protein says nothing of
            // this one.
            let on_transcript = match c {
                Coding::One(one) => SequenceVariant::Coding(one.clone()),
                Coding::Cis(cis) => SequenceVariant::CisPhased(cis.clone()),
            };
            match self
                .mapper
                .check_stated(&on_transcript)
                .and_then(|()| self.encodes(transcript, ac))
            {
                Ok(true) => {}
                Ok(false) => continue,
                Err(e) => {
                    failed.get_or_insert(e);
                    continue;
                }
            }
            let predicted = match c {
                Coding::One(one) => self
                    .mapper
                    .c_to_p(one, Some(ac))
                    .map(SequenceVariant::Protein),
                Coding::Cis(cis) => self.mapper.c_to_p_cis(cis, Some(ac)),
            };
            let level = match (predicted, described) {
                (Ok(SequenceVariant::Protein(predicted)), SequenceVariant::Protein(vp)) => {
                    self.protein_level(&predicted, || self.left_by(c, &predicted), vp)
                }
                (Ok(predicted), _) => self.equivalent_level_single(&predicted, described),
                // An identity of no whole codon says nothing of the protein,
                // so it matches no protein description; c_to_p refuses c.=
                // named against another protein, of which it says nothing.
                (Err(HgvsError::NoProteinStatement(_)), _) => continue,
                (Err(HgvsError::ValidationError(_)), _) if matches!(c, Coding::One(one) if matches!(one.posedit.edit, NaEdit::WholeIdentity)) => {
                    continue
                }
                (Err(e), _) => Err(e),
            };
            match level {
                Ok(EquivalenceLevel::Identity) => return Ok(EquivalenceLevel::Identity),
                Ok(EquivalenceLevel::Analogous) => best = EquivalenceLevel::Analogous,
                Ok(_) => {}
                Err(e) => {
                    failed.get_or_insert(e);
                }
            }
        }
        if best.is_equivalent() {
            return Ok(best);
        }
        // A description that cannot be read on its protein is the error,
        // whatever the transcripts predict.
        match described {
            SequenceVariant::Protein(vp) => {
                self.protein_outcome_on(vp, &vp.ac)?;
            }
            SequenceVariant::CisPhased(allele) => {
                for member in &allele.members {
                    if let SequenceVariant::Protein(vp) = member {
                        self.protein_outcome_on(vp, &vp.ac)?;
                    }
                }
            }
            _ => {}
        }
        failed.map_or(Ok(best), Err)
    }

    /// Whether `transcript` encodes `protein`: the provider pairs the two, in
    /// any version of the protein; or the transcript's CDS translates to the
    /// protein's sequence and the provider names the protein among those of
    /// the transcript's gene, as for another publisher's transcript, another
    /// isoform, or one paired with nothing. A paralogue that encodes the same
    /// sequence is of another gene, whose proteins the protein is not among.
    ///
    /// # Errors
    ///
    /// What the provider raises for a transcript or a protein sequence it
    /// does not hold; [`HgvsError::ValidationError`] for a transcript with no
    /// CDS to translate, and for one whose translation is the protein's when
    /// the transcript names no gene or the provider names no protein of its
    /// gene, so that whether it is the protein's or a paralogue's is unknown.
    fn encodes(&self, transcript: &str, protein: &str) -> Result<bool, HgvsError> {
        let provider = self.mapper.provider();
        let names = |symbol: &str, source: IdentifierKind| -> Result<Vec<String>, HgvsError> {
            Ok(provider
                .get_symbol_accessions(symbol, source, IdentifierKind::Protein)?
                .into_iter()
                .filter(|(kind, _)| *kind == crate::data::IdentifierType::ProteinAccession)
                .map(|(_, p)| p)
                .collect())
        };
        let named = |proteins: &[String]| {
            proteins
                .iter()
                .any(|p| base_accession(p) == base_accession(protein))
        };
        if named(&names(transcript, IdentifierKind::Transcript)?) {
            return Ok(true);
        }
        if self.mapper.encoded_protein(transcript)? != self.mapper.reference_protein(protein)? {
            return Ok(false);
        }
        let gene = provider.get_transcript(transcript, None)?.gene;
        let unknown = |why: String| {
            Err(HgvsError::ValidationError(format!(
                "{transcript} translates to {protein}, but {why}, so whether it encodes \
                 {protein} or a paralogue's copy is unknown"
            )))
        };
        if gene.is_empty() {
            return unknown(format!("{transcript} names no gene"));
        }
        let of_gene = names(&gene, IdentifierKind::Gene)?;
        if of_gene.is_empty() {
            return unknown(format!("the provider names no protein of its gene {gene}"));
        }
        Ok(named(&of_gene))
    }

    /// The coding forms of a nucleotide variant, a single one or an allele in
    /// cis: its own for a `c.` variant, its projection onto every coding
    /// transcript it lies on for a `g.` or `n.` one, its `c.` spelling for an
    /// `r.` one on a coding transcript.
    fn coding_forms(&self, var: &SequenceVariant) -> Result<Vec<Coding>, HgvsError> {
        let reference_of = |ac: &str| -> Result<String, HgvsError> {
            Ok(self
                .mapper
                .transcript((ac, crate::mapper::ReadAs::Description("n")), None)?
                .reference_accession)
        };
        Ok(match var {
            SequenceVariant::Coding(c) => vec![Coding::One(c.clone())],
            SequenceVariant::Genomic(g) => self.coding_of_genomic(g)?,
            SequenceVariant::NonCoding(n) => {
                let g = self.mapper.tx_to_g(n, Some(&reference_of(&n.ac)?))?;
                self.coding_of_genomic(&g)?
            }
            SequenceVariant::CisPhased(cis) => match cis.coordinate_type() {
                "c" => vec![Coding::Cis(cis.clone())],
                "g" => self.coding_of_genomic_cis(cis)?,
                "n" => {
                    let g = self.mapper.n_to_g_cis(cis, Some(&reference_of(&cis.ac)?))?;
                    self.coding_of_genomic_cis(&g)?
                }
                "r" => {
                    let tx = self
                        .mapper
                        .transcript((&cis.ac, crate::mapper::ReadAs::Description("r")), None)?;
                    if tx.cds_start_index.is_some() {
                        vec![Coding::Cis(self.mapper.r_to_c_cis(cis)?)]
                    } else {
                        let n = self.mapper.r_to_n_cis(cis)?;
                        let g = self.mapper.n_to_g_cis(&n, Some(&tx.reference_accession))?;
                        self.coding_of_genomic_cis(&g)?
                    }
                }
                _ => vec![],
            },
            _ => vec![],
        })
    }

    fn coding_of_genomic(&self, g: &crate::structs::GVariant) -> Result<Vec<Coding>, HgvsError> {
        Ok(self
            .mapper
            .g_to_c_all(g, self.searcher)?
            .into_iter()
            .map(Coding::One)
            .collect())
    }

    fn coding_of_genomic_cis(&self, g: &CisPhasedVariant) -> Result<Vec<Coding>, HgvsError> {
        Ok(self
            .mapper
            .g_to_c_all_cis(g, self.searcher)?
            .into_iter()
            .map(Coding::Cis)
            .collect())
    }

    /// The protein the coding variant `c` leaves, whose prediction is
    /// `predicted`: the residues the edited transcript encodes, exactly; for a
    /// prediction that leaves the protein as it is, the residues it states;
    /// and, where the prediction names no sequence (`p.(Met1?)`, `p.?`), what
    /// the prediction says.
    fn left_by(&self, c: &Coding, predicted: &PVariant) -> Result<Option<Outcome>, HgvsError> {
        if predicted.posedit.edit.is_identity() {
            let reference = self.mapper.reference_protein(&predicted.ac)?;
            return Ok(Some(silent_prediction(predicted, reference.len())));
        }
        let left = match c {
            Coding::One(one) => self.mapper.predicted_protein(one, Some(&predicted.ac))?,
            Coding::Cis(cis) => self
                .mapper
                .predicted_protein_cis(cis, Some(&predicted.ac))?,
        };
        match left {
            Some(left) => {
                let reference = self.mapper.reference_protein(&predicted.ac)?;
                Ok(Some(Outcome::Leaves {
                    changed_from: first_difference(&reference, &left.residues),
                    residues: left.residues,
                    open: !left.stops,
                }))
            }
            None => self.protein_outcome_on(predicted, &predicted.ac),
        }
    }

    /// The protein a p. description leaves: its residues in 1-letter code up
    /// to the stop, `X` where the description does not say (a frameshift's
    /// unnamed residues), and `open` when more residues follow whose number
    /// the description does not give (a stop loss written `Ter#Xxx`, an
    /// extension or frameshift of unknown length); where the change starts,
    /// for `p.Met1?`; or, for a description that leaves the protein as it
    /// is, the residues it says that of. `None` when it says nothing: `p.?`.
    fn protein_outcome_on(&self, vp: &PVariant, ac: &str) -> Result<Option<Outcome>, HgvsError> {
        match self.mapper.protein_read(vp, ac)? {
            ProteinRead::Said(outcome) => Ok(outcome),
            ProteinRead::Edit(edit) => Ok(Some(edit.outcome(&self.mapper.reference_protein(ac)?))),
        }
    }
}

/// The text two descriptions are one spelling by: as printed, without the
/// parentheses of a prediction, and with a protein's residues in one-letter
/// code. A prediction is taken as its confirmation, so `p.(Val600Glu)` is
/// `p.Val600Glu` and `r.(=)` is `r.=`; every other parenthesis is part of
/// what is said and stays, as in `g.(3_4)del`, `p.(Lys2_Leu3)del` and `r.(?)`.
/// A cis allele's members are unordered, so `g.[(3_4)del;10G>A]` is
/// `g.[10G>A;(3_4)del]`; a member written twice is still two members here.
fn spelling(var: &SequenceVariant) -> String {
    let mut var = var.clone();
    as_confirmed(&mut var);
    if let SequenceVariant::CisPhased(cis) = &mut var {
        cis.members.sort_by_cached_key(spelling);
    }
    let text = var.to_string();
    if var.coordinate_type() == "p" {
        one_letter(&text)
    } else {
        text
    }
}

/// A description's reading on `ac`, `None` where `ac` is another version
/// of its accession that contradicts what it states or lacks a residue it
/// names: checked on its own version, it names a change `ac` does not have.
fn read_across_versions<T>(
    read: Result<T, HgvsError>,
    other_version: bool,
) -> Result<Option<T>, HgvsError> {
    match read {
        Ok(read) => Ok(Some(read)),
        Err(HgvsError::ReferenceMismatch { .. } | HgvsError::PositionOutsideSequence { .. })
            if other_version =>
        {
            Ok(None)
        }
        Err(e) => Err(e),
    }
}

/// An r. allele of one member as that member.
fn rna_alone(var: std::borrow::Cow<'_, SequenceVariant>) -> std::borrow::Cow<'_, SequenceVariant> {
    match &*var {
        SequenceVariant::CisPhased(cis)
            if cis.coordinate_type() == "r" && cis.members.len() == 1 =>
        {
            std::borrow::Cow::Owned(cis.members[0].clone())
        }
        _ => var,
    }
}

/// Whether `e` refuses a description outright, alone or as a member of an
/// allele, rather than failing to find what it names.
fn is_a_refusal(e: &HgvsError) -> bool {
    match e {
        HgvsError::UnsupportedOperation(_) => true,
        HgvsError::CisMember { source, .. } => is_a_refusal(source),
        _ => false,
    }
}

/// [`spelling`] of a protein description.
fn protein_spelling(vp: &PVariant) -> String {
    spelling(&SequenceVariant::Protein(vp.clone()))
}

/// [`spelling`] of a protein posedit alone, for descriptions on two versions
/// of one accession.
fn posedit_spelling(posedit: &PosEdit<AaInterval, AaEdit>) -> String {
    let mut posedit = posedit.clone();
    posedit.predicted = false;
    one_letter(&posedit.to_string())
}

/// Clears the prediction flag of `var`, or of a cis allele's members. Only
/// `p.` and `r.` carry one.
fn as_confirmed(var: &mut SequenceVariant) {
    match var {
        SequenceVariant::Protein(p) => p.posedit.predicted = false,
        SequenceVariant::Rna(r) => r.posedit.predicted = false,
        SequenceVariant::CisPhased(cis) => cis.members.iter_mut().for_each(as_confirmed),
        _ => {}
    }
}

/// Protein text with its three-letter residue codes in one letter. A `?` is
/// kept: `p.Met1?` is an unknown change from Met1, not `p.Met1Xaa`.
fn one_letter(text: &str) -> String {
    const CODES: [(&str, &str); 24] = [
        ("Ala", "A"),
        ("Arg", "R"),
        ("Asn", "N"),
        ("Asp", "D"),
        ("Cys", "C"),
        ("Gln", "Q"),
        ("Glu", "E"),
        ("Gly", "G"),
        ("His", "H"),
        ("Ile", "I"),
        ("Leu", "L"),
        ("Lys", "K"),
        ("Met", "M"),
        ("Phe", "F"),
        ("Pro", "P"),
        ("Ser", "S"),
        ("Thr", "T"),
        ("Trp", "W"),
        ("Tyr", "Y"),
        ("Val", "V"),
        ("Asx", "B"),
        ("Glx", "Z"),
        ("Xaa", "X"),
        ("Ter", "*"),
    ];
    CODES
        .iter()
        .fold(text.to_string(), |s, (from, to)| s.replace(from, to))
}

/// The value two single nucleotide variants are the same change by, on the
/// genome. A variant that changes the sequence states no base unchanged, and
/// one that makes the reference changes nothing, so neither is ever the other.
#[derive(Debug, PartialEq, Eq)]
enum Comparand {
    /// The change, fully justified: one value for every place it could be
    /// written.
    Change(CanonicalAllele),
    /// The bases a variant that makes the reference names, as written: no
    /// widening over a run, so `c.10=` and `c.20=` differ.
    Statement(StatedBases),
    /// A copy number, which has no bases to apply: its accession and edit as
    /// projected.
    Projection(String, String),
}

/// The value a member of a cis allele is compared by, where the allele is
/// compared member by member or the member is not applied to a sequence.
#[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
enum MemberKey {
    /// An exact member's canonical allele, as SPDI.
    Allele(String),
    /// A copy number's accession and edit, projected to the genome.
    Projection(String, String),
    /// A member compared by its text: an uncertain one by the [`spelling`] of
    /// its `g.`, `c.` or `n.` form; one that is applied to no sequence, a
    /// conversion or unstated inserted bases, by that of its genomic form.
    Spelling(String),
}

/// Whether a member of a nucleotide allele names a position. A statement,
/// the whole-sequence `=`, `?`, `spl` or `0`, names none.
fn names_a_position(m: &SequenceVariant) -> bool {
    match m {
        SequenceVariant::Genomic(v) => v.posedit.pos.is_some(),
        SequenceVariant::Mitochondrial(v) => v.posedit.pos.is_some(),
        SequenceVariant::Coding(v) => v.posedit.pos.is_some(),
        SequenceVariant::NonCoding(v) => v.posedit.pos.is_some(),
        SequenceVariant::Rna(v) => v.posedit.pos.is_some(),
        SequenceVariant::Protein(_) | SequenceVariant::CisPhased(_) => false,
    }
}

/// Whether a member of a nucleotide allele has an exact position. A
/// statement, which has no position, makes no sequence.
fn is_exact_member(m: &SequenceVariant) -> bool {
    m.has_exact_position() && names_a_position(m)
}

/// Whether every member of a nucleotide allele has an exact position, so
/// that the allele can be read into the one sequence its members make.
fn reads_as_a_sequence(cis: &CisPhasedVariant) -> bool {
    cis.coordinate_type() != "p" && cis.members.iter().all(is_exact_member)
}

/// The value an allele of nucleotide members is the same change by, on the
/// genome: the sequence its members make, the bases it states unchanged, and
/// the members no sequence is built from, as a set.
struct AlleleComparand {
    reading: CisReading,
    unapplied: Vec<MemberKey>,
}

impl AlleleComparand {
    fn same(&self, other: &AlleleComparand, eq: &VariantEquivalence) -> Result<bool, HgvsError> {
        let (a, b) = (&self.reading, &other.reading);
        if (a.accession(), a.kind()) != (b.accession(), b.kind())
            || a.stated() != b.stated()
            || self.unapplied != other.unapplied
        {
            return Ok(false);
        }
        let genome = eq.mapper.refs.reference(a.accession(), a.kind());
        a.same_sequence(b, &genome)
    }
}

/// The statement a silent prediction of `c_to_p` makes, `p.(=)` or
/// `p.(Leu3=)`, read from its positions on a protein `len` residues long
/// before its stop: the prediction was translated from the transcript, so
/// its residues are not read again on the protein.
fn silent_prediction(p: &PVariant, len: usize) -> Outcome {
    Outcome::Unchanged(match &p.posedit.pos {
        None => Stated::whole(),
        Some(pos) => {
            let index = |a: &crate::structs::AAPosition| a.base.to_index().0 as usize;
            let start = index(&pos.start);
            Stated::residues(start..pos.end.as_ref().map_or(start, index) + 1, len)
        }
    })
}

/// The accession without its version: `NP_000050` of `NP_000050.3`.
fn base_accession(ac: &str) -> &str {
    ac.split('.').next().unwrap_or(ac)
}
