use crate::allele::{CanonicalAllele, StatedBases};
use crate::data::{IdentifierKind, TranscriptSearch};
use crate::edits::{AaEdit, NaEdit, StopDistance};
use crate::error::HgvsError;
use crate::mapper::VariantMapper;
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

pub struct VariantEquivalence<'a> {
    /// The mapper whose provider, cache and refget lookup the comparison uses.
    pub mapper: &'a VariantMapper<'a>,
    pub searcher: &'a dyn TranscriptSearch,
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
        // An r. variant is its c. or n. spelling in RNA letters; compare it as that.
        let var1 = self.rna_as_transcript(var1)?;
        let var2 = self.rna_as_transcript(var2)?;
        // Expand gene symbols if present
        let vars1 = self.expand_if_gene_symbol(&var1)?;
        let vars2 = self.expand_if_gene_symbol(&var2)?;

        for v1 in &vars1 {
            for v2 in &vars2 {
                let lvl = self.equivalent_level_single(v1, v2)?;
                if lvl.is_equivalent() {
                    return Ok(lvl);
                }
            }
        }
        Ok(EquivalenceLevel::Different)
    }

    fn rna_as_transcript<'v>(
        &self,
        var: std::borrow::Cow<'v, SequenceVariant>,
    ) -> Result<std::borrow::Cow<'v, SequenceVariant>, HgvsError> {
        match &*var {
            SequenceVariant::Rna(r) if !r.posedit.edit.is_whole_transcript() => {
                Ok(std::borrow::Cow::Owned(self.mapper.r_as_transcript(r)?))
            }
            _ => Ok(var),
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
        // An allele whose exact members are contradictory or ambiguous is
        // refused before any spelling can match it; one whose readings are not
        // listed still matches its own spelling.
        for var in [var1, var2] {
            if let SequenceVariant::CisPhased(cis) = var {
                if let Some(exact) = exact_members(cis)? {
                    if let Err(e @ (HgvsError::ValidationError(_) | HgvsError::Ambiguous { .. })) =
                        self.mapper.read_cis(&exact)
                    {
                        return Err(e);
                    }
                }
            }
        }
        if spelling(var1) == spelling(var2) {
            return Ok(EquivalenceLevel::Identity);
        }
        if let Some(level) = self.cis_phased_level(var1, var2)? {
            return Ok(level);
        }
        let same = match (var1, var2) {
            (SequenceVariant::Protein(p1), SequenceVariant::Protein(p2)) => {
                return self.protein_level(p1, || self.protein_outcome_on(p1, &p1.ac), p2);
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
                let m = self.rna_as_transcript(mito_as_genomic(m))?;
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
    /// read on `p1`'s accession. `left` is asked for only when the spellings
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
        // Versions of one protein accession are compared on ours.
        let same = base_accession(&p1.ac) == base_accession(&p2.ac)
            && (posedit_spelling(&p1.posedit) == posedit_spelling(&p2.posedit)
                || same_protein(left()?, self.protein_outcome_on(p2, &p1.ac)?));
        Ok(if same {
            EquivalenceLevel::Analogous
        } else {
            EquivalenceLevel::Different
        })
    }

    /// A nucleotide variant against a protein description: through every
    /// coding transcript the variant lies on, its prediction compared with the
    /// description as two protein descriptions are, the prediction leaving
    /// the protein the edited transcript encodes. The level through the best
    /// transcript; a transcript whose prediction fails is the error only when
    /// no transcript makes the two equivalent. A nucleotide allele predicts
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
        if let SequenceVariant::Protein(vp) = described {
            self.protein_outcome_on(vp, &vp.ac)?;
        }
        failed.map_or(Ok(best), Err)
    }

    /// The coding forms of a nucleotide variant, a single one or an allele in
    /// cis: its own for a `c.` variant, its projection onto every coding
    /// transcript it lies on for a `g.` or `n.` one, its `c.` spelling for an
    /// `r.` one on a coding transcript.
    fn coding_forms(&self, var: &SequenceVariant) -> Result<Vec<Coding>, HgvsError> {
        let reference_of = |ac: &str| -> Result<String, HgvsError> {
            Ok(self
                .mapper
                .provider()
                .get_transcript(ac, None)?
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
                    let tx = self.mapper.provider().get_transcript(&cis.ac, None)?;
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
            let reference = self.reference_protein(&predicted.ac)?;
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
                let reference = self.reference_protein(&predicted.ac)?;
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
        use crate::edits::{AaEdit, ProteinStatement};
        let reference = self
            .mapper
            .refs
            .reference(ac, crate::data::IdentifierType::ProteinAccession);
        let whole = reference.whole()?;
        let protein = whole.trim_end_matches('*');
        let len = protein.len();
        if vp.posedit.pos.is_none() {
            return Ok(match &vp.posedit.edit {
                AaEdit::Statement(ProteinStatement::Unchanged) => {
                    Some(Outcome::Unchanged(Stated::whole()))
                }
                AaEdit::Statement(ProteinStatement::NoProtein { .. }) => Some(Outcome::Leaves {
                    residues: String::new(),
                    open: false,
                    changed_from: 0,
                }),
                _ => None,
            });
        }
        // An uncertain range, `p.(Lys2_Leu3)del`, does not say which residues change.
        if vp.posedit.pos.as_ref().is_some_and(|p| p.uncertain) {
            return Ok(None);
        }
        let (start, end) = crate::mapper::residue_range(vp, &reference)?;
        // `p.Met1?`: something happens from this residue on; what, it does not say.
        if let AaEdit::Statement(statement) = &vp.posedit.edit {
            return Ok((*statement == ProteinStatement::Unknown)
                .then(|| Outcome::ChangesFrom(start.min(len))));
        }
        // The residues an edit that resolves to a sequence names: those its
        // position names (an insertion's two flanks) and those it resolves
        // over (a repeat's whole run).
        let mut written_over = None;
        let (start, end, alt, open) = match &vp.posedit.edit {
            // `Xxx#Yyyfs*N`: Yyy, then N - 2 residues it does not name, then a stop.
            AaEdit::Fs { alt, stop } => {
                let mut alt = crate::utils::residues_1(alt)?;
                let count = match stop {
                    Some(StopDistance::Residues { count }) => Some(count.get() as usize),
                    _ => None,
                };
                match count {
                    Some(n) => {
                        while alt.len() < n - 1 {
                            alt.push('X');
                        }
                        alt.push('*');
                        (start, end, alt, false)
                    }
                    _ if alt == "*" => (start, end, alt, false),
                    _ => (start, end, alt, true),
                }
            }
            // `Ter#Xxxext*N`: the stop becomes Xxx, then N - 1 residues, then a stop.
            AaEdit::Ext { alt, end } => {
                let mut alt = crate::utils::residues_1(alt)?;
                match end.residues_to_stop().map(|n| n.get() as usize) {
                    Some(n) => {
                        while alt.len() < n {
                            alt.push('X');
                        }
                        (len, len, alt, false)
                    }
                    _ => (len, len, alt, true),
                }
            }
            AaEdit::Statement(_) | AaEdit::None => return Ok(None),
            edit => {
                // Read with its stop, so that `Ter7=` one past a protein
                // served without it resolves to the stop it names.
                let with_stop = format!("{protein}*");
                let r = edit.resolve_with(reference.accession(), start, end, |s, e| {
                    let s = s.min(with_stop.len());
                    Ok(with_stop[s..e.clamp(s, with_stop.len())].to_string())
                })?;
                // A change that replaces the stop and does not stop again reads
                // on for an unknown distance: ClinVar's `Ter#Xxx` for a stop
                // loss. An insertion before the stop keeps it.
                let open = r.end > len && !r.alt.contains('*');
                written_over = Some(start.min(r.start)..end.max(r.end));
                (r.start, r.end, r.alt, open)
            }
        };
        // Nothing of the reference follows a stop, nor an open end.
        let (end, alt) = match alt.find('*') {
            Some(k) => (len, alt[..k].to_string()),
            None if open => (len, alt),
            None => (end, alt),
        };
        let start = start.min(len);
        let residues = format!(
            "{}{}{}",
            &protein[..start],
            alt,
            &protein[end.clamp(start, len)..]
        );
        // An edit that leaves the protein as it is states the residues it
        // names: `p.Leu3=`, `p.Leu3Leu`, `p.Ter7Ter`, `p.Arg6_Ter7insTer`.
        if let Some(span) = written_over.filter(|_| !open && residues == protein) {
            return Ok(Some(Outcome::Unchanged(Stated::residues(span, len))));
        }
        Ok(Some(Outcome::Leaves {
            changed_from: first_difference(protein, &residues),
            residues,
            open,
        }))
    }

    /// The protein sequence of `ac`, without a trailing stop.
    fn reference_protein(&self, ac: &str) -> Result<String, HgvsError> {
        let whole = self
            .mapper
            .refs
            .reference(ac, crate::data::IdentifierType::ProteinAccession)
            .whole()?;
        Ok(whole.trim_end_matches('*').to_string())
    }

    fn expand_if_gene_symbol(
        &self,
        var: &SequenceVariant,
    ) -> Result<Vec<SequenceVariant>, HgvsError> {
        let ac = var.ac();

        // Use DataProvider to determine if this is a symbol or an accession.
        let id_type = self.mapper.provider().get_identifier_type(ac)?;

        if id_type == crate::data::IdentifierType::GeneSymbol {
            let target_kind = match var {
                SequenceVariant::Protein(_) => IdentifierKind::Protein,
                SequenceVariant::Coding(_)
                | SequenceVariant::NonCoding(_)
                | SequenceVariant::Rna(_) => IdentifierKind::Transcript,
                _ => IdentifierKind::Genomic,
            };

            // Try symbol expansion.
            let accessions = self.mapper.provider().get_symbol_accessions(
                ac,
                IdentifierKind::Genomic,
                target_kind,
            )?;

            if !accessions.is_empty() {
                let mut expanded = Vec::new();
                for (ac_type, new_ac) in accessions {
                    // Only expand into compatible types.
                    let is_compatible = match (var, ac_type) {
                        (
                            SequenceVariant::Protein(_),
                            crate::data::IdentifierType::ProteinAccession,
                        ) => true,
                        (
                            SequenceVariant::Coding(_)
                            | SequenceVariant::NonCoding(_)
                            | SequenceVariant::Rna(_),
                            crate::data::IdentifierType::TranscriptAccession,
                        ) => true,
                        (
                            SequenceVariant::Genomic(_) | SequenceVariant::Mitochondrial(_),
                            crate::data::IdentifierType::GenomicAccession,
                        ) => true,
                        // Allow g. on transcripts if specifically provided
                        (
                            SequenceVariant::Genomic(_),
                            crate::data::IdentifierType::TranscriptAccession,
                        ) => true,
                        _ => false,
                    };

                    if is_compatible {
                        let mut v = var.clone();
                        v.set_ac(new_ac);
                        expanded.push(v);
                    }
                }
                if !expanded.is_empty() {
                    return Ok(expanded);
                }
            }
        }

        // Default: return as-is.
        Ok(vec![var.clone()])
    }
}

/// The text two descriptions are one spelling by: as printed, without the
/// parentheses of a prediction, and with a protein's residues in one-letter
/// code. A prediction is taken as its confirmation, so `p.(Val600Glu)` is
/// `p.Val600Glu` and `r.(=)` is `r.=`; every other parenthesis is part of
/// what is said and stays, as in `g.(3_4)del`, `p.(Lys2_Leu3)del` and `r.(?)`.
/// A cis allele's members are unordered, so `g.[(3_4)del;10G>A]` is
/// `g.[10G>A;(3_4)del]`; a member written twice is still two members here.
/// A coding variant, or an allele in cis of coding members, whose prediction
/// is compared with a protein description.
enum Coding {
    One(CVariant),
    Cis(CisPhasedVariant),
}

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

/// Whether a member of a nucleotide allele has an exact position. A
/// statement, which has no position, makes no sequence.
fn is_exact_member(m: &SequenceVariant) -> bool {
    m.has_exact_position()
        && match m {
            SequenceVariant::Genomic(v) => v.posedit.pos.is_some(),
            SequenceVariant::Mitochondrial(v) => v.posedit.pos.is_some(),
            SequenceVariant::Coding(v) => v.posedit.pos.is_some(),
            SequenceVariant::NonCoding(v) => v.posedit.pos.is_some(),
            SequenceVariant::Rna(v) => v.posedit.pos.is_some(),
            SequenceVariant::Protein(_) | SequenceVariant::CisPhased(_) => false,
        }
}

/// Whether every member of a nucleotide allele has an exact position, so
/// that the allele can be read into the one sequence its members make.
fn reads_as_a_sequence(cis: &CisPhasedVariant) -> bool {
    cis.coordinate_type() != "p" && cis.members.iter().all(is_exact_member)
}

/// The allele of a nucleotide allele's members at exact positions; `None`
/// for a protein allele or one with no such member.
fn exact_members(cis: &CisPhasedVariant) -> Result<Option<CisPhasedVariant>, HgvsError> {
    if cis.coordinate_type() == "p" {
        return Ok(None);
    }
    let members: Vec<SequenceVariant> = cis
        .members
        .iter()
        .filter(|m| is_exact_member(m))
        .cloned()
        .collect();
    if members.is_empty() {
        return Ok(None);
    }
    CisPhasedVariant::new(cis.ac.clone(), cis.gene.clone(), members).map(Some)
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

/// What a protein description says of the protein; see `protein_outcome_on`.
#[derive(Clone, Debug, PartialEq, Eq)]
enum Outcome {
    /// The protein left, which differs from the reference.
    Leaves {
        residues: String,
        /// More residues follow, their number unknown.
        open: bool,
        /// The first residue that differs from the reference (its length if none).
        changed_from: usize,
    },
    /// `p.Met1?`: a change from this residue on, of which nothing more is
    /// said; the other side must start its change there.
    ChangesFrom(usize),
    /// The protein as it is, stated of these residues.
    Unchanged(Stated),
}

use stated::Stated;

mod stated {
    /// The residues a description that leaves the protein unchanged speaks
    /// of: the whole protein, `p.(=)`, or some of its residues, `p.(Leu3=)`.
    /// Every residue through the stop is the whole protein, so the two
    /// spellings of it are one value.
    #[derive(Clone, Debug, PartialEq, Eq)]
    pub(super) struct Stated {
        /// 0-based and half-open; `None` for the whole protein.
        residues: Option<std::ops::Range<usize>>,
    }

    impl Stated {
        pub(super) fn whole() -> Self {
            Stated { residues: None }
        }

        /// `residues` of a protein `len` residues long before its stop.
        pub(super) fn residues(residues: std::ops::Range<usize>, len: usize) -> Self {
            if residues.start == 0 && residues.end > len {
                return Self::whole();
            }
            Stated {
                residues: Some(residues),
            }
        }
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

/// The first index at which two proteins differ, or the shorter length.
fn first_difference(a: &str, b: &str) -> usize {
    a.bytes()
        .zip(b.bytes())
        .position(|(x, y)| x != y)
        .unwrap_or(a.len().min(b.len()))
}

/// Whether two outcomes can be the same protein: residue for residue with
/// `X` standing for any, over the whole of both when both are closed, over
/// the shorter when one is open (the other must then be at least as long).
/// A change from a residue matches a change starting there; a protein left
/// unchanged matches only the same statement about the same residues.
/// `None` (nothing said) matches nothing.
fn same_protein(a: Option<Outcome>, b: Option<Outcome>) -> bool {
    match (a, b) {
        (Some(Outcome::Unchanged(a)), Some(Outcome::Unchanged(b))) => a == b,
        (Some(Outcome::ChangesFrom(a)), Some(Outcome::ChangesFrom(b))) => a == b,
        (Some(Outcome::ChangesFrom(at)), Some(Outcome::Leaves { changed_from, .. }))
        | (Some(Outcome::Leaves { changed_from, .. }), Some(Outcome::ChangesFrom(at))) => {
            at == changed_from
        }
        (
            Some(Outcome::Leaves {
                residues: a,
                open: a_open,
                ..
            }),
            Some(Outcome::Leaves {
                residues: b,
                open: b_open,
                ..
            }),
        ) => {
            let (na, nb) = (a.len(), b.len());
            let fits = match (a_open, b_open) {
                (false, false) => na == nb,
                (true, false) => nb >= na,
                (false, true) => na >= nb,
                (true, true) => true,
            };
            fits && a
                .chars()
                .zip(b.chars())
                .all(|(x, y)| x == y || x == 'X' || y == 'X')
        }
        _ => false,
    }
}
