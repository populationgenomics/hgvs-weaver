use crate::allele::{CanonicalAllele, StatedBases};
use crate::data::{IdentifierKind, TranscriptSearch};
use crate::edits::{AaEdit, NaEdit, StopDistance};
use crate::error::HgvsError;
use crate::mapper::VariantMapper;
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
        if spelling(var1) == spelling(var2) {
            return Ok(EquivalenceLevel::Identity);
        }
        if let Some(level) = self.cis_phased_level(var1, var2)? {
            return Ok(level);
        }
        let same = match (var1, var2) {
            (SequenceVariant::Protein(p1), SequenceVariant::Protein(p2)) => {
                // Versions of one protein accession are compared on ours: the
                // same description in another spelling, or the same protein left.
                base_accession(&p1.ac) == base_accession(&p2.ac)
                    && (posedit_spelling(&p1.posedit) == posedit_spelling(&p2.posedit)
                        || same_protein(
                            self.protein_outcome_on(p1, &p1.ac)?,
                            self.protein_outcome_on(p2, &p1.ac)?,
                        ))
            }
            (SequenceVariant::Protein(p), nucleotide)
            | (nucleotide, SequenceVariant::Protein(p)) => {
                return self.nucleotide_vs_protein(nucleotide, p);
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
    /// is. Two cis alleles are `Analogous` when the sets of their members'
    /// canonical alleles are equal, whatever the order written, and
    /// `Different` otherwise. A cis allele of one member is that member; one
    /// of several is `Different` from any plain variant. Nothing is projected
    /// to protein: what several changes on one molecule do to the protein is
    /// not the sum of what each does.
    fn cis_phased_level(
        &self,
        var1: &SequenceVariant,
        var2: &SequenceVariant,
    ) -> Result<Option<EquivalenceLevel>, HgvsError> {
        Ok(Some(match (var1, var2) {
            (SequenceVariant::CisPhased(a), SequenceVariant::CisPhased(b)) => {
                match (self.member_alleles(a)?, self.member_alleles(b)?) {
                    (Some(x), Some(y)) if x == y => EquivalenceLevel::Analogous,
                    _ => EquivalenceLevel::Different,
                }
            }
            (SequenceVariant::CisPhased(cis), other) | (other, SequenceVariant::CisPhased(cis)) => {
                match cis.members.as_slice() {
                    [member] => return self.equivalent_level_single(member, other).map(Some),
                    _ => EquivalenceLevel::Different,
                }
            }
            _ => return Ok(None),
        }))
    }

    /// The canonical alleles of a cis allele's members, as a sorted set of
    /// SPDI strings; `None` when a member has no canonical allele.
    fn member_alleles(&self, cis: &CisPhasedVariant) -> Result<Option<Vec<String>>, HgvsError> {
        let mut alleles = Vec::with_capacity(cis.members.len());
        for m in &cis.members {
            match self.nucleotide_allele(m)? {
                Some(a) => alleles.push(a.spdi()),
                None => return Ok(None),
            }
        }
        alleles.sort_unstable();
        alleles.dedup();
        Ok(Some(alleles))
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

    /// A nucleotide variant against a protein description: through every
    /// transcript the variant lies on, `Identity` if the predicted description
    /// is the same text, `Analogous` if the protein left is the same.
    fn nucleotide_vs_protein(
        &self,
        nucleotide: &SequenceVariant,
        vp: &PVariant,
    ) -> Result<EquivalenceLevel, HgvsError> {
        // An uncertain position names no codons, so it predicts no protein.
        if !nucleotide.has_exact_position() {
            return Ok(EquivalenceLevel::Different);
        }
        let coding: Vec<CVariant> = match nucleotide {
            SequenceVariant::Coding(c) => vec![c.clone()],
            SequenceVariant::Genomic(g) => self.mapper.g_to_c_all(g, self.searcher)?,
            SequenceVariant::NonCoding(n) => {
                let tx = self.mapper.provider().get_transcript(&n.ac, None)?;
                let g = self
                    .mapper
                    .tx_to_g(n, Some(tx.reference_accession.as_str()))?;
                self.mapper.g_to_c_all(&g, self.searcher)?
            }
            _ => vec![],
        };
        // First by description: the prediction written exactly as given is
        // Identity (c. implies p.(...) exactly); the same description in
        // another spelling is Analogous. Neither needs the protein sequence.
        let described = protein_spelling(vp);
        let mut analogous = false;
        let mut predictions = Vec::with_capacity(coding.len());
        for c in &coding {
            let predicted = self.mapper.c_to_p(c, Some(&vp.ac));
            if let Ok(predicted) = &predicted {
                if predicted.to_string() == vp.to_string() {
                    return Ok(EquivalenceLevel::Identity);
                }
                analogous |= protein_spelling(predicted) == described;
            }
            predictions.push(predicted);
        }
        // Then by the protein left behind, or, for a prediction that leaves
        // the protein unchanged, by the residues it states.
        if !analogous {
            let outcome = self.protein_outcome_on(vp, &vp.ac)?;
            let reference = self.reference_protein(&vp.ac)?;
            for (c, predicted) in coding.iter().zip(&predictions) {
                let candidate = match predicted {
                    Ok(p) if p.posedit.edit.is_identity() => {
                        Some(silent_prediction(p, reference.len()))
                    }
                    // An identity of no whole codon says nothing of the
                    // protein, so it matches no protein description.
                    Err(HgvsError::NoProteinStatement(_)) => continue,
                    // c_to_p refuses c.= named against another protein, of
                    // which it says nothing.
                    _ if matches!(c.posedit.edit, NaEdit::WholeIdentity) => continue,
                    _ => self
                        .mapper
                        .predicted_protein(c, Some(&vp.ac))?
                        .map(|residues| Outcome::Leaves {
                            changed_from: first_difference(&reference, &residues),
                            residues,
                            open: false,
                        }),
                };
                if same_protein(candidate, outcome.clone()) {
                    analogous = true;
                }
            }
        }
        Ok(if analogous {
            EquivalenceLevel::Analogous
        } else {
            EquivalenceLevel::Different
        })
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
fn spelling(var: &SequenceVariant) -> String {
    let mut var = var.clone();
    as_confirmed(&mut var);
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
