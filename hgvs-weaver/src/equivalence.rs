use crate::allele::{CanonicalAllele, StatedBases};
use crate::data::{IdentifierKind, TranscriptSearch};
use crate::edits::{NaEdit, RnaStatement, StopDistance};
use crate::error::HgvsError;
use crate::mapper::VariantMapper;
use crate::structs::{
    CVariant, CisPhasedVariant, LinearVariant, PVariant, SequenceVariant, Variant,
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
        // The text comparison drops parentheses, and those of r.(?) are part
        // of what it states.
        if rna_statements(var1) != rna_statements(var2) {
            return Ok(EquivalenceLevel::Different);
        }
        if self.normalize_format(&var1.to_string()) == self.normalize_format(&var2.to_string()) {
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
                    && (self.normalize_format(&p1.posedit.to_string())
                        == self.normalize_format(&p2.posedit.to_string())
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
            let edit = self.normalize_format(&g.posedit.to_string());
            return Ok(Some(Comparand::Projection(g.ac, edit)));
        }
        Ok(self
            .nucleotide_allele(&SequenceVariant::Genomic(g))?
            .map(|a| {
                if a.is_reference() {
                    Comparand::Statement(StatedBases::new(&a.accession, [(a.start, a.end)]))
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
        let described = self.normalize_format(&vp.to_string());
        let mut analogous = false;
        for c in &coding {
            if let Ok(predicted) = self.mapper.c_to_p(c, Some(&vp.ac)) {
                if predicted.to_string() == vp.to_string() {
                    return Ok(EquivalenceLevel::Identity);
                }
                analogous |= self.normalize_format(&predicted.to_string()) == described;
            }
        }
        // Then by the protein left behind.
        if !analogous {
            let outcome = self.protein_outcome_on(vp, &vp.ac)?;
            let reference = self.reference_protein(&vp.ac)?;
            for c in &coding {
                // predicted_protein raises for c.= named against another protein;
                // on its own, p.(=) names no residues to compare.
                if matches!(c.posedit.edit, crate::edits::NaEdit::WholeIdentity) {
                    continue;
                }
                let predicted = self
                    .mapper
                    .predicted_protein(c, Some(&vp.ac))?
                    .map(|residues| Outcome {
                        changed_from: first_difference(&reference, &residues),
                        residues,
                        open: false,
                        anchored: false,
                    });
                if same_protein(predicted, outcome.clone()) {
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
    /// extension or frameshift of unknown length). `None` when it says
    /// nothing about the sequence: `p.?`, `p.Met1?`.
    fn protein_outcome_on(&self, vp: &PVariant, ac: &str) -> Result<Option<Outcome>, HgvsError> {
        use crate::edits::{AaEdit, ProteinStatement};
        let reference = self
            .mapper
            .refs
            .reference(ac, crate::data::IdentifierType::ProteinAccession);
        let whole = reference.whole()?;
        let protein = whole.trim_end_matches('*');
        let len = protein.len();
        let closed = |s: String, changed_from: usize| Outcome {
            residues: s,
            open: false,
            changed_from,
            anchored: false,
        };
        if vp.posedit.pos.is_none() {
            return Ok(match &vp.posedit.edit {
                AaEdit::Statement(ProteinStatement::Unchanged) => {
                    Some(closed(protein.to_string(), len))
                }
                AaEdit::Statement(ProteinStatement::NoProtein { .. }) => {
                    Some(closed(String::new(), 0))
                }
                _ => None,
            });
        }
        let (start, end) = crate::mapper::residue_range(vp, &reference)?;
        // `p.Met1?`: something happens from this residue on; what, it does not say.
        if let AaEdit::Statement(statement) = &vp.posedit.edit {
            return Ok((*statement == ProteinStatement::Unknown).then(|| Outcome {
                residues: protein[..start.min(len)].to_string(),
                open: true,
                changed_from: start.min(len),
                anchored: true,
            }));
        }
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
                let r = edit.resolve(&reference, start, end)?;
                // A change at the stop itself that does not stop again reads on
                // for an unknown distance: ClinVar's `Ter#Xxx` for a stop loss.
                let open = r.start >= len && !r.alt.contains('*');
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
        Ok(Some(Outcome {
            changed_from: first_difference(protein, &residues),
            residues,
            open,
            anchored: false,
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

    fn normalize_format(&self, s: &str) -> String {
        // A ? is kept: p.Met1? is an unknown change from Met1, not p.Met1Xaa.
        let mut s = s.replace(['(', ')'], "");
        // Normalize 3-letter AA codes to 1-letter
        let replacements = [
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
        for (from, to) in replacements {
            s = s.replace(from, to);
        }
        s
    }
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

/// What a protein description leaves behind; see `protein_outcome`.
#[derive(Clone, Debug, PartialEq, Eq)]
struct Outcome {
    residues: String,
    /// More residues follow, their number unknown.
    open: bool,
    /// The first residue that differs from the reference (its length if none).
    changed_from: usize,
    /// A statement (`p.Met1?`): only where the change starts is known, and
    /// the other side must start its change there.
    anchored: bool,
}

/// What `var` states about the whole RNA (`r.0`, `r.(?)`), itself or through
/// the members of a cis allele, in the order written.
fn rna_statements(var: &SequenceVariant) -> Vec<RnaStatement> {
    match var {
        SequenceVariant::Rna(r) => match r.posedit.edit {
            NaEdit::Statement(statement) => vec![statement],
            _ => vec![],
        },
        SequenceVariant::CisPhased(cis) => cis.members.iter().flat_map(rna_statements).collect(),
        _ => vec![],
    }
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
/// `None` (nothing said) matches nothing.
fn same_protein(a: Option<Outcome>, b: Option<Outcome>) -> bool {
    let (Some(a), Some(b)) = (a, b) else {
        return false;
    };
    // A statement says only where the change starts: the other side must
    // start its change at the same residue.
    if a.anchored || b.anchored {
        return a.changed_from == b.changed_from;
    }
    let (na, nb) = (a.residues.len(), b.residues.len());
    let fits = match (a.open, b.open) {
        (false, false) => na == nb,
        (true, false) => nb >= na,
        (false, true) => na >= nb,
        (true, true) => true,
    };
    fits && a
        .residues
        .chars()
        .zip(b.residues.chars())
        .all(|(x, y)| x == y || x == 'X' || y == 'X')
}
