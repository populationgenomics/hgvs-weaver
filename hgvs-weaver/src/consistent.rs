//! A description on a name that is not one sequence, a gene symbol or an
//! unversioned accession, as the variants it can be on the versioned
//! accessions a store holds under that name; choices.md, "A name that is not
//! one sequence is read on each the store holds".

use crate::coords::SequenceVariant;
use crate::data::{AccessionCatalogue, AccessionInfo, IdentifierKind, IdentifierType};
use crate::error::HgvsError;
use crate::mapper::{ReadAs, VariantMapper};
use crate::structs::Variant;
use crate::structs::{Anchor, BaseOffsetPosition};

/// A description on an accession the store holds as one sequence,
/// consistent with that sequence.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct ConsistentVariant {
    /// The description on one accession.
    pub variant: SequenceVariant,
    /// Whether the description states bases or residues the sequence was
    /// found to hold. Without any, it is consistent with every sequence its
    /// positions lie on, which is no evidence for this one.
    pub states_bases: bool,
    /// What the catalogue records about the accession; `None` without a
    /// catalogue, or where it records nothing.
    pub info: Option<AccessionInfo>,
}

impl<'a> VariantMapper<'a> {
    /// Every variant `var` can be on an accession the store holds as one
    /// sequence, that is consistent with that sequence as
    /// [`validate`](Self::validate) finds it. A gene symbol gives the
    /// description on each of the gene's accessions that can carry it:
    /// coding transcripts for `c.`, any transcript for `n.`, and for `r.`
    /// unless it counts a position from the CDS, proteins for `p.`, genomic
    /// sequences for `g.`. An unversioned accession gives it on each
    /// version `catalogue` lists, newest first, or on the accession itself
    /// where the catalogue lists that alone, for a record with no version of
    /// its own. A versioned accession gives the description itself. Empty
    /// when the store holds sequences under the name and none is consistent
    /// with the description. Each variant it returns, read again with the
    /// same catalogue, gives itself alone.
    ///
    /// # Errors
    ///
    /// [`HgvsError::DataProviderError`] when the store holds nothing under
    /// the name that can carry the description, when the name is an
    /// unversioned accession and there is no `catalogue`, or when the
    /// catalogue lists something other than the name's versions; any error
    /// the provider or the catalogue raises, and any `validate` raises
    /// rather than answering.
    pub fn consistent_variants(
        &self,
        var: &SequenceVariant,
        catalogue: Option<&dyn AccessionCatalogue>,
    ) -> Result<Vec<ConsistentVariant>, HgvsError> {
        let mut consistent = Vec::new();
        for candidate in self.on_each_accession(var, catalogue)? {
            if !self.validate(&candidate)? {
                continue;
            }
            let info = match catalogue {
                Some(c) => c.get_accession_info(candidate.ac())?,
                None => None,
            };
            consistent.push(ConsistentVariant {
                states_bases: states_bases(&candidate),
                variant: candidate,
                info,
            });
        }
        Ok(consistent)
    }

    /// `var` on each versioned accession its name stands for, unchecked:
    /// itself on a versioned accession, on each of a gene's accessions of its
    /// kind, or on each version of an unversioned accession.
    fn on_each_accession(
        &self,
        var: &SequenceVariant,
        catalogue: Option<&dyn AccessionCatalogue>,
    ) -> Result<Vec<SequenceVariant>, HgvsError> {
        let name = var.ac();
        if self.refs.identifier_type(name)? == IdentifierType::GeneSymbol {
            return self.on_each_of_the_genes_accessions(var);
        }
        if is_versioned(name) {
            return Ok(vec![var.clone()]);
        }
        let catalogue = catalogue.ok_or_else(|| {
            HgvsError::DataProviderError(format!(
                "{name} is not a versioned accession, and there is no accession catalogue to read it"
            ))
        })?;
        let versions = catalogue.get_accession_versions(name)?;
        if versions.is_empty() {
            return Err(HgvsError::DataProviderError(format!(
                "no version of {name} is held"
            )));
        }
        let bare = versions.len() == 1 && versions[0] == name;
        if !bare {
            if let Some(other) = versions.iter().find(|ac| !is_a_version_of(ac, name)) {
                return Err(HgvsError::DataProviderError(format!(
                    "the accession catalogue lists {other} as a version of {name}"
                )));
            }
            if let Some(twice) = versions
                .iter()
                .enumerate()
                .find_map(|(i, ac)| versions[..i].contains(ac).then_some(ac))
            {
                return Err(HgvsError::DataProviderError(format!(
                    "the accession catalogue lists {twice} twice as a version of {name}"
                )));
            }
        }
        Ok(versions.into_iter().map(|ac| on(var, ac)).collect())
    }

    /// `var`, on a gene symbol, on each of the gene's accessions of its kind:
    /// transcripts for `c.`, `n.` and `r.`, proteins for `p.`, genomic
    /// sequences for `g.`.
    pub(crate) fn on_each_of_the_genes_accessions(
        &self,
        var: &SequenceVariant,
    ) -> Result<Vec<SequenceVariant>, HgvsError> {
        use IdentifierType as T;
        let (target, compatible): (IdentifierKind, &[T]) = match var.coordinate_type() {
            "p" => (IdentifierKind::Protein, &[T::ProteinAccession]),
            "c" | "n" | "r" => (IdentifierKind::Transcript, &[T::TranscriptAccession]),
            _ => (IdentifierKind::Genomic, &[T::GenomicAccession]),
        };
        let mut on_each = Vec::new();
        for (kind, ac) in
            self.provider()
                .get_symbol_accessions(var.ac(), IdentifierKind::Gene, target)?
        {
            if compatible.contains(&kind) && self.carries(&ac, var)? {
                on_each.push(on(var, ac));
            }
        }
        if on_each.is_empty() {
            return Err(HgvsError::DataProviderError(format!(
                "no accession of {} for a {}. description is held",
                var.ac(),
                var.coordinate_type()
            )));
        }
        Ok(on_each)
    }
}

impl VariantMapper<'_> {
    /// Whether the transcript `ac` can carry `var` as written: a `c.`
    /// description needs a CDS, and an `r.` one a CDS where it names a
    /// position counted from one; an `n.` one numbers any transcript.
    fn carries(&self, ac: &str, var: &SequenceVariant) -> Result<bool, HgvsError> {
        let kind = var.coordinate_type();
        if !matches!(kind, "c" | "r") {
            return Ok(true);
        }
        let coding = self.has_cds(ac, ReadAs::Description(kind))?;
        Ok(coding || (kind == "r" && !counts_from_the_cds(var)))
    }
}

/// Whether an `r.` description, or a member of one, names a position counted
/// from the CDS: before its first base, `r.-2`, or after its stop, `r.*5`.
fn counts_from_the_cds(var: &SequenceVariant) -> bool {
    let from_the_cds = |p: &BaseOffsetPosition| p.anchor == Anchor::CdsEnd || p.base.0 < 1;
    match var {
        SequenceVariant::Rna(r) => r.posedit.pos.as_ref().is_some_and(|pos| {
            from_the_cds(&pos.start) || pos.end.as_ref().is_some_and(from_the_cds)
        }),
        SequenceVariant::CisPhased(cis) => cis.members.iter().any(counts_from_the_cds),
        _ => false,
    }
}

/// `var` on `ac`.
fn on(var: &SequenceVariant, ac: String) -> SequenceVariant {
    let mut var = var.clone();
    var.set_ac(ac);
    var
}

/// Whether `ac` is `name` followed by `.` and a version number.
fn is_a_version_of(ac: &str, name: &str) -> bool {
    ac.strip_prefix(name)
        .and_then(|v| v.strip_prefix('.'))
        .is_some_and(|v| is_digits(v) && !v.starts_with('0'))
}

fn is_digits(s: &str) -> bool {
    !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit())
}

/// Whether `ac` ends in a version, `.` and digits.
fn is_versioned(ac: &str) -> bool {
    ac.rsplit_once('.').is_some_and(|(_, v)| is_digits(v))
}

/// Whether `var` states bases or residues its sequence must hold: a
/// nucleotide edit's reference bases or repeat unit, a protein position's
/// residue, or any member's.
fn states_bases(var: &SequenceVariant) -> bool {
    match var {
        SequenceVariant::Genomic(v) => v.posedit.edit.states_bases(),
        SequenceVariant::Mitochondrial(v) => v.posedit.edit.states_bases(),
        SequenceVariant::Coding(v) => v.posedit.edit.states_bases(),
        SequenceVariant::NonCoding(v) => v.posedit.edit.states_bases(),
        SequenceVariant::Rna(v) => v.posedit.edit.states_bases(),
        SequenceVariant::Protein(v) => v.posedit.pos.is_some(),
        SequenceVariant::CisPhased(cis) => cis.members.iter().any(states_bases),
    }
}
