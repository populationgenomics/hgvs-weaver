use crate::structs::{LinearVariant, TranscriptVariant};
use pest::Parser;
use pest_derive::Parser as PestParser;

#[derive(PestParser)]
#[grammar = "grammar.pest"]
pub struct HgvsParser;

/// How [`parse_hgvs_variant_with`] reads an `X` where a stop and any residue
/// could both stand, in a `p.` description of either code. An
/// `X` with one possible reading is read as it, whatever the treatment.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum XTreatment {
    /// Refuse the description as [`HgvsError::AmbiguousResidue`].
    #[default]
    Error,
    /// Read the `X` as a stop, `*`: the convention before `X` meant any residue.
    Stop,
    /// Read the `X` as any residue, `Xaa`.
    Unknown,
}

/// How [`parse_hgvs_variant_with`] reads what a description leaves open.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub struct ParseOptions {
    pub x_treatment: XTreatment,
}

impl ParseOptions {
    pub fn with_x_treatment(self, x_treatment: XTreatment) -> Self {
        Self { x_treatment }
    }
}

/// Parses an HGVS string into a `SequenceVariant`, refusing an ambiguous `X`;
/// [`parse_hgvs_variant_with`] with the default options.
pub fn parse_hgvs_variant(hgvs_str: &str) -> Result<SequenceVariant, HgvsError> {
    parse_hgvs_variant_with(hgvs_str, ParseOptions::default())
}

/// Parses an HGVS string into a `SequenceVariant`.
///
/// A `p.` description writes every residue in one amino acid code; `*`,
/// `Xaa` and `X` stand in either. One that mixes the codes is
/// [`HgvsError::MixedResidueCodes`]. An `X` where a stop and any residue
/// could both stand is read as `options.x_treatment` says: the description
/// with `*` or `Xaa` written in its place, or [`HgvsError::AmbiguousResidue`].
pub fn parse_hgvs_variant_with(
    hgvs_str: &str,
    options: ParseOptions,
) -> Result<SequenceVariant, HgvsError> {
    let mut pairs = HgvsParser::parse(Rule::hgvs_variant, hgvs_str)
        .map_err(|e| HgvsError::PestError(e.to_string()))?;

    let pair = pairs
        .next()
        .ok_or_else(|| HgvsError::PestError("Empty input".into()))?;
    parser::one_residue_code(&pair)?;
    let ambiguous = parser::ambiguous_x(&pair);
    if !ambiguous.is_empty() {
        let written = match options.x_treatment {
            XTreatment::Error => {
                return Err(HgvsError::AmbiguousResidue {
                    description: hgvs_str.to_string(),
                    ambiguous,
                });
            }
            XTreatment::Stop => "*",
            XTreatment::Unknown => "Xaa",
        };
        let mut resolved = String::with_capacity(hgvs_str.len() + 2 * ambiguous.len());
        let mut at = 0;
        for x in &ambiguous {
            resolved.push_str(&hgvs_str[at..x.start]);
            resolved.push_str(written);
            at = x.end;
        }
        resolved.push_str(&hgvs_str[at..]);
        return parse_hgvs_variant(&resolved);
    }
    let inner = pair
        .into_inner()
        .next()
        .ok_or_else(|| HgvsError::PestError("Missing inner variant".into()))?;

    let system = inner.as_rule();
    let (ac, gene, posedit, posedits) = variant_parts(inner)?;
    Ok(match system {
        Rule::g_variant => SequenceVariant::Genomic(GVariant::from_parts(
            ac,
            gene,
            parser::parse_g_posedit(posedit)?,
        )),
        Rule::m_variant => SequenceVariant::Mitochondrial(MVariant::from_parts(
            ac,
            gene,
            parser::parse_g_posedit(posedit)?,
        )),
        Rule::c_variant => SequenceVariant::Coding(CVariant::from_parts(
            ac,
            gene,
            parser::parse_tx_posedit(posedit, CVariant::DEFAULT_ANCHOR)?,
        )),
        Rule::n_variant => SequenceVariant::NonCoding(NVariant::from_parts(
            ac,
            gene,
            parser::parse_tx_posedit(posedit, NVariant::DEFAULT_ANCHOR)?,
        )),
        Rule::r_variant => SequenceVariant::Rna(RVariant {
            ac,
            gene,
            posedit: parser::parse_tx_posedit(posedit, coords::Anchor::TranscriptStart)?,
        }),
        Rule::p_variant => SequenceVariant::Protein(PVariant {
            ac,
            gene,
            posedit: parser::parse_p_posedit(posedit)?,
        }),
        Rule::g_cis_variant
        | Rule::m_cis_variant
        | Rule::c_cis_variant
        | Rule::n_cis_variant
        | Rule::r_cis_variant
        | Rule::p_cis_variant
        | Rule::p_cis_predicted => {
            let mut members = Vec::with_capacity(1 + posedits.len());
            for p in std::iter::once(posedit).chain(posedits) {
                members.push(cis_member(system, ac.clone(), gene.clone(), p)?);
            }
            SequenceVariant::CisPhased(CisPhasedVariant::new(ac, gene, members)?)
        }
        Rule::trans_alleles => {
            return Err(HgvsError::UnsupportedOperation(format!(
                "{hgvs_str} describes two molecules (alleles in trans, `[..];[..]`), \
                 which is not one variant; parse each `[..]` as a cis allele instead"
            )))
        }
        other => {
            return Err(HgvsError::PestError(format!(
                "Unsupported variant type: {:?}",
                other
            )))
        }
    })
}

/// A member of a cis allele: the posedit `p` as a plain variant of the system
/// the `<x>_cis_variant` rule `cis` names, on `ac`.
fn cis_member(
    cis: Rule,
    ac: String,
    gene: Option<String>,
    p: pest::iterators::Pair<Rule>,
) -> Result<SequenceVariant, HgvsError> {
    Ok(match cis {
        Rule::g_cis_variant => {
            SequenceVariant::Genomic(GVariant::from_parts(ac, gene, parser::parse_g_posedit(p)?))
        }
        Rule::m_cis_variant => SequenceVariant::Mitochondrial(MVariant::from_parts(
            ac,
            gene,
            parser::parse_g_posedit(p)?,
        )),
        Rule::c_cis_variant => SequenceVariant::Coding(CVariant::from_parts(
            ac,
            gene,
            parser::parse_tx_posedit(p, CVariant::DEFAULT_ANCHOR)?,
        )),
        Rule::n_cis_variant => SequenceVariant::NonCoding(NVariant::from_parts(
            ac,
            gene,
            parser::parse_tx_posedit(p, NVariant::DEFAULT_ANCHOR)?,
        )),
        Rule::r_cis_variant => SequenceVariant::Rna(RVariant {
            ac,
            gene,
            posedit: parser::parse_tx_posedit(p, coords::Anchor::TranscriptStart)?,
        }),
        Rule::p_cis_variant => SequenceVariant::Protein(PVariant {
            ac,
            gene,
            posedit: parser::parse_p_posedit(p)?,
        }),
        Rule::p_cis_predicted => {
            let mut posedit = parser::parse_p_posedit(p)?;
            posedit.predicted = true;
            SequenceVariant::Protein(PVariant { ac, gene, posedit })
        }
        other => {
            return Err(HgvsError::PestError(format!(
                "Not a cis allele rule: {:?}",
                other
            )))
        }
    })
}

/// Splits a `<x>_variant` or `<x>_cis_variant` pair into accession, optional
/// gene symbol, the first posedit pair and any further ones (a cis allele's
/// other members); every coordinate system is written the same way up to there.
#[allow(clippy::type_complexity)]
fn variant_parts(
    pair: pest::iterators::Pair<Rule>,
) -> Result<
    (
        String,
        Option<String>,
        pest::iterators::Pair<Rule>,
        Vec<pest::iterators::Pair<Rule>>,
    ),
    HgvsError,
> {
    let mut inner = pair.into_inner();
    let ac = inner
        .next()
        .ok_or_else(|| HgvsError::PestError("Missing accession".into()))?
        .as_str()
        .to_string();
    let gene = parse_gene_expr(
        inner
            .next()
            .ok_or_else(|| HgvsError::PestError("Missing gene expr".into()))?,
    );
    let posedit = inner
        .next()
        .ok_or_else(|| HgvsError::PestError("Missing posedit".into()))?;
    Ok((ac, gene, posedit, inner.collect()))
}

fn parse_gene_expr(pair: pest::iterators::Pair<Rule>) -> Option<String> {
    let s = pair.as_str();
    if s.is_empty() {
        return None;
    }
    Some(s.replace(['(', ')'], ""))
}

pub mod allele;
pub mod cigar;
pub mod cis;
mod cis_protein;
mod cis_vrs;
pub mod consistent;
pub mod coords;
pub mod data;
pub mod edits;
pub mod equivalence;
pub mod error;
pub mod fmt;
pub mod mapper;
pub mod normalize;
pub mod parser;
pub mod protein;
mod protein_reading;
pub mod reading;
pub mod reference;
pub mod refget;
pub mod structs;
pub mod transcript_mapper;
pub mod utils;
pub mod vrs;

// Re-exports for public usage
pub use coords::SequenceVariant;
pub use data::{DataProvider, IdentifierKind, TranscriptData, TranscriptSearch};
pub use equivalence::VariantEquivalence;
pub use error::HgvsError;
pub use mapper::VariantMapper;
pub use structs::{
    CVariant, CisPhasedVariant, GVariant, MVariant, NVariant, PVariant, RVariant, Variant,
};
