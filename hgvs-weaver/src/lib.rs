use crate::structs::{LinearVariant, TranscriptVariant};
use pest::Parser;
use pest_derive::Parser as PestParser;

#[derive(PestParser)]
#[grammar = "grammar.pest"]
pub struct HgvsParser;

/// Parses an HGVS string into a `SequenceVariant`.
pub fn parse_hgvs_variant(hgvs_str: &str) -> Result<SequenceVariant, HgvsError> {
    let mut pairs = HgvsParser::parse(Rule::hgvs_variant, hgvs_str)
        .map_err(|e| HgvsError::PestError(e.to_string()))?;

    let pair = pairs
        .next()
        .ok_or_else(|| HgvsError::PestError("Empty input".into()))?;
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
        | Rule::p_cis_variant => {
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
pub mod reference;
pub mod refget;
pub mod structs;
pub mod transcript_mapper;
pub mod transform;
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
pub use transform::{transform_variant, StartCodonConvention, VariantTransformSettings};
