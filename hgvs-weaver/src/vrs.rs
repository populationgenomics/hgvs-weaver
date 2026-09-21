//! GA4GH VRS 2.0 rendering of a [`CanonicalAllele`] and of a copy-number
//! count, with computed identifiers.
//!
//! A VRS `Allele` is a `SequenceLocation` on a refget-identified sequence plus
//! a `state`: a literal sequence, a `ReferenceLengthExpression` when the
//! alternate is the reference repeated, or a `LengthExpression` when only the
//! number of inserted bases is known. A `CopyNumberCount` is a
//! `SequenceLocation` plus the number of copies of it; a `CopyNumberChange`
//! is a `SequenceLocation` plus the direction of a change in copies (a gain
//! or a loss, as an EFO copy-number term). Identifiers are
//! `sha512t24u` digests of an RFC 8785 canonical JSON serialisation of each
//! object's inherent properties, nested identifiable objects replaced by their
//! digests.

use crate::allele::CanonicalAllele;
use crate::error::HgvsError;
use base64::Engine;
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use serde_json::{json, Value};
use sha2::{Digest, Sha512};

/// The GA4GH truncated digest: SHA-512, first 24 bytes, base64url without padding.
pub fn sha512t24u(data: &[u8]) -> String {
    let digest = Sha512::digest(data);
    base64::engine::general_purpose::URL_SAFE_NO_PAD.encode(&digest[..24])
}

/// The refget accession of a sequence, `SQ.` plus the digest of its residues.
pub fn refget_accession(sequence: &str) -> String {
    format!("SQ.{}", sha512t24u(sequence.as_bytes()))
}

/// Canonical JSON per RFC 8785 for the objects VRS digests: keys sorted, no
/// whitespace. `serde_json::Value` maps are already key-sorted.
fn canonical_json(value: &Value) -> String {
    serde_json::to_string(value).expect("serialising a JSON value cannot fail")
}

/// What the sequence an allele sits on is, for VRS's `SequenceReference`.
/// Not part of any computed identifier.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VrsMolecule {
    Genomic,
    Protein,
}

impl VrsMolecule {
    pub fn residue_alphabet(self) -> &'static str {
        match self {
            VrsMolecule::Genomic => "na",
            VrsMolecule::Protein => "aa",
        }
    }

    pub fn molecule_type(self) -> &'static str {
        match self {
            VrsMolecule::Genomic => "genomic",
            VrsMolecule::Protein => "protein",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct VrsSequenceReference {
    #[serde(rename = "type")]
    pub type_: String,
    #[serde(rename = "refgetAccession")]
    pub refget_accession: String,
    #[serde(rename = "residueAlphabet", default)]
    pub residue_alphabet: String,
    #[serde(rename = "moleculeType", default)]
    pub molecule_type: String,
}

/// A VRS integer-or-`Range` value: an exact number, or `[min, max]` when it is
/// only known to lie within a range (`None` for an unbounded side). One end of
/// a `SequenceLocation`, or the `copies` of a `CopyNumberCount`. Serialises as
/// a number or a two-element array with `null`s, which is also its form in
/// computed identifiers.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VrsBound {
    Exact(usize),
    Range(Option<usize>, Option<usize>),
}

impl<'de> Deserialize<'de> for VrsBound {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        #[derive(Deserialize)]
        #[serde(untagged)]
        enum Raw {
            Exact(usize),
            Range(Option<usize>, Option<usize>),
        }
        Ok(match Raw::deserialize(deserializer)? {
            Raw::Exact(n) => VrsBound::Exact(n),
            Raw::Range(min, max) => VrsBound::Range(min, max),
        })
    }
}

impl Serialize for VrsBound {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        match self {
            VrsBound::Exact(n) => serializer.serialize_u64(*n as u64),
            VrsBound::Range(min, max) => (min, max).serialize(serializer),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct VrsSequenceLocation {
    #[serde(default)]
    pub id: String,
    #[serde(rename = "type")]
    pub type_: String,
    #[serde(default)]
    pub digest: String,
    #[serde(rename = "sequenceReference")]
    pub sequence_reference: VrsSequenceReference,
    pub start: VrsBound,
    pub end: VrsBound,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(untagged)]
pub enum VrsState {
    ReferenceLength {
        #[serde(rename = "type")]
        type_: String,
        length: usize,
        #[serde(rename = "repeatSubunitLength")]
        repeat_subunit_length: usize,
        #[serde(default)]
        sequence: String,
    },
    Literal {
        #[serde(rename = "type")]
        type_: String,
        sequence: String,
    },
    /// A sequence known only by its length, HGVS `insN[20]`. Last of the
    /// untagged variants: it needs `length` alone, which the others also have.
    Length {
        #[serde(rename = "type")]
        type_: String,
        length: VrsBound,
    },
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct VrsExpression {
    pub syntax: String,
    pub value: String,
}

/// The `SequenceLocation` over interbase `[start, end)` of the sequence
/// `refget`, with its computed identifier. The inherent properties are the
/// type, the sequence reference (not itself identifiable, so serialised in
/// full without its descriptive fields) and the two bounds.
fn sequence_location(
    refget: &str,
    start: VrsBound,
    end: VrsBound,
    molecule: VrsMolecule,
) -> VrsSequenceLocation {
    let digest = sha512t24u(
        canonical_json(&json!({
            "type": "SequenceLocation",
            "sequenceReference": {"type": "SequenceReference", "refgetAccession": refget},
            "start": start,
            "end": end,
        }))
        .as_bytes(),
    );
    VrsSequenceLocation {
        id: format!("ga4gh:SL.{digest}"),
        type_: "SequenceLocation".into(),
        digest,
        sequence_reference: VrsSequenceReference {
            type_: "SequenceReference".into(),
            refget_accession: refget.to_string(),
            residue_alphabet: molecule.residue_alphabet().into(),
            molecule_type: molecule.molecule_type().into(),
        },
        start,
        end,
    }
}

/// The `expressions` list carrying an HGVS string, empty when there is none.
fn expressions(hgvs: Option<(&str, &str)>) -> Vec<VrsExpression> {
    hgvs.map(|(syntax, value)| {
        vec![VrsExpression {
            syntax: syntax.into(),
            value: value.into(),
        }]
    })
    .unwrap_or_default()
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct VrsAllele {
    #[serde(default)]
    pub id: String,
    #[serde(rename = "type")]
    pub type_: String,
    #[serde(default)]
    pub digest: String,
    pub location: VrsSequenceLocation,
    pub state: VrsState,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub expressions: Vec<VrsExpression>,
}

impl VrsAllele {
    /// Builds the VRS Allele for `allele` on the sequence identified by
    /// `refget`, carrying `hgvs` (syntax such as `hgvs.g` and the string) as an
    /// expression when given.
    pub fn new(
        allele: &CanonicalAllele,
        refget: &str,
        molecule: VrsMolecule,
        hgvs: Option<(&str, &str)>,
    ) -> Self {
        let state = match allele.repeat_subunit {
            Some(unit) => VrsState::ReferenceLength {
                type_: "ReferenceLengthExpression".into(),
                length: allele.alternate.len(),
                repeat_subunit_length: unit,
                sequence: allele.alternate.clone(),
            },
            None => VrsState::Literal {
                type_: "LiteralSequenceExpression".into(),
                sequence: allele.alternate.clone(),
            },
        };
        Self::build(
            refget,
            VrsBound::Exact(allele.start),
            VrsBound::Exact(allele.end),
            state,
            molecule,
            hgvs,
        )
    }

    /// An insertion of bases known only by number, HGVS `insN[20]` or
    /// `delinsN[20]`: interbase `[start, end)` (empty for a pure insertion)
    /// becomes `length` unspecified bases, a `LengthExpression` state. Unknown
    /// bases cannot slide, so the allele is rendered where it was written.
    pub fn length_expression(
        refget: &str,
        start: VrsBound,
        end: VrsBound,
        length: VrsBound,
        molecule: VrsMolecule,
        hgvs: Option<(&str, &str)>,
    ) -> Self {
        let state = VrsState::Length {
            type_: "LengthExpression".into(),
            length,
        };
        Self::build(refget, start, end, state, molecule, hgvs)
    }

    /// A deletion whose breakpoints are only known to lie within ranges, HGVS
    /// `g.(a_b)_(c_d)del`. The location carries the ranges and the state is
    /// the empty literal sequence. Such an allele cannot be normalised, so it
    /// is rendered as given.
    pub fn imprecise_deletion(
        refget: &str,
        start: VrsBound,
        end: VrsBound,
        molecule: VrsMolecule,
        hgvs: Option<(&str, &str)>,
    ) -> Self {
        let state = VrsState::Literal {
            type_: "LiteralSequenceExpression".into(),
            sequence: String::new(),
        };
        Self::build(refget, start, end, state, molecule, hgvs)
    }

    fn build(
        refget: &str,
        start: VrsBound,
        end: VrsBound,
        state: VrsState,
        molecule: VrsMolecule,
        hgvs: Option<(&str, &str)>,
    ) -> Self {
        let location = sequence_location(refget, start, end, molecule);
        let state_inherent = match &state {
            VrsState::Literal { sequence, .. } => {
                json!({"type": "LiteralSequenceExpression", "sequence": sequence})
            }
            VrsState::ReferenceLength {
                length,
                repeat_subunit_length,
                ..
            } => json!({
                "type": "ReferenceLengthExpression",
                "length": length,
                "repeatSubunitLength": repeat_subunit_length,
            }),
            // The inherent property of LengthExpression is `length` (VRS 2.0
            // vrs-source.yaml, `LengthExpression.ga4gh.inherent`), plus `type`
            // like every digested object.
            VrsState::Length { length, .. } => {
                json!({"type": "LengthExpression", "length": length})
            }
        };
        let allele_digest = sha512t24u(
            canonical_json(&json!({
                "type": "Allele",
                "location": location.digest,
                "state": state_inherent,
            }))
            .as_bytes(),
        );
        VrsAllele {
            id: format!("ga4gh:VA.{allele_digest}"),
            type_: "Allele".into(),
            digest: allele_digest,
            location,
            state,
            expressions: expressions(hgvs),
        }
    }

    /// Parses a VRS 2.0 Allele from JSON. Properties this module does not
    /// model are ignored; the `type` must be `Allele`.
    pub fn from_json(json: &str) -> Result<VrsAllele, HgvsError> {
        let allele: VrsAllele = serde_json::from_str(json)
            .map_err(|e| HgvsError::ValidationError(format!("Not a VRS Allele: {e}")))?;
        if allele.type_ != "Allele" {
            return Err(HgvsError::ValidationError(format!(
                "Expected a VRS Allele, got a {}",
                allele.type_
            )));
        }
        Ok(allele)
    }

    pub fn to_json(&self) -> String {
        serde_json::to_string(self).expect("serialising a VRS allele cannot fail")
    }
}

/// A VRS 2.0 `CopyNumberCount`: the number of copies of a location in a
/// genome, HGVS `g.1000_2000copy3`.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct VrsCopyNumberCount {
    #[serde(default)]
    pub id: String,
    #[serde(rename = "type")]
    pub type_: String,
    #[serde(default)]
    pub digest: String,
    pub location: VrsSequenceLocation,
    /// An exact count from HGVS; a `Range` is accepted when reading, as the
    /// schema allows it.
    pub copies: VrsBound,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub expressions: Vec<VrsExpression>,
}

impl VrsCopyNumberCount {
    /// `copies` copies of interbase `[start, end)` of the sequence `refget`,
    /// carrying `hgvs` (syntax such as `hgvs.g` and the string) as an
    /// expression when given.
    pub fn new(
        refget: &str,
        start: VrsBound,
        end: VrsBound,
        copies: VrsBound,
        molecule: VrsMolecule,
        hgvs: Option<(&str, &str)>,
    ) -> Self {
        let location = sequence_location(refget, start, end, molecule);
        // The inherent properties of CopyNumberCount are `location` and
        // `copies` (VRS 2.0 vrs-source.yaml, `CopyNumberCount.ga4gh.inherent`,
        // prefix `CN`); `type` goes into every digest by the computed
        // identifier convention.
        let digest = sha512t24u(
            canonical_json(&json!({
                "type": "CopyNumberCount",
                "location": location.digest,
                "copies": copies,
            }))
            .as_bytes(),
        );
        VrsCopyNumberCount {
            id: format!("ga4gh:CN.{digest}"),
            type_: "CopyNumberCount".into(),
            digest,
            location,
            copies,
            expressions: expressions(hgvs),
        }
    }

    /// Parses a VRS 2.0 CopyNumberCount from JSON. Properties this module does
    /// not model are ignored; the `type` must be `CopyNumberCount`.
    pub fn from_json(json: &str) -> Result<VrsCopyNumberCount, HgvsError> {
        let count: VrsCopyNumberCount = serde_json::from_str(json)
            .map_err(|e| HgvsError::ValidationError(format!("Not a VRS CopyNumberCount: {e}")))?;
        if count.type_ != "CopyNumberCount" {
            return Err(HgvsError::ValidationError(format!(
                "Expected a VRS CopyNumberCount, got a {}",
                count.type_
            )));
        }
        Ok(count)
    }

    pub fn to_json(&self) -> String {
        serde_json::to_string(self).expect("serialising a VRS copy number count cannot fail")
    }
}

/// The direction and degree of a `CopyNumberChange`: the EFO copy-number
/// terms VRS 2.0.1 enumerates by label (vrs-source.yaml,
/// `CopyNumberChange.properties.copyChange.enum`). VRS 2.0.0 gave the same
/// terms as EFO codes in a `MappableConcept`; both forms are read, the label
/// is written. `Other` keeps a term this module does not know.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum VrsCopyChange {
    /// EFO:0030069, complete genomic deletion.
    CompleteGenomicLoss,
    /// EFO:0020073, high-level copy number loss.
    HighLevelLoss,
    /// EFO:0030068, low-level copy number loss.
    LowLevelLoss,
    /// EFO:0030067, copy number loss.
    Loss,
    /// EFO:0030064, regional base ploidy.
    RegionalBasePloidy,
    /// EFO:0030070, copy number gain.
    Gain,
    /// EFO:0030071, low-level copy number gain.
    LowLevelGain,
    /// EFO:0030072, high-level copy number gain.
    HighLevelGain,
    Other(String),
}

impl VrsCopyChange {
    const KNOWN: [VrsCopyChange; 8] = [
        VrsCopyChange::CompleteGenomicLoss,
        VrsCopyChange::HighLevelLoss,
        VrsCopyChange::LowLevelLoss,
        VrsCopyChange::Loss,
        VrsCopyChange::RegionalBasePloidy,
        VrsCopyChange::Gain,
        VrsCopyChange::LowLevelGain,
        VrsCopyChange::HighLevelGain,
    ];

    /// The label VRS 2.0.1 serialises, or the unknown term as given.
    pub fn label(&self) -> &str {
        match self {
            VrsCopyChange::CompleteGenomicLoss => "complete genomic loss",
            VrsCopyChange::HighLevelLoss => "high-level loss",
            VrsCopyChange::LowLevelLoss => "low-level loss",
            VrsCopyChange::Loss => "loss",
            VrsCopyChange::RegionalBasePloidy => "regional base ploidy",
            VrsCopyChange::Gain => "gain",
            VrsCopyChange::LowLevelGain => "low-level gain",
            VrsCopyChange::HighLevelGain => "high-level gain",
            VrsCopyChange::Other(term) => term,
        }
    }

    /// The EFO code, `None` for an unknown term.
    pub fn efo(&self) -> Option<&'static str> {
        Some(match self {
            VrsCopyChange::CompleteGenomicLoss => "EFO:0030069",
            VrsCopyChange::HighLevelLoss => "EFO:0020073",
            VrsCopyChange::LowLevelLoss => "EFO:0030068",
            VrsCopyChange::Loss => "EFO:0030067",
            VrsCopyChange::RegionalBasePloidy => "EFO:0030064",
            VrsCopyChange::Gain => "EFO:0030070",
            VrsCopyChange::LowLevelGain => "EFO:0030071",
            VrsCopyChange::HighLevelGain => "EFO:0030072",
            VrsCopyChange::Other(_) => return None,
        })
    }

    /// Whether the term is a gain of copies: a duplication in HGVS.
    pub fn is_gain(&self) -> bool {
        matches!(
            self,
            VrsCopyChange::Gain | VrsCopyChange::LowLevelGain | VrsCopyChange::HighLevelGain
        )
    }

    /// Whether the term is a loss of copies: a deletion in HGVS.
    pub fn is_loss(&self) -> bool {
        matches!(
            self,
            VrsCopyChange::Loss
                | VrsCopyChange::LowLevelLoss
                | VrsCopyChange::HighLevelLoss
                | VrsCopyChange::CompleteGenomicLoss
        )
    }

    /// The term whose label or EFO code is `text`; `Other` when neither.
    pub fn parse(text: &str) -> VrsCopyChange {
        Self::KNOWN
            .iter()
            .find(|c| c.label() == text || c.efo() == Some(text))
            .cloned()
            .unwrap_or_else(|| VrsCopyChange::Other(text.to_string()))
    }
}

impl Serialize for VrsCopyChange {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.serialize_str(self.label())
    }
}

impl<'de> Deserialize<'de> for VrsCopyChange {
    /// A label or EFO code as a string (VRS 2.0.1), or a `MappableConcept`
    /// whose `primaryCoding.code` (or `primaryCode`) is the EFO code (VRS
    /// 2.0.0).
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let value = Value::deserialize(deserializer)?;
        let text = match &value {
            Value::String(s) => Some(s.as_str()),
            Value::Object(concept) => concept
                .get("primaryCoding")
                .and_then(|coding| coding.get("code"))
                .or_else(|| concept.get("primaryCode"))
                .and_then(Value::as_str),
            _ => None,
        };
        text.map(VrsCopyChange::parse).ok_or_else(|| {
            serde::de::Error::custom(format!(
                "copyChange must be a copy-number term or a MappableConcept, not {value}"
            ))
        })
    }
}

/// A VRS 2.0 `CopyNumberChange`: a gain or loss of copies of a location
/// relative to the baseline ploidy, without a count. HGVS
/// `g.(a_b)_(c_d)dup` (a gain) or `del` (a loss).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct VrsCopyNumberChange {
    #[serde(default)]
    pub id: String,
    #[serde(rename = "type")]
    pub type_: String,
    #[serde(default)]
    pub digest: String,
    pub location: VrsSequenceLocation,
    #[serde(rename = "copyChange")]
    pub copy_change: VrsCopyChange,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub expressions: Vec<VrsExpression>,
}

impl VrsCopyNumberChange {
    /// A `copy_change` of interbase `[start, end)` of the sequence `refget`,
    /// carrying `hgvs` (syntax such as `hgvs.g` and the string) as an
    /// expression when given.
    pub fn new(
        refget: &str,
        start: VrsBound,
        end: VrsBound,
        copy_change: VrsCopyChange,
        molecule: VrsMolecule,
        hgvs: Option<(&str, &str)>,
    ) -> Self {
        let location = sequence_location(refget, start, end, molecule);
        // The inherent properties of CopyNumberChange are `location` and
        // `copyChange` (VRS 2.0.1 vrs-source.yaml,
        // `CopyNumberChange.ga4gh.inherent`, prefix `CX`); `type` goes into
        // every digest by the computed identifier convention.
        let digest = sha512t24u(
            canonical_json(&json!({
                "type": "CopyNumberChange",
                "location": location.digest,
                "copyChange": copy_change,
            }))
            .as_bytes(),
        );
        VrsCopyNumberChange {
            id: format!("ga4gh:CX.{digest}"),
            type_: "CopyNumberChange".into(),
            digest,
            location,
            copy_change,
            expressions: expressions(hgvs),
        }
    }

    /// Parses a VRS 2.0 CopyNumberChange from JSON. Properties this module
    /// does not model are ignored; the `type` must be `CopyNumberChange`.
    pub fn from_json(json: &str) -> Result<VrsCopyNumberChange, HgvsError> {
        let change: VrsCopyNumberChange = serde_json::from_str(json)
            .map_err(|e| HgvsError::ValidationError(format!("Not a VRS CopyNumberChange: {e}")))?;
        if change.type_ != "CopyNumberChange" {
            return Err(HgvsError::ValidationError(format!(
                "Expected a VRS CopyNumberChange, got a {}",
                change.type_
            )));
        }
        Ok(change)
    }

    pub fn to_json(&self) -> String {
        serde_json::to_string(self).expect("serialising a VRS copy number change cannot fail")
    }
}

/// The VRS object a variant renders as: an `Allele` for a sequence change, a
/// `CopyNumberCount` for a copy-number edit, a `CopyNumberChange` for a gain
/// or loss without a count. The JSON `type` tells them apart.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum VrsVariation {
    Allele(VrsAllele),
    CopyNumberCount(VrsCopyNumberCount),
    CopyNumberChange(VrsCopyNumberChange),
}

impl VrsVariation {
    /// The computed identifier, `ga4gh:VA.`, `ga4gh:CN.` or `ga4gh:CX.` plus
    /// the digest.
    pub fn id(&self) -> &str {
        match self {
            VrsVariation::Allele(a) => &a.id,
            VrsVariation::CopyNumberCount(c) => &c.id,
            VrsVariation::CopyNumberChange(c) => &c.id,
        }
    }

    pub fn to_json(&self) -> String {
        match self {
            VrsVariation::Allele(a) => a.to_json(),
            VrsVariation::CopyNumberCount(c) => c.to_json(),
            VrsVariation::CopyNumberChange(c) => c.to_json(),
        }
    }
}

/// The `type` of a VRS object given as JSON, so a reader can pick the class
/// before parsing it in full. Empty when there is none.
pub fn vrs_type(json: &str) -> Result<String, HgvsError> {
    let value: Value = serde_json::from_str(json)
        .map_err(|e| HgvsError::ValidationError(format!("Not a VRS object: {e}")))?;
    Ok(value
        .get("type")
        .and_then(Value::as_str)
        .unwrap_or_default()
        .to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn digest_and_identifier_match_the_spec_example() {
        // https://vrs.ga4gh.org/en/stable/conventions/computed_identifiers.html
        let allele = CanonicalAllele {
            accession: "NC_000019.10".into(),
            start: 44908821,
            end: 44908822,
            reference: "C".into(),
            alternate: "T".into(),
            repeat_subunit: None,
        };
        let vrs = VrsAllele::new(
            &allele,
            "SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl",
            VrsMolecule::Genomic,
            None,
        );
        assert_eq!(vrs.location.digest, "wIlaGykfwHIpPY2Fcxtbx4TINbbODFVz");
        assert_eq!(vrs.id, "ga4gh:VA.0AePZIWZUNsUlQTamyLrjm2HWUw2opLt");
    }

    #[test]
    fn copy_number_count_matches_the_spec_example() {
        // https://vrs.ga4gh.org/en/stable/concepts/SystemicVariation/CopyNumberCount.html
        // gives ga4gh:CN.ezEUXykQvIhX8jHADILwC9f8k-jp8tZC for three or more
        // copies of [44905795, 44909393) of SQ.jdEWLvLvT8827O59m1Agh5H3n6kTzBsJ.
        let count = VrsCopyNumberCount::new(
            "SQ.jdEWLvLvT8827O59m1Agh5H3n6kTzBsJ",
            VrsBound::Exact(44905795),
            VrsBound::Exact(44909393),
            VrsBound::Range(Some(3), None),
            VrsMolecule::Genomic,
            None,
        );
        assert_eq!(count.id, "ga4gh:CN.ezEUXykQvIhX8jHADILwC9f8k-jp8tZC");
        assert_eq!(vrs_type(&count.to_json()).unwrap(), "CopyNumberCount");
        assert_eq!(
            VrsCopyNumberCount::from_json(&count.to_json()).unwrap(),
            count
        );
    }

    #[test]
    fn copy_number_change_digests_the_label_over_the_location() {
        // https://vrs.ga4gh.org/en/stable/concepts/SystemicVariation/CopyNumberChange.html
        // shows a low-level gain of [44905795, 44909393) of
        // SQ.jdEWLvLvT8827O59m1Agh5H3n6kTzBsJ with `"copyChange": "low-level
        // gain"`, the VRS 2.0.1 label enum (vrs-source.yaml, vrs-python
        // `CopyChange`), which is what is digested here.
        let change = VrsCopyNumberChange::new(
            "SQ.jdEWLvLvT8827O59m1Agh5H3n6kTzBsJ",
            VrsBound::Exact(44905795),
            VrsBound::Exact(44909393),
            VrsCopyChange::LowLevelGain,
            VrsMolecule::Genomic,
            None,
        );
        assert_eq!(change.location.digest, "d9h3FkfTWFkJSH56L1A26y-N2oq_SSuB");
        assert_eq!(change.id, "ga4gh:CX._rPTdFeOE9elAozZsakJGTqCvlaiEyr6");
        assert!(
            change
                .to_json()
                .contains(r#""copyChange":"low-level gain""#),
            "{}",
            change.to_json()
        );
        assert_eq!(vrs_type(&change.to_json()).unwrap(), "CopyNumberChange");
        assert_eq!(
            VrsCopyNumberChange::from_json(&change.to_json()).unwrap(),
            change
        );
        // The identifier the page prints, ga4gh:CX.2_fT_6-IpUm5aS0wp8ZAkJ01MCE569L2,
        // predates the enum: it is the digest with the bare EFO code as the
        // `copyChange` string. Reproduced here so that the digest machinery
        // is pinned to a published value.
        let pre_release = VrsCopyNumberChange::new(
            "SQ.jdEWLvLvT8827O59m1Agh5H3n6kTzBsJ",
            VrsBound::Exact(44905795),
            VrsBound::Exact(44909393),
            VrsCopyChange::Other("EFO:0030071".into()),
            VrsMolecule::Genomic,
            None,
        );
        assert_eq!(pre_release.id, "ga4gh:CX.2_fT_6-IpUm5aS0wp8ZAkJ01MCE569L2");
    }

    #[test]
    fn copy_change_terms_are_read_by_label_or_efo_code() {
        for term in VrsCopyChange::KNOWN {
            assert_eq!(VrsCopyChange::parse(term.label()), term);
            assert_eq!(VrsCopyChange::parse(term.efo().unwrap()), term);
        }
        assert!(VrsCopyChange::Gain.is_gain() && !VrsCopyChange::Gain.is_loss());
        assert!(VrsCopyChange::CompleteGenomicLoss.is_loss());
        assert!(!VrsCopyChange::RegionalBasePloidy.is_gain());
        assert!(!VrsCopyChange::RegionalBasePloidy.is_loss());
        assert_eq!(
            VrsCopyChange::parse("EFO:0000001"),
            VrsCopyChange::Other("EFO:0000001".into())
        );
        // The VRS 2.0.0 MappableConcept form is read too.
        let json = r#"{"type":"CopyNumberChange","copyChange":{"primaryCoding":{"code":"EFO:0030067","system":"https://www.ebi.ac.uk/efo/"}},"location":{"type":"SequenceLocation","sequenceReference":{"type":"SequenceReference","refgetAccession":"SQ.x"},"start":1,"end":2}}"#;
        let change = VrsCopyNumberChange::from_json(json).unwrap();
        assert_eq!(change.copy_change, VrsCopyChange::Loss);
        let json = json.replace(
            r#"{"primaryCoding":{"code":"EFO:0030067","system":"https://www.ebi.ac.uk/efo/"}}"#,
            r#"{"primaryCode":"EFO:0030072"}"#,
        );
        assert_eq!(
            VrsCopyNumberChange::from_json(&json).unwrap().copy_change,
            VrsCopyChange::HighLevelGain
        );
        assert!(VrsCopyNumberChange::from_json(
            &json.replace(r#"{"primaryCode":"EFO:0030072"}"#, "7")
        )
        .is_err());
    }

    #[test]
    fn a_length_expression_state_carries_the_number_of_bases() {
        let exact = VrsAllele::length_expression(
            "SQ.test",
            VrsBound::Exact(10),
            VrsBound::Exact(10),
            VrsBound::Exact(20),
            VrsMolecule::Genomic,
            Some(("hgvs.g", "X:g.10_11insN[20]")),
        );
        assert!(exact.id.starts_with("ga4gh:VA."), "{}", exact.id);
        assert!(
            exact
                .to_json()
                .contains(r#""state":{"type":"LengthExpression","length":20}"#),
            "{}",
            exact.to_json()
        );
        assert_eq!(VrsAllele::from_json(&exact.to_json()).unwrap(), exact);
        assert!(matches!(
            VrsAllele::from_json(&exact.to_json()).unwrap().state,
            VrsState::Length {
                length: VrsBound::Exact(20),
                ..
            }
        ));
        // The digest is over the inherent properties only: the same allele
        // without the expression has the same identifier.
        let bare = VrsAllele::length_expression(
            "SQ.test",
            VrsBound::Exact(10),
            VrsBound::Exact(10),
            VrsBound::Exact(20),
            VrsMolecule::Genomic,
            None,
        );
        assert_eq!(bare.id, exact.id);
        // A range of lengths serialises as [min, max] and changes the digest.
        let ranged = VrsAllele::length_expression(
            "SQ.test",
            VrsBound::Exact(10),
            VrsBound::Exact(10),
            VrsBound::Range(Some(20), Some(30)),
            VrsMolecule::Genomic,
            None,
        );
        assert!(
            ranged.to_json().contains(r#""length":[20,30]"#),
            "{}",
            ranged.to_json()
        );
        assert_ne!(ranged.id, exact.id);
        // The other states are still read as themselves.
        let literal = VrsAllele::imprecise_deletion(
            "SQ.test",
            VrsBound::Exact(10),
            VrsBound::Exact(12),
            VrsMolecule::Genomic,
            None,
        );
        assert!(matches!(
            VrsAllele::from_json(&literal.to_json()).unwrap().state,
            VrsState::Literal { .. }
        ));
    }

    #[test]
    fn sha512t24u_of_empty_input_is_the_known_value() {
        // Documented refget example: the digest of the empty sequence.
        assert_eq!(sha512t24u(b""), "z4PhNX7vuL3xVChQ1m2AB9Yg5AULVxXc");
    }

    #[test]
    fn reference_derived_state_is_a_reference_length_expression() {
        let allele = CanonicalAllele {
            accession: "X".into(),
            start: 2,
            end: 8,
            reference: "CAGCAG".into(),
            alternate: "CAGCAGCAG".into(),
            repeat_subunit: Some(3),
        };
        let vrs = VrsAllele::new(
            &allele,
            "SQ.test",
            VrsMolecule::Genomic,
            Some(("hgvs.g", "X:g.3_8dup")),
        );
        match &vrs.state {
            VrsState::ReferenceLength {
                length,
                repeat_subunit_length,
                ..
            } => {
                assert_eq!((*length, *repeat_subunit_length), (9, 3));
            }
            other => panic!("expected ReferenceLengthExpression, got {other:?}"),
        }
        assert!(vrs
            .to_json()
            .contains("\"expressions\":[{\"syntax\":\"hgvs.g\",\"value\":\"X:g.3_8dup\"}]"));
    }
}
