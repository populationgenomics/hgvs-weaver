//! GA4GH VRS 2.0 rendering of a [`CanonicalAllele`] and of a copy-number
//! count, with computed identifiers.
//!
//! A VRS `Allele` is a `SequenceLocation` on a refget-identified sequence plus
//! a `state`: a literal sequence, or a `ReferenceLengthExpression` when the
//! alternate is the reference repeated. A `CopyNumberCount` is a
//! `SequenceLocation` plus the number of copies of it. Identifiers are
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

/// The VRS object a variant renders as: an `Allele` for a sequence change, a
/// `CopyNumberCount` for a copy-number edit. The JSON `type` tells them apart.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum VrsVariation {
    Allele(VrsAllele),
    CopyNumberCount(VrsCopyNumberCount),
}

impl VrsVariation {
    /// The computed identifier, `ga4gh:VA.` or `ga4gh:CN.` plus the digest.
    pub fn id(&self) -> &str {
        match self {
            VrsVariation::Allele(a) => &a.id,
            VrsVariation::CopyNumberCount(c) => &c.id,
        }
    }

    pub fn to_json(&self) -> String {
        match self {
            VrsVariation::Allele(a) => a.to_json(),
            VrsVariation::CopyNumberCount(c) => c.to_json(),
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
