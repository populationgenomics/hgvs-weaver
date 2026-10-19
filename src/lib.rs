use ::hgvs_weaver::{
    DataProvider, HgvsError, IdentifierKind, SequenceVariant, TranscriptData, TranscriptSearch,
    Variant as VariantTrait, VariantMapper,
};
use pyo3::prelude::*;
use pyo3::types::PyDict;
use pyo3::{Bound, PyErr};

pyo3::create_exception!(_weaver, HGVSError, pyo3::exceptions::PyException);
pyo3::create_exception!(_weaver, ParseError, HGVSError);
pyo3::create_exception!(_weaver, ValidationError, HGVSError);
pyo3::create_exception!(_weaver, DataProviderError, HGVSError);
pyo3::create_exception!(_weaver, UnsupportedOperationError, HGVSError);
pyo3::create_exception!(_weaver, CigarError, HGVSError);
pyo3::create_exception!(
    _weaver,
    ReferenceMismatchError,
    HGVSError,
    "A base or residue a description states is not what its reference holds: the transcript record for an exonic c., n. or r. base, the genome (on the transcript's strand) for an intronic one, the genome for g. and m., the protein for p. Raised by every operation that reads the stated bases; Variant.validate answers False instead."
);

pyo3::create_exception!(
    _weaver,
    PositionOutsideSequenceError,
    HGVSError,
    "A description names a position its reference does not have: past the end of the genome for g. and m., of the transcript record for an exonic c., n. or r. position, of the protein for p., where a protein served without its stop has the stop one past its last residue. An insertion needs both its flanking bases, so one after the last base is outside too. Raised by every operation that reads the positions; Variant.validate answers False instead."
);

pyo3::create_exception!(
    _weaver,
    AmbiguousDescriptionError,
    HGVSError,
    "A description with several readings that make different sequences or state different bases. weaver does not pick one: readings holds each, as a Variant with one reading. A transcript allele with an intronic member is read on the genome, and its readings are written in g."
);

pyo3::create_exception!(
    _weaver,
    NoProteinStatementError,
    HGVSError,
    "A description makes no statement about the protein: an identity that covers no whole codon, which states bases rather than residues. HGVS has no p. form for it: p.(=) would claim the whole protein, and p.? an effect."
);

fn map_hgvs_error(e: HgvsError) -> PyErr {
    py_error(e, String::new())
}

/// The Python exception for `e`, its message after `context`. A member of an
/// allele in cis fails with its own kind of error, named in the message.
fn py_error(e: HgvsError, context: String) -> PyErr {
    match e {
        HgvsError::PestError(msg) => ParseError::new_err(context + &msg),
        HgvsError::ValidationError(msg) => ValidationError::new_err(context + &msg),
        HgvsError::DataProviderError(msg) => DataProviderError::new_err(context + &msg),
        HgvsError::UnsupportedOperation(msg) => UnsupportedOperationError::new_err(context + &msg),
        HgvsError::CigarError(msg) => CigarError::new_err(context + &msg),
        HgvsError::ReferenceMismatch {
            accession,
            expected,
            found,
            start,
            end,
        } => ReferenceMismatchError::new_err(format!(
            "{context}expected {expected}, found {found} on {accession} at indices {start}..{end}"
        )),
        HgvsError::PositionOutsideSequence {
            accession,
            start,
            end,
            length,
        } => PositionOutsideSequenceError::new_err(format!(
            "{context}{accession} has {length} positions, and the 0-based range [{start}, {end}) runs past its end"
        )),
        HgvsError::NoProteinStatement(msg) => NoProteinStatementError::new_err(context + &msg),
        HgvsError::Other(msg) => HGVSError::new_err(context + &msg),
        HgvsError::Ambiguous {
            description,
            readings,
        } => ambiguity_error(context, description, readings),
        HgvsError::CisMember {
            allele,
            index,
            member,
            source,
        } => py_error(
            *source,
            format!("{context}Member {} of {allele}, {member}: ", index + 1),
        ),
        // HgvsError is non_exhaustive; a variant added without its own arm is
        // still raised, as the base class, with its message.
        other => HGVSError::new_err(context + &other.to_string()),
    }
}

/// An `AmbiguousDescriptionError` whose `readings` are the readings, as
/// variants.
fn ambiguity_error(context: String, description: String, readings: Vec<SequenceVariant>) -> PyErr {
    let listed: Vec<String> = readings.iter().map(ToString::to_string).collect();
    let err = AmbiguousDescriptionError::new_err(format!(
        "{context}{description} reads as {}",
        listed.join(" | ")
    ));
    Python::attach(|py| {
        let variants = readings
            .into_iter()
            .map(|inner| Py::new(py, PyVariant { inner }))
            .collect::<PyResult<Vec<_>>>()?;
        let list = pyo3::types::PyList::new(py, variants)?;
        err.value(py).setattr("readings", list)?;
        Ok::<(), PyErr>(())
    })
    .err()
    .unwrap_or(err)
}

#[pyclass(name = "IdentifierType", module = "weaver._weaver")]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PyIdentifierType {
    GenomicAccession,
    TranscriptAccession,
    ProteinAccession,
    GeneSymbol,
    Unknown,
}

impl From<PyIdentifierType> for ::hgvs_weaver::data::IdentifierType {
    fn from(it: PyIdentifierType) -> Self {
        match it {
            PyIdentifierType::GenomicAccession => {
                ::hgvs_weaver::data::IdentifierType::GenomicAccession
            }
            PyIdentifierType::TranscriptAccession => {
                ::hgvs_weaver::data::IdentifierType::TranscriptAccession
            }
            PyIdentifierType::ProteinAccession => {
                ::hgvs_weaver::data::IdentifierType::ProteinAccession
            }
            PyIdentifierType::GeneSymbol => ::hgvs_weaver::data::IdentifierType::GeneSymbol,
            PyIdentifierType::Unknown => ::hgvs_weaver::data::IdentifierType::Unknown,
        }
    }
}

impl From<::hgvs_weaver::data::IdentifierType> for PyIdentifierType {
    fn from(it: ::hgvs_weaver::data::IdentifierType) -> Self {
        match it {
            ::hgvs_weaver::data::IdentifierType::GenomicAccession => {
                PyIdentifierType::GenomicAccession
            }
            ::hgvs_weaver::data::IdentifierType::TranscriptAccession => {
                PyIdentifierType::TranscriptAccession
            }
            ::hgvs_weaver::data::IdentifierType::ProteinAccession => {
                PyIdentifierType::ProteinAccession
            }
            ::hgvs_weaver::data::IdentifierType::GeneSymbol => PyIdentifierType::GeneSymbol,
            ::hgvs_weaver::data::IdentifierType::Unknown => PyIdentifierType::Unknown,
        }
    }
}

#[pymethods]
impl PyIdentifierType {
    fn __repr__(&self) -> String {
        format!("IdentifierType.{:?}", self)
    }
    fn __eq__(&self, other: &Self) -> bool {
        self == other
    }
    fn __hash__(&self) -> u64 {
        let mut s = std::collections::hash_map::DefaultHasher::new();
        std::hash::Hash::hash(self, &mut s);
        std::hash::Hasher::finish(&s)
    }
}

#[pyclass(name = "EquivalenceLevel", module = "weaver._weaver")]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PyEquivalenceLevel {
    Identity,
    Analogous,
    Different,
    Unknown,
}

impl From<::hgvs_weaver::equivalence::EquivalenceLevel> for PyEquivalenceLevel {
    fn from(el: ::hgvs_weaver::equivalence::EquivalenceLevel) -> Self {
        match el {
            ::hgvs_weaver::equivalence::EquivalenceLevel::Identity => Self::Identity,
            ::hgvs_weaver::equivalence::EquivalenceLevel::Analogous => Self::Analogous,
            ::hgvs_weaver::equivalence::EquivalenceLevel::Different => Self::Different,
            ::hgvs_weaver::equivalence::EquivalenceLevel::Unknown => Self::Unknown,
        }
    }
}

#[pymethods]
impl PyEquivalenceLevel {
    fn __repr__(&self) -> String {
        format!("EquivalenceLevel.{:?}", self)
    }
    fn __eq__(&self, other: &Self) -> bool {
        self == other
    }
    fn __hash__(&self) -> u64 {
        let mut s = std::collections::hash_map::DefaultHasher::new();
        std::hash::Hash::hash(self, &mut s);
        std::hash::Hasher::finish(&s)
    }
}

#[pyclass(name = "Variant", module = "weaver._weaver")]
#[doc = "Represents a parsed HGVS variant.\n\nProvides access to the variant's accession, gene symbol, and coordinate type.\nVariants can be formatted back to HGVS strings or converted to JSON/dict representations.\n\nAn allele in cis, NM_004006.2:c.[145C>T;147C>G], is one Variant whose\nmembers are the plain variants written between the brackets."]
#[derive(Clone)]
pub struct PyVariant {
    pub inner: SequenceVariant,
}

#[pymethods]
impl PyVariant {
    #[getter]
    #[doc = "The primary accession of the variant (e.g., 'NM_000051.3')."]
    fn ac(&self) -> String {
        self.inner.ac().to_string()
    }

    #[getter]
    #[doc = "The gene symbol associated with the variant, if available."]
    fn gene(&self) -> Option<String> {
        self.inner.gene().map(|s| s.to_string())
    }

    #[getter]
    #[doc = "The coordinate type of the variant ('g', 'c', 'p', etc.); for an allele in cis, its members'."]
    fn coordinate_type(&self) -> String {
        self.inner.coordinate_type().to_string()
    }

    #[getter]
    #[doc = "The members of an allele in cis, c.[145C>T;147C>G], as plain Variants in the\norder written; None for a plain variant.\n\nThe projections between g., c., n. and r. take the whole allele and\nproject it member by member."]
    fn members(&self) -> Option<Vec<PyVariant>> {
        match &self.inner {
            SequenceVariant::CisPhased(cis) => Some(
                cis.members
                    .iter()
                    .map(|m| PyVariant { inner: m.clone() })
                    .collect(),
            ),
            _ => None,
        }
    }

    #[doc = "Formats the variant back into a standard HGVS string."]
    fn format(&self) -> String {
        self.inner.to_string()
    }

    #[doc = "Returns a JSON string representation of the internal variant structure."]
    fn to_json(&self) -> PyResult<String> {
        serde_json::to_string(&self.inner)
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))
    }

    #[doc = "Returns a dictionary representation of the internal variant structure."]
    fn to_dict(&self, py: Python) -> PyResult<Py<PyAny>> {
        let json_str = self.to_json()?;
        let json_mod = py.import("json")?;
        let dict = json_mod.call_method1("loads", (json_str,))?;
        Ok(dict.unbind())
    }

    #[staticmethod]
    #[doc = "Constructs a Variant from a dictionary produced by to_dict.\n\nArgs:\n    d: A dict with the same structure as returned by to_dict.\n\nReturns:\n    A Variant object.\n\nRaises:\n    ValueError: If the dict cannot be deserialised into a valid variant."]
    fn from_dict(py: Python, d: Py<PyAny>) -> PyResult<PyVariant> {
        let json_mod = py.import("json")?;
        let json_str: String = json_mod.call_method1("dumps", (d,))?.extract::<String>()?;
        let inner: SequenceVariant = serde_json::from_str(&json_str)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;
        Ok(PyVariant { inner })
    }

    fn __str__(&self) -> String {
        self.format()
    }

    fn __repr__(&self) -> String {
        format!("<weaver.Variant {}>", self.format())
    }

    fn __reduce__(&self, py: Python) -> PyResult<Py<PyAny>> {
        let module = py.import("weaver._weaver")?;
        let parse_fn: pyo3::Bound<'_, pyo3::PyAny> = module.getattr("parse")?;
        let args: pyo3::Bound<'_, pyo3::PyAny> =
            pyo3::types::PyTuple::new(py, [self.format()])?.into_any();
        Ok(pyo3::types::PyTuple::new(py, [parse_fn, args])?
            .into_any()
            .unbind())
    }

    #[doc = "Validates the positions the variant names, and the bases and residues it states, against a DataProvider.\n\nAn exonic c., n. or r. base is read on the transcript record, an intronic one\non the genome (on the transcript's strand), a g. or m. base on the genome and a\np. residue on the protein. An allele in cis of nucleotide members must also have\none reading: members that reach for the same base are read every way the\ndescription allows.\n\nArgs:\n    provider: The data provider instance for sequence retrieval.\n\nReturns:\n    True if every position the variant names and every stated base is there,\n    False where an operation reading the variant would raise\n    ReferenceMismatchError or PositionOutsideSequenceError. c.= names no bases\n    and is True once its transcript is found to have a CDS.\n\nRaises:\n    ValidationError: If the variant names no position or a range that runs\n        backwards, or c.= is on a transcript without a CDS. An intronic c. or n. variant\n        is projected to the genome to be read, so a projection's error is raised too.\n        An allele in cis whose members leave it no reading is a ValidationError naming two.\n    AmbiguousDescriptionError: If an allele in cis has several readings; its readings\n        are each of them.\n    UnsupportedOperationError: If an allele in cis has readings weaver does not list.\n    DataProviderError: If sequence data cannot be retrieved."]
    fn validate(&self, _py: Python, provider: Py<PyAny>) -> PyResult<bool> {
        let bridge = PyDataProviderBridge { provider };
        VariantMapper::new(&bridge)
            .validate(&self.inner)
            .map_err(map_hgvs_error)
    }
}

#[pyfunction]
#[doc = "Parses an HGVS string into a Variant object.\n\nSupported types include genomic (g.), coding cDNA (c.), non-coding (n.),\nmitochondrial (m.), RNA (r.) and protein (p.) variants, and alleles in cis\nin any of them, c.[145C>T;147C>G]: several changes on one molecule, one\nVariant whose members are the changes. The trans form c.[145C>T];[147C>G]\ndescribes two molecules and is refused.\n\nArgs:\n    input: The HGVS string to parse.\n\nReturns:\n    A Variant object.\n\nRaises:\n    ParseError: If the HGVS string is malformed or unsupported.\n    UnsupportedOperationError: If the string describes alleles in trans."]
fn parse(input: &str) -> PyResult<PyVariant> {
    match ::hgvs_weaver::parse_hgvs_variant(input) {
        Ok(inner) => Ok(PyVariant { inner }),
        Err(e) => Err(map_hgvs_error(e)),
    }
}

/// The typed variants a mapper method accepts, or the ValueError Python callers expect.
macro_rules! expect_variant {
    ($name:ident, $variant:ident, $ty:ty, $what:literal) => {
        fn $name(var: &SequenceVariant) -> PyResult<&$ty> {
            match var {
                SequenceVariant::$variant(v) => Ok(v),
                SequenceVariant::CisPhased(_) => Err(pyo3::exceptions::PyValueError::new_err(
                    concat!("Expected a ", $what, ", not an allele in cis"),
                )),
                _ => Err(pyo3::exceptions::PyValueError::new_err(concat!(
                    "Expected a ",
                    $what
                ))),
            }
        }
    };
}

/// What a projection takes: a plain variant, or an allele in cis of members
/// in the same system, which projects member by member.
enum OneOrCis<'v, T> {
    One(&'v T),
    Cis(&'v ::hgvs_weaver::CisPhasedVariant),
}

/// The plain variant or allele in cis a projection accepts, or the ValueError
/// Python callers expect.
macro_rules! expect_one_or_cis {
    ($name:ident, $variant:ident, $ty:ty, $system:literal, $what:literal) => {
        fn $name(var: &SequenceVariant) -> PyResult<OneOrCis<'_, $ty>> {
            match var {
                SequenceVariant::$variant(v) => Ok(OneOrCis::One(v)),
                SequenceVariant::CisPhased(c) if c.coordinate_type() == $system => {
                    Ok(OneOrCis::Cis(c))
                }
                _ => Err(pyo3::exceptions::PyValueError::new_err(concat!(
                    "Expected a ",
                    $what,
                    " or an allele in cis of them"
                ))),
            }
        }
    };
}
expect_one_or_cis!(
    expect_genomic,
    Genomic,
    ::hgvs_weaver::GVariant,
    "g",
    "genomic variant (g.)"
);
expect_variant!(
    expect_coding,
    Coding,
    ::hgvs_weaver::CVariant,
    "coding variant (c.)"
);
expect_one_or_cis!(
    expect_noncoding,
    NonCoding,
    ::hgvs_weaver::structs::NVariant,
    "n",
    "non-coding variant (n.)"
);
expect_variant!(
    expect_protein,
    Protein,
    ::hgvs_weaver::PVariant,
    "protein variant (p.)"
);
expect_one_or_cis!(
    expect_rna,
    Rna,
    ::hgvs_weaver::RVariant,
    "r",
    "RNA variant (r.)"
);
expect_one_or_cis!(
    expect_coding_or_cis,
    Coding,
    ::hgvs_weaver::CVariant,
    "c",
    "coding variant (c.)"
);

// --- Mapper and DataProvider Bridge ---

pub struct PyDataProviderBridge {
    provider: Py<PyAny>,
}

/// A Python object with `get_refget_accession(ac)` and
/// `get_accession_for_refget(refget)`, each returning `str | None`.
pub struct PyRefgetBridge {
    refget: Py<PyAny>,
}

impl PyRefgetBridge {
    fn lookup(&self, method: &str, arg: &str) -> Result<Option<String>, HgvsError> {
        Python::attach(|py| {
            self.refget
                .bind(py)
                .call_method1(method, (arg,))
                .and_then(|r| r.extract::<Option<String>>())
                .map_err(|e: PyErr| HgvsError::DataProviderError(e.to_string()))
        })
    }
}

impl ::hgvs_weaver::refget::Refget for PyRefgetBridge {
    fn refget_accession(&self, ac: &str) -> Result<Option<String>, HgvsError> {
        self.lookup("get_refget_accession", ac)
    }

    fn accession_for_refget(&self, refget: &str) -> Result<Option<String>, HgvsError> {
        self.lookup("get_accession_for_refget", refget)
    }
}

impl DataProvider for PyDataProviderBridge {
    fn get_transcript(
        &self,
        transcript_ac: &str,
        reference_ac: Option<&str>,
    ) -> Result<TranscriptData, HgvsError> {
        Python::attach(|py| {
            let res = self
                .provider
                .bind(py)
                .call_method1("get_transcript", (transcript_ac, reference_ac))
                .map_err(|e| HgvsError::DataProviderError(e.to_string()))?;
            let dict = res.cast_into::<PyDict>().map_err(|e| {
                HgvsError::DataProviderError(format!("Failed to cast to PyDict: {}", e))
            })?;
            let json_mod = py
                .import("json")
                .map_err(|e| HgvsError::DataProviderError(e.to_string()))?;
            let json_str: String = json_mod
                .call_method1("dumps", (dict,))
                .map_err(|e| HgvsError::DataProviderError(e.to_string()))?
                .extract::<String>()
                .map_err(|e| HgvsError::DataProviderError(e.to_string()))?;
            let data: ::hgvs_weaver::data::TranscriptData = serde_json::from_str(&json_str)
                .map_err(|e| HgvsError::DataProviderError(e.to_string()))?;
            Ok(data)
        })
    }

    fn get_seq(
        &self,
        ac: &str,
        start: i32,
        end: Option<i32>,
        kind: ::hgvs_weaver::data::IdentifierType,
    ) -> Result<String, HgvsError> {
        Python::attach(|py| {
            let py_kind: PyIdentifierType = kind.into();
            let res = self
                .provider
                .bind(py)
                .call_method1("get_seq", (ac, start, end, py_kind))
                .map_err(|e: PyErr| HgvsError::DataProviderError(e.to_string()))?;
            res.extract::<String>()
                .map_err(|e: PyErr| HgvsError::DataProviderError(e.to_string()))
        })
    }

    fn get_symbol_accessions(
        &self,
        symbol: &str,
        source_kind: IdentifierKind,
        target_kind: IdentifierKind,
    ) -> Result<Vec<(::hgvs_weaver::data::IdentifierType, String)>, HgvsError> {
        Python::attach(|py| {
            let sk = match source_kind {
                IdentifierKind::Genomic => "g",
                IdentifierKind::Transcript => "c",
                IdentifierKind::Protein => "p",
            };
            let tk = match target_kind {
                IdentifierKind::Genomic => "g",
                IdentifierKind::Transcript => "c",
                IdentifierKind::Protein => "p",
            };
            let res = self
                .provider
                .bind(py)
                .call_method1("get_symbol_accessions", (symbol, sk, tk))
                .map_err(|e: PyErr| HgvsError::DataProviderError(e.to_string()))?;
            let raw_list: Vec<(Bound<'_, PyAny>, String)> = res
                .extract::<Vec<(Bound<'_, PyAny>, String)>>()
                .map_err(|e: PyErr| {
                    HgvsError::DataProviderError(format!(
                        "Failed to extract symbol accessions: {}",
                        e
                    ))
                })?;

            let mut result = Vec::new();
            for (type_any, ac) in raw_list {
                let it = if let Ok(s) = type_any.extract::<String>() {
                    match s.as_str() {
                        "genomic_accession" => {
                            ::hgvs_weaver::data::IdentifierType::GenomicAccession
                        }
                        "transcript_accession" => {
                            ::hgvs_weaver::data::IdentifierType::TranscriptAccession
                        }
                        "protein_accession" => {
                            ::hgvs_weaver::data::IdentifierType::ProteinAccession
                        }
                        "gene_symbol" => ::hgvs_weaver::data::IdentifierType::GeneSymbol,
                        _ => ::hgvs_weaver::data::IdentifierType::Unknown,
                    }
                } else if let Ok(py_it) = type_any.extract::<PyIdentifierType>() {
                    py_it.into()
                } else {
                    ::hgvs_weaver::data::IdentifierType::Unknown
                };
                result.push((it, ac));
            }
            Ok(result)
        })
    }

    fn get_identifier_type(
        &self,
        identifier: &str,
    ) -> Result<::hgvs_weaver::data::IdentifierType, HgvsError> {
        Python::attach(|py| {
            let res = self
                .provider
                .bind(py)
                .call_method1("get_identifier_type", (identifier,))
                .map_err(|e| HgvsError::DataProviderError(e.to_string()))?;

            // Try to extract as String first (for backward compatibility or simpler mocks)
            if let Ok(s) = res.extract::<String>() {
                return Ok(match s.as_str() {
                    "genomic_accession" => ::hgvs_weaver::data::IdentifierType::GenomicAccession,
                    "transcript_accession" => {
                        ::hgvs_weaver::data::IdentifierType::TranscriptAccession
                    }
                    "protein_accession" => ::hgvs_weaver::data::IdentifierType::ProteinAccession,
                    "gene_symbol" => ::hgvs_weaver::data::IdentifierType::GeneSymbol,
                    _ => ::hgvs_weaver::data::IdentifierType::Unknown,
                });
            }

            // Otherwise try to extract as the enum type
            let py_it: PyIdentifierType = res.extract::<PyIdentifierType>().map_err(|e| {
                HgvsError::DataProviderError(format!("Failed to extract IdentifierType: {}", e))
            })?;
            Ok(py_it.into())
        })
    }
}

pub struct PyTranscriptSearchBridge {
    searcher: Py<PyAny>,
}

impl TranscriptSearch for PyTranscriptSearchBridge {
    fn get_transcripts_for_region(
        &self,
        chrom: &str,
        start: i32,
        end: i32,
    ) -> Result<Vec<String>, HgvsError> {
        Python::attach(|py| {
            let res = self
                .searcher
                .bind(py)
                .call_method1("get_transcripts_for_region", (chrom, start, end))
                .map_err(|e| HgvsError::DataProviderError(e.to_string()))?;
            res.extract::<Vec<String>>()
                .map_err(|e| HgvsError::DataProviderError(e.to_string()))
        })
    }
}

#[pyclass(name = "VariantMapper", module = "weaver._weaver")]
#[doc = "High-level variant mapping engine.\n\nCoordinates mapping between different reference sequences (e.g., g. to c.)\nand projects cDNA variants onto protein sequences (c. to p.).\nRequires a DataProvider to retrieve transcript and sequence information."]
pub struct PyVariantMapper {
    pub bridge: PyDataProviderBridge,
    /// Sequence blocks and refget accessions, kept for the mapper's lifetime.
    cache: std::sync::Arc<::hgvs_weaver::reference::SequenceCache>,
    refget: Option<PyRefgetBridge>,
}

impl PyVariantMapper {
    /// A mapper over this object's provider, cache and refget lookup.
    fn mapper(&self) -> VariantMapper<'_> {
        let refget = self
            .refget
            .as_ref()
            .map(|r| r as &dyn ::hgvs_weaver::refget::Refget);
        VariantMapper::from_store(::hgvs_weaver::reference::ReferenceStore::shared(
            &self.bridge,
            std::sync::Arc::clone(&self.cache),
            refget,
        ))
    }
}

#[pymethods]
impl PyVariantMapper {
    #[new]
    #[pyo3(signature = (provider, refget=None))]
    #[doc = "Creates a new VariantMapper with the given DataProvider.\n\nThe mapper caches the sequence blocks and refget accessions it fetches for\nas long as it lives, so keep one and reuse it.\n\nArgs:\n    provider: The DataProvider for transcripts and sequences.\n    refget: Optional Refget lookup: an object with get_refget_accession(ac)\n        and get_accession_for_refget(refget), each returning str | None (a\n        weaver.refget.RefgetProvider is one). Without it refget accessions\n        are computed from the whole sequence and from_vrs needs the\n        accession passed."]
    fn new(provider: Py<PyAny>, refget: Option<Py<PyAny>>) -> Self {
        PyVariantMapper {
            bridge: PyDataProviderBridge { provider },
            cache: std::sync::Arc::new(::hgvs_weaver::reference::SequenceCache::new()),
            refget: refget.map(|refget| PyRefgetBridge { refget }),
        }
    }
    #[pyo3(signature = (var_g, transcript_ac))]
    #[doc = "Maps a genomic variant (g.) to a coding cDNA variant (c.) for a specific transcript.\n\nThe edit is re-read against the transcript: a stated reference the genome has\nbut the record does not becomes the record's base, and a change the record\nalready carries is written =; no change on the genome (g.123=) is the record's\nbase becoming the genome's where the two differ, and a duplication or\ninversion carries the genome's bases. An intronic position has no transcript base and\nis carried as given. A range whose ends lie in two exons is written over\nthe c. range that numbers the same bases, intron included: g.19_31dup\nover c.4, c.5, the intron and c.6 is c.4_6dup. Where the record's bases\nthere differ from the genome's, only a statement projects, as the\nrecord's correction.\n\nAn allele in cis of g. members, g.[a;b], projects member by member to an\nallele in cis on the transcript, its members written in reverse on a\nminus-strand transcript. A member that fails raises its own error,\nnaming the member.\n\nArgs:\n    var_g: The genomic Variant to map.\n    transcript_ac: The accession of the target transcript.\n\nReturns:\n    A new Variant object in 'c.' coordinates.\n\nRaises:\n    ValueError: If var_g is not a genomic variant or an allele in cis of them.\n    UnsupportedOperationError: If a range whose ends lie in two exons\n        changes bases where the record differs from the genome, or ends\n        on bases only one of them holds.\n    HGVSError: If mapping fails due to data or alignment issues."]
    fn g_to_c(&self, _py: Python, var_g: &PyVariant, transcript_ac: String) -> PyResult<PyVariant> {
        let mapper = self.mapper();
        let inner = match expect_genomic(&var_g.inner)? {
            OneOrCis::One(v) => mapper
                .g_to_c(v, &transcript_ac)
                .map(SequenceVariant::Coding),
            OneOrCis::Cis(c) => mapper
                .g_to_c_cis(c, &transcript_ac)
                .map(SequenceVariant::CisPhased),
        }
        .map_err(map_hgvs_error)?;
        Ok(PyVariant { inner })
    }

    #[pyo3(signature = (var_g, searcher))]
    #[doc = "Maps a genomic variant (g.) to all overlapping transcripts discovered via the searcher.\n\nAn allele in cis of g. members gives an allele in cis on each transcript\nevery member projects to, as g_to_c does.\n\nArgs:\n    var_g: The genomic Variant to map.\n    searcher: An object implementing the TranscriptSearch protocol.\n\nReturns:\n    A list of Variant objects in 'c.' coordinates.\n\nRaises:\n    ValueError: If var_g is not a genomic variant or an allele in cis of them.\n    HGVSError: If mapping fails due to data or alignment issues."]
    fn g_to_c_all(
        &self,
        _py: Python,
        var_g: &PyVariant,
        searcher: Py<PyAny>,
    ) -> PyResult<Vec<PyVariant>> {
        let mapper = self.mapper();
        let bridge_searcher = PyTranscriptSearchBridge { searcher };
        let res: Vec<SequenceVariant> = match expect_genomic(&var_g.inner)? {
            OneOrCis::One(v) => mapper
                .g_to_c_all(v, &bridge_searcher)
                .map(|all| all.into_iter().map(SequenceVariant::Coding).collect()),
            OneOrCis::Cis(c) => mapper
                .g_to_c_all_cis(c, &bridge_searcher)
                .map(|all| all.into_iter().map(SequenceVariant::CisPhased).collect()),
        }
        .map_err(map_hgvs_error)?;
        Ok(res.into_iter().map(|inner| PyVariant { inner }).collect())
    }

    #[pyo3(signature = (var_c, reference_ac = None))]
    #[doc = "Maps a coding cDNA variant (c.) to a genomic variant (g.).\n\nA range across an exon junction is the genomic range it numbers, the\nintron included: c.4_6dup across the junction after c.5 duplicates\nc.4, c.5, the intron and c.6.\n\nThe edit is re-read against the genome: a stated reference the record has\nbut the genome does not becomes the genome's base, and a change the genome\nalready carries is written =; no change on the record (c.123=) is the genome's\nbase becoming the record's where the two differ, and a duplication or\ninversion carries the record's bases (c.123dup over a differing genomic base\nis g.…delinsCC, as VariantValidator writes it).\n\nAn allele in cis of c. members, c.[a;b], projects member by member to an\nallele in cis on the genome, its members written in reverse from a\nminus-strand transcript. A member that fails raises its own error,\nnaming the member.\n\nArgs:\n    var_c: The coding Variant to map.\n    reference_ac: Optional chromosomal accession. If not provided, the primary chromosome for the transcript will be used.\n\nReturns:\n    A new Variant object in 'g.' coordinates.\n\nRaises:\n    ValueError: If var_c is not a coding variant or an allele in cis of them.\n    ReferenceMismatchError: If a base the variant states is not the record's\n        (exonic) or the genome's on the transcript's strand (intronic).\n    PositionOutsideSequenceError: If an intronic position lies past the end\n        of the genome.\n    UnsupportedOperationError: If the variant inserts at an exon\n        junction (c.5_6insT where c.5 and c.6 lie in different\n        exons), which is before or after the intron on the genome.\n        The message names the spelling at each end of the intron,\n        c.5_5+1insT and c.6-1_6insT, either of which projects.\n    HGVSError: If mapping fails due to data or alignment issues."]
    fn c_to_g(
        &self,
        _py: Python,
        var_c: &PyVariant,
        reference_ac: Option<String>,
    ) -> PyResult<PyVariant> {
        let mapper = self.mapper();
        let reference_ac = reference_ac.as_deref();
        let inner = match expect_coding_or_cis(&var_c.inner)? {
            OneOrCis::One(v) => mapper.c_to_g(v, reference_ac).map(SequenceVariant::Genomic),
            OneOrCis::Cis(c) => mapper
                .c_to_g_cis(c, reference_ac)
                .map(SequenceVariant::CisPhased),
        }
        .map_err(map_hgvs_error)?;
        Ok(PyVariant { inner })
    }

    #[pyo3(signature = (var_n, reference_ac = None))]
    #[doc = "Maps a non-coding cDNA variant (n.) to a genomic variant (g.).\n\nAn allele in cis of n. members projects member by member, as c_to_g\nprojects one of c. members.\n\nArgs:\n    var_n: The non-coding Variant to map.\n    reference_ac: Optional chromosomal accession.\n\nReturns:\n    A new Variant object in 'g.' coordinates.\n\nRaises:\n    ValueError: If var_n is not a non-coding variant or an allele in cis of them.\n    UnsupportedOperationError: If the variant inserts at an exon\n        junction, as c_to_g refuses.\n    HGVSError: If mapping fails due to data or alignment issues."]
    fn n_to_g(
        &self,
        _py: Python,
        var_n: &PyVariant,
        reference_ac: Option<String>,
    ) -> PyResult<PyVariant> {
        let mapper = self.mapper();
        let reference_ac = reference_ac.as_deref();
        let inner = match expect_noncoding(&var_n.inner)? {
            OneOrCis::One(v) => mapper.n_to_g(v, reference_ac).map(SequenceVariant::Genomic),
            OneOrCis::Cis(c) => mapper
                .n_to_g_cis(c, reference_ac)
                .map(SequenceVariant::CisPhased),
        }
        .map_err(map_hgvs_error)?;
        Ok(PyVariant { inner })
    }

    #[pyo3(signature = (var_c, protein_ac=None))]
    #[doc = "Projects a coding cDNA variant (c.) or an RNA variant (r.) to its protein consequence (p.).\n\nAn r. variant is predicted from its c. spelling; a statement about the\ntranscript (r.0, r.spl, r.?, r.=) becomes the matching statement about the\nprotein (p.0?, p.?, p.(=)). c.=, the whole transcript unchanged, predicts\np.(=) on the protein the transcript encodes; a protein_ac naming another\nprotein is an error. Every form needs a transcript with a CDS. r.(?)\nsays the RNA follows a DNA change it does not carry, so it has no\nprotein consequence of its own and raises UnsupportedOperationError.\n\nAn allele in cis of c. members, c.[a;b], is predicted from the sequence\nits members make together, as one change: c.[4A>T;6A>T] predicts what\nc.4_6delinsTAT does. Changes apart are written as the members of a\npredicted protein allele, c.[4A>C;13T>G] as p.[(Lys2Gln;Tyr5Asp)], and\nresidues whose whole codon the allele states unchanged, and the\nchanges leave as they were, as members beside them, p.[(Ala4=;Tyr5Asp)].\n\nArgs:\n    var_c: The coding or RNA Variant to project, or an allele in cis of\n        c. members.\n    protein_ac: Optional protein accession. If not provided, it will be retrieved from the DataProvider.\n\nReturns:\n    A new Variant object in 'p.' coordinates: a single variant, or an\n    allele in cis of p. members.\n\nRaises:\n    ValueError: If var_c is none of a coding variant, an RNA variant or an\n        allele in cis of c. members.\n    UnsupportedOperationError: If var_c is r.(?), which has no protein consequence of its own,\n        or an allele with a member at an uncertain position or whose bases\n        are not stated, two members whose 3'-most spellings overlap, or a\n        member past the stop that shifts the frame of the read-through.\n    NoProteinStatementError: If var_c is an identity that covers no whole codon,\n        which states bases and no residue, or an allele that changes\n        nothing and states no whole codon.\n    AmbiguousDescriptionError: If an allele in cis has several readings.\n    ValidationError: If an allele in cis has no reading, or the\n        transcript's CDS is open at the 5' end.\n    HGVSError: If projection fails due to data retrieval or out-of-bounds coordinates."]
    fn c_to_p(
        &self,
        _py: Python,
        var_c: &PyVariant,
        protein_ac: Option<String>,
    ) -> PyResult<PyVariant> {
        let mapper = self.mapper();
        let protein_ac = protein_ac.as_deref();
        let inner = match &var_c.inner {
            SequenceVariant::Coding(v) => mapper.c_to_p(v, protein_ac).map(SequenceVariant::Protein),
            SequenceVariant::Rna(r) => mapper.r_to_p(r, protein_ac).map(SequenceVariant::Protein),
            SequenceVariant::CisPhased(c) if c.coordinate_type() == "c" => {
                mapper.c_to_p_cis(c, protein_ac)
            }
            _ => {
                return Err(pyo3::exceptions::PyValueError::new_err(
                    "Expected a coding variant (c.), an RNA variant (r.) or an allele in cis of c. members",
                ))
            }
        }
        .map_err(map_hgvs_error)?;
        Ok(PyVariant { inner })
    }

    #[pyo3(signature = (var_c, protein_ac=None))]
    #[doc = "Returns the GA4GH VRS 2.0 Allele of a coding variant's protein consequence, on the protein sequence.\n\nThe residues from the first change to the end of the protein become the\nresidues the edited transcript encodes, up to its new stop, then the allele\nis canonicalised. Unlike to_vrs on a p. variant this covers frameshifts,\nextensions and stop losses. The predicted p. description is carried as an\nhgvs.p expression.\n\nArgs:\n    var_c: The coding Variant.\n    protein_ac: Optional protein accession; else the provider's mapping for the transcript.\n\nRaises:\n    ValueError: If var_c is not a coding variant.\n    HGVSError: If the consequence is a statement (p.?, p.(Met1?)), or the translated CDS is not the protein the provider serves."]
    fn protein_vrs(
        &self,
        py: Python,
        var_c: &PyVariant,
        protein_ac: Option<String>,
    ) -> PyResult<Py<PyAny>> {
        let v = expect_coding(&var_c.inner)?;
        let mapper = self.mapper();
        let json_str = mapper
            .protein_vrs(v, protein_ac.as_deref())
            .map_err(map_hgvs_error)?
            .to_json();
        let json_mod = py.import("json")?;
        Ok(json_mod.call_method1("loads", (json_str,))?.unbind())
    }

    #[pyo3(signature = (var_c, protein_ac=None))]
    #[doc = "Returns the SPDI of a coding variant's protein consequence, on the protein sequence.\n\nSee protein_vrs for what the allele is.\n\nArgs:\n    var_c: The coding Variant.\n    protein_ac: Optional protein accession; else the provider's mapping for the transcript."]
    fn protein_spdi(
        &self,
        _py: Python,
        var_c: &PyVariant,
        protein_ac: Option<String>,
    ) -> PyResult<String> {
        let v = expect_coding(&var_c.inner)?;
        let mapper = self.mapper();
        Ok(mapper
            .protein_allele(v, protein_ac.as_deref())
            .map_err(map_hgvs_error)?
            .spdi())
    }

    #[pyo3(signature = (var_r, reference_ac = None))]
    #[doc = "Maps an RNA variant (r.) to a genomic variant (g.).\n\nOnly a change within one exon has a genomic form; one spanning a splice\njunction describes the spliced RNA and raises UnsupportedOperationError.\n\nAn allele in cis of r. members projects member by member, as c_to_g\nprojects one of c. members: a member spanning a splice junction fails the\nallele.\n\nArgs:\n    var_r: The RNA Variant to map.\n    reference_ac: Optional chromosomal accession.\n\nReturns:\n    A new Variant object in 'g.' coordinates.\n\nRaises:\n    ValueError: If var_r is not an RNA variant or an allele in cis of them.\n    HGVSError: If mapping fails."]
    fn r_to_g(
        &self,
        _py: Python,
        var_r: &PyVariant,
        reference_ac: Option<String>,
    ) -> PyResult<PyVariant> {
        let mapper = self.mapper();
        let reference_ac = reference_ac.as_deref();
        let inner = match expect_rna(&var_r.inner)? {
            OneOrCis::One(v) => mapper.r_to_g(v, reference_ac).map(SequenceVariant::Genomic),
            OneOrCis::Cis(c) => mapper
                .r_to_g_cis(c, reference_ac)
                .map(SequenceVariant::CisPhased),
        }
        .map_err(map_hgvs_error)?;
        Ok(PyVariant { inner })
    }

    #[pyo3(signature = (var_r))]
    #[doc = "Respells an RNA variant (r.) as the coding variant (c.) it is numbered from.\n\nr. positions on a coding transcript are c. positions; the bases become\nuppercase DNA letters (u to T).\n\nAn allele in cis of r. members is respelled member by member.\n\nRaises:\n    ValueError: If var_r is not an RNA variant or an allele in cis of them.\n    HGVSError: If the transcript has no CDS, or the variant is a statement about the transcript (r.0, r.spl)."]
    fn r_to_c(&self, _py: Python, var_r: &PyVariant) -> PyResult<PyVariant> {
        let mapper = self.mapper();
        let inner = match expect_rna(&var_r.inner)? {
            OneOrCis::One(v) => mapper.r_to_c(v).map(SequenceVariant::Coding),
            OneOrCis::Cis(c) => mapper.r_to_c_cis(c).map(SequenceVariant::CisPhased),
        }
        .map_err(map_hgvs_error)?;
        Ok(PyVariant { inner })
    }

    #[pyo3(signature = (var_r))]
    #[doc = "Respells an RNA variant (r.) on a non-coding transcript as the n. variant it is numbered from.\n\nAn allele in cis of r. members is respelled member by member.\n\nRaises:\n    ValueError: If var_r is not an RNA variant or an allele in cis of them.\n    HGVSError: If the transcript has a CDS (use r_to_c), or the variant is a statement about the transcript."]
    fn r_to_n(&self, _py: Python, var_r: &PyVariant) -> PyResult<PyVariant> {
        let mapper = self.mapper();
        let inner = match expect_rna(&var_r.inner)? {
            OneOrCis::One(v) => mapper.r_to_n(v).map(SequenceVariant::NonCoding),
            OneOrCis::Cis(c) => mapper.r_to_n_cis(c).map(SequenceVariant::CisPhased),
        }
        .map_err(map_hgvs_error)?;
        Ok(PyVariant { inner })
    }

    #[pyo3(signature = (var_c))]
    #[doc = "Respells a coding variant (c.) as an RNA variant (r.): the same positions, bases in lowercase RNA letters (T to u).\n\nAn allele in cis of c. members is respelled member by member.\n\nRaises:\n    ValueError: If var_c is not a coding variant or an allele in cis of them."]
    fn c_to_r(&self, _py: Python, var_c: &PyVariant) -> PyResult<PyVariant> {
        let mapper = self.mapper();
        let inner = match expect_coding_or_cis(&var_c.inner)? {
            OneOrCis::One(v) => mapper.c_to_r(v).map(SequenceVariant::Rna),
            OneOrCis::Cis(c) => mapper.c_to_r_cis(c).map(SequenceVariant::CisPhased),
        }
        .map_err(map_hgvs_error)?;
        Ok(PyVariant { inner })
    }

    #[pyo3(signature = (var_n))]
    #[doc = "Respells a non-coding variant (n.) as an RNA variant (r.): the same positions, bases in lowercase RNA letters (T to u).\n\nAn allele in cis of n. members is respelled member by member.\n\nRaises:\n    ValueError: If var_n is not a non-coding variant or an allele in cis of them."]
    fn n_to_r(&self, _py: Python, var_n: &PyVariant) -> PyResult<PyVariant> {
        let mapper = self.mapper();
        let inner = match expect_noncoding(&var_n.inner)? {
            OneOrCis::One(v) => mapper.n_to_r(v).map(SequenceVariant::Rna),
            OneOrCis::Cis(c) => mapper.n_to_r_cis(c).map(SequenceVariant::CisPhased),
        }
        .map_err(map_hgvs_error)?;
        Ok(PyVariant { inner })
    }

    #[pyo3(signature = (var_p, transcript_ac=None))]
    #[doc = "Back-converts a protein substitution (p.) to a coding variant (c.).\n\nCurrently handles single amino acid substitutions only. When multiple codons\ncould produce the target amino acid, the one requiring the fewest nucleotide\nchanges is chosen. The residue var_p names is read on the protein, so the\nprovider must serve the protein's sequence.\n\nArgs:\n    var_p: The protein Variant to back-convert.\n    transcript_ac: Optional transcript accession (NM_). Required if the DataProvider cannot resolve NP to NM.\n\nReturns:\n    A tuple of (Variant in 'c.' coordinates, is_unique: bool).\n    is_unique is True if the back-conversion is unambiguous.\n\nRaises:\n    ValueError: If var_p is not a protein variant.\n    ReferenceMismatchError: If the residue var_p names is not the protein's.\n    PositionOutsideSequenceError: If var_p names a position past the protein's stop.\n    HGVSError: If back-conversion fails due to unsupported consequence or missing transcript."]
    fn p_to_c(
        &self,
        _py: Python,
        var_p: &PyVariant,
        transcript_ac: Option<String>,
    ) -> PyResult<(PyVariant, bool)> {
        let v = expect_protein(&var_p.inner)?;
        let mapper = self.mapper();
        let (res, is_unique) = mapper
            .p_to_c(v, transcript_ac.as_deref())
            .map_err(map_hgvs_error)?;
        Ok((
            PyVariant {
                inner: SequenceVariant::Coding(res),
            },
            is_unique,
        ))
    }

    #[pyo3(signature = (var))]
    #[doc = "Returns the GA4GH VRS 2.0 object for a variant as a dict: an Allele, a\nCopyNumberCount for a copy-number edit, a CopyNumberChange for a\nduplication with uncertain breakpoints, or a CisPhasedBlock for an allele\nin cis.\n\nA nucleotide variant is projected to its genomic reference; a protein\nvariant stays on its protein. Either is canonicalised (fully justified\nover its region of ambiguity) and rendered with computed identifiers. The\nsequence is identified by its refget accession, from the Refget given at\nconstruction or computed from the whole sequence.\n\nA g. or m. copy-number edit, g.1000_2000copy3, becomes a CopyNumberCount\nover the range with the count as ``copies``; it is not normalised.\n\nA g. or m. duplication with uncertain breakpoints,\ng.(100_200)_(300_400)dup, has no Allele (its bases are unknown) and\nbecomes a CopyNumberChange over the Range bounds with ``copyChange``\n\"gain\" (EFO:0030070, copy number gain); a deletion with uncertain\nbreakpoints stays an Allele, as below.\n\nAn insertion of bases known only by number, c.123_124insN[20] (also the\nolder ins(20); insN[(20_30)] for a range of lengths), becomes an Allele\nover the insertion point whose state is a LengthExpression,\n``{\"type\": \"LengthExpression\", \"length\": 20}``; delinsN[20] does the\nsame over the deleted range. Unknown bases cannot slide, so neither is\nnormalised.\n\nAn allele in cis of nucleotide members is projected to the genome, where the\nmembers that overlap, touch or would shift into each other become one Allele;\nno codon groups them there, so a c. allele and its g. projection get one\nidentifier. One that is a single change there is that change's object;\notherwise it becomes a CisPhasedBlock whose ``members`` are the Alleles, in\norder along the genome, with the shared sequence as ``sequenceReference``. Its\nidentifier (ga4gh:CPB.) does not depend on the order of the members. An allele\nof protein members becomes a block of its members as written.\n\nThe dict's ``type`` says which object was returned.\n\nArgs:\n    var: A g., m., c., n., r. or p. Variant, or an allele in cis of any of\n        them. Protein variants must describe a sequence: frameshifts,\n        extensions and p.? have no allele. A g. or m. deletion with\n        uncertain breakpoints, g.(?_100)_(200_?)del, is rendered as given,\n        with Range bounds ([min, max], null when unbounded) and an empty\n        literal state; it is not normalised.\n\nReturns:\n    A dict in the VRS 2.0 Allele, CopyNumberCount, CopyNumberChange or\n    CisPhasedBlock schema.\n\nRaises:\n    HGVSError: If the variant cannot be resolved against the reference.\n    AmbiguousDescriptionError: If an allele in cis has several readings.\n    ValidationError: If an allele in cis has no reading.\n    UnsupportedOperationError: If an allele in cis cannot be written on the\n        genome as the Alleles it makes, or a member has no Allele."]
    fn to_vrs(&self, py: Python, var: &PyVariant) -> PyResult<Py<PyAny>> {
        let mapper = self.mapper();
        let json_str = mapper
            .to_vrs_variation(&var.inner)
            .map_err(map_hgvs_error)?
            .to_json();
        let json_mod = py.import("json")?;
        Ok(json_mod.call_method1("loads", (json_str,))?.unbind())
    }

    #[pyo3(signature = (var))]
    #[doc = "Returns the GA4GH VRS computed identifier of a variant: ga4gh:VA.<digest>\nfor an Allele, ga4gh:CN.<digest> for a copy-number edit, ga4gh:CX.<digest>\nfor a duplication with uncertain breakpoints, ga4gh:CPB.<digest> for an\nallele in cis.\n\nTwo variants describing the same change on the same sequence have the same\nidentifier; so do two cis alleles of the same members in either order.\n\nArgs:\n    var: A g., m., c., n., r. or p. Variant or allele in cis, as for to_vrs.\n\nRaises:\n    HGVSError: If the variant cannot be resolved against the reference."]
    fn vrs_id(&self, _py: Python, var: &PyVariant) -> PyResult<String> {
        let mapper = self.mapper();
        let variation = mapper
            .to_vrs_variation(&var.inner)
            .map_err(map_hgvs_error)?;
        Ok(variation.id().to_string())
    }

    #[pyo3(signature = (allele, accession=None))]
    #[doc = "Returns the Variant a GA4GH VRS 2.0 Allele, CopyNumberCount,\nCopyNumberChange or CisPhasedBlock names.\n\nAn Allele is written in HGVS on its own sequence, trimmed to the change and\nnormalised (3'-shifted): g. for a nucleotide sequence, p. for a protein.\nLiteral, ReferenceLengthExpression and LengthExpression states are read;\na LengthExpression comes back as written, g.<a>_<a+1>insN[n] over an\nempty location and g.<start+1>_<end>delinsN[n] over a non-empty one\n(N[(min_max)] for a range of lengths). Range bounds are accepted for a\ndeletion, which comes back as g.(a_b)_(c_d)del. A\nCopyNumberCount comes back as g.<start+1>_<end>copyN; its ``copies`` must\nbe an exact count, as HGVS has no syntax for a range of counts. A\nCopyNumberChange comes back as g.<start+1>_<end>dup for a gain\n(``copyChange`` gain, low-level gain or high-level gain, or the EFO codes\nEFO:0030070-72) and del for a loss (loss, low-level loss, high-level loss\nor complete genomic loss, or EFO:0030067-69 and EFO:0020073), Range\nbounds as g.(a_b)_(c_d); other terms have no HGVS form. A\nCisPhasedBlock comes back as the allele in cis of its members, g.[a;b] (or\np.[a;b]), each member read as an Allele; the members must all lie on one\nsequence, which may be stated once on the block as the spec allows.\n\nThe sequence behind the object's refget accession is named by ``accession``\nwhen given, else looked up through the Refget given at construction; the\ndigest is checked against the sequence either way.\n\nArgs:\n    allele: The Allele, CopyNumberCount, CopyNumberChange or CisPhasedBlock as\n        a dict (as to_vrs returns) or a JSON string.\n    accession: The accession of the sequence, when the provider cannot look\n        it up from the refget accession.\n\nRaises:\n    HGVSError: If the object is malformed, unsupported, or does not match the\n        sequence."]
    // Public Python API: a constructor-style name on the mapper, kept as is.
    #[allow(clippy::wrong_self_convention)]
    fn from_vrs(
        &self,
        py: Python,
        allele: &Bound<'_, PyAny>,
        accession: Option<String>,
    ) -> PyResult<PyVariant> {
        let json: String = if let Ok(s) = allele.extract::<String>() {
            s
        } else {
            py.import("json")?
                .call_method1("dumps", (allele,))?
                .extract()?
        };
        let mapper = self.mapper();
        let inner = mapper
            .from_vrs(&json, accession.as_deref())
            .map_err(map_hgvs_error)?;
        Ok(PyVariant { inner })
    }

    #[pyo3(signature = (spdi))]
    #[doc = "Returns the Variant an SPDI string names.\n\n``accession:position:deletion:insertion`` with an interbase position and the\ndeletion given as bases or as a length. The variant is written in HGVS on\nthe accession's own sequence, trimmed to the change and normalised\n(3'-shifted): g. for a nucleotide sequence, p. for a protein.\n\nArgs:\n    spdi: The SPDI string.\n\nRaises:\n    HGVSError: If the string is malformed or the deletion does not match the\n        sequence."]
    // Public Python API: a constructor-style name on the mapper, kept as is.
    #[allow(clippy::wrong_self_convention)]
    fn from_spdi(&self, _py: Python, spdi: &str) -> PyResult<PyVariant> {
        let mapper = self.mapper();
        let inner = mapper.from_spdi(spdi).map_err(map_hgvs_error)?;
        Ok(PyVariant { inner })
    }

    #[pyo3(signature = (var))]
    #[doc = "Normalizes a variant by shifting it to its 3'-most position.\n\nNormalization is performed in the coordinate space of the input variant.\n\nAn allele in cis of nucleotide members, c.[a;b], has its members that meet\nwritten as the one change they make, statements as identities, and each\nmember that meets no other normalised as a single variant is, in order of\nposition; when one member is left, the result is that plain variant.\n\nArgs:\n    var: The Variant object to normalize.\n\nReturns:\n    A new normalized Variant object.\n\nRaises:\n    HGVSError: If normalization fails due to reference sequence boundaries.\n    AmbiguousDescriptionError: If an allele in cis has several readings.\n    ValidationError: If an allele in cis has no reading.\n    UnsupportedOperationError: If an allele in cis has readings weaver does not\n        list, members that share bases across an exon junction, a form it\n        would be written in that does not read as it does, or members that\n        all project to the genome but have no reading there."]
    fn normalize_variant(&self, _py: Python, var: &PyVariant) -> PyResult<PyVariant> {
        let mapper = self.mapper();
        let res = mapper
            .normalize_variant(var.inner.clone())
            .map_err(map_hgvs_error)?;
        Ok(PyVariant { inner: res })
    }

    #[pyo3(signature = (var1, var2, searcher))]
    #[doc = "Determines if two variants are biologically equivalent.\n\nHandles normalization, cross-coordinate mapping (g. vs c.), and gene symbol expansion.\n\nArgs:\n    var1: The first Variant object.\n    var2: The second Variant object.\n    searcher: An object implementing the TranscriptSearch protocol.\n\nReturns:\n    True if the variants are equivalent, False otherwise.\n\nRaises:\n    HGVSError: If mapping, normalization, or data retrieval fails during equivalence checks."]
    fn equivalent(
        &self,
        _py: Python,
        var1: &PyVariant,
        var2: &PyVariant,
        searcher: Py<PyAny>,
    ) -> PyResult<bool> {
        let bridge_searcher = PyTranscriptSearchBridge { searcher };
        let mapper = self.mapper();
        let equiv = ::hgvs_weaver::equivalence::VariantEquivalence::new(&mapper, &bridge_searcher);
        equiv
            .equivalent(&var1.inner, &var2.inner)
            .map_err(map_hgvs_error)
    }

    #[pyo3(signature = (var1, var2, searcher))]
    #[doc = "Determines the granular equivalence level of two variants.\n\nArgs:\n    var1: The first Variant object.\n    var2: The second Variant object.\n    searcher: An object implementing the TranscriptSearch protocol.\n\nReturns:\n    An EquivalenceLevel enum value.\n\nRaises:\n    HGVSError: If mapping, normalization, or data retrieval fails during equivalence checks."]
    fn equivalent_level(
        &self,
        _py: Python,
        var1: &PyVariant,
        var2: &PyVariant,
        searcher: Py<PyAny>,
    ) -> PyResult<PyEquivalenceLevel> {
        let bridge_searcher = PyTranscriptSearchBridge { searcher };
        let mapper = self.mapper();
        let equiv = ::hgvs_weaver::equivalence::VariantEquivalence::new(&mapper, &bridge_searcher);
        let res = equiv
            .equivalent_level(&var1.inner, &var2.inner)
            .map_err(map_hgvs_error)?;
        Ok(res.into())
    }

    #[pyo3(signature = (var, unambiguous = false))]
    #[doc = "Converts a variant to a SPDI string format.\n\nArgs:\n    var: The Variant object to convert.\n    unambiguous: If True, expands the variant range to cover the entire ambiguous region of a repeat or homopolymer. Default is False.\n\nReturns:\n    A string representing the variant in SPDI format.\n\nRaises:\n    HGVSError: If the variant type cannot be converted to SPDI or sequence data is unavailable."]
    fn to_spdi(&self, _py: Python, var: &PyVariant, unambiguous: bool) -> PyResult<String> {
        let mapper = self.mapper();
        mapper
            .to_spdi(&var.inner, unambiguous)
            .map_err(map_hgvs_error)
    }

    #[pyo3(signature = (var))]
    #[doc = "Converts a variant to an unambiguous SPDI string format.\n\nThis format is independent of specific shifting conventions (like 3' or 5' shifting)\nby expanding the variant range to cover the entire ambiguous region of a repeat or homopolymer.\n\nArgs:\n    var: The Variant object to convert.\n\nReturns:\n    A string representing the variant in unambiguous SPDI format.\n\nRaises:\n    HGVSError: If the variant type cannot be converted to SPDI or sequence data is unavailable."]
    fn to_spdi_unambiguous(&self, _py: Python, var: &PyVariant) -> PyResult<String> {
        let mapper = self.mapper();
        mapper.to_spdi(&var.inner, true).map_err(map_hgvs_error)
    }
}

#[pymodule]
fn _weaver(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(parse, m)?)?;
    m.add_class::<PyVariant>()?;
    m.add_class::<PyVariantMapper>()?;
    m.add_class::<PyIdentifierType>()?;
    m.add_class::<PyEquivalenceLevel>()?;
    m.add("HGVSError", m.py().get_type::<HGVSError>())?;
    m.add("ParseError", m.py().get_type::<ParseError>())?;
    m.add("ValidationError", m.py().get_type::<ValidationError>())?;
    m.add("DataProviderError", m.py().get_type::<DataProviderError>())?;
    m.add(
        "UnsupportedOperationError",
        m.py().get_type::<UnsupportedOperationError>(),
    )?;
    m.add("CigarError", m.py().get_type::<CigarError>())?;
    m.add(
        "ReferenceMismatchError",
        m.py().get_type::<ReferenceMismatchError>(),
    )?;
    m.add(
        "PositionOutsideSequenceError",
        m.py().get_type::<PositionOutsideSequenceError>(),
    )?;
    m.add(
        "AmbiguousDescriptionError",
        m.py().get_type::<AmbiguousDescriptionError>(),
    )?;
    m.add(
        "NoProteinStatementError",
        m.py().get_type::<NoProteinStatementError>(),
    )?;
    Ok(())
}
