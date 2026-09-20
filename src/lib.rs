use ::hgvs_weaver::transform::{transform_variant, StartCodonConvention, VariantTransformSettings};
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
pyo3::create_exception!(_weaver, TranscriptMismatchError, HGVSError);

fn map_hgvs_error(e: HgvsError) -> PyErr {
    match e {
        HgvsError::PestError(msg) => ParseError::new_err(msg),
        HgvsError::ValidationError(msg) => ValidationError::new_err(msg),
        HgvsError::DataProviderError(msg) => DataProviderError::new_err(msg),
        HgvsError::UnsupportedOperation(msg) => UnsupportedOperationError::new_err(msg),
        HgvsError::CigarError(msg) => CigarError::new_err(msg),
        HgvsError::TranscriptMismatch {
            expected,
            found,
            start,
            end,
        } => TranscriptMismatchError::new_err(format!(
            "expected {}, found {} at transcript indices {}..{}",
            expected, found, start, end
        )),
        HgvsError::Other(msg) => HGVSError::new_err(msg),
    }
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

#[pyclass(name = "StartCodonConvention", module = "weaver._weaver")]
#[doc = "Controls how start-codon protein variants are represented.\n\nUsed in VariantTransformSettings to select between keeping the specific\npredicted amino acid change or using the HGVS p.Met1? notation."]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PyStartCodonConvention {
    /// Keep the specific predicted amino acid change (e.g., `p.(Met1Val)`). Default.
    Specific,
    /// Use the HGVS `p.Met1?` notation for any non-silent change at the first codon.
    HgvsQuestion,
}

#[pymethods]
impl PyStartCodonConvention {
    fn __repr__(&self) -> String {
        format!("StartCodonConvention.{:?}", self)
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

impl From<PyStartCodonConvention> for StartCodonConvention {
    fn from(c: PyStartCodonConvention) -> Self {
        match c {
            PyStartCodonConvention::Specific => StartCodonConvention::Specific,
            PyStartCodonConvention::HgvsQuestion => StartCodonConvention::HgvsQuestion,
        }
    }
}

#[pyclass(name = "VariantTransformSettings", module = "weaver._weaver")]
#[doc = "Settings that control how a variant is transformed before formatting or comparison.\n\nCreate with keyword arguments:\n    settings = VariantTransformSettings(start_codon=StartCodonConvention.HgvsQuestion)"]
#[derive(Clone)]
pub struct PyVariantTransformSettings {
    pub inner: VariantTransformSettings,
}

#[pymethods]
impl PyVariantTransformSettings {
    #[new]
    #[pyo3(signature = (start_codon = PyStartCodonConvention::Specific))]
    #[doc = "Creates a new VariantTransformSettings.\n\nArgs:\n    start_codon: Convention for start-codon protein variants. Defaults to Specific."]
    fn new(start_codon: PyStartCodonConvention) -> Self {
        PyVariantTransformSettings {
            inner: VariantTransformSettings {
                start_codon: start_codon.into(),
            },
        }
    }

    #[getter]
    fn start_codon(&self) -> PyStartCodonConvention {
        match self.inner.start_codon {
            StartCodonConvention::Specific => PyStartCodonConvention::Specific,
            StartCodonConvention::HgvsQuestion => PyStartCodonConvention::HgvsQuestion,
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "VariantTransformSettings(start_codon={:?})",
            self.inner.start_codon
        )
    }
}

#[pyclass(name = "Variant", module = "weaver._weaver")]
#[doc = "Represents a parsed HGVS variant.\n\nProvides access to the variant's accession, gene symbol, and coordinate type.\nVariants can be formatted back to HGVS strings or converted to JSON/dict representations."]
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
    #[doc = "The coordinate type of the variant ('g', 'c', 'p', etc.)."]
    fn coordinate_type(&self) -> String {
        self.inner.coordinate_type().to_string()
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

    #[doc = "Validates the variant's reference sequence against the provided DataProvider.\n\nArgs:\n    provider: The data provider instance for sequence retrieval.\n\nReturns:\n    True if the reference sequence matches, False otherwise.\n\nRaises:\n    ValidationError: If transcript sequence is too short or coordinates are out of bounds.\n    DataProviderError: If sequence data cannot be retrieved."]
    fn validate(&self, _py: Python, provider: Py<PyAny>) -> PyResult<bool> {
        let bridge = PyDataProviderBridge { provider };
        VariantMapper::new(&bridge)
            .validate(&self.inner)
            .map_err(map_hgvs_error)
    }

    #[doc = "Returns a new variant with the given transform settings applied.\n\nCurrently transforms protein variants according to the start_codon convention.\nAll other variant types are returned unchanged.\n\nArgs:\n    settings: A VariantTransformSettings object.\n\nReturns:\n    A new Variant with the settings applied."]
    fn transform(&self, settings: &PyVariantTransformSettings) -> PyVariant {
        PyVariant {
            inner: transform_variant(&self.inner, &settings.inner),
        }
    }
}

#[pyfunction]
#[doc = "Parses an HGVS string into a Variant object.\n\nSupported types include genomic (g.), coding cDNA (c.), non-coding (n.),\nmitochondrial (m.), and protein (p.) variants.\n\nArgs:\n    input: The HGVS string to parse.\n\nReturns:\n    A Variant object.\n\nRaises:\n    ParseError: If the HGVS string is malformed or unsupported."]
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
                _ => Err(pyo3::exceptions::PyValueError::new_err(concat!(
                    "Expected a ",
                    $what
                ))),
            }
        }
    };
}
expect_variant!(
    expect_genomic,
    Genomic,
    ::hgvs_weaver::GVariant,
    "genomic variant (g.)"
);
expect_variant!(
    expect_coding,
    Coding,
    ::hgvs_weaver::CVariant,
    "coding variant (c.)"
);
expect_variant!(
    expect_noncoding,
    NonCoding,
    ::hgvs_weaver::structs::NVariant,
    "non-coding variant (n.)"
);
expect_variant!(
    expect_protein,
    Protein,
    ::hgvs_weaver::PVariant,
    "protein variant (p.)"
);
expect_variant!(expect_rna, Rna, ::hgvs_weaver::RVariant, "RNA variant (r.)");

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
    #[doc = "Maps a genomic variant (g.) to a coding cDNA variant (c.) for a specific transcript.\n\nArgs:\n    var_g: The genomic Variant to map.\n    transcript_ac: The accession of the target transcript.\n\nReturns:\n    A new Variant object in 'c.' coordinates.\n\nRaises:\n    ValueError: If var_g is not a genomic variant.\n    HGVSError: If mapping fails due to data or alignment issues."]
    fn g_to_c(&self, _py: Python, var_g: &PyVariant, transcript_ac: String) -> PyResult<PyVariant> {
        let v = expect_genomic(&var_g.inner)?;
        let mapper = self.mapper();
        let res = mapper.g_to_c(v, &transcript_ac).map_err(map_hgvs_error)?;
        Ok(PyVariant {
            inner: SequenceVariant::Coding(res),
        })
    }

    #[pyo3(signature = (var_g, searcher))]
    #[doc = "Maps a genomic variant (g.) to all overlapping transcripts discovered via the searcher.\n\nArgs:\n    var_g: The genomic Variant to map.\n    searcher: An object implementing the TranscriptSearch protocol.\n\nReturns:\n    A list of Variant objects in 'c.' coordinates.\n\nRaises:\n    ValueError: If var_g is not a genomic variant.\n    HGVSError: If mapping fails due to data or alignment issues."]
    fn g_to_c_all(
        &self,
        _py: Python,
        var_g: &PyVariant,
        searcher: Py<PyAny>,
    ) -> PyResult<Vec<PyVariant>> {
        let v = expect_genomic(&var_g.inner)?;
        let mapper = self.mapper();
        let bridge_searcher = PyTranscriptSearchBridge { searcher };
        let res = mapper
            .g_to_c_all(v, &bridge_searcher)
            .map_err(map_hgvs_error)?;
        Ok(res
            .into_iter()
            .map(|v| PyVariant {
                inner: SequenceVariant::Coding(v),
            })
            .collect())
    }

    #[pyo3(signature = (var_c, reference_ac = None))]
    #[doc = "Maps a coding cDNA variant (c.) to a genomic variant (g.).\n\nArgs:\n    var_c: The coding Variant to map.\n    reference_ac: Optional chromosomal accession. If not provided, the primary chromosome for the transcript will be used.\n\nReturns:\n    A new Variant object in 'g.' coordinates.\n\nRaises:\n    ValueError: If var_c is not a coding variant.\n    HGVSError: If mapping fails due to data or alignment issues."]
    fn c_to_g(
        &self,
        _py: Python,
        var_c: &PyVariant,
        reference_ac: Option<String>,
    ) -> PyResult<PyVariant> {
        let v = expect_coding(&var_c.inner)?;
        let mapper = self.mapper();
        let res = mapper
            .c_to_g(v, reference_ac.as_deref())
            .map_err(map_hgvs_error)?;
        Ok(PyVariant {
            inner: SequenceVariant::Genomic(res),
        })
    }

    #[pyo3(signature = (var_n, reference_ac = None))]
    #[doc = "Maps a non-coding cDNA variant (n.) to a genomic variant (g.).\n\nArgs:\n    var_n: The non-coding Variant to map.\n    reference_ac: Optional chromosomal accession.\n\nReturns:\n    A new Variant object in 'g.' coordinates.\n\nRaises:\n    ValueError: If var_n is not a non-coding variant.\n    HGVSError: If mapping fails due to data or alignment issues."]
    fn n_to_g(
        &self,
        _py: Python,
        var_n: &PyVariant,
        reference_ac: Option<String>,
    ) -> PyResult<PyVariant> {
        let v = expect_noncoding(&var_n.inner)?;
        let mapper = self.mapper();
        let res = mapper
            .n_to_g(v, reference_ac.as_deref())
            .map_err(map_hgvs_error)?;
        Ok(PyVariant {
            inner: SequenceVariant::Genomic(res),
        })
    }

    #[pyo3(signature = (var_c, protein_ac=None))]
    #[doc = "Projects a coding cDNA variant (c.) or an RNA variant (r.) to its protein consequence (p.).\n\nAn r. variant is predicted from its c. spelling; a statement about the\ntranscript (r.0, r.spl, r.?, r.=) becomes the matching statement about the\nprotein (p.0, p.?, p.(=)).\n\nArgs:\n    var_c: The coding or RNA Variant to project.\n    protein_ac: Optional protein accession. If not provided, it will be retrieved from the DataProvider.\n\nReturns:\n    A new Variant object in 'p.' coordinates.\n\nRaises:\n    ValueError: If var_c is not a coding or RNA variant.\n    HGVSError: If projection fails due to data retrieval or out-of-bounds coordinates."]
    fn c_to_p(
        &self,
        _py: Python,
        var_c: &PyVariant,
        protein_ac: Option<String>,
    ) -> PyResult<PyVariant> {
        let mapper = self.mapper();
        let res = match &var_c.inner {
            SequenceVariant::Coding(v) => mapper.c_to_p(v, protein_ac.as_deref()),
            SequenceVariant::Rna(r) => mapper.r_to_p(r, protein_ac.as_deref()),
            _ => {
                return Err(pyo3::exceptions::PyValueError::new_err(
                    "Expected a coding variant (c.) or an RNA variant (r.)",
                ))
            }
        }
        .map_err(map_hgvs_error)?;
        Ok(PyVariant {
            inner: SequenceVariant::Protein(res),
        })
    }

    #[pyo3(signature = (var_c, protein_ac=None))]
    #[doc = "Returns the GA4GH VRS 2.0 Allele of a coding variant's protein consequence, on the protein sequence.\n\nThe residues from the first change to the end of the protein become the\nresidues the edited transcript encodes, up to its new stop, then the allele\nis canonicalised. Unlike to_vrs on a p. variant this covers frameshifts,\nextensions and stop losses. The predicted p. description is carried as an\nhgvs.p expression.\n\nArgs:\n    var_c: The coding Variant.\n    protein_ac: Optional protein accession; else the provider's mapping for the transcript.\n\nRaises:\n    ValueError: If var_c is not a coding variant.\n    HGVSError: If the consequence is a statement (p.?, p.Met1?), or the translated CDS is not the protein the provider serves."]
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
    #[doc = "Maps an RNA variant (r.) to a genomic variant (g.).\n\nOnly a change within one exon has a genomic form; one spanning a splice\njunction describes the spliced RNA and raises UnsupportedOperationError.\n\nArgs:\n    var_r: The RNA Variant to map.\n    reference_ac: Optional chromosomal accession.\n\nReturns:\n    A new Variant object in 'g.' coordinates.\n\nRaises:\n    ValueError: If var_r is not an RNA variant.\n    HGVSError: If mapping fails."]
    fn r_to_g(
        &self,
        _py: Python,
        var_r: &PyVariant,
        reference_ac: Option<String>,
    ) -> PyResult<PyVariant> {
        let v = expect_rna(&var_r.inner)?;
        let mapper = self.mapper();
        let res = mapper
            .r_to_g(v, reference_ac.as_deref())
            .map_err(map_hgvs_error)?;
        Ok(PyVariant {
            inner: SequenceVariant::Genomic(res),
        })
    }

    #[pyo3(signature = (var_r))]
    #[doc = "Respells an RNA variant (r.) as the coding variant (c.) it is numbered from.\n\nr. positions on a coding transcript are c. positions; the bases become\nuppercase DNA letters (u to T).\n\nRaises:\n    ValueError: If var_r is not an RNA variant.\n    HGVSError: If the transcript has no CDS, or the variant is a statement about the transcript (r.0, r.spl)."]
    fn r_to_c(&self, _py: Python, var_r: &PyVariant) -> PyResult<PyVariant> {
        let v = expect_rna(&var_r.inner)?;
        let mapper = self.mapper();
        let res = mapper.r_to_c(v).map_err(map_hgvs_error)?;
        Ok(PyVariant {
            inner: SequenceVariant::Coding(res),
        })
    }

    #[pyo3(signature = (var_r))]
    #[doc = "Respells an RNA variant (r.) on a non-coding transcript as the n. variant it is numbered from.\n\nRaises:\n    ValueError: If var_r is not an RNA variant.\n    HGVSError: If the transcript has a CDS (use r_to_c), or the variant is a statement about the transcript."]
    fn r_to_n(&self, _py: Python, var_r: &PyVariant) -> PyResult<PyVariant> {
        let v = expect_rna(&var_r.inner)?;
        let mapper = self.mapper();
        let res = mapper.r_to_n(v).map_err(map_hgvs_error)?;
        Ok(PyVariant {
            inner: SequenceVariant::NonCoding(res),
        })
    }

    #[pyo3(signature = (var_c))]
    #[doc = "Respells a coding variant (c.) as an RNA variant (r.): the same positions, bases in lowercase RNA letters (T to u).\n\nRaises:\n    ValueError: If var_c is not a coding variant."]
    fn c_to_r(&self, _py: Python, var_c: &PyVariant) -> PyResult<PyVariant> {
        let v = expect_coding(&var_c.inner)?;
        let mapper = self.mapper();
        let res = mapper.c_to_r(v).map_err(map_hgvs_error)?;
        Ok(PyVariant {
            inner: SequenceVariant::Rna(res),
        })
    }

    #[pyo3(signature = (var_n))]
    #[doc = "Respells a non-coding variant (n.) as an RNA variant (r.): the same positions, bases in lowercase RNA letters (T to u).\n\nRaises:\n    ValueError: If var_n is not a non-coding variant."]
    fn n_to_r(&self, _py: Python, var_n: &PyVariant) -> PyResult<PyVariant> {
        let v = expect_noncoding(&var_n.inner)?;
        let mapper = self.mapper();
        let res = mapper.n_to_r(v).map_err(map_hgvs_error)?;
        Ok(PyVariant {
            inner: SequenceVariant::Rna(res),
        })
    }

    #[pyo3(signature = (var_p, transcript_ac=None))]
    #[doc = "Back-converts a protein substitution (p.) to a coding variant (c.).\n\nCurrently handles single amino acid substitutions only. When multiple codons\ncould produce the target amino acid, the one requiring the fewest nucleotide\nchanges is chosen.\n\nArgs:\n    var_p: The protein Variant to back-convert.\n    transcript_ac: Optional transcript accession (NM_). Required if the DataProvider cannot resolve NP to NM.\n\nReturns:\n    A tuple of (Variant in 'c.' coordinates, is_unique: bool).\n    is_unique is True if the back-conversion is unambiguous.\n\nRaises:\n    ValueError: If var_p is not a protein variant.\n    HGVSError: If back-conversion fails due to unsupported consequence or missing transcript."]
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
    #[doc = "Returns the GA4GH VRS 2.0 object for a variant as a dict: an Allele, or a\nCopyNumberCount for a copy-number edit.\n\nA nucleotide variant is projected to its genomic reference; a protein\nvariant stays on its protein. Either is canonicalised (fully justified\nover its region of ambiguity) and rendered with computed identifiers. The\nsequence is identified by its refget accession, from the Refget given at\nconstruction or computed from the whole sequence.\n\nA g. or m. copy-number edit, g.1000_2000copy3, becomes a CopyNumberCount\nover the range with the count as ``copies``; it is not normalised. The\ndict's ``type`` says which object was returned.\n\nArgs:\n    var: A g., m., c., n. or p. Variant. Protein variants must describe a\n        sequence: frameshifts, extensions and p.? have no allele. A g. or\n        m. deletion with uncertain breakpoints, g.(?_100)_(200_?)del, is\n        rendered as given, with Range bounds ([min, max], null when\n        unbounded) and an empty literal state; it is not normalised.\n\nReturns:\n    A dict in the VRS 2.0 Allele or CopyNumberCount schema.\n\nRaises:\n    HGVSError: If the variant cannot be resolved against the reference."]
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
    #[doc = "Returns the GA4GH VRS computed identifier of a variant: ga4gh:VA.<digest>\nfor an Allele, ga4gh:CN.<digest> for a copy-number edit.\n\nTwo variants describing the same change on the same sequence have the same\nidentifier.\n\nArgs:\n    var: A g., m., c., n. or p. Variant, as for to_vrs.\n\nRaises:\n    HGVSError: If the variant cannot be resolved against the reference."]
    fn vrs_id(&self, _py: Python, var: &PyVariant) -> PyResult<String> {
        let mapper = self.mapper();
        let variation = mapper
            .to_vrs_variation(&var.inner)
            .map_err(map_hgvs_error)?;
        Ok(variation.id().to_string())
    }

    #[pyo3(signature = (allele, accession=None))]
    #[doc = "Returns the Variant a GA4GH VRS 2.0 Allele or CopyNumberCount names.\n\nAn Allele is written in HGVS on its own sequence, trimmed to the change and\nnormalised (3'-shifted): g. for a nucleotide sequence, p. for a protein.\nLiteral and ReferenceLengthExpression states are read; Range bounds are\naccepted for a deletion, which comes back as g.(a_b)_(c_d)del. A\nCopyNumberCount comes back as g.<start+1>_<end>copyN; its ``copies`` must\nbe an exact count, as HGVS has no syntax for a range of counts.\n\nThe sequence behind the object's refget accession is named by ``accession``\nwhen given, else looked up through the Refget given at construction; the\ndigest is checked against the sequence either way.\n\nArgs:\n    allele: The Allele or CopyNumberCount as a dict (as to_vrs returns) or a\n        JSON string.\n    accession: The accession of the sequence, when the provider cannot look\n        it up from the refget accession.\n\nRaises:\n    HGVSError: If the object is malformed, unsupported, or does not match the\n        sequence."]
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
    #[doc = "Normalizes a variant by shifting it to its 3'-most position.\n\nNormalization is performed in the coordinate space of the input variant.\n\nArgs:\n    var: The Variant object to normalize.\n\nReturns:\n    A new normalized Variant object.\n\nRaises:\n    HGVSError: If normalization fails due to reference sequence boundaries."]
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
    m.add_class::<PyStartCodonConvention>()?;
    m.add_class::<PyVariantTransformSettings>()?;
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
        "TranscriptMismatchError",
        m.py().get_type::<TranscriptMismatchError>(),
    )?;
    Ok(())
}
