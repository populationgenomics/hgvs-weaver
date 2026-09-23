pub use crate::coords::{
    Anchor, GenomicPos, HgvsGenomicPos, HgvsProteinPos, HgvsTranscriptPos, IntronicOffset,
    ProteinPos, SequenceVariant, TranscriptPos,
};
pub use crate::data::{IdentifierKind, IdentifierType};
pub use crate::edits::{AaEdit, NaEdit};
use crate::error::HgvsError;
use serde::{Deserialize, Serialize};

/// Common trait for all HGVS variants.
pub trait Variant {
    /// Returns the primary accession (e.g., "NM_000051.3").
    fn ac(&self) -> &str;
    /// Returns the optional gene symbol (e.g., "ATM").
    fn gene(&self) -> Option<&str>;
    /// Returns the coordinate type code ("g", "c", "p", etc.).
    fn coordinate_type(&self) -> &str;
    /// Replaces the accession (used when a gene symbol resolves to accessions).
    fn set_ac(&mut self, ac: String);
}

macro_rules! impl_variant {
    ($struct_name:ident, $type_code:expr) => {
        impl Variant for $struct_name {
            fn ac(&self) -> &str {
                &self.ac
            }
            fn gene(&self) -> Option<&str> {
                self.gene.as_deref()
            }
            fn coordinate_type(&self) -> &str {
                $type_code
            }
            fn set_ac(&mut self, ac: String) {
                self.ac = ac;
            }
        }
    };
}

/// A variant on a linear reference (`g.`, `m.`): positions are plain indices
/// on the accession itself. Everything that maps, normalises or renders a
/// genomic variant works on any implementor.
pub trait LinearVariant: Variant + Clone {
    fn posedit(&self) -> &PosEdit<SimpleInterval, NaEdit>;
    fn posedit_mut(&mut self) -> &mut PosEdit<SimpleInterval, NaEdit>;
    fn from_parts(
        ac: String,
        gene: Option<String>,
        posedit: PosEdit<SimpleInterval, NaEdit>,
    ) -> Self;

    /// The same variant as a `g.` variant. Mitochondrial and genomic variants
    /// differ only in the letter they are written with.
    fn to_genomic(&self) -> GVariant {
        GVariant {
            ac: self.ac().to_string(),
            gene: self.gene().map(str::to_string),
            posedit: self.posedit().clone(),
        }
    }
}

macro_rules! impl_linear_variant {
    ($struct_name:ident) => {
        impl LinearVariant for $struct_name {
            fn posedit(&self) -> &PosEdit<SimpleInterval, NaEdit> {
                &self.posedit
            }
            fn posedit_mut(&mut self) -> &mut PosEdit<SimpleInterval, NaEdit> {
                &mut self.posedit
            }
            fn from_parts(
                ac: String,
                gene: Option<String>,
                posedit: PosEdit<SimpleInterval, NaEdit>,
            ) -> Self {
                $struct_name { ac, gene, posedit }
            }
        }
    };
}

/// A variant in transcript space (`c.`, `n.`): positions carry an anchor and
/// may carry an intronic offset, and resolve through a `TranscriptMapper`.
/// The two systems differ only in how an index is written back as a position.
pub trait TranscriptVariant: Variant + Clone {
    /// The anchor a position has when it states none.
    const DEFAULT_ANCHOR: Anchor;

    fn posedit(&self) -> &PosEdit<BaseOffsetInterval, NaEdit>;
    fn posedit_mut(&mut self) -> &mut PosEdit<BaseOffsetInterval, NaEdit>;
    fn from_parts(
        ac: String,
        gene: Option<String>,
        posedit: PosEdit<BaseOffsetInterval, NaEdit>,
    ) -> Self;

    /// The position, in this system's numbering, of a 0-based transcript index.
    fn position_from_index(
        am: &crate::transcript_mapper::TranscriptMapper,
        index: i32,
    ) -> Result<BaseOffsetPosition, HgvsError>;
}

macro_rules! impl_transcript_variant {
    ($struct_name:ident, $anchor:expr, $position_from_index:expr) => {
        impl TranscriptVariant for $struct_name {
            const DEFAULT_ANCHOR: Anchor = $anchor;
            fn posedit(&self) -> &PosEdit<BaseOffsetInterval, NaEdit> {
                &self.posedit
            }
            fn posedit_mut(&mut self) -> &mut PosEdit<BaseOffsetInterval, NaEdit> {
                &mut self.posedit
            }
            fn from_parts(
                ac: String,
                gene: Option<String>,
                posedit: PosEdit<BaseOffsetInterval, NaEdit>,
            ) -> Self {
                $struct_name { ac, gene, posedit }
            }
            fn position_from_index(
                am: &crate::transcript_mapper::TranscriptMapper,
                index: i32,
            ) -> Result<BaseOffsetPosition, HgvsError> {
                let f: fn(
                    &crate::transcript_mapper::TranscriptMapper,
                    i32,
                ) -> Result<BaseOffsetPosition, HgvsError> = $position_from_index;
                f(am, index)
            }
        }
    };
}

/// c. numbering: the CDS anchors (`c.-5`, `c.*3`) come from the transcript model.
fn coding_position_from_index(
    am: &crate::transcript_mapper::TranscriptMapper,
    index: i32,
) -> Result<BaseOffsetPosition, HgvsError> {
    let (c_pos, offset, anchor) = am.n_to_c(TranscriptPos(index))?;
    Ok(BaseOffsetPosition {
        base: c_pos.to_hgvs(),
        offset: (offset.0 != 0).then_some(offset),
        anchor,
        uncertain: false,
    })
}

/// n. numbering: index plus one, always from the transcript start.
fn noncoding_position_from_index(
    _am: &crate::transcript_mapper::TranscriptMapper,
    index: i32,
) -> Result<BaseOffsetPosition, HgvsError> {
    Ok(BaseOffsetPosition {
        base: TranscriptPos(index).to_hgvs(),
        offset: None,
        anchor: Anchor::TranscriptStart,
        uncertain: false,
    })
}

/// Represents a genomic variant (g.).
#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub struct GVariant {
    pub ac: String,
    pub gene: Option<String>,
    pub posedit: PosEdit<SimpleInterval, NaEdit>,
}
impl_variant!(GVariant, "g");
impl_linear_variant!(GVariant);

/// Represents a coding cDNA variant (c.).
#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub struct CVariant {
    pub ac: String,
    pub gene: Option<String>,
    pub posedit: PosEdit<BaseOffsetInterval, NaEdit>,
}
impl_variant!(CVariant, "c");
impl_transcript_variant!(CVariant, Anchor::CdsStart, coding_position_from_index);

/// Represents a protein variant (p.).
#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub struct PVariant {
    pub ac: String,
    pub gene: Option<String>,
    pub posedit: PosEdit<AaInterval, AaEdit>,
}
impl_variant!(PVariant, "p");

/// Represents a mitochondrial variant (m.).
#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub struct MVariant {
    pub ac: String,
    pub gene: Option<String>,
    pub posedit: PosEdit<SimpleInterval, NaEdit>,
}
impl_variant!(MVariant, "m");
impl_linear_variant!(MVariant);

/// Represents a non-coding transcript variant (n.).
#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub struct NVariant {
    pub ac: String,
    pub gene: Option<String>,
    pub posedit: PosEdit<BaseOffsetInterval, NaEdit>,
}
impl_variant!(NVariant, "n");
impl_transcript_variant!(
    NVariant,
    Anchor::TranscriptStart,
    noncoding_position_from_index
);

/// Represents an RNA variant (r.).
#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub struct RVariant {
    pub ac: String,
    pub gene: Option<String>,
    pub posedit: PosEdit<BaseOffsetInterval, NaEdit>,
}
impl_variant!(RVariant, "r");

/// Changes in cis, on one molecule: HGVS `NM_004006.2:c.[145C>T;147C>G]`.
/// Every member is a plain variant in the same coordinate system on the same
/// accession; the members are kept in the order written. The trans form,
/// `c.[145C>T];[147C>G]`, describes two molecules and is not a variant.
#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub struct CisPhasedVariant {
    pub ac: String,
    pub gene: Option<String>,
    pub members: Vec<SequenceVariant>,
}

impl CisPhasedVariant {
    /// A cis allele of `members`, which must be one or more plain variants of
    /// one coordinate system on the accession `ac`.
    pub fn new(
        ac: String,
        gene: Option<String>,
        members: Vec<SequenceVariant>,
    ) -> Result<Self, HgvsError> {
        let Some(first) = members.first() else {
            return Err(HgvsError::ValidationError(
                "A cis allele needs at least one member".into(),
            ));
        };
        for m in &members {
            if matches!(m, SequenceVariant::CisPhased(_)) {
                return Err(HgvsError::ValidationError(
                    "A cis allele's members are plain variants, not cis alleles".into(),
                ));
            }
            if m.coordinate_type() != first.coordinate_type() {
                return Err(HgvsError::ValidationError(format!(
                    "Cis allele members are in one coordinate system, not {} and {}",
                    first.coordinate_type(),
                    m.coordinate_type()
                )));
            }
            if m.ac() != ac {
                return Err(HgvsError::ValidationError(format!(
                    "Cis allele members are on {ac}, not {}",
                    m.ac()
                )));
            }
        }
        Ok(CisPhasedVariant { ac, gene, members })
    }
}

impl Variant for CisPhasedVariant {
    fn ac(&self) -> &str {
        &self.ac
    }
    fn gene(&self) -> Option<&str> {
        self.gene.as_deref()
    }
    /// The members' coordinate system letter; empty with no members.
    fn coordinate_type(&self) -> &str {
        self.members.first().map_or("", |m| m.coordinate_type())
    }
    fn set_ac(&mut self, ac: String) {
        for m in &mut self.members {
            m.set_ac(ac.clone());
        }
        self.ac = ac;
    }
}

/// Combines an interval and an edit (e.g., `123A>G`).
#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub struct PosEdit<I, E> {
    /// The location of the variant.
    pub pos: Option<I>,
    /// The type of change (substitution, deletion, etc.).
    pub edit: E,
    /// Whether the variant is uncertain (indicated by `?`).
    pub uncertain: bool,
    /// Whether the variant is predicted (indicated by `()`).
    pub predicted: bool,
}

impl<I, E> PosEdit<I, E>
where
    I: IntervalSpdi,
    E: EditSpdi,
{
    pub fn to_spdi(
        &self,
        ac: &str,
        refs: &crate::reference::ReferenceStore<'_>,
    ) -> Result<String, HgvsError> {
        let (start, end, spdi_ac) = if let Some(pos) = &self.pos {
            pos.spdi_interval(ac, refs.provider())?
        } else {
            return Err(HgvsError::ValidationError(
                "SPDI requires a position".into(),
            ));
        };

        // SPDI is 0-based. HGVS is 1-based (mostly).
        // IntervalSpdi trait handles the coordinate conversion to 0-based.

        self.edit.to_spdi(&spdi_ac, start, end, refs)
    }
}

pub trait IntervalSpdi {
    /// Returns (start, end, ac) as 0-based integer coordinates and the accession to use for SPDI.
    fn spdi_interval(
        &self,
        ac: &str,
        data_provider: &dyn crate::data::DataProvider,
    ) -> Result<(i32, i32, String), HgvsError>;
}

impl IntervalSpdi for SimpleInterval {
    fn spdi_interval(
        &self,
        ac: &str,
        _data_provider: &dyn crate::data::DataProvider,
    ) -> Result<(i32, i32, String), HgvsError> {
        let start = self.start.base.to_index().0;
        let end = self
            .end
            .as_ref()
            .map_or(start + 1, |e| e.base.to_index().0 + 1);
        Ok((start, end, ac.to_string()))
    }
}

impl IntervalSpdi for BaseOffsetInterval {
    fn spdi_interval(
        &self,
        ac: &str,
        data_provider: &dyn crate::data::DataProvider,
    ) -> Result<(i32, i32, String), HgvsError> {
        // SPDI is expressed on the chromosomal accession, so resolve the
        // transcript interval to genomic coordinates via the transcript model.
        let transcript = data_provider.get_transcript(ac, None)?;
        let reference_ac = transcript.reference_accession.clone();
        let am = crate::transcript_mapper::TranscriptMapper::new(transcript)?;
        let (start, end) = am.interval_to_g(self)?;
        Ok((start.0, end.0, reference_ac))
    }
}

pub trait EditSpdi {
    /// Renders the edit as SPDI over the 0-based half-open genomic range
    /// `[start, end)` of `ac`, reading reference bases from `refs` as needed.
    fn to_spdi(
        &self,
        ac: &str,
        start: i32,
        end: i32,
        refs: &crate::reference::ReferenceStore<'_>,
    ) -> Result<String, HgvsError>;
}

impl EditSpdi for NaEdit {
    fn to_spdi(
        &self,
        ac: &str,
        start: i32,
        end: i32,
        refs: &crate::reference::ReferenceStore<'_>,
    ) -> Result<String, HgvsError> {
        if matches!(
            self,
            NaEdit::None
                | NaEdit::Con { .. }
                | NaEdit::NACopy { .. }
                | NaEdit::Special { .. }
                | NaEdit::InsLength { .. }
                | NaEdit::DelInsLength { .. }
        ) {
            return Err(HgvsError::UnsupportedOperation(format!(
                "Edit type {:?} not yet supported for SPDI",
                self
            )));
        }
        let range = |v: i32, what: &str| {
            usize::try_from(v)
                .map_err(|_| HgvsError::ValidationError(format!("Negative SPDI {} {}", what, v)))
        };
        // SPDI is always on the chromosomal accession.
        let reference = refs.reference(ac, IdentifierType::GenomicAccession);
        let resolved = self.resolve(&reference, range(start, "start")?, range(end, "end")?)?;
        let (p_start, r_strip, a_strip) =
            strip_common_prefix_suffix(resolved.start as i32, &resolved.ref_, &resolved.alt);
        Ok(format!("{}:{}:{}:{}", ac, p_start, r_strip, a_strip))
    }
}

pub fn strip_common_prefix_suffix(
    start: i32,
    ref_seq: &str,
    alt_seq: &str,
) -> (i32, String, String) {
    let mut r_bytes = ref_seq.as_bytes();
    let mut a_bytes = alt_seq.as_bytes();
    let mut p_start = start;

    // Strip prefix
    let mut prefix_len = 0;
    while prefix_len < r_bytes.len()
        && prefix_len < a_bytes.len()
        && r_bytes[prefix_len] == a_bytes[prefix_len]
    {
        prefix_len += 1;
    }
    r_bytes = &r_bytes[prefix_len..];
    a_bytes = &a_bytes[prefix_len..];
    p_start += prefix_len as i32;

    // Strip suffix
    let mut suffix_len = 0;
    while suffix_len < r_bytes.len() && suffix_len < a_bytes.len() {
        let r_idx = r_bytes.len() - 1 - suffix_len;
        let a_idx = a_bytes.len() - 1 - suffix_len;
        if r_bytes[r_idx] == a_bytes[a_idx] {
            suffix_len += 1;
        } else {
            break;
        }
    }
    r_bytes = &r_bytes[..r_bytes.len() - suffix_len];
    a_bytes = &a_bytes[..a_bytes.len() - suffix_len];

    (
        p_start,
        String::from_utf8_lossy(r_bytes).to_string(),
        String::from_utf8_lossy(a_bytes).to_string(),
    )
}

/// An interval spanning simple genomic or mitochondrial coordinates.
#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub struct SimpleInterval {
    pub start: SimplePosition,
    pub end: Option<SimplePosition>,
    pub uncertain: bool,
}

impl SimpleInterval {
    pub fn length(&self) -> Result<i32, HgvsError> {
        match &self.end {
            Some(end) => Ok(end.base.0 - self.start.base.0 + 1),
            None => Ok(1),
        }
    }
}

/// A simple position in genomic or mitochondrial coordinates.
#[derive(Debug, PartialEq, Clone, Copy, Serialize, Deserialize)]
pub struct SimplePosition {
    pub base: HgvsGenomicPos,
    pub end: Option<HgvsGenomicPos>,
    pub uncertain: bool,
}

/// An interval spanning cDNA, n. or r. coordinates.
#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub struct BaseOffsetInterval {
    pub start: BaseOffsetPosition,
    pub end: Option<BaseOffsetPosition>,
    pub uncertain: bool,
}

impl BaseOffsetInterval {
    pub fn length(&self) -> Result<i32, HgvsError> {
        match &self.end {
            Some(end) => {
                if self.start.anchor != end.anchor
                    || self.start.offset.is_some()
                    || end.offset.is_some()
                {
                    return Err(HgvsError::UnsupportedOperation(
                        "Complex interval length calculation not implemented".into(),
                    ));
                }
                Ok(end.base.0 - self.start.base.0 + 1)
            }
            None => Ok(1),
        }
    }
}

/// A position in cDNA, n. or r. coordinates.
#[derive(Debug, PartialEq, Clone, Copy, Serialize, Deserialize)]
pub struct BaseOffsetPosition {
    pub base: HgvsTranscriptPos,
    pub offset: Option<IntronicOffset>,
    pub anchor: Anchor,
    pub uncertain: bool,
}

impl BaseOffsetPosition {
    /// Where the position falls along the transcript, for ordering two
    /// positions of one system: the 5'UTR and CDS (negative and positive
    /// `c.` bases) come before the 3'UTR (`c.*`), then the base, then the
    /// intronic offset (`c.88-1` before `c.88` before `c.88+1`).
    pub fn order_key(&self) -> (u8, i32, i32) {
        let region = match self.anchor {
            Anchor::TranscriptStart => 0,
            Anchor::CdsStart => 1,
            Anchor::CdsEnd => 2,
        };
        (region, self.base.0, self.offset.map_or(0, |o| o.0))
    }
}

/// An interval spanning amino acid positions.
#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub struct AaInterval {
    pub start: AAPosition,
    pub end: Option<AAPosition>,
    pub uncertain: bool,
}

impl AaInterval {
    pub fn length(&self) -> Result<i32, HgvsError> {
        match &self.end {
            Some(end) => Ok(end.base.0 - self.start.base.0 + 1),
            None => Ok(1),
        }
    }
}

/// A position in a protein sequence.
#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub struct AAPosition {
    pub base: HgvsProteinPos,
    pub aa: String,
    pub uncertain: bool,
}
