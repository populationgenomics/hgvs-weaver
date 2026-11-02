//! A `p.` description, or an allele of them, read on its protein: the
//! residues it replaces and what it puts there, so that descriptions are
//! compared by the protein they leave and an allele whose members reach for
//! one residue is refused wherever it is read.

use crate::edits::StopDistance;
use crate::error::HgvsError;
use crate::mapper::VariantMapper;
use crate::structs::{PVariant, SequenceVariant};

impl VariantMapper<'_> {
    /// What a p. description does to the protein `ac`: the edit it makes,
    /// or, for one that names no stretch of residues to replace (`p.?`, an
    /// uncertain range, `p.Met1?`, `p.=`, `p.0`), what it says.
    pub(crate) fn protein_read(&self, vp: &PVariant, ac: &str) -> Result<ProteinRead, HgvsError> {
        use crate::edits::{AaEdit, ExtensionEnd, ProteinStatement};
        let reference = self
            .refs
            .reference(ac, crate::data::IdentifierType::ProteinAccession);
        let whole = reference.whole()?;
        let protein = whole.trim_end_matches('*');
        let len = protein.len();
        if vp.posedit.pos.is_none() {
            return Ok(ProteinRead::Said(match &vp.posedit.edit {
                AaEdit::Statement(ProteinStatement::Unchanged) => {
                    Some(Outcome::Unchanged(Stated::whole()))
                }
                AaEdit::Statement(ProteinStatement::NoProtein { .. }) => Some(Outcome::Leaves {
                    residues: String::new(),
                    open: false,
                    changed_from: 0,
                }),
                _ => None,
            }));
        }
        // An uncertain range, `p.(Lys2_Leu3)del`, does not say which residues change.
        if vp.posedit.pos.as_ref().is_some_and(|p| p.uncertain) {
            return Ok(ProteinRead::Said(None));
        }
        let (start, end) = crate::mapper::residue_range(vp, &reference)?;
        // `p.Met1?`: something happens from this residue on; what, it does not say.
        if let AaEdit::Statement(statement) = &vp.posedit.edit {
            return Ok(ProteinRead::Said(
                (*statement == ProteinStatement::Unknown)
                    .then(|| Outcome::ChangesFrom(start.min(len))),
            ));
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
            // `Met1ext-N`: translation starts N residues upstream, at a Met
            // unless another residue is named, and Met1 becomes the residue
            // written before `ext`, if one is.
            AaEdit::Ext {
                alt,
                end: ExtensionEnd::Upstream { residue, distance },
            } => {
                if start != 0 {
                    return Err(HgvsError::ValidationError(format!(
                        "{vp} extends the protein upstream of a residue other than Met1, \
                         so it has no reading"
                    )));
                }
                let mut added = match residue {
                    Some(r) => crate::utils::residues_1(r)?,
                    None => "M".to_string(),
                };
                while added.len() < distance.get() as usize {
                    added.push('X');
                }
                let met1 = crate::utils::residues_1(alt)?;
                if met1.contains('*') || added.contains('*') {
                    return Err(HgvsError::ValidationError(format!(
                        "{vp} writes a stop in an N-terminal extension, which then ends before \
                         the protein it extends, so it has no reading"
                    )));
                }
                let end = if met1.is_empty() { 0 } else { 1 };
                added.push_str(&met1);
                (0, end, added, false)
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
            AaEdit::Statement(_) | AaEdit::None => return Ok(ProteinRead::Said(None)),
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
        // A lost stop that writes nothing in place of the reference's own
        // residues reads on by at least one: a stop straight after them
        // would leave the protein as it is.
        let alt = match open && protein.get(start..) == Some(alt.as_str()) {
            true => alt + "X",
            false => alt,
        };
        Ok(ProteinRead::Edit(ProteinEdit {
            start,
            end: end.clamp(start, len),
            alt,
            open,
            written_over,
        }))
    }

    /// What a protein allele, or a single p. description, is compared by
    /// when read as the protein it leaves: `None` when a member names no
    /// stretch of residues to replace, so that the allele is compared member
    /// by member.
    ///
    /// # Errors
    ///
    /// [`HgvsError::ValidationError`] when two members reach for the same
    /// residue, a member after a frameshift, a stop or an open end among
    /// them; and what reading a member on the protein raises.
    pub(crate) fn protein_comparand(
        &self,
        var: &SequenceVariant,
        ac: &str,
    ) -> Result<Option<ProteinComparand>, HgvsError> {
        let members: Vec<&PVariant> = match var {
            SequenceVariant::Protein(p) => vec![p],
            SequenceVariant::CisPhased(cis) => cis
                .members
                .iter()
                .filter_map(|m| match m {
                    SequenceVariant::Protein(p) => Some(p),
                    _ => None,
                })
                .collect(),
            _ => return Ok(None),
        };
        let protein = self.reference_protein(ac)?;
        let mut changes: Vec<ProteinEdit> = Vec::new();
        let mut stated: Vec<std::ops::Range<usize>> = Vec::new();
        for vp in members {
            match self.protein_read(vp, ac)? {
                ProteinRead::Said(_) => return Ok(None),
                ProteinRead::Edit(edit) => match edit.statement(&protein) {
                    Some(span) => stated.push(span),
                    None => changes.push(edit),
                },
            }
        }
        changes.sort_by_key(|e| (e.start, e.end));
        let reaches = |e: &ProteinEdit| if e.open { protein.len() } else { e.end };
        for pair in changes.windows(2) {
            if pair[1].start < reaches(&pair[0]) || pair[1].start == pair[0].start {
                let why = match reaches(&pair[0]) == protein.len() {
                    true => "a member after a frameshift, a stop or an open end",
                    false => "two members that reach for one residue",
                };
                return Err(HgvsError::ValidationError(format!(
                    "{var} has {why}, so it has no reading"
                )));
            }
        }
        for span in &stated {
            if changes
                .iter()
                .any(|e| span.start < reaches(e).max(e.start + 1) && e.start < span.end)
            {
                return Err(HgvsError::ValidationError(format!(
                    "{var} states a residue unchanged that another member changes, so it has no reading"
                )));
            }
        }
        Ok(Some(ProteinComparand {
            left: (!changes.is_empty()).then(|| ProteinEdit::applied(&changes, &protein)),
            stated: merged_ranges(stated),
        }))
    }

    /// The protein sequence of `ac`, without a trailing stop.
    pub(crate) fn reference_protein(&self, ac: &str) -> Result<String, HgvsError> {
        let whole = self
            .refs
            .reference(ac, crate::data::IdentifierType::ProteinAccession)
            .whole()?;
        Ok(whole.trim_end_matches('*').to_string())
    }
}

/// What a protein description says of the protein; see `protein_outcome_on`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum Outcome {
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

/// What a p. description does to its protein; see `protein_read`.
pub(crate) enum ProteinRead {
    /// It names no stretch of residues to replace; what it says, if anything.
    Said(Option<Outcome>),
    Edit(ProteinEdit),
}

/// A p. description's replacement of residues `[start, end)` of the
/// protein, before its stop, by `alt` in 1-letter code, `X` where the
/// description does not say; a stop among the new residues has already cut
/// the protein there.
#[derive(Clone, Debug)]
pub(crate) struct ProteinEdit {
    start: usize,
    end: usize,
    alt: String,
    /// More residues follow, their number unknown.
    open: bool,
    /// The residues the description names, where it resolves to a sequence:
    /// its position's (an insertion's two flanks) and those it resolves over.
    written_over: Option<std::ops::Range<usize>>,
}

impl ProteinEdit {
    /// The residues the edit names unchanged, when it leaves `protein` as it is:
    /// `p.Leu3=`, `p.Leu3Leu`, `p.Ter7Ter`, `p.Arg6_Ter7insTer`.
    fn statement(&self, protein: &str) -> Option<std::ops::Range<usize>> {
        let leaves = !self.open && protein.get(self.start..self.end) == Some(self.alt.as_str());
        self.written_over.clone().filter(|_| leaves)
    }

    /// What the edit says of `protein`.
    pub(crate) fn outcome(&self, protein: &str) -> Outcome {
        if let Some(span) = self.statement(protein) {
            return Outcome::Unchanged(Stated::residues(span, protein.len()));
        }
        Self::applied(std::slice::from_ref(self), protein)
    }

    /// The protein `edits`, sorted and apart, leave of `protein`.
    fn applied(edits: &[ProteinEdit], protein: &str) -> Outcome {
        let mut residues = String::with_capacity(protein.len());
        let mut at = 0;
        for e in edits {
            residues.push_str(&protein[at..e.start]);
            residues.push_str(&e.alt);
            at = e.end;
        }
        residues.push_str(&protein[at..]);
        Outcome::Leaves {
            changed_from: first_difference(protein, &residues),
            residues,
            open: edits.iter().any(|e| e.open),
        }
    }
}

/// What a protein allele is compared by when read as a protein: the protein
/// its changes leave (`None` when it changes nothing), and the residues its
/// statements name unchanged, sorted, those that meet made one.
#[derive(Debug)]
pub(crate) struct ProteinComparand {
    left: Option<Outcome>,
    stated: Vec<std::ops::Range<usize>>,
}

impl ProteinComparand {
    /// Whether the two can leave the same protein and state the same residues.
    pub(crate) fn same(self, other: ProteinComparand) -> bool {
        self.stated == other.stated
            && match (self.left, other.left) {
                (None, None) => true,
                (Some(a), Some(b)) => same_protein(Some(a), Some(b)),
                _ => false,
            }
    }
}

/// `ranges`, sorted, with those that overlap or touch made one.
fn merged_ranges(mut ranges: Vec<std::ops::Range<usize>>) -> Vec<std::ops::Range<usize>> {
    ranges.sort_by_key(|r| (r.start, r.end));
    let mut merged: Vec<std::ops::Range<usize>> = Vec::with_capacity(ranges.len());
    for r in ranges {
        match merged.last_mut() {
            Some(last) if r.start <= last.end => last.end = last.end.max(r.end),
            _ => merged.push(r),
        }
    }
    merged
}

pub(crate) use stated::Stated;

mod stated {
    /// The residues a description that leaves the protein unchanged speaks
    /// of: the whole protein, `p.(=)`, or some of its residues, `p.(Leu3=)`.
    /// Every residue through the stop is the whole protein, so the two
    /// spellings of it are one value.
    #[derive(Clone, Debug, PartialEq, Eq)]
    pub(crate) struct Stated {
        /// 0-based and half-open; `None` for the whole protein.
        residues: Option<std::ops::Range<usize>>,
    }

    impl Stated {
        pub(crate) fn whole() -> Self {
            Stated { residues: None }
        }

        /// `residues` of a protein `len` residues long before its stop.
        pub(crate) fn residues(residues: std::ops::Range<usize>, len: usize) -> Self {
            if residues.start == 0 && residues.end > len {
                return Self::whole();
            }
            Stated {
                residues: Some(residues),
            }
        }
    }
}

/// The first index at which two proteins differ, or the shorter length.
pub(crate) fn first_difference(a: &str, b: &str) -> usize {
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
pub(crate) fn same_protein(a: Option<Outcome>, b: Option<Outcome>) -> bool {
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
