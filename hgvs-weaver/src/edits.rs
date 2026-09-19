use crate::coords::SequenceVariant;
use crate::error::HgvsError;
use crate::reference::Reference;
use serde::{Deserialize, Serialize};

/// Nucleic acid edits (substitutions, deletions, insertions, etc.).
#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
#[serde(tag = "type")]
pub enum NaEdit {
    /// Substitution, deletion, or insertion represented by reference and alternate sequences.
    RefAlt {
        ref_: Option<String>,
        alt: Option<String>,
        uncertain: bool,
    },
    /// Deletion of a sequence.
    Del {
        ref_: Option<String>,
        uncertain: bool,
    },
    /// Insertion of a sequence.
    Ins {
        alt: Option<String>,
        uncertain: bool,
    },
    /// Duplication of a sequence.
    Dup {
        ref_: Option<String>,
        uncertain: bool,
    },
    /// Inversion of a sequence.
    Inv {
        ref_: Option<String>,
        uncertain: bool,
    },
    /// Conversion to another variant sequence.
    Con {
        con: Box<SequenceVariant>,
        uncertain: bool,
    },
    /// Repeat sequence (e.g., `[10]`).
    Repeat {
        ref_: Option<String>,
        min: i32,
        max: i32,
        uncertain: bool,
    },
    /// Copy number change.
    NACopy { copy: i32, uncertain: bool },
    /// A statement about the whole molecule, with no position: `r.0` (no
    /// transcript), `r.?`, `r.spl` (splicing affected), `r.=`.
    Special { value: String, uncertain: bool },
    /// No change (identity).
    None,
}

/// Amino acid edits (substitutions, frameshifts, extensions, etc.).
#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
#[serde(tag = "type")]
pub enum AaEdit {
    /// Simple substitution.
    Subst {
        ref_: String,
        alt: String,
        uncertain: bool,
    },
    /// Deletion of amino acids.
    Del { ref_: String, uncertain: bool },
    /// Insertion of amino acids.
    Ins { alt: String, uncertain: bool },
    /// Deletion-insertion.
    DelIns {
        ref_: String,
        alt: String,
        uncertain: bool,
    },
    /// Substitution, deletion, or insertion represented by reference and alternate sequences.
    RefAlt {
        ref_: Option<String>,
        alt: Option<String>,
        uncertain: bool,
    },
    /// Frameshift.
    Fs {
        ref_: String,
        alt: String,
        term: Option<String>,
        length: Option<String>,
        uncertain: bool,
    },
    /// Extension of the protein (stop codon loss).
    Ext {
        ref_: String,
        alt: String,
        aaterm: Option<String>,
        length: Option<String>,
        uncertain: bool,
    },
    /// Repeat sequence.
    Repeat {
        ref_: Option<String>,
        min: i32,
        max: i32,
        uncertain: bool,
    },
    /// Duplication.
    Dup {
        ref_: Option<String>,
        uncertain: bool,
    },
    /// Silent variant (no change).
    Identity { uncertain: bool },
    /// Special cases (e.g., `p.0`, `p.?`).
    Special { value: String, uncertain: bool },
    /// Placeholder for no edit.
    None,
}

impl AaEdit {
    pub fn is_identity(&self) -> bool {
        match self {
            AaEdit::Identity { .. } => true,
            AaEdit::Special { value, .. } if value == "=" => true,
            _ => false,
        }
    }
}

/// True for the strings HGVS allows in place of bases: a length (`del3`) or nothing.
pub(crate) fn is_length(s: &str) -> bool {
    s.is_empty() || s.chars().all(|c| c.is_ascii_digit())
}

/// A nucleotide edit resolved against a reference: concrete bases over a
/// concrete 0-based half-open range.
///
/// An insertion has an empty range: `start == end` is the index of the base
/// the new bases go in front of. Every other edit replaces `ref_`, the bases
/// actually at `[start, end)`, with `alt`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ResolvedEdit {
    pub start: usize,
    pub end: usize,
    pub ref_: String,
    pub alt: String,
}

impl ResolvedEdit {
    /// Change in sequence length the edit causes.
    pub fn length_change(&self) -> i64 {
        self.alt.len() as i64 - self.ref_.len() as i64
    }

    /// Whether the edit shifts the reading frame.
    pub fn is_frameshift(&self) -> bool {
        self.length_change() % 3 != 0
    }
}

impl NaEdit {
    /// The reference bases the edit spells out, if it spells any.
    ///
    /// A length (`del3`), an empty string, or an edit whose `ref_` is not the
    /// bases of its own range (a repeat's unit) all count as unstated.
    pub fn stated_ref(&self) -> Option<&str> {
        match self {
            NaEdit::RefAlt { ref_: Some(r), .. }
            | NaEdit::Del { ref_: Some(r), .. }
            | NaEdit::Dup { ref_: Some(r), .. }
                if !is_length(r) =>
            {
                Some(r)
            }
            _ => None,
        }
    }

    /// Resolves the edit over `[start, end)` of `reference`.
    ///
    /// `[start, end)` is the range HGVS writes: the two flanking bases for an
    /// insertion, the edited bases for everything else. Implied bases (an
    /// unstated deletion, a duplication, an inversion) are read from the
    /// reference here, once.
    pub fn resolve(
        &self,
        reference: &Reference<'_, '_>,
        start: usize,
        end: usize,
    ) -> Result<ResolvedEdit, HgvsError> {
        self.resolve_with(start, end, |s, e| reference.slice(s, e))
    }

    /// As [`resolve`](Self::resolve), reading reference bases through `fetch`.
    pub fn resolve_with(
        &self,
        start: usize,
        end: usize,
        fetch: impl Fn(usize, usize) -> Result<String, HgvsError>,
    ) -> Result<ResolvedEdit, HgvsError> {
        let stated_or_fetched = |stated: &Option<String>| -> Result<String, HgvsError> {
            match stated {
                Some(r) if !is_length(r) => Ok(r.clone()),
                _ => fetch(start, end),
            }
        };
        let (ref_, alt) = match self {
            NaEdit::RefAlt { ref_, alt, .. } => {
                let r = stated_or_fetched(ref_)?;
                let a = if ref_.is_none() && alt.is_none() {
                    r.clone()
                } else {
                    alt.clone().unwrap_or_default()
                };
                (r, a)
            }
            NaEdit::Del { ref_, .. } => (stated_or_fetched(ref_)?, String::new()),
            NaEdit::Ins { alt, .. } => {
                let anchor = if end > start { end - 1 } else { start };
                return Ok(ResolvedEdit {
                    start: anchor,
                    end: anchor,
                    ref_: String::new(),
                    alt: alt.clone().unwrap_or_default(),
                });
            }
            NaEdit::Dup { ref_, .. } => {
                let r = stated_or_fetched(ref_)?;
                let a = format!("{r}{r}");
                (r, a)
            }
            NaEdit::Inv { .. } => {
                let r = fetch(start, end)?;
                let a = crate::utils::reverse_complement(&r);
                (r, a)
            }
            NaEdit::Repeat { ref_, max, .. } => {
                // `unit[n]`: every existing copy of the unit, starting at `start`,
                // becomes n copies. The reference is the whole run.
                let unit = match ref_ {
                    Some(u) if !is_length(u) => u.clone(),
                    _ => fetch(start, end)?,
                };
                let mut run_end = start;
                if !unit.is_empty() {
                    while fetch(run_end, run_end + unit.len())? == unit {
                        run_end += unit.len();
                    }
                }
                if run_end == start {
                    // The unit is not there at all; read the stated range.
                    run_end = end;
                }
                return Ok(ResolvedEdit {
                    start,
                    end: run_end,
                    ref_: fetch(start, run_end)?,
                    alt: unit.repeat((*max).max(0) as usize),
                });
            }
            NaEdit::None => {
                let r = fetch(start, end)?;
                (r.clone(), r)
            }
            NaEdit::Con { .. } | NaEdit::NACopy { .. } | NaEdit::Special { .. } => {
                return Err(HgvsError::UnsupportedOperation(format!(
                    "Edit type {:?} cannot be resolved to reference and alternate bases",
                    self
                )))
            }
        };
        Ok(ResolvedEdit {
            start,
            end,
            ref_,
            alt,
        })
    }

    /// Applies a function to all sequence strings within the edit.
    pub fn map_sequence<F>(self, f: F) -> NaEdit
    where
        F: Fn(&str) -> String + Copy,
    {
        match self {
            NaEdit::RefAlt {
                ref_,
                alt,
                uncertain,
            } => NaEdit::RefAlt {
                ref_: ref_.map(|s| f(&s)),
                alt: alt.map(|s| f(&s)),
                uncertain,
            },
            NaEdit::Del { ref_, uncertain } => NaEdit::Del {
                ref_: ref_.map(|s| f(&s)),
                uncertain,
            },
            NaEdit::Ins { alt, uncertain } => NaEdit::Ins {
                alt: alt.map(|s| f(&s)),
                uncertain,
            },
            NaEdit::Dup { ref_, uncertain } => NaEdit::Dup {
                ref_: ref_.map(|s| f(&s)),
                uncertain,
            },
            NaEdit::Inv { ref_, uncertain } => NaEdit::Inv {
                ref_: ref_.map(|s| f(&s)),
                uncertain,
            },
            NaEdit::Repeat {
                ref_,
                min,
                max,
                uncertain,
            } => NaEdit::Repeat {
                ref_: ref_.map(|s| f(&s)),
                min,
                max,
                uncertain,
            },
            _ => self,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn resolve(edit: NaEdit, start: usize, end: usize) -> ResolvedEdit {
        //          0123456789
        let seq = "TTCAGCAGTT";
        edit.resolve_with(start, end, |s, e| {
            Ok(seq[s.min(seq.len())..e.min(seq.len())].to_string())
        })
        .unwrap()
    }
    fn r(start: usize, end: usize, ref_: &str, alt: &str) -> ResolvedEdit {
        ResolvedEdit {
            start,
            end,
            ref_: ref_.into(),
            alt: alt.into(),
        }
    }

    #[test]
    fn implied_bases_are_read_from_the_reference_once() {
        assert_eq!(
            resolve(
                NaEdit::Del {
                    ref_: None,
                    uncertain: false
                },
                2,
                5
            ),
            r(2, 5, "CAG", "")
        );
        assert_eq!(
            resolve(
                NaEdit::Del {
                    ref_: Some("3".into()),
                    uncertain: false
                },
                2,
                5
            ),
            r(2, 5, "CAG", "")
        );
        assert_eq!(
            resolve(
                NaEdit::Del {
                    ref_: Some("CAG".into()),
                    uncertain: false
                },
                2,
                5
            ),
            r(2, 5, "CAG", "")
        );
        assert_eq!(
            resolve(
                NaEdit::Dup {
                    ref_: None,
                    uncertain: false
                },
                2,
                5
            ),
            r(2, 5, "CAG", "CAGCAG")
        );
        assert_eq!(
            resolve(
                NaEdit::Inv {
                    ref_: None,
                    uncertain: false
                },
                2,
                5
            ),
            r(2, 5, "CAG", "CTG")
        );
        // A repeat covers every existing copy of its unit: CAG twice at 2..8.
        let cag4 = NaEdit::Repeat {
            ref_: Some("CAG".into()),
            min: 4,
            max: 4,
            uncertain: false,
        };
        assert_eq!(resolve(cag4, 2, 3), r(2, 8, "CAGCAG", "CAGCAGCAGCAG"));
        let unstated1 = NaEdit::Repeat {
            ref_: None,
            min: 1,
            max: 1,
            uncertain: false,
        };
        assert_eq!(resolve(unstated1, 2, 5), r(2, 8, "CAGCAG", "CAG"));
        assert_eq!(resolve(NaEdit::None, 2, 5), r(2, 5, "CAG", "CAG"));
    }

    #[test]
    fn substitutions_and_delins_keep_their_stated_bases() {
        let sub = NaEdit::RefAlt {
            ref_: Some("C".into()),
            alt: Some("T".into()),
            uncertain: false,
        };
        assert_eq!(resolve(sub, 2, 3), r(2, 3, "C", "T"));
        // An unstated delins reference ("" from the parser) is fetched.
        let delins = NaEdit::RefAlt {
            ref_: Some("".into()),
            alt: Some("A".into()),
            uncertain: false,
        };
        assert_eq!(resolve(delins, 2, 5), r(2, 5, "CAG", "A"));
        let identity = NaEdit::RefAlt {
            ref_: None,
            alt: None,
            uncertain: false,
        };
        assert_eq!(resolve(identity, 2, 5), r(2, 5, "CAG", "CAG"));
    }

    #[test]
    fn insertion_resolves_to_the_empty_range_at_its_second_flank() {
        let ins = NaEdit::Ins {
            alt: Some("GG".into()),
            uncertain: false,
        };
        assert_eq!(resolve(ins.clone(), 4, 6), r(5, 5, "", "GG"));
        // Already-placed (empty) ranges are left where they are.
        assert_eq!(resolve(ins, 5, 5), r(5, 5, "", "GG"));
        assert!(r(5, 5, "", "GG").is_frameshift());
        assert!(!r(2, 5, "CAG", "").is_frameshift());
    }

    #[test]
    fn stated_ref_ignores_lengths_units_and_empties() {
        assert_eq!(
            NaEdit::Del {
                ref_: Some("CAG".into()),
                uncertain: false
            }
            .stated_ref(),
            Some("CAG")
        );
        assert_eq!(
            NaEdit::Del {
                ref_: Some("3".into()),
                uncertain: false
            }
            .stated_ref(),
            None
        );
        assert_eq!(
            NaEdit::RefAlt {
                ref_: Some("".into()),
                alt: Some("A".into()),
                uncertain: false
            }
            .stated_ref(),
            None
        );
        assert_eq!(
            NaEdit::Repeat {
                ref_: Some("CAG".into()),
                min: 1,
                max: 1,
                uncertain: false
            }
            .stated_ref(),
            None
        );
    }
}

impl AaEdit {
    /// Resolves the edit over `[start, end)` of a protein `reference` (0-based
    /// residue indices) to the residues it removes and the residues it puts
    /// there, in 1-letter code. The protein counterpart of [`NaEdit::resolve`].
    pub fn resolve(
        &self,
        reference: &Reference<'_, '_>,
        start: usize,
        end: usize,
    ) -> Result<ResolvedEdit, HgvsError> {
        self.resolve_with(start, end, |s, e| reference.slice(s, e))
    }

    /// [`AaEdit::resolve`] with `fetch(start, end)` supplying reference residues.
    ///
    /// Edits that describe a consequence rather than a sequence (frameshift,
    /// extension, `p.?`, `p.0`, `p.Xxx1?`) cannot be resolved and are
    /// unsupported. Stated residues are not consulted; the reference is what
    /// the sequence holds.
    pub fn resolve_with(
        &self,
        start: usize,
        end: usize,
        fetch: impl Fn(usize, usize) -> Result<String, HgvsError>,
    ) -> Result<ResolvedEdit, HgvsError> {
        use crate::utils::residues_1;
        let unsupported = || {
            Err(HgvsError::UnsupportedOperation(format!(
                "Protein edit {self:?} describes a consequence, not a sequence"
            )))
        };
        let (ref_, alt) = match self {
            AaEdit::Subst { alt, .. } => {
                if alt == "?" {
                    return unsupported();
                }
                (fetch(start, end)?, residues_1(alt)?)
            }
            AaEdit::Del { .. } => (fetch(start, end)?, String::new()),
            AaEdit::Ins { alt, .. } => {
                let anchor = if end > start { end - 1 } else { start };
                return Ok(ResolvedEdit {
                    start: anchor,
                    end: anchor,
                    ref_: String::new(),
                    alt: residues_1(alt)?,
                });
            }
            AaEdit::DelIns { alt, .. } => (fetch(start, end)?, residues_1(alt)?),
            AaEdit::RefAlt { alt, .. } => {
                let r = fetch(start, end)?;
                let a = match alt {
                    Some(a) => residues_1(a)?,
                    None => r.clone(),
                };
                (r, a)
            }
            AaEdit::Dup { .. } => {
                let r = fetch(start, end)?;
                let a = format!("{r}{r}");
                (r, a)
            }
            AaEdit::Repeat { ref_, max, .. } => {
                // `unit[n]`: every existing copy of the unit, starting at
                // `start`, becomes n copies. The reference is the whole run.
                let unit = match ref_ {
                    Some(u) if !u.is_empty() && !u.chars().all(|c| c.is_ascii_digit()) => {
                        residues_1(u)?
                    }
                    _ => fetch(start, end)?,
                };
                let mut run_end = start;
                if !unit.is_empty() {
                    while fetch(run_end, run_end + unit.len())? == unit {
                        run_end += unit.len();
                    }
                }
                if run_end == start {
                    run_end = end;
                }
                return Ok(ResolvedEdit {
                    start,
                    end: run_end,
                    ref_: fetch(start, run_end)?,
                    alt: unit.repeat((*max).max(0) as usize),
                });
            }
            AaEdit::Identity { .. } => {
                let r = fetch(start, end)?;
                (r.clone(), r)
            }
            AaEdit::Fs { .. } | AaEdit::Ext { .. } | AaEdit::Special { .. } | AaEdit::None => {
                return unsupported()
            }
        };
        Ok(ResolvedEdit {
            start,
            end,
            ref_,
            alt,
        })
    }
}
