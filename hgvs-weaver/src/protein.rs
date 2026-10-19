//! Protein consequence of an edit to a coding sequence.
//!
//! The input is the reference coding sequence, the edit resolved against it
//! (concrete bases over a concrete range), and how long the CDS is. Everything
//! else is derived here: both translations, the first residue that changes,
//! and, crucially, where every stop codon in the alternate protein came from.
//! Because the edit is known at the nucleotide level, a stop codon is the
//! reference's own stop exactly when the edit is in frame, lies wholly before
//! the stop codon, and the stop lands at its shifted index. Every other stop
//! was created by the edit, so nothing has to guess from the edit's kind.

use crate::edits::{AaEdit, ExtensionEnd, ProteinStatement, ResolvedEdit, StopDistance};
use crate::error::HgvsError;
use crate::structs::{AAPosition, AaInterval, PVariant, PosEdit, ProteinPos};
use crate::utils::{aa1_to_aa3, translate};

/// A nucleotide edit placed on a coding sequence.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct CodingChange {
    /// Reference transcript bases from the first base of the CDS to the end
    /// of the transcript, so that read-through past the stop can be translated.
    pub coding: String,

    /// Length of the CDS in bases, stop codon included, as the data source
    /// declares it. This, not the first stop codon in the translation, says
    /// where the protein ends: a selenoprotein has an in-frame TGA that codes
    /// selenocysteine and translates here as `*`.
    pub cds_len: usize,
    /// The CDS runs off the 3' end of the record: `cds_len` bases are all it
    /// carries and no stop codon is declared, so the protein ends at the
    /// first in-frame TAA or TAG, or where the bases do.
    pub end_open: bool,
    /// The edit, with `start` and `end` relative to `coding`.
    pub edit: ResolvedEdit,
    pub protein_ac: String,
    /// Whether translation initiates at the first codon. Residue 1 is then
    /// Met whatever the codon, and a change to the codon is `p.(Met1?)`;
    /// otherwise the codon is read by the standard table, as any other is.
    pub initiates: bool,
}

impl CodingChange {
    /// `edit` on a CDS of `cds_len` bases that starts at the first base of
    /// `coding`, initiates translation there and ends at a stop codon the
    /// record carries. [`with_end_open`](Self::with_end_open) and
    /// [`with_initiates`](Self::with_initiates) set the other cases.
    pub fn new(
        coding: impl Into<String>,
        cds_len: usize,
        edit: ResolvedEdit,
        protein_ac: impl Into<String>,
    ) -> Self {
        CodingChange {
            coding: coding.into(),
            cds_len,
            end_open: false,
            edit,
            protein_ac: protein_ac.into(),
            initiates: true,
        }
    }

    /// The same change on a CDS that runs off the record's 3' end, or not.
    pub fn with_end_open(self, end_open: bool) -> Self {
        CodingChange { end_open, ..self }
    }

    /// The same change with translation initiating at the first codon, or not.
    pub fn with_initiates(self, initiates: bool) -> Self {
        CodingChange { initiates, ..self }
    }
}

/// `seq` with `[start, end)` replaced by `insert`. Indices past the end clamp.
fn splice(seq: &str, start: usize, end: usize, insert: &str) -> String {
    let start = start.min(seq.len());
    let end = end.min(seq.len()).max(start);
    let mut out = String::with_capacity(seq.len() + insert.len());
    out.push_str(&seq[..start]);
    out.push_str(insert);
    out.push_str(&seq[end..]);
    out
}

/// A coding change read at the protein level: both translations, and the
/// indices every rule about stops and frames refers to. Building one is two
/// translations, so a caller that asks whether translation
/// [runs off](Self::runs_off) and then describes the change translates once.
pub(crate) struct Reading {
    change: CodingChange,
    ref_aa: Vec<char>,
    alt_aa: Vec<char>,
    /// The reference stop: the codon the CDS end declares, when that codon
    /// really is a stop; otherwise (a loosely marked CDS end) the first stop
    /// in the translation; otherwise the length.
    stop: usize,
    /// Bases added minus bases removed.
    net: i64,
    in_frame: bool,
    /// The codons the edit touches, for describing a silent change.
    first_codon: usize,
    last_codon: usize,
    /// The first residue that differs, at or after the first touched codon.
    first_changed: usize,
    /// Whether the edit left the reference's own stop codon intact.
    stop_intact: bool,
    /// Whether translation initiates at the first codon and the edited
    /// sequence's first codon differs from the reference's.
    start_changed: bool,
}

impl Reading {
    pub(crate) fn new(mut change: CodingChange) -> Self {
        if change.end_open {
            change.coding.truncate(change.cds_len);
        }
        let ResolvedEdit {
            start, end, alt, ..
        } = &change.edit;
        let (start, end) = (*start, *end);
        let alt_nt = splice(&change.coding, start, end, alt);
        let mut ref_aa: Vec<char> = translate(&change.coding).chars().collect();
        let mut alt_aa: Vec<char> = translate(&alt_nt).chars().collect();
        // An initiator is Met whatever its codon.
        let start_changed = change.initiates
            && match (alt_nt.get(..3), change.coding.get(..3)) {
                (Some(a), Some(r)) => !a.eq_ignore_ascii_case(r),
                (a, r) => a != r,
            };
        if let Some(first) = ref_aa.first_mut().filter(|_| change.initiates) {
            *first = 'M';
            if !start_changed {
                alt_aa[0] = 'M';
            }
        }
        let net = alt.len() as i64 - (end - start) as i64;
        // An open CDS has no declared stop: a TGA in frame is selenocysteine,
        // and any other stop codon ends the protein.
        let stop = match change.end_open {
            true => (0..ref_aa.len())
                .find(|&j| {
                    ref_aa[j] == '*' && !change.coding[3 * j..3 * j + 3].eq_ignore_ascii_case("TGA")
                })
                .unwrap_or(ref_aa.len()),
            false => reference_stop(&ref_aa, change.cds_len),
        };
        // A selenocysteine TGA is U in the reference before its stop, and in
        // the alternate where a TGA the two sequences share reads one in
        // frame: in their common prefix, or in their common suffix after an
        // in-frame change. Read from the sequences, not the edit's
        // coordinates, every spelling of one change reads alike.
        let is_tga = |nt: &str, j: usize| {
            nt.get(3 * j..3 * j + 3)
                .is_some_and(|c| c.eq_ignore_ascii_case("TGA"))
        };
        let (ref_nt, edited) = (change.coding.as_bytes(), alt_nt.as_bytes());
        let prefix = ref_nt
            .iter()
            .zip(edited)
            .take_while(|(a, b)| a == b)
            .count();
        let suffix = ref_nt
            .iter()
            .rev()
            .zip(edited.iter().rev())
            .take_while(|(a, b)| a == b)
            .count();
        let ref_sec = |r: usize| r < stop && ref_aa.get(r) == Some(&'*');
        let alt_sec: Vec<usize> = (0..alt_aa.len())
            .filter(|&j| {
                let r = if (j + 1) * 3 <= prefix {
                    Some(j)
                } else if net % 3 == 0 && j * 3 >= edited.len() - suffix {
                    usize::try_from(j as i64 - net / 3).ok()
                } else {
                    None
                };
                alt_aa[j] == '*' && is_tga(&alt_nt, j) && r.is_some_and(ref_sec)
            })
            .collect();
        let before_stop = stop.min(ref_aa.len());
        for c in ref_aa[..before_stop].iter_mut().filter(|c| **c == '*') {
            *c = 'U';
        }
        for j in alt_sec {
            alt_aa[j] = 'U';
        }
        let first_codon = start / 3;
        let last_codon = if end > start {
            (end - 1) / 3
        } else {
            start / 3
        };
        let mut first_changed = first_codon;
        while first_changed < ref_aa.len()
            && first_changed < alt_aa.len()
            && ref_aa[first_changed] == alt_aa[first_changed]
        {
            first_changed += 1;
        }
        Reading {
            change,
            ref_aa,
            alt_aa,
            stop,
            net,
            in_frame: net % 3 == 0,
            first_codon,
            last_codon,
            first_changed,
            stop_intact: end <= stop * 3,
            start_changed,
        }
    }

    /// The edit changes the first codon: whether translation starts, and
    /// where, is not predicted, so the change is the unknown change from
    /// Met1. Read from the sequences, so an edit that rewrites the codon as
    /// itself, or inserts after it, is not one.
    pub(crate) fn changes_start(&self) -> bool {
        self.start_changed
    }

    /// The change this is a reading of.
    pub(crate) fn change(&self) -> &CodingChange {
        &self.change
    }

    /// Where the edit starts in the coding sequence.
    fn start(&self) -> usize {
        self.change.edit.start
    }

    /// Where the bases the edit inserts end in the edited sequence.
    fn inserted_end(&self) -> usize {
        self.change.edit.start + self.change.edit.alt.len()
    }

    /// Refuses a change whose protein cannot be predicted: one that
    /// [runs off](Self::runs_off) the edited transcript, or one that
    /// [changes the start](Self::changes_start).
    fn check_predictable(&self) -> Result<(), HgvsError> {
        self.check_reads_on()?;
        if self.start_changed {
            return Err(HgvsError::UnsupportedOperation(
                "the edit changes the first codon, so whether and where translation starts, and the \
                 protein it makes, is not predicted"
                    .into(),
            ));
        }
        Ok(())
    }

    /// On a CDS open at its 3' end, refuses a change whose protein needs a
    /// base past the last whole codon the record carries: an edit reaching
    /// past it; an out-of-frame edit, or one losing a stop the record
    /// carries, whose new frame meets no stop within the record; an in-frame
    /// change whose first difference lies past the end
    /// of the reference or the edited protein, where the residue that
    /// follows decides how far HGVS shifts it 3'.
    pub(crate) fn check_within_record(&self, ac: &str) -> Result<(), HgvsError> {
        if !self.change.end_open {
            return Ok(());
        }
        let carried = self.change.cds_len / 3 * 3;
        let edit = &self.change.edit;
        let needs_more = if edit.start >= carried || edit.end > carried {
            true
        } else if self.start_changed || self.first_changed > self.stop {
            false
        } else if !self.in_frame || (self.stop < self.ref_aa.len() && !self.stop_intact) {
            !self.alt_aa[self.first_changed.min(self.alt_aa.len())..].contains(&'*')
        } else {
            self.ref_aa != self.alt_aa
                && (self.first_changed >= self.ref_aa.len()
                    || self.first_changed >= self.alt_aa.len())
        };
        match needs_more {
            false => Ok(()),
            true => Err(HgvsError::ValidationError(format!(
                "{ac}: its CDS is open at the 3' end, so the record carries no stop codon; a protein \
                 consequence is predicted only where it needs no base past the last whole codon the \
                 record carries"
            ))),
        }
    }

    /// Refuses a change that [runs off](Self::runs_off) the edited transcript.
    fn check_reads_on(&self) -> Result<(), HgvsError> {
        if self.runs_off() {
            return Err(HgvsError::UnsupportedOperation(format!(
                "the edited transcript ends before a codon at residue {}: translation runs off its end \
                 with no residue to name, so what the protein becomes is not predicted",
                self.first_changed + 1
            )));
        }
        Ok(())
    }

    /// The edited transcript ends before a codon at the first changed
    /// residue, and the edit lost the reference's stop or shifted the frame:
    /// there is no residue there to name. An in-frame edit of a reference
    /// with no stop (a CDS open at its 3' end) loses none, and an alternate
    /// that is a prefix of it is a deletion.
    pub(crate) fn runs_off(&self) -> bool {
        !self.unchanged()
            && (self.stop < self.ref_aa.len() || !self.in_frame)
            && self.first_changed >= self.alt_aa.len()
    }

    /// The protein is the same: nothing differs, or the difference lies past the stop.
    fn unchanged(&self) -> bool {
        self.ref_aa == self.alt_aa || self.first_changed > self.stop
    }

    /// Where the reference's own stop sits in the alternate, for an in-frame
    /// edit that left it intact.
    fn original_stop_in_alt(&self) -> Option<usize> {
        (self.in_frame && self.stop_intact)
            .then(|| (self.stop as i64 + self.net / 3).max(0) as usize)
    }

    /// Where the alternate protein ends: its first stop.
    fn alt_stop(&self) -> usize {
        self.alt_aa
            .iter()
            .position(|&c| c == '*')
            .unwrap_or(self.alt_aa.len())
    }

    /// The last reference codon an edit's description covers when a new stop
    /// ends it: the touched codons, or at least the first changed one.
    fn ref_end_at_new_stop(&self) -> usize {
        (self.last_codon + 1)
            .max(self.first_changed + 1)
            .min(self.ref_aa.len())
    }

    /// Whether the alternate protein ends at a stop rather than at the end of
    /// the edited transcript.
    pub(crate) fn alternate_stops(&self) -> bool {
        self.alt_stop() < self.alt_aa.len()
    }

    /// Whether the edit adds or removes a multiple of three bases.
    pub(crate) fn in_frame(&self) -> bool {
        self.in_frame
    }

    /// Whether the first residue that differs lies before the reference's stop.
    pub(crate) fn changes_before_stop(&self) -> bool {
        self.first_changed < self.stop
    }

    /// Where translation of the edited sequence ends, as an index into the
    /// reference coding sequence: one past the alternate stop codon, mapped
    /// back over the bases the edit adds or removes; `None` when translation
    /// runs off the end of the edited transcript.
    pub(crate) fn translation_end(&self) -> Option<i64> {
        self.alternate_stops()
            .then(|| 3 * (self.alt_stop() as i64 + 1) - self.net)
    }

    /// [`proteins`](crate::protein::proteins) of the change read.
    pub(crate) fn proteins(&self) -> Result<(String, String), HgvsError> {
        self.check_predictable()?;
        let reference: String = self.ref_aa[..self.stop.min(self.ref_aa.len())]
            .iter()
            .collect();
        let alternate: String = self.alt_aa[..self.alt_stop()].iter().collect();
        Ok((reference, alternate))
    }

    /// [`describe`](crate::protein::describe) of the change read.
    pub(crate) fn describe(&self) -> Result<PVariant, HgvsError> {
        let out = Writer {
            ref_aa: &self.ref_aa,
            alt_aa: &self.alt_aa,
            protein_ac: &self.change.protein_ac,
        };
        self.check_reads_on()?;
        if self.start_changed {
            return out.variant(Some(met1()), AaEdit::Statement(ProteinStatement::Unknown));
        }
        if self.unchanged() {
            // A difference entirely past the stop is no change to the protein;
            // it is not written against a residue the protein does not have.
            // Nor is one over every codon through the stop written residue by
            // residue: it is the whole protein.
            if self.first_codon > self.stop
                || (self.first_codon == 0 && self.last_codon >= self.stop)
            {
                return out.no_change();
            }
            return out.identity(self.first_codon, self.last_codon.min(self.stop));
        }
        if self.in_frame {
            describe_in_frame(self, &out)
        } else {
            describe_frameshift(self, &out)
        }
    }
}

/// Residue 1, the initiator, which reads as Met whatever its codon.
pub(crate) fn met1() -> AaInterval {
    AaInterval {
        start: AAPosition {
            base: ProteinPos(0).to_hgvs(),
            aa: "Met".to_string(),
        },
        end: None,
        uncertain: false,
    }
}

/// The index of the reference's stop in `ref_aa`, the translation of a CDS
/// `cds_len` bases long and what follows it: the codon the CDS end declares,
/// when that codon really is a stop; otherwise (a loosely marked CDS end)
/// the first stop in the translation; otherwise the length.
pub(crate) fn reference_stop(ref_aa: &[char], cds_len: usize) -> usize {
    let declared = cds_len.saturating_sub(1) / 3;
    if ref_aa.get(declared) == Some(&'*') {
        declared
    } else {
        ref_aa
            .iter()
            .position(|&c| c == '*')
            .unwrap_or(ref_aa.len())
    }
}

/// The reference protein and the protein `change` produces, in 1-letter code
/// without their stops, as a protein record writes them (the initiator as Met,
/// a selenocysteine as `U`): what a protein allele is made of.
pub fn proteins(change: &CodingChange) -> Result<(String, String), HgvsError> {
    Reading::new(change.clone()).proteins()
}

/// Describes the protein consequence of `change` in HGVS p. terms.
pub fn describe(change: &CodingChange) -> Result<PVariant, HgvsError> {
    Reading::new(change.clone()).describe()
}

/// An out-of-frame edit: a frameshift, unless a stop formed entirely from the
/// inserted bases ends translation before the new frame reads a reference
/// base (then the touched codons are replaced up to that stop), or the first
/// changed residue is the stop itself (then it reads through: an extension).
fn describe_frameshift(t: &Reading, out: &Writer<'_>) -> Result<PVariant, HgvsError> {
    let i = t.first_changed;
    let inserted_end = t.inserted_end();
    if let Some(j) = (i..t.alt_aa.len()).find(|&j| t.alt_aa[j] == '*') {
        if j * 3 >= t.start() && (j + 1) * 3 <= inserted_end {
            if j == i {
                return out.substitution(i, t.ref_aa[i], '*');
            }
            return out.delins(i, t.ref_end_at_new_stop() - 1, &t.alt_aa[i..=j]);
        }
    }
    if i == t.stop && t.alt_aa.get(i).is_some_and(|&c| c != '*') {
        return out.extension(i);
    }
    out.frameshift(i)
}

/// An in-frame edit: written at its 3'-most equivalent residues (the shared
/// tail is trimmed), cut at a stop the edit creates, and classified as a
/// stop loss, nonsense, duplication, insertion, substitution, deletion or
/// delins.
fn describe_in_frame(t: &Reading, out: &Writer<'_>) -> Result<PVariant, HgvsError> {
    let i = t.first_changed;
    // Stop lost: the edit reaches into the stop codon and the residue there
    // is no longer a stop. An in-frame insertion just before the stop can also
    // make the stop position the first differing residue, but then the stop
    // is intact further along, and the change is an insertion, not a loss.
    let original_stop_in_alt = t.original_stop_in_alt();
    if original_stop_in_alt.is_none() && t.ref_aa.get(i) == Some(&'*') {
        return out.extension(i);
    }

    let mut ref_end = t.ref_aa.len();
    let mut alt_end = t.alt_aa.len();
    while ref_end > i && alt_end > i && t.ref_aa[ref_end - 1] == t.alt_aa[alt_end - 1] {
        ref_end -= 1;
        alt_end -= 1;
    }

    // A stop created by the edit, before the original stop: the change is
    // written up to and including it, over the codons the edit touched.
    let search_end = original_stop_in_alt
        .unwrap_or(t.alt_aa.len())
        .min(t.alt_aa.len());
    if let Some(j) = (i..search_end).find(|&j| t.alt_aa[j] == '*') {
        alt_end = j + 1;
        ref_end = t.ref_end_at_new_stop();
    }

    // Nonsense: the first changed residue is a stop.
    if t.alt_aa.get(i) == Some(&'*') && i < alt_end {
        return out.substitution(i, t.ref_aa[i], '*');
    }

    let del: &[char] = &t.ref_aa[i..ref_end.max(i)];
    let ins: &[char] = &t.alt_aa[i..alt_end.max(i)];
    if del.is_empty() && !ins.is_empty() {
        // Duplication: the inserted residues repeat those just before them.
        if i >= ins.len() && t.ref_aa[i - ins.len()..i] == *ins {
            return out.duplication(i - ins.len(), i - 1);
        }
        // No residue precedes residue 1 to insert after, so residue 1 is
        // replaced by the inserted residues and itself.
        let Some(before) = i.checked_sub(1) else {
            let with_first: Vec<char> = ins.iter().chain(&t.ref_aa[..1]).copied().collect();
            return out.delins(0, 0, &with_first);
        };
        return out.insertion(before, i, ins);
    }
    if del.len() == 1 && ins.len() == 1 {
        return out.substitution(i, del[0], ins[0]);
    }
    if ins.is_empty() {
        return out.deletion(i, ref_end - 1);
    }
    out.delins(i, ref_end - 1, ins)
}

fn aa3(residues: &[char]) -> String {
    residues.iter().map(|c| aa1_to_aa3(*c)).collect()
}

struct Writer<'a> {
    ref_aa: &'a [char],
    alt_aa: &'a [char],
    protein_ac: &'a str,
}

impl Writer<'_> {
    fn residue(&self, i: usize) -> AAPosition {
        AAPosition {
            base: ProteinPos(i as i32).to_hgvs(),
            aa: aa1_to_aa3(self.ref_aa.get(i).copied().unwrap_or('*')).to_string(),
        }
    }

    fn variant(&self, pos: Option<AaInterval>, edit: AaEdit) -> Result<PVariant, HgvsError> {
        Ok(PVariant {
            ac: self.protein_ac.to_string(),
            gene: None,
            posedit: PosEdit {
                pos,
                edit,
                uncertain: false,
                predicted: false,
            },
        })
    }

    fn span(&self, first: usize, last: Option<usize>) -> AaInterval {
        AaInterval {
            start: self.residue(first),
            end: last.map(|l| self.residue(l)),
            uncertain: false,
        }
    }

    /// `p.(=)`: the protein unchanged, the statement the parser reads.
    fn no_change(&self) -> Result<PVariant, HgvsError> {
        self.variant(None, AaEdit::Statement(ProteinStatement::Unchanged))
    }

    /// `p.(Xxx1=)` or `p.(Xxx1_Yyy2=)` over the codons the edit touched, or a
    /// bare `p.(=)` when they lie past the protein.
    fn identity(&self, first_codon: usize, last_codon: usize) -> Result<PVariant, HgvsError> {
        if first_codon >= self.ref_aa.len() {
            return self.no_change();
        }
        let last = last_codon.min(self.ref_aa.len() - 1);
        let end = (last > first_codon).then_some(last);
        self.variant(
            Some(self.span(first_codon, end)),
            AaEdit::Identity { uncertain: false },
        )
    }

    fn substitution(&self, i: usize, from: char, to: char) -> Result<PVariant, HgvsError> {
        let alt = if to == '*' {
            "Ter".to_string()
        } else {
            aa1_to_aa3(to).to_string()
        };
        self.variant(
            Some(self.span(i, None)),
            AaEdit::Subst {
                ref_: aa1_to_aa3(from).to_string(),
                alt,
                uncertain: false,
            },
        )
    }

    fn deletion(&self, first: usize, last: usize) -> Result<PVariant, HgvsError> {
        self.variant(
            Some(self.span(first, (last > first).then_some(last))),
            AaEdit::Del {
                ref_: aa3(&self.ref_aa[first..=last]),
                uncertain: false,
            },
        )
    }

    fn delins(&self, first: usize, last: usize, ins: &[char]) -> Result<PVariant, HgvsError> {
        self.variant(
            Some(self.span(first, (last > first).then_some(last))),
            AaEdit::DelIns {
                ref_: aa3(&self.ref_aa[first..=last]),
                alt: aa3(ins),
                uncertain: false,
            },
        )
    }

    fn insertion(&self, before: usize, after: usize, ins: &[char]) -> Result<PVariant, HgvsError> {
        self.variant(
            Some(self.span(before, Some(after))),
            AaEdit::Ins {
                alt: aa3(ins),
                uncertain: false,
            },
        )
    }

    fn duplication(&self, first: usize, last: usize) -> Result<PVariant, HgvsError> {
        self.variant(
            Some(self.span(first, (last > first).then_some(last))),
            AaEdit::Dup {
                ref_: Some(aa3(&self.ref_aa[first..=last])),
                uncertain: false,
            },
        )
    }

    /// `p.(Xxx#Yyyfs Ter N)`: the first changed residue and the distance to the
    /// first stop the new frame reaches. An immediate stop is a plain nonsense
    /// (or silent) substitution, as HGVS prefers.
    fn frameshift(&self, i: usize) -> Result<PVariant, HgvsError> {
        let ref_c = self.ref_aa.get(i).copied().unwrap_or('*');
        let alt_c = self.alt_aa.get(i).copied().unwrap_or('*');
        if alt_c == '*' {
            if ref_c == '*' {
                return self.variant(
                    Some(self.span(i, None)),
                    AaEdit::Identity { uncertain: false },
                );
            }
            return self.substitution(i, ref_c, '*');
        }
        let stop = match self.alt_aa[i..].iter().position(|&c| c == '*') {
            Some(k) => StopDistance::Residues {
                count: residues_on(k)?,
            },
            None => StopDistance::Unknown,
        };
        self.variant(
            Some(self.span(i, None)),
            AaEdit::Fs {
                alt: aa1_to_aa3(alt_c).to_string(),
                stop: Some(stop),
            },
        )
    }

    /// `p.(Ter#Xxxext*N)`: the stop is read through; N is the distance to the
    /// next stop, or `?` if translation runs off the end.
    fn extension(&self, i: usize) -> Result<PVariant, HgvsError> {
        let alt_c = self.alt_aa.get(i).copied().unwrap_or('*');
        let distance = match self.alt_aa[i + 1..].iter().position(|&c| c == '*') {
            Some(k) => StopDistance::Residues {
                count: residues_on(k)?,
            },
            None => StopDistance::Unknown,
        };
        self.variant(
            Some(AaInterval {
                start: AAPosition {
                    base: ProteinPos(i as i32).to_hgvs(),
                    aa: "Ter".to_string(),
                },
                end: None,
                uncertain: false,
            }),
            AaEdit::Ext {
                alt: aa1_to_aa3(alt_c).to_string(),
                end: ExtensionEnd::Stop { distance },
            },
        )
    }
}

/// The distance to a stop `k` residues after the first one counted.
fn residues_on(k: usize) -> Result<std::num::NonZeroU32, HgvsError> {
    u32::try_from(k)
        .ok()
        .and_then(|k| std::num::NonZeroU32::MIN.checked_add(k))
        .ok_or_else(|| HgvsError::ValidationError(format!("a stop {} residues on", k + 1)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::edits::NaEdit;

    /// Describe an HGVS-range edit against a coding sequence whose CDS is the
    /// whole of `cds` and whose 3' UTR is `utr`.
    fn p(cds: &str, utr: &str, edit: NaEdit, start: usize, end: usize) -> String {
        let coding = format!("{cds}{utr}");
        let resolved = edit
            .resolve_with("NM", start, end, |s, e| {
                Ok(coding[s.min(coding.len())..e.min(coding.len())].to_string())
            })
            .unwrap();
        describe(&CodingChange {
            coding: coding.clone(),
            cds_len: cds.len(),
            end_open: false,
            edit: resolved,
            protein_ac: "NP".into(),
            initiates: true,
        })
        .map(|v| v.to_string())
        .unwrap_or_else(|e| format!("ERR:{e}"))
    }
    fn sub(r: &str, a: &str) -> NaEdit {
        NaEdit::RefAlt {
            ref_: Some(r.into()),
            alt: Some(a.into()),
            uncertain: false,
        }
    }
    fn del() -> NaEdit {
        NaEdit::Del {
            ref_: None,
            uncertain: false,
        }
    }
    fn ins(a: &str) -> NaEdit {
        NaEdit::Ins {
            alt: Some(a.into()),
            uncertain: false,
        }
    }

    // M  K  L  A  Y  R  *     then a UTR that translates to "P" and a stop.
    const CDS: &str = "ATGAAACTGGCCTATCGCTAA";
    const UTR: &str = "CCGTAG";

    #[test]
    fn substitutions_nonsense_and_silence() {
        assert_eq!(p(CDS, UTR, sub("A", "C"), 3, 4), "NP:p.Lys2Gln");
        assert_eq!(p(CDS, UTR, sub("A", "T"), 3, 4), "NP:p.Lys2Ter");
        assert_eq!(p(CDS, UTR, sub("A", "G"), 5, 6), "NP:p.Lys2=");
    }

    #[test]
    fn in_frame_deletion_and_duplication_take_the_3_prime_position() {
        // Deleting the first of two identical codons is written at the second.
        //          M  K  K  L  *
        let two_k = "ATGAAAAAACTGTAA";
        assert_eq!(p(two_k, UTR, del(), 3, 6), "NP:p.Lys3del");
        assert_eq!(
            p(
                two_k,
                UTR,
                NaEdit::Dup {
                    ref_: None,
                    uncertain: false
                },
                3,
                6
            ),
            "NP:p.Lys3dup"
        );
        assert_eq!(p(CDS, UTR, del(), 6, 12), "NP:p.Leu3_Ala4del");
    }

    #[test]
    fn a_stop_inside_the_inserted_bases_is_a_delins_not_a_frameshift() {
        // Replace codons 2..4 (K L A) with L, then TAA: net +? no, 9 out, 6 in
        // is in frame; but 9 out and 8 in is not, and still ends in a new stop.
        assert_eq!(
            p(CDS, UTR, sub("AAACTGGCC", "CTGTAA"), 3, 12),
            "NP:p.Lys2_Ala4delinsLeuTer"
        );
        assert_eq!(
            p(CDS, UTR, sub("AAACTGGCC", "CTGTAAGG"), 3, 12),
            "NP:p.Lys2_Ala4delinsLeuTer"
        );
    }

    #[test]
    fn an_in_frame_deletion_reaching_the_stop_keeps_it_original() {
        // Delete R (codon 5): the stop shifts left but is the same stop.
        assert_eq!(p(CDS, UTR, del(), 15, 18), "NP:p.Arg6del");
    }

    #[test]
    fn frameshift_reports_first_changed_residue_and_distance_to_stop() {
        // Delete one base in codon 2: K L A Y R * P * becomes ... frame shifts.
        let out = p(CDS, UTR, del(), 3, 4);
        assert!(out.starts_with("NP:p.Lys2"), "{out}");
        assert!(out.contains("fsTer"), "{out}");
    }

    #[test]
    fn stop_loss_is_an_extension_to_the_next_stop() {
        // TAA -> CAA (Gln), read through P then the UTR stop: ext*2.
        assert_eq!(p(CDS, UTR, sub("T", "C"), 18, 19), "NP:p.Ter7GlnextTer2");
    }

    #[test]
    fn a_frameshift_in_the_stop_codon_is_an_extension() {
        // Deleting the AA of TAA: the new frame reads TCC GTA TAA (Ser Val *),
        // so the stop becomes Ser and a new stop follows two codons on.
        assert_eq!(p(CDS, "CCGTATAA", del(), 19, 21), "NP:p.Ter7SerextTer2");
        // Two bases earlier the first changed residue is Arg6, so it is a
        // frameshift: CTA ACC GTA TAA.
        assert_eq!(p(CDS, "CCGTATAA", del(), 16, 18), "NP:p.Arg6LeufsTer4");
    }

    #[test]
    fn a_selenocysteine_codon_is_not_the_stop() {
        // M K U(TGA) L *  : the declared CDS runs to the real stop.
        let sec = "ATGAAATGACTGTAA";
        assert_eq!(p(sec, UTR, sub("A", "C"), 3, 4), "NP:p.Lys2Gln");
        assert_eq!(p(sec, UTR, sub("C", "G"), 9, 10), "NP:p.Leu4Val");
        assert_eq!(p(sec, UTR, del(), 9, 12), "NP:p.Leu4del");
    }

    #[test]
    fn an_insertion_before_the_stop_that_repeats_the_last_residue_is_not_an_extension() {
        // M L F V L C R L *  ; insert TTG TCT (Leu Ser) before the last CTT (Leu).
        // The first differing residue is at the stop's index (Leu == Leu), but
        // the stop itself is intact two codons on.
        let cds = "ATGCTGTTTGTATTGTGTCGTCTTTAA";
        let out = p(cds, "AGTGCTTTAAG", ins("TTGTCT"), 20, 22);
        assert_eq!(out, "NP:p.Leu8_Ter9insSerLeu");
    }

    #[test]
    fn a_lowercase_start_codon_written_back_as_itself_is_unchanged() {
        let delins = NaEdit::RefAlt {
            ref_: None,
            alt: Some("ATG".into()),
            uncertain: false,
        };
        assert_eq!(
            p(&CDS.to_ascii_lowercase(), UTR, delins, 0, 3),
            "NP:p.Met1="
        );
    }

    #[test]
    fn an_insertion_inside_the_first_codon_is_an_unknown_change_from_met1() {
        // ATG becomes AGTGTG.
        assert_eq!(p(CDS, UTR, ins("GTG"), 0, 1), "NP:p.Met1?");
    }

    #[test]
    fn a_change_past_the_stop_is_no_change_to_the_protein() {
        // c.22C>G lies in the 3'UTR: the protein has no residue 8 to name.
        assert_eq!(p(CDS, UTR, sub("C", "G"), 21, 22), "NP:p.=");
        // A silent change in the stop codon itself is named at the stop.
        assert_eq!(p(CDS, UTR, sub("A", "G"), 20, 21), "NP:p.Ter7=");
    }
}
