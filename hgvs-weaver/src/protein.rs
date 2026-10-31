//! Protein consequence of an edit to a coding sequence.
//!
//! The input is the reference coding sequence, the edit resolved against it
//! (concrete bases over a concrete range), and how long the CDS is. Everything
//! else is derived here: both translations, the first residue that changes,
//! and, crucially, where every stop codon in the alternate protein came from.
//! Because the edit is known at the nucleotide level, a stop codon is the
//! reference's own stop exactly when the edit is in frame and a stop is read
//! at the reference stop's shifted index. Every other stop was created by the
//! edit or read on from the 3'UTR, so nothing has to guess from the edit's
//! kind.

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
    /// The bases that differ, the edit with the prefix and tail the two
    /// sequences share taken off, so however many shared bases an edit
    /// spells reads alike: where they start, where they end in the
    /// reference, and where they end in the edited sequence.
    changed_start: usize,
    changed_ref_end: usize,
    changed_end: usize,
    /// Where the reference's own stop is read in the alternate, when the
    /// change keeps it, in frame whatever it does past the stop.
    original_stop: Option<usize>,
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
        let tail = suffix.min(ref_nt.len().min(edited.len()) - prefix);
        let original_stop = kept_stop(ref_nt, edited, &alt_aa, stop, prefix, suffix);
        let changed_end = edited.len() - tail;
        let changed_ref_end = ref_nt.len() - tail;
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
            changed_start: prefix,
            changed_ref_end,
            changed_end,
            original_stop,
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
        } else if !self.in_frame || self.loses_stop() {
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

    /// The reference has a stop, and the edit does not keep it.
    fn loses_stop(&self) -> bool {
        self.stop < self.ref_aa.len() && self.original_stop.is_none()
    }

    /// Where the alternate protein ends: its first stop.
    fn alt_stop(&self) -> usize {
        self.alt_aa
            .iter()
            .position(|&c| c == '*')
            .unwrap_or(self.alt_aa.len())
    }

    /// One past the last reference residue an edit's description covers when
    /// a new stop ends it: the codons the bases that differ touch, or at least
    /// the first changed one, and never the stop, which the protein record
    /// does not carry as a residue.
    fn ref_end_at_new_stop(&self) -> usize {
        let last = self.changed_ref_end.max(self.changed_start + 1) - 1;
        (last / 3 + 1)
            .min(self.stop)
            .max(self.first_changed + 1)
            .min(self.ref_aa.len())
    }

    /// The first stop of the alternate at or after residue `i`, when its codon
    /// ends within the bases that differ from the reference: a stop the edit
    /// made, not one read on from the 3'UTR.
    fn stop_made_by_edit(&self, i: usize) -> Option<usize> {
        (i..self.alt_aa.len())
            .find(|&j| self.alt_aa[j] == '*')
            .filter(|&j| (j + 1) * 3 <= self.changed_end)
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
        if self.in_frame || self.original_stop.is_some() {
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

/// Where the reference's stop codon, residue `stop`, is read in the edited
/// sequence `alt_nt`, when the change keeps it. The two sequences share a
/// `prefix` and a `suffix`; between them, the bases that differ are aligned
/// with the fewest substitutions, insertions and deletions, and the stop is
/// kept when an alignment that costs no more aligns its codon, base for base,
/// with a stop codon read in frame. The first such stop is where it is read.
/// What the change does past it is no change to the protein, whatever frame
/// it leaves the 3'UTR in. Read from the two sequences, every spelling of one
/// change reads alike. A change whose differing runs are too long to align,
/// more than [`MAX_ALIGNED`] pairs of bases, is read as losing the stop.
fn kept_stop(
    ref_nt: &[u8],
    alt_nt: &[u8],
    alt_aa: &[char],
    stop: usize,
    prefix: usize,
    suffix: usize,
) -> Option<usize> {
    let codon = 3 * stop..3 * stop + 3;
    if codon.end > ref_nt.len() {
        return None;
    }
    let net = alt_nt.len() as i64 - ref_nt.len() as i64;
    let reads_stop_at = |base: i64| {
        (base % 3 == 0 && base >= 0)
            .then_some((base / 3) as usize)
            .filter(|&k| alt_aa.get(k) == Some(&'*'))
    };
    let tail = suffix.min(ref_nt.len().min(alt_nt.len()) - prefix);
    let ref_end = ref_nt.len() - tail;
    if codon.end <= prefix {
        return reads_stop_at(codon.start as i64);
    }
    if codon.start >= ref_end {
        return reads_stop_at(codon.start as i64 + net);
    }
    // Shared bases on either side, as many as the change adds or removes, let
    // the alignment place the change wholly before or after the codon.
    let shift = net.unsigned_abs() as usize;
    let start = prefix.min(codon.start.saturating_sub(shift));
    let end = ref_end
        .max(codon.end + if net < 0 { shift } else { 0 })
        .min(ref_nt.len());
    let upper = |s: &[u8]| s.to_ascii_uppercase();
    let r = upper(&ref_nt[start..end]);
    let a = upper(&alt_nt[start..(end as i64 + net) as usize]);
    if r.len().saturating_mul(a.len()) > MAX_ALIGNED {
        return None;
    }
    let at = codon.start - start;
    // The cost of aligning r[..at] with each prefix of a, and of r[at + 3..]
    // with each suffix, and of the whole.
    let before = edit_rows(&r[..at], &a);
    let after: Vec<usize> = {
        let (rr, ra): (Vec<u8>, Vec<u8>) = (
            r[at + 3..].iter().rev().copied().collect(),
            a.iter().rev().copied().collect(),
        );
        edit_rows(&rr, &ra).into_iter().rev().collect()
    };
    let whole = *edit_rows(&r, &a)
        .last()
        .expect("a row has a column per base and one more");
    (0..(a.len() + 1).saturating_sub(3)).find_map(|j| {
        let substitutions = r[at..at + 3]
            .iter()
            .zip(&a[j..j + 3])
            .filter(|(x, y)| x != y)
            .count();
        (before[j] + substitutions + after[j + 3] == whole)
            .then(|| reads_stop_at((start + j) as i64))
            .flatten()
    })
}

/// The most pairs of bases [`kept_stop`] aligns.
const MAX_ALIGNED: usize = 1 << 24;

/// The last row of the edit-distance table of `r` against `a`: entry `j` is
/// the fewest substitutions, insertions and deletions that make `r` into
/// `a[..j]`.
fn edit_rows(r: &[u8], a: &[u8]) -> Vec<usize> {
    let mut row: Vec<usize> = (0..=a.len()).collect();
    for (i, x) in r.iter().enumerate() {
        let mut diagonal = row[0];
        row[0] = i + 1;
        for (j, y) in a.iter().enumerate() {
            let next = (diagonal + usize::from(x != y))
                .min(row[j] + 1)
                .min(row[j + 1] + 1);
            diagonal = row[j + 1];
            row[j + 1] = next;
        }
    }
    row
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

/// An out-of-frame edit that does not keep the stop: a frameshift, unless the
/// first changed residue is the stop itself (then it reads through: an
/// extension), or a stop formed entirely from the bases that differ ends
/// translation before the new frame reads a reference base (then the codons
/// they touch are replaced up to that stop).
fn describe_frameshift(t: &Reading, out: &Writer<'_>) -> Result<PVariant, HgvsError> {
    let i = t.first_changed;
    if i == t.stop && t.alt_aa.get(i).is_some_and(|&c| c != '*') {
        return out.extension(i);
    }
    if let Some(j) = (i..t.alt_aa.len()).find(|&j| t.alt_aa[j] == '*') {
        if j * 3 >= t.changed_start && (j + 1) * 3 <= t.changed_end {
            if j == i {
                return out.substitution(i, t.ref_aa[i], '*');
            }
            return out.delins(i, t.ref_end_at_new_stop() - 1, &t.alt_aa[i..=j]);
        }
    }
    out.frameshift(i)
}

/// An in-frame edit, or one that keeps the stop and shifts the frame only past
/// it: an extension when it first differs at the stop, whether the stop is
/// lost or read again after residues added before it; a
/// [stop loss](describe_stop_lost) when the stop is not kept; otherwise the
/// change before the stop, written at its 3'-most equivalent residues (the
/// shared tail is trimmed), cut at a stop the edit creates, and classified as
/// nonsense, duplication, insertion, substitution, deletion or delins.
fn describe_in_frame(t: &Reading, out: &Writer<'_>) -> Result<PVariant, HgvsError> {
    let i = t.first_changed;
    if i == t.stop {
        return out.extension(i);
    }
    if t.loses_stop() {
        return describe_stop_lost(t, out);
    }

    // Both proteins end at the reference's stop, or, with none, at the end of
    // the sequence: what follows is read-through, and no part of the change.
    let mut ref_end = t.stop.min(t.ref_aa.len());
    let protein_end = t.original_stop.unwrap_or(t.alt_aa.len());
    let mut alt_end = protein_end;
    while ref_end > i && alt_end > i && t.ref_aa[ref_end - 1] == t.alt_aa[alt_end - 1] {
        ref_end -= 1;
        alt_end -= 1;
    }

    // A stop created by the edit, before the original stop: the change is
    // written up to and including it, over the codons the edit touched.
    if let Some(j) = (i..protein_end).find(|&j| t.alt_aa[j] == '*') {
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

/// An in-frame edit that changes a residue before the stop, after which no
/// stop is read where the reference's was: nonsense, a delins ending at a stop
/// the edit made, or, where the protein reads on into the 3'UTR, a frameshift
/// from the first changed residue.
fn describe_stop_lost(t: &Reading, out: &Writer<'_>) -> Result<PVariant, HgvsError> {
    let i = t.first_changed;
    if t.alt_aa.get(i) == Some(&'*') {
        return out.substitution(i, t.ref_aa[i], '*');
    }
    if let Some(j) = t.stop_made_by_edit(i) {
        return out.delins(i, t.ref_end_at_new_stop() - 1, &t.alt_aa[i..=j]);
    }
    out.frameshift(i)
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

    fn delins(alt: &str) -> NaEdit {
        NaEdit::RefAlt {
            ref_: None,
            alt: Some(alt.into()),
            uncertain: false,
        }
    }
    fn inv() -> NaEdit {
        NaEdit::Inv {
            ref_: None,
            uncertain: false,
        }
    }

    #[test]
    fn an_in_frame_change_that_loses_only_the_stop_is_an_extension_wherever_it_starts() {
        // CGC TAA becomes CGC CAA: Arg6 is rewritten as itself, the stop as Gln.
        assert_eq!(p(CDS, UTR, delins("CGCCAA"), 15, 21), "NP:p.Ter7GlnextTer2");
    }

    #[test]
    fn an_in_frame_change_from_before_the_stop_that_loses_it_is_a_frameshift() {
        // CGC TAA becomes TGG CAA: Arg6 is Trp, the stop Gln, and the 3'UTR is
        // read on to its stop: Trp Gln Pro Ter.
        assert_eq!(p(CDS, UTR, delins("TGGCAA"), 15, 21), "NP:p.Arg6TrpfsTer4");
        // Deleting Arg6 and the stop reads Pro from the 3'UTR, then its stop.
        assert_eq!(p(CDS, UTR, del(), 15, 21), "NP:p.Arg6ProfsTer2");
    }

    #[test]
    fn an_in_frame_change_that_loses_the_stop_and_finds_no_other_is_a_frameshift_of_unknown_length()
    {
        //          M  A  E  Q  P  L  S  E  A  *   then GGC CGG CCG GC: no stop.
        let cds = "ATGGCTGAACAACCACTTTCTGAAGCTTAA";
        // c.26_28delinsTTG: GCT TAA becomes GTT GAA, Val Glu.
        assert_eq!(
            p(cds, "GGCCGGCCGGC", delins("TTG"), 25, 28),
            "NP:p.Ala9ValfsTer?"
        );
        assert_eq!(
            p(cds, "GGCCGGCCGTAA", delins("TTG"), 25, 28),
            "NP:p.Ala9ValfsTer6"
        );
    }

    #[test]
    fn an_in_frame_change_over_the_stop_names_no_residue_past_it() {
        // c.4_*2inv makes GGC TAA from codon 2 on, a stop among the inverted
        // bases; the codons it touches run two bases into the 3'UTR. The
        // deleted range ends at the last residue the protein has.
        assert_eq!(
            p(CDS, "CCGTATAAGTAAGG", inv(), 3, 23),
            "NP:p.Lys2_Arg6delinsGlyTer"
        );
        // c.16_*3delinsTGGCAATAA: CGC TAA CCG becomes TGG CAA TAA.
        assert_eq!(
            p(CDS, UTR, delins("TGGCAATAA"), 15, 24),
            "NP:p.Arg6delinsTrpGlnTer"
        );
    }

    #[test]
    fn an_in_frame_change_that_rewrites_the_stop_as_a_stop_keeps_it() {
        // CGC TAA becomes TGG TAA: only Arg6 changes.
        assert_eq!(p(CDS, UTR, delins("TGGTAA"), 15, 21), "NP:p.Arg6Trp");
        // The same, with the 3'UTR's CCG read-through codon made GGG.
        assert_eq!(p(CDS, UTR, delins("TGGTAAGGG"), 15, 24), "NP:p.Arg6Trp");
        // The same, with the 3'UTR's CCG deleted.
        assert_eq!(p(CDS, UTR, delins("TGGTAA"), 15, 24), "NP:p.Arg6Trp");
    }

    #[test]
    fn residues_duplicated_before_the_stop_are_an_extension() {
        // TAT CGC TAA becomes TAT CGC TAT CGC TAA: Tyr Arg follow Arg6.
        assert_eq!(
            p(
                CDS,
                UTR,
                NaEdit::Dup {
                    ref_: None,
                    uncertain: false
                },
                12,
                18
            ),
            "NP:p.Ter7TyrextTer2"
        );
    }

    #[test]
    fn a_change_on_both_sides_of_a_kept_stop_is_the_change_before_it() {
        // c.[16C>G;*2dup]: CGC TAA CCG becomes GGC TAA CCCG, the frame shifted
        // only past the stop.
        assert_eq!(p(CDS, UTR, delins("GGCTAAC"), 15, 21), "NP:p.Arg6Gly");
        assert_eq!(p(CDS, UTR, delins("GGCTAACC"), 15, 22), "NP:p.Arg6Gly");
    }

    #[test]
    fn residues_inserted_before_a_kept_stop_are_an_extension_however_the_edit_ends() {
        // CAA CCG inserted before TAA, and a C after it in the last two.
        assert_eq!(
            p(CDS, UTR, delins("CAACCGTAA"), 18, 21),
            "NP:p.Ter7GlnextTer2"
        );
        assert_eq!(
            p(CDS, UTR, delins("CAACCGTAAC"), 18, 21),
            "NP:p.Ter7GlnextTer2"
        );
        assert_eq!(
            p(CDS, UTR, delins("CAACCGTAACC"), 18, 22),
            "NP:p.Ter7GlnextTer2"
        );
        // c.18_19insCAA, and CAA before TAA with a G after it.
        assert_eq!(p(CDS, UTR, ins("CAA"), 17, 19), "NP:p.Ter7GlnextTer1");
        assert_eq!(
            p(CDS, UTR, delins("CAATAAG"), 18, 21),
            "NP:p.Ter7GlnextTer1"
        );
    }

    #[test]
    fn a_delins_written_through_a_stop_it_writes_back_is_the_change_before_it() {
        // CGC TAA becomes GGC GCT AAC TAA, also c.16_17delinsGGCGCTAA.
        assert_eq!(
            p(CDS, UTR, delins("GGCGCTAACTAA"), 15, 21),
            "NP:p.Arg6delinsGlyAlaAsn"
        );
        assert_eq!(
            p(CDS, UTR, delins("GGCGCTAA"), 15, 17),
            "NP:p.Arg6delinsGlyAlaAsn"
        );
        assert_eq!(p(CDS, UTR, delins("AGCTAA"), 15, 21), "NP:p.Arg6Ser");
        assert_eq!(p(CDS, UTR, delins("TATTAA"), 12, 21), "NP:p.Arg6del");
        assert_eq!(p(CDS, UTR, delins("GCCTAA"), 9, 21), "NP:p.Tyr5_Arg6del");
    }

    #[test]
    fn a_stop_the_edit_makes_ends_a_delins_over_the_codons_that_differ_however_it_is_spelled() {
        // c.4_12delinsCTGTAA, spelled through the TAT that follows.
        assert_eq!(
            p(CDS, UTR, delins("CTGTAATAT"), 3, 15),
            "NP:p.Lys2_Ala4delinsLeuTer"
        );
        // Out of frame: c.4_12delinsCTGTAAGG, spelled through the TAT.
        assert_eq!(
            p(CDS, UTR, delins("CTGTAAGGTAT"), 3, 15),
            "NP:p.Lys2_Ala4delinsLeuTer"
        );
    }

    #[test]
    fn a_deletion_before_the_stop_that_could_cut_into_its_codon_keeps_it() {
        // M K L A Y *: TAT TAA becomes TAA, which also reads as the deletion
        // of the stop codon's TA and the T before it.
        assert_eq!(p("ATGAAACTGGCCTATTAA", UTR, del(), 12, 15), "NP:p.Tyr5del");
    }

    #[test]
    fn a_lost_stop_reads_alike_however_many_3_prime_utr_bases_the_edit_spells() {
        // c.16_21delinsTGGCAA, spelled through the 3'UTR's own CCG TAG.
        assert_eq!(
            p(CDS, UTR, delins("TGGCAACCGTAG"), 15, 27),
            "NP:p.Arg6TrpfsTer4"
        );
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
    fn an_insertion_before_the_last_residue_that_repeats_it_is_an_extension() {
        // M L F V L C R L *  ; insert TTG TCT (Leu Ser) before the last CTT (Leu).
        // Written 3'-most, Ser Leu follow Leu8: the protein is extended.
        let cds = "ATGCTGTTTGTATTGTGTCGTCTTTAA";
        let out = p(cds, "AGTGCTTTAAG", ins("TTGTCT"), 20, 22);
        assert_eq!(out, "NP:p.Ter9SerextTer2");
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
