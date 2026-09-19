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

use crate::edits::{AaEdit, ResolvedEdit};
use crate::error::HgvsError;
use crate::structs::{AAPosition, AaInterval, PVariant, PosEdit, ProteinPos};
use crate::utils::{aa1_to_aa3, translate};

/// A nucleotide edit placed on a coding sequence.
pub struct CodingChange<'a> {
    /// Reference transcript bases from the first base of the CDS to the end
    /// of the transcript, so that read-through past the stop can be translated.
    pub coding: &'a str,

    /// Length of the CDS in bases, stop codon included, as the data source
    /// declares it. This, not the first stop codon in the translation, says
    /// where the protein ends: a selenoprotein has an in-frame TGA that codes
    /// selenocysteine and translates here as `*`.
    pub cds_len: usize,
    /// The edit, with `start` and `end` relative to `coding`.
    pub edit: ResolvedEdit,
    pub protein_ac: String,
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

/// The reference protein and the protein `change` produces, in 1-letter code
/// without their stops: what a protein allele is made of. The stops are read
/// as `describe` reads them: the declared CDS end for the reference; for the
/// alternate the first stop the edit creates or the read-through reaches,
/// where a `*` the reference already has in frame (a selenocysteine TGA) is
/// not a stop.
pub fn proteins(change: &CodingChange<'_>) -> Result<(String, String), HgvsError> {
    let ResolvedEdit {
        start, end, alt, ..
    } = &change.edit;
    let (start, end, alt) = (*start, *end, alt.as_str());
    let alt_nt = splice(change.coding, start, end, alt);
    let ref_aa: Vec<char> = translate(change.coding).chars().collect();
    let alt_aa: Vec<char> = translate(&alt_nt).chars().collect();
    let net = alt.len() as i64 - (end - start) as i64;
    let in_frame = net % 3 == 0;
    let declared = change.cds_len.saturating_sub(1) / 3;
    let stop = if ref_aa.get(declared) == Some(&'*') {
        declared
    } else {
        ref_aa
            .iter()
            .position(|&c| c == '*')
            .unwrap_or(ref_aa.len())
    };
    let reference: String = ref_aa[..stop.min(ref_aa.len())].iter().collect();

    // The first residue that differs; everything before it is the reference's.
    let mut i = start / 3;
    while i < ref_aa.len() && i < alt_aa.len() && ref_aa[i] == alt_aa[i] {
        i += 1;
    }
    // The reference residue an alternate position corresponds to, if the
    // frame is kept; a `*` there before the stop is a selenocysteine.
    let edit_alt_end = (start + alt.len()).div_ceil(3);
    let is_sec = |j: usize| {
        let ref_j = if j < i {
            j
        } else if !in_frame {
            return false;
        } else if j >= edit_alt_end {
            (j as i64 - net / 3).max(0) as usize
        } else {
            j
        };
        ref_j < stop && ref_aa.get(ref_j) == Some(&'*')
    };
    let alt_stop = (0..alt_aa.len())
        .find(|&j| alt_aa[j] == '*' && !is_sec(j))
        .unwrap_or(alt_aa.len());
    let alternate: String = alt_aa[..alt_stop].iter().collect();
    Ok((reference, alternate))
}

/// Describes the protein consequence of `change` in HGVS p. terms.
pub fn describe(change: &CodingChange<'_>) -> Result<PVariant, HgvsError> {
    let ResolvedEdit {
        start, end, alt, ..
    } = &change.edit;
    let (start, end, alt) = (*start, *end, alt.as_str());
    let alt_nt = splice(change.coding, start, end, alt);
    let ref_aa: Vec<char> = translate(change.coding).chars().collect();
    let alt_aa: Vec<char> = translate(&alt_nt).chars().collect();
    let net = alt.len() as i64 - (end - start) as i64;
    let in_frame = net % 3 == 0;
    // The reference stop is the codon the CDS end declares, when that codon
    // really is a stop; otherwise (a loosely marked CDS end) the first stop in
    // the translation.
    let declared = change.cds_len.saturating_sub(1) / 3;
    let stop = if ref_aa.get(declared) == Some(&'*') {
        declared
    } else {
        ref_aa
            .iter()
            .position(|&c| c == '*')
            .unwrap_or(ref_aa.len())
    };
    let first_codon = start / 3;
    // Codons the edit touches, for describing a silent change.
    let last_codon = if end > start {
        (end - 1) / 3
    } else {
        start / 3
    };

    let out = Writer {
        ref_aa: &ref_aa,
        alt_aa: &alt_aa,
        protein_ac: &change.protein_ac,
    };

    if ref_aa == alt_aa {
        return out.identity(first_codon, last_codon);
    }

    // The first residue that differs, at or after the first codon the edit touches.
    let mut i = first_codon;
    while i < ref_aa.len() && i < alt_aa.len() && ref_aa[i] == alt_aa[i] {
        i += 1;
    }
    if i > stop {
        // The protein is unchanged; the difference lies past its stop.
        return out.identity(first_codon, last_codon);
    }

    if !in_frame {
        // A stop codon formed entirely from inserted bases ends translation
        // before the new frame reads a single reference base, so there is no
        // frameshift to report: the touched codons are replaced up to that stop.
        let inserted_end = start + alt.len();
        if let Some(j) = (i..alt_aa.len()).find(|&j| alt_aa[j] == '*') {
            if j * 3 >= start && (j + 1) * 3 <= inserted_end {
                let ref_end = (last_codon + 1).max(i + 1).min(ref_aa.len());
                if j == i {
                    return out.substitution(i, ref_aa[i], '*');
                }
                return out.delins(i, ref_end - 1, &alt_aa[i..=j]);
            }
        }
        // A frameshift whose first changed residue is the stop itself reads
        // through it: HGVS writes that as an extension, not a frameshift.
        if i == stop && alt_aa.get(i).is_some_and(|&c| c != '*') {
            return out.extension(i);
        }
        return out.frameshift(i);
    }

    // --- In frame ---

    // Where the reference's own stop codon sits in the alternate protein, if
    // the edit left it intact.
    let original_stop_in_alt = (end <= stop * 3).then(|| (stop as i64 + net / 3).max(0) as usize);

    // Stop lost: the edit reaches into the stop codon and the residue there
    // is no longer a stop. An in-frame insertion just before the stop can also
    // make the stop position the first differing residue, but then the stop
    // is intact further along, and the change is an insertion, not a loss.
    if original_stop_in_alt.is_none() && ref_aa.get(i) == Some(&'*') {
        return out.extension(i);
    }

    // Trim the shared tail. This is also the protein-level 3' rule: an
    // in-frame change is written at its 3'-most equivalent residues.
    let mut ref_end = ref_aa.len();
    let mut alt_end = alt_aa.len();
    while ref_end > i && alt_end > i && ref_aa[ref_end - 1] == alt_aa[alt_end - 1] {
        ref_end -= 1;
        alt_end -= 1;
    }

    // A stop created by the edit, before the original stop: the change is
    // written up to and including it, over the codons the edit touched. A
    // `*` the reference already has at the corresponding codon (a
    // selenocysteine TGA) was not created by the edit.
    let search_end = original_stop_in_alt
        .unwrap_or(alt_aa.len())
        .min(alt_aa.len());
    let edit_alt_end = (start + alt.len()).div_ceil(3);
    let already_in_ref = |j: usize| {
        let ref_j = if j >= edit_alt_end {
            (j as i64 - net / 3).max(0) as usize
        } else {
            j
        };
        ref_aa.get(ref_j) == Some(&'*')
    };
    if let Some(j) = (i..search_end).find(|&j| alt_aa[j] == '*' && !already_in_ref(j)) {
        alt_end = j + 1;
        ref_end = (last_codon + 1).max(i + 1).min(ref_aa.len());
    }

    // Nonsense: the first changed residue is a stop.
    if alt_aa.get(i) == Some(&'*') && i < alt_end {
        return out.substitution(i, ref_aa[i], '*');
    }

    let del: &[char] = &ref_aa[i..ref_end.max(i)];
    let ins: &[char] = &alt_aa[i..alt_end.max(i)];

    if del.is_empty() && !ins.is_empty() {
        // Duplication: the inserted residues repeat those just before them.
        if i >= ins.len() && ref_aa[i - ins.len()..i] == *ins {
            return out.duplication(i - ins.len(), i - 1);
        }
        if i == 0 {
            return Err(HgvsError::UnsupportedOperation(
                "N-terminal protein insertions are not supported".into(),
            ));
        }
        return out.insertion(i - 1, i, ins);
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
            uncertain: false,
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

    /// `p.(Xxx1=)` or `p.(Xxx1_Yyy2=)` over the codons the edit touched, or a
    /// bare `p.(=)` when they lie past the protein.
    fn identity(&self, first_codon: usize, last_codon: usize) -> Result<PVariant, HgvsError> {
        if first_codon >= self.ref_aa.len() {
            return self.variant(None, AaEdit::Identity { uncertain: false });
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
        let length = self.alt_aa[i..]
            .iter()
            .position(|&c| c == '*')
            .map_or("?".to_string(), |k| (k + 1).to_string());
        self.variant(
            Some(self.span(i, None)),
            AaEdit::Fs {
                ref_: String::new(),
                alt: aa1_to_aa3(alt_c).to_string(),
                term: Some("Ter".to_string()),
                length: Some(length),
                uncertain: false,
            },
        )
    }

    /// `p.(Ter#Xxxext*N)`: the stop is read through; N is the distance to the
    /// next stop, or `?` if translation runs off the end.
    fn extension(&self, i: usize) -> Result<PVariant, HgvsError> {
        let alt_c = self.alt_aa.get(i).copied().unwrap_or('*');
        let next_stop = self.alt_aa[i + 1..].iter().position(|&c| c == '*');
        self.variant(
            Some(AaInterval {
                start: AAPosition {
                    base: ProteinPos(i as i32).to_hgvs(),
                    aa: "Ter".to_string(),
                    uncertain: false,
                },
                end: None,
                uncertain: false,
            }),
            AaEdit::Ext {
                ref_: "Ter".into(),
                alt: aa1_to_aa3(alt_c).to_string(),
                aaterm: Some("*".to_string()),
                length: Some(next_stop.map_or("?".to_string(), |k| (k + 1).to_string())),
                uncertain: next_stop.is_none(),
            },
        )
    }
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
            .resolve_with(start, end, |s, e| {
                Ok(coding[s.min(coding.len())..e.min(coding.len())].to_string())
            })
            .unwrap();
        describe(&CodingChange {
            coding: &coding,
            cds_len: cds.len(),
            edit: resolved,
            protein_ac: "NP".into(),
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
    fn n_terminal_insertion_is_unsupported() {
        assert!(p(CDS, UTR, ins("GTG"), 0, 1).contains("N-terminal"));
    }

    #[test]
    fn a_change_past_the_stop_is_silent() {
        // Written against the codon it touches, as the previous implementation did.
        assert_eq!(p(CDS, UTR, sub("C", "G"), 21, 22), "NP:p.Pro8=");
    }
}
