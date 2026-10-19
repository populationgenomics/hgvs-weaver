//! Shared generators for the property tests: sequences with runs, edits placed
//! on them, and transcripts with random exon structure on both strands, all
//! served by one in-memory provider.

use hgvs_weaver::data::{ExonData, Strand, TranscriptData};
use hgvs_weaver::edits::NaEdit;
use proptest::prelude::*;

pub use crate::support::Provider;

pub const BASES: [char; 4] = ['A', 'C', 'G', 'T'];

pub fn revcomp(s: &str) -> String {
    s.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            other => other,
        })
        .collect()
}

/// A DNA string biased towards repeats and homopolymers, so that shifting and
/// ambiguity have something to bite on.
pub fn dna(min: usize, max: usize) -> impl Strategy<Value = String> {
    prop_oneof![
        3 => prop::collection::vec(prop::sample::select(BASES.to_vec()), min..=max)
            .prop_map(|v| v.into_iter().collect()),
        2 => (prop::collection::vec(prop::sample::select(BASES.to_vec()), 1..=4), min..=max)
            .prop_map(|(unit, len)| unit.iter().cycle().take(len).collect()),
        1 => (
            prop::collection::vec(prop::sample::select(BASES.to_vec()), min..=max),
            prop::collection::vec(prop::sample::select(BASES.to_vec()), 1..=3),
            1usize..=4,
        )
            .prop_map(|(mut v, unit, copies)| {
                // splice a run of the unit into a random sequence
                let at = v.len() / 2;
                let run: Vec<char> = unit.iter().cycle().take(unit.len() * copies).cloned().collect();
                v.splice(at..at, run);
                v.into_iter().collect()
            }),
    ]
}

/// An edit written over an HGVS range on a sequence of `len` bases:
/// `(start, end, edit)` with `[start, end)` the range HGVS names (the two
/// flanking bases for an insertion). Insertions are always well inside.
#[derive(Debug, Clone)]
pub struct PlacedHgvs {
    pub start: usize,
    pub end: usize,
    pub edit: NaEdit,
}

/// A sequence together with an edit placed on it, so nothing is rejected.
pub fn seq_and_edit(del_ins_only: bool) -> impl Strategy<Value = (String, PlacedHgvs)> {
    dna(8, 80).prop_flat_map(move |seq| {
        let len = seq.len();
        edit_on_kinds(len, del_ins_only).prop_map(move |p| (seq.clone(), p))
    })
}

pub fn edit_on_kinds(len: usize, del_ins_only: bool) -> impl Strategy<Value = PlacedHgvs> {
    let span = (1usize..=len.saturating_sub(2).max(1))
        .prop_flat_map(move |l| (0..=len.saturating_sub(l + 1)).prop_map(move |s| (s, s + l)));
    let bases = |n: usize| {
        prop::collection::vec(prop::sample::select(BASES.to_vec()), 1..=n)
            .prop_map(|v| v.into_iter().collect::<String>())
    };
    let del = span.clone().prop_map(|(s, e)| PlacedHgvs {
        start: s,
        end: e,
        edit: NaEdit::Del {
            ref_: None,
            uncertain: false,
        },
    });
    let dup = span.clone().prop_map(|(s, e)| PlacedHgvs {
        start: s,
        end: e,
        edit: NaEdit::Dup {
            ref_: None,
            uncertain: false,
        },
    });
    let ins = (1usize..len.saturating_sub(1).max(2), bases(6)).prop_map(|(b, alt)| PlacedHgvs {
        start: b,
        end: b + 2,
        edit: NaEdit::Ins {
            alt: Some(alt),
            uncertain: false,
        },
    });
    let delins = (span, bases(6)).prop_map(|((s, e), alt)| PlacedHgvs {
        start: s,
        end: e,
        // `delinsX`, as the parser reads it: the deleted bases unstated
        edit: NaEdit::RefAlt {
            ref_: Some(String::new()),
            alt: Some(alt),
            uncertain: false,
        },
    });
    if del_ins_only {
        prop_oneof![del, ins].boxed()
    } else {
        prop_oneof![del, dup, ins, delins].boxed()
    }
}

/// What the edit does to the sequence, computed the dumb way.
pub fn apply(seq: &str, p: &PlacedHgvs) -> String {
    let (s, e) = (p.start.min(seq.len()), p.end.min(seq.len()));
    match &p.edit {
        NaEdit::Del { .. } => format!("{}{}", &seq[..s], &seq[e..]),
        NaEdit::Dup { .. } => format!("{}{}{}", &seq[..e], &seq[s..e], &seq[e..]),
        NaEdit::Ins { alt, .. } => {
            // Written between p.start and p.end - 1; the insert goes in front of
            // the second flank, which may be one past the last base.
            let anchor = (p.end - 1).min(seq.len());
            format!(
                "{}{}{}",
                &seq[..anchor],
                alt.as_deref().unwrap_or(""),
                &seq[anchor..]
            )
        }
        // `=`, with no bases stated: the range as it is.
        NaEdit::RefAlt {
            ref_: None,
            alt: None,
            ..
        } => seq.to_string(),
        NaEdit::RefAlt { alt, .. } => {
            format!("{}{}{}", &seq[..s], alt.as_deref().unwrap_or(""), &seq[e..])
        }
        NaEdit::Inv { .. } => {
            let inverted: String = seq[s..e]
                .bytes()
                .rev()
                .map(|b| match b {
                    b'A' => 'T',
                    b'C' => 'G',
                    b'G' => 'C',
                    b'T' => 'A',
                    other => panic!("apply: not a base: {}", other as char),
                })
                .collect();
            format!("{}{inverted}{}", &seq[..s], &seq[e..])
        }
        other => panic!("apply: unsupported edit {other:?}"),
    }
}

/// A transcript with random exon structure on a random strand, and the
/// genome it sits on. Transcript sequence is `utr5 + cds + utr3`; the CDS
/// starts with ATG, has no internal stop, and ends with a stop codon.
#[derive(Debug, Clone)]
pub struct Gene {
    pub ac: String,
    pub reference_ac: String,
    pub strand: Strand,
    pub transcript_seq: String,
    pub genome: String,
    pub exons: Vec<ExonData>,
    pub cds_start: usize,
    /// Index of the last base of the stop codon.
    pub cds_end: usize,
    /// Exon boundaries in transcript indices, `[start, end)`, transcript order.
    pub tx_exons: Vec<(usize, usize)>,
    /// Intron lengths between consecutive exons, transcript order.
    pub introns: Vec<usize>,
}

const CODONS_NO_STOP: [&str; 61] = [
    "TTT", "TTC", "TTA", "TTG", "CTT", "CTC", "CTA", "CTG", "ATT", "ATC", "ATA", "ATG", "GTT",
    "GTC", "GTA", "GTG", "TCT", "TCC", "TCA", "TCG", "CCT", "CCC", "CCA", "CCG", "ACT", "ACC",
    "ACA", "ACG", "GCT", "GCC", "GCA", "GCG", "TAT", "TAC", "CAT", "CAC", "CAA", "CAG", "AAT",
    "AAC", "AAA", "AAG", "GAT", "GAC", "GAA", "GAG", "TGT", "TGC", "TGG", "CGT", "CGC", "CGA",
    "CGG", "AGT", "AGC", "AGA", "AGG", "GGT", "GGC", "GGA", "GGG",
];

pub fn gene() -> impl Strategy<Value = Gene> {
    let cds = prop::collection::vec(prop::sample::select(CODONS_NO_STOP.to_vec()), 3..=20)
        .prop_map(|codons| format!("ATG{}TAA", codons.concat()));
    (
        dna(0, 12),                                // utr5
        cds,                                       // cds
        dna(3, 15),                                // utr3
        prop::collection::vec(1usize..=40, 0..=3), // intron lengths
        prop::bool::ANY,                           // minus strand?
        dna(2, 10),                                // genomic flank
    )
        .prop_flat_map(|(utr5, cds, utr3, introns, minus, flank)| {
            let transcript_seq = format!("{utr5}{cds}{utr3}");
            let tx_len = transcript_seq.len();
            let n_exons = introns.len() + 1;
            // n_exons - 1 cut points strictly inside the transcript
            prop::collection::vec(1usize..tx_len, n_exons - 1).prop_map(move |mut cuts| {
                cuts.sort_unstable();
                cuts.dedup();
                let mut bounds = vec![0usize];
                bounds.extend(cuts.iter().cloned());
                bounds.push(tx_len);
                let tx_exons: Vec<(usize, usize)> =
                    bounds.windows(2).map(|w| (w[0], w[1])).collect();
                let introns: Vec<usize> =
                    introns.iter().cloned().take(tx_exons.len() - 1).collect();

                // Genome on the plus strand. Transcript order == genomic order
                // for plus; reversed and complemented for minus.
                let strand = if minus { Strand::Minus } else { Strand::Plus };
                let intron_seq = |i: usize| -> String {
                    // deterministic intron bases distinct-ish from exons
                    "T".repeat(i)
                };
                let mut genome = flank.clone();
                let mut exons: Vec<ExonData> = Vec::new();
                let order: Vec<usize> = if minus {
                    (0..tx_exons.len()).rev().collect()
                } else {
                    (0..tx_exons.len()).collect()
                };
                for (k, &ei) in order.iter().enumerate() {
                    let (ts, te) = tx_exons[ei];
                    let piece = &transcript_seq[ts..te];
                    let g_start = genome.len();
                    if minus {
                        genome.push_str(&revcomp(piece));
                    } else {
                        genome.push_str(piece);
                    }
                    let g_end = genome.len() - 1;
                    exons.push(ExonData::new(
                        hgvs_weaver::coords::TranscriptPos(ts as i32),
                        hgvs_weaver::coords::TranscriptPos(te as i32),
                        hgvs_weaver::coords::GenomicPos(g_start as i32),
                        hgvs_weaver::coords::GenomicPos(g_end as i32),
                        strand,
                        format!("{}M", te - ts),
                    ));
                    if k + 1 < order.len() {
                        // the intron following this exon in genomic order
                        let intron_index = if minus { order[k + 1] } else { ei };
                        genome.push_str(&intron_seq(introns[intron_index]));
                    }
                }
                genome.push_str(&flank);
                Gene {
                    ac: "NM_PROP.1".into(),
                    reference_ac: "NC_PROP.1".into(),
                    strand,
                    cds_start: utr5.len(),
                    cds_end: utr5.len() + cds.len() - 1,
                    transcript_seq: transcript_seq.clone(),
                    genome,
                    exons,
                    tx_exons,
                    introns,
                }
            })
        })
}

impl Gene {
    pub fn transcript_data(&self) -> TranscriptData {
        TranscriptData::new(
            self.ac.clone(),
            "PROP",
            Some((
                hgvs_weaver::coords::TranscriptPos(self.cds_start as i32),
                hgvs_weaver::coords::TranscriptPos(self.cds_end as i32),
            )),
            self.strand,
            self.reference_ac.clone(),
            self.exons.clone(),
        )
    }
    /// The gene's transcript, genome and protein, served by one provider.
    pub fn provider(&self) -> Provider {
        let protein =
            hgvs_weaver::utils::translate(&self.transcript_seq[self.cds_start..=self.cds_end]);
        Provider::new()
            .sequence(&self.ac, &self.transcript_seq)
            .sequence(&self.reference_ac, &self.genome)
            .sequence("NP_PROP.1", &protein)
            .transcript(self.transcript_data())
            .protein_for(&self.ac, "NP_PROP.1")
    }

    /// Bases of genomic flank on either side of the transcript.
    pub fn flank(&self) -> usize {
        self.exons
            .iter()
            .map(|e| e.reference_start.0 as usize)
            .min()
            .unwrap()
    }

    /// Every base of the gene's genome in transcript order, named the way
    /// HGVS names it: exonic bases by index, intron bases from the nearer
    /// exon (the earlier on a tie), flanking bases as offsets from the
    /// transcript's first and last base.
    pub fn points(&self) -> Vec<TxPoint> {
        let flank = self.flank() as i32;
        let tx_len = self.transcript_seq.len() as i32;
        let mut points: Vec<TxPoint> = (1..=flank)
            .rev()
            .map(|k| TxPoint {
                index: 0,
                offset: -k,
            })
            .collect();
        for (k, &(ts, te)) in self.tx_exons.iter().enumerate() {
            points.extend((ts..te).map(|i| TxPoint {
                index: i as i32,
                offset: 0,
            }));
            if let Some(&len) = self.introns.get(k) {
                let len = len as i32;
                points.extend((1..=len).map(|i| {
                    if i <= len + 1 - i {
                        TxPoint {
                            index: te as i32 - 1,
                            offset: i,
                        }
                    } else {
                        TxPoint {
                            index: te as i32,
                            offset: -(len + 1 - i),
                        }
                    }
                }));
            }
        }
        points.extend((1..=flank).map(|k| TxPoint {
            index: tx_len - 1,
            offset: k,
        }));
        points
    }

    /// Whether every exonic base of the genome is the transcript's base:
    /// false once [`gene_and_edit`] has made them disagree at one.
    pub fn genome_agrees_with_transcript(&self) -> bool {
        let exonic: Vec<TxPoint> = self
            .points()
            .into_iter()
            .filter(|p| p.offset == 0)
            .collect();
        self.bases_at(&exonic) == self.transcript_seq
    }

    /// The genome index a point names.
    pub fn genome_index(&self, p: TxPoint) -> usize {
        let exon = self
            .exons
            .iter()
            .find(|e| e.transcript_start.0 <= p.index && p.index < e.transcript_end.0)
            .unwrap();
        let into = p.index - exon.transcript_start.0;
        let g = match self.strand {
            Strand::Plus => exon.reference_start.0 + into + p.offset,
            Strand::Minus => exon.reference_end.0 - into - p.offset,
        };
        g as usize
    }

    /// The bases at `points`, read in the transcript's sense.
    pub fn bases_at(&self, points: &[TxPoint]) -> String {
        let genomic: String = points
            .iter()
            .map(|&p| self.genome.as_bytes()[self.genome_index(p)] as char)
            .collect();
        match self.strand {
            Strand::Plus => genomic,
            // complemented base by base; `points` is already in transcript order
            Strand::Minus => revcomp(&genomic).chars().rev().collect(),
        }
    }
}

/// A base on the genome named from the transcript: a transcript index and an
/// offset from it into an intron or a flank.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TxPoint {
    pub index: i32,
    pub offset: i32,
}

/// An edit over the transcript points `start..=end`; an insertion lies
/// between them, so `end` is the point after `start`.
#[derive(Debug, Clone)]
pub struct TxEdit {
    pub start: TxPoint,
    pub end: TxPoint,
    pub edit: NaEdit,
}

/// The genome with the base at `p` changed, so the transcript and the genome
/// disagree there.
fn discordant_at(g: &Gene, p: TxPoint, shift: usize) -> String {
    let at = g.genome_index(p);
    let mut genome = g.genome.clone().into_bytes();
    let k = BASES.iter().position(|&b| b as u8 == genome[at]).unwrap();
    genome[at] = BASES[(k + shift) % 4] as u8;
    String::from_utf8(genome).unwrap()
}

/// A gene and an edit anywhere on it: substitution, identity (bare or with
/// its bases stated), deletion, duplication, insertion, deletion-insertion,
/// inversion or repeat, over exonic, intronic, UTR and flanking bases, the
/// start and stop codons included, on either strand. Bases an edit states are
/// the transcript's; sometimes one exonic genome base under the edit then
/// differs, so a projection has to rewrite the edit against the genome's.
pub fn gene_and_edit() -> impl Strategy<Value = (Gene, TxEdit)> {
    gene().prop_flat_map(edit_on)
}

/// A gene and two edits on it, each as [`gene_and_edit`] makes one; the
/// genome may differ from the transcript under either.
pub fn gene_and_two_edits() -> impl Strategy<Value = (Gene, TxEdit, TxEdit)> {
    gene_and_edit().prop_flat_map(|(g, first)| {
        edit_on(g).prop_map(move |(g, second)| (g, first.clone(), second))
    })
}

/// An edit on `g` as [`gene_and_edit`] makes one, with `g`'s genome changed
/// under it when the edit makes the two disagree.
fn edit_on(g: Gene) -> impl Strategy<Value = (Gene, TxEdit)> {
    let n = g.points().len();
    (
        Just(g),
        0..n,
        // mostly short; sometimes long enough to cover the whole CDS
        prop_oneof![4 => 0usize..=8, 1 => 0..n],
        0u8..8,
        dna(1, 6),
        1usize..=3,
        1i32..=6,
        prop::option::weighted(0.25, (0usize..10_000, 1usize..=3)),
    )
        .prop_map(|(mut g, s, span, kind, alt, shift, copies, discord)| {
            let points = g.points();
            let last = points.len() - 1;
            let (start, end, edit) = match kind {
                0 => {
                    let r = g.bases_at(&points[s..=s]);
                    let k = BASES.iter().position(|&b| r == b.to_string()).unwrap();
                    let a = BASES[(k + shift) % 4].to_string();
                    let edit = NaEdit::RefAlt {
                        ref_: Some(r),
                        alt: Some(a),
                        uncertain: false,
                    };
                    (s, s, edit)
                }
                1 => {
                    // `=`, or `CAG=` with the unchanged bases stated
                    let e = (s + span % 4).min(last);
                    let stated = (shift == 1).then(|| g.bases_at(&points[s..=e]));
                    let edit = NaEdit::RefAlt {
                        ref_: stated.clone(),
                        alt: stated,
                        uncertain: false,
                    };
                    (s, e, edit)
                }
                2 => (
                    s,
                    (s + span).min(last),
                    NaEdit::Del {
                        ref_: None,
                        uncertain: false,
                    },
                ),
                3 => (
                    s,
                    (s + span).min(last),
                    NaEdit::Dup {
                        ref_: None,
                        uncertain: false,
                    },
                ),
                4 => (
                    s.min(last - 1),
                    s.min(last - 1) + 1,
                    NaEdit::Ins {
                        alt: Some(alt),
                        uncertain: false,
                    },
                ),
                5 => (
                    s,
                    (s + span).min(last),
                    // `delinsX`, as the parser reads it: the deleted bases unstated
                    NaEdit::RefAlt {
                        ref_: Some(String::new()),
                        alt: Some(alt),
                        uncertain: false,
                    },
                ),
                6 => (
                    s.min(last - 1),
                    (s.min(last - 1) + span.max(1)).min(last),
                    NaEdit::Inv {
                        ref_: None,
                        uncertain: false,
                    },
                ),
                _ => {
                    let unit = g.bases_at(&points[s..=(s + span % 3).min(last)]);
                    let edit = NaEdit::Repeat {
                        ref_: Some(unit),
                        min: copies,
                        max: copies,
                        uncertain: false,
                    };
                    (s, s, edit)
                }
            };
            let exonic_under_edit: Vec<TxPoint> = points[start..=end]
                .iter()
                .copied()
                .filter(|p| p.offset == 0)
                .collect();
            if let (Some((pick, genome_shift)), false) = (discord, exonic_under_edit.is_empty()) {
                let p = exonic_under_edit[pick % exonic_under_edit.len()];
                g.genome = discordant_at(&g, p, genome_shift);
            }
            let edit = TxEdit {
                start: points[start],
                end: points[end],
                edit,
            };
            (g, edit)
        })
}

/// A gene and an allele of one to four edits on it, in transcript order, over
/// exonic and intronic bases (not the flanks) on either strand: substitution,
/// deletion, duplication, insertion, deletion-insertion or inversion, none
/// over two exons. At least one unchanged base lies between two members, so
/// no member's change meets another's.
pub fn gene_and_cis_edits() -> impl Strategy<Value = (Gene, Vec<TxEdit>)> {
    gene()
        .prop_flat_map(|g| {
            let member = (0usize..10_000, 0usize..=3, 0u8..6, dna(1, 4), 1usize..=3);
            (Just(g), prop::collection::vec(member, 1..=4))
        })
        .prop_map(|(g, picks)| {
            let points = g.points();
            let flank = g.flank();
            // the points between the flanks, in transcript order
            let inner = &points[flank..points.len() - flank];
            let mut starts: Vec<usize> = picks.iter().map(|p| p.0 % inner.len()).collect();
            starts.sort_unstable();
            let mut edits = Vec::new();
            let mut free = 0;
            for (start, (_, span, kind, alt, shift)) in starts.into_iter().zip(&picks) {
                // the first member always fits; a later one may not
                let s = if edits.is_empty() {
                    start.min(inner.len() - 2)
                } else {
                    start.max(free)
                };
                if s + 1 >= inner.len() {
                    break;
                }
                let mut e = match kind {
                    0 => s,
                    3 => s + 1,
                    _ => (s + span).min(inner.len() - 1),
                };
                // A range over two exons is read on the spliced transcript,
                // which has no genomic counterpart to come back from.
                let exon_of = |p: TxPoint| {
                    g.tx_exons
                        .iter()
                        .position(|&(a, b)| a as i32 <= p.index && p.index < b as i32)
                };
                let exons: Vec<_> = inner[s..=e]
                    .iter()
                    .filter(|p| p.offset == 0)
                    .map(|&p| exon_of(p))
                    .collect();
                if *kind != 3 && exons.windows(2).any(|w| w[0] != w[1]) {
                    e = s;
                }
                let edit = match kind {
                    0 => {
                        let r = g.bases_at(&inner[s..=s]);
                        let k = BASES.iter().position(|&b| r == b.to_string()).unwrap();
                        NaEdit::RefAlt {
                            ref_: Some(r),
                            alt: Some(BASES[(k + shift) % 4].to_string()),
                            uncertain: false,
                        }
                    }
                    1 => NaEdit::Del {
                        ref_: None,
                        uncertain: false,
                    },
                    2 => NaEdit::Dup {
                        ref_: None,
                        uncertain: false,
                    },
                    3 => NaEdit::Ins {
                        alt: Some(alt.clone()),
                        uncertain: false,
                    },
                    4 => NaEdit::RefAlt {
                        ref_: Some(String::new()),
                        alt: Some(alt.clone()),
                        uncertain: false,
                    },
                    _ => NaEdit::Inv {
                        ref_: None,
                        uncertain: false,
                    },
                };
                edits.push(TxEdit {
                    start: inner[s],
                    end: inner[e],
                    edit,
                });
                free = e + 2;
            }
            (g, edits)
        })
}
