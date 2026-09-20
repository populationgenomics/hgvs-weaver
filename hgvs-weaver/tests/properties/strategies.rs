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
        edit: NaEdit::RefAlt {
            ref_: None,
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
        NaEdit::RefAlt { alt, .. } => {
            format!("{}{}{}", &seq[..s], alt.as_deref().unwrap_or(""), &seq[e..])
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
                    exons.push(ExonData {
                        transcript_start: hgvs_weaver::coords::TranscriptPos(ts as i32),
                        transcript_end: hgvs_weaver::coords::TranscriptPos(te as i32),
                        reference_start: hgvs_weaver::coords::GenomicPos(g_start as i32),
                        reference_end: hgvs_weaver::coords::GenomicPos(g_end as i32),
                        alt_strand: strand,
                        cigar: format!("{}M", te - ts),
                    });
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
        TranscriptData {
            ac: self.ac.clone(),
            gene: "PROP".into(),
            cds_start_index: Some(hgvs_weaver::coords::TranscriptPos(self.cds_start as i32)),
            cds_end_index: Some(hgvs_weaver::coords::TranscriptPos(self.cds_end as i32)),
            strand: self.strand,
            reference_accession: self.reference_ac.clone(),
            exons: self.exons.clone(),
        }
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
}
