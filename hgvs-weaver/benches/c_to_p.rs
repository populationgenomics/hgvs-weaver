//! `c_to_p` on a synthetic three-exon transcript of typical size: a 1,500-base
//! CDS between a 100-base 5' UTR and a 1,400-base 3' UTR. Each case is one kind
//! of edit, since the protein reading differs by kind.

#[path = "../tests/support/mod.rs"]
mod support;

use criterion::{criterion_group, criterion_main, Criterion};
use hgvs_weaver::data::Strand;
use hgvs_weaver::mapper::VariantMapper;
use hgvs_weaver::{parse_hgvs_variant, CVariant, SequenceVariant};
use std::hint::black_box;
use support::{exon, transcript, Provider};

const UTR5: usize = 100;
const CODONS: usize = 499;
const UTR3: usize = 1_400;
/// Transcript-index boundaries of the three exons.
const EXONS: [(usize, usize); 3] = [
    (0, 700),
    (700, 1_900),
    (1_900, UTR5 + (CODONS + 1) * 3 + UTR3),
];
/// Bases of intron between consecutive exons.
const INTRON: usize = 500;
const GENOMIC_START: usize = 1_000;

/// A deterministic sequence of `n` bases.
fn bases(n: usize, seed: u64) -> String {
    let mut state = seed;
    (0..n)
        .map(|_| {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1);
            b"ACGT"[(state >> 62) as usize] as char
        })
        .collect()
}

/// ATG, `CODONS - 1` sense codons, then TAA.
fn cds() -> String {
    let mut out = String::from("ATG");
    let mut i = 0u64;
    while out.len() < CODONS * 3 {
        let codon = bases(3, 0x9e37_79b9 ^ i);
        if !matches!(codon.as_str(), "TAA" | "TAG" | "TGA") {
            out.push_str(&codon);
        }
        i += 1;
    }
    out.push_str("TAA");
    out
}

fn provider() -> (Provider, String) {
    let tx = format!("{}{}{}", bases(UTR5, 1), cds(), bases(UTR3, 2));
    let mut genome = bases(GENOMIC_START, 3);
    let mut exons = Vec::new();
    for (i, &(start, end)) in EXONS.iter().enumerate() {
        if i > 0 {
            genome.push_str(&bases(INTRON, 4 + i as u64));
        }
        let g = genome.len() as i32;
        genome.push_str(&tx[start..end]);
        exons.push(exon(
            (start as i32, end as i32),
            (g, g + (end - start) as i32 - 1),
            Strand::Plus,
        ));
    }
    genome.push_str(&bases(GENOMIC_START, 9));
    let cds_end = (UTR5 + (CODONS + 1) * 3 - 1) as i32;
    let provider = Provider::new()
        .sequence("NC_000099.1", &genome)
        .sequence("NM_000099.1", &tx)
        .transcript(transcript(
            "NM_000099.1",
            "NC_000099.1",
            Strand::Plus,
            Some((UTR5 as i32, cds_end)),
            exons,
        ))
        .protein_for("NM_000099.1", "NP_000099.1");
    (provider, tx)
}

fn coding(s: &str) -> CVariant {
    match parse_hgvs_variant(s).unwrap() {
        SequenceVariant::Coding(c) => c,
        other => panic!("{other} is not c."),
    }
}

fn bench_c_to_p(c: &mut Criterion) {
    let (provider, tx) = provider();
    let mapper = VariantMapper::new(&provider);
    let base = |c_pos: usize| &tx[UTR5 + c_pos - 1..UTR5 + c_pos];
    let other = |b: &str| if b == "A" { "G" } else { "A" };
    let missense_pos = 700;
    let cases = [
        (
            "substitution",
            format!(
                "NM_000099.1:c.{missense_pos}{}>{}",
                base(missense_pos),
                other(base(missense_pos))
            ),
        ),
        ("frameshift", "NM_000099.1:c.701del".to_string()),
        ("in_frame_deletion", "NM_000099.1:c.701_703del".to_string()),
    ];
    let mut group = c.benchmark_group("c_to_p");
    for (name, hgvs) in &cases {
        let variant = coding(hgvs);
        mapper
            .c_to_p(&variant, None)
            .unwrap_or_else(|e| panic!("{hgvs}: {e}"));
        group.bench_function(*name, |b| {
            b.iter(|| mapper.c_to_p(black_box(&variant), None).unwrap())
        });
    }
    group.finish();
}

criterion_group!(benches, bench_c_to_p);
criterion_main!(benches);
