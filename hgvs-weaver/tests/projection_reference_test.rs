//! A projection names the target's bases. The transcript record and the
//! genome can differ at a base; carrying the stated edit across then writes a
//! reference the target does not have. After projecting, the reference is
//! re-read from the target and a change the target already carries is `=`.

mod support;

use hgvs_weaver::data::Strand;
use hgvs_weaver::mapper::VariantMapper;
use hgvs_weaver::parse_hgvs_variant;
use hgvs_weaver::SequenceVariant;
use support::{exon, transcript, Provider};

const GENOME_AC: &str = "NC_D.1";

/// 120 bases of ACGT, so index i holds "ACGT"[i % 4].
fn genome() -> String {
    "ACGT".repeat(30)
}

fn revcomp(s: &str) -> String {
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

/// Three transcripts: plus and minus strand over genome[0..100) with one
/// base that disagrees with the genome at transcript index 20 (c.11, with
/// the CDS at 10..=39), and a two-exon plus-strand one whose intron is
/// genome[50..60) and whose CDS reaches into exon 2.
fn provider() -> Provider {
    let g = genome();
    let mut plus: Vec<char> = g[..100].chars().collect();
    assert_eq!(plus[20], 'A');
    plus[20] = 'C'; // the record says C where the genome says A
    let mut minus: Vec<char> = revcomp(&g[..100]).chars().collect();
    assert_eq!(minus[20], 'A'); // genome index 79 is T; complemented, A
    minus[20] = 'G';
    let spliced = format!("{}{}", &g[..50], &g[60..110]);
    Provider::new()
        .sequence(GENOME_AC, &g)
        .sequence("NM_PLUS.1", &plus.into_iter().collect::<String>())
        .sequence("NM_MINUS.1", &minus.into_iter().collect::<String>())
        .sequence("NM_SPLICED.1", &spliced)
        .transcript(transcript(
            "NM_PLUS.1",
            GENOME_AC,
            Strand::Plus,
            Some((10, 39)),
            vec![exon((0, 100), (0, 99), Strand::Plus)],
        ))
        .transcript(transcript(
            "NM_MINUS.1",
            GENOME_AC,
            Strand::Minus,
            Some((10, 39)),
            vec![exon((0, 100), (0, 99), Strand::Minus)],
        ))
        .transcript(transcript(
            "NM_SPLICED.1",
            GENOME_AC,
            Strand::Plus,
            Some((10, 69)),
            vec![
                exon((0, 50), (0, 49), Strand::Plus),
                exon((50, 100), (60, 109), Strand::Plus),
            ],
        ))
}

fn c_to_g(mapper: &VariantMapper, s: &str) -> String {
    match parse_hgvs_variant(s).unwrap() {
        SequenceVariant::Coding(c) => mapper.c_to_g(&c, None).unwrap().to_string(),
        other => panic!("{other} is not c."),
    }
}

fn g_to_c(mapper: &VariantMapper, s: &str, tx: &str) -> String {
    match parse_hgvs_variant(s).unwrap() {
        SequenceVariant::Genomic(g) => mapper.g_to_c(&g, tx).unwrap().to_string(),
        other => panic!("{other} is not g."),
    }
}

#[test]
fn a_projection_states_the_genomes_bases_not_the_records() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    // c.11 is transcript index 20: the record says C, the genome (g.21) says A.
    assert_eq!(c_to_g(&mapper, "NM_PLUS.1:c.11C>G"), "NC_D.1:g.21A>G");
    // The genome already holds the alternate: nothing changes on the genome.
    assert_eq!(c_to_g(&mapper, "NM_PLUS.1:c.11C>A"), "NC_D.1:g.21=");
    // Stated bases on a deletion, delins and inversion are re-read too.
    assert_eq!(c_to_g(&mapper, "NM_PLUS.1:c.11delC"), "NC_D.1:g.21delA");
    assert_eq!(
        c_to_g(&mapper, "NM_PLUS.1:c.11_12delCCinsTT"),
        "NC_D.1:g.21_22delACinsTT"
    );
    // No change on the record is a change on the genome where they differ,
    // and no change where they agree.
    assert_eq!(c_to_g(&mapper, "NM_PLUS.1:c.11="), "NC_D.1:g.21A>C");
    assert_eq!(c_to_g(&mapper, "NM_PLUS.1:c.12="), "NC_D.1:g.22=");
    assert_eq!(
        g_to_c(&mapper, "NC_D.1:g.21=", "NM_PLUS.1"),
        "NM_PLUS.1:c.11C>A"
    );
    assert_eq!(c_to_g(&mapper, "NM_MINUS.1:c.11="), "NC_D.1:g.80T>C");
    // Where record and genome agree nothing changes, and unstated bases are
    // never filled in.
    assert_eq!(c_to_g(&mapper, "NM_PLUS.1:c.12C>T"), "NC_D.1:g.22C>T");
    assert_eq!(c_to_g(&mapper, "NM_PLUS.1:c.11del"), "NC_D.1:g.21del");
    assert_eq!(
        c_to_g(&mapper, "NM_PLUS.1:c.11_12insTT"),
        "NC_D.1:g.21_22insTT"
    );
}

#[test]
fn a_projection_states_the_transcripts_bases_not_the_genomes() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_eq!(
        g_to_c(&mapper, "NC_D.1:g.21A>G", "NM_PLUS.1"),
        "NM_PLUS.1:c.11C>G"
    );
    // The record already carries the alternate: the honest answer is =.
    assert_eq!(
        g_to_c(&mapper, "NC_D.1:g.21A>C", "NM_PLUS.1"),
        "NM_PLUS.1:c.11="
    );
    assert_eq!(
        g_to_c(&mapper, "NC_D.1:g.22C>T", "NM_PLUS.1"),
        "NM_PLUS.1:c.12C>T"
    );
}

#[test]
fn the_re_read_is_in_the_targets_orientation_on_the_minus_strand() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    // Transcript index 20 is genome index 79 (g.80), a T; the record says G.
    // The edit is complemented to the genome, then the reference re-read.
    assert_eq!(c_to_g(&mapper, "NM_MINUS.1:c.11G>C"), "NC_D.1:g.80T>G");
    assert_eq!(c_to_g(&mapper, "NM_MINUS.1:c.11G>A"), "NC_D.1:g.80=");
    // On the transcript the genome's T>C reads A>G, and the record holds G.
    assert_eq!(
        g_to_c(&mapper, "NC_D.1:g.80T>C", "NM_MINUS.1"),
        "NM_MINUS.1:c.11="
    );
    assert_eq!(
        g_to_c(&mapper, "NC_D.1:g.80T>A", "NM_MINUS.1"),
        "NM_MINUS.1:c.11G>T"
    );
}

#[test]
fn an_intronic_position_has_no_transcript_base_to_re_read() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    // Genome index 52 is in the intron, an A; the edit is carried as given.
    assert_eq!(
        g_to_c(&mapper, "NC_D.1:g.53A>G", "NM_SPLICED.1"),
        "NM_SPLICED.1:c.40+3A>G"
    );
    // And back to the genome it is re-read against the genome, which has it.
    assert_eq!(c_to_g(&mapper, "NM_SPLICED.1:c.40+3A>G"), "NC_D.1:g.53A>G");
    assert_eq!(c_to_g(&mapper, "NM_SPLICED.1:c.40+3C>G"), "NC_D.1:g.53A>G");
}
