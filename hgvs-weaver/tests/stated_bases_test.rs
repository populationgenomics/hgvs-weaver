//! A base a description states, which the reference it is read against does
//! not carry, is an error from every operation that reads the description: a
//! `c.` exonic base against the transcript record, an intronic one against the
//! genome on the transcript's strand, a `g.` base against the genome and a
//! `p.` residue against the protein. Each test names one operation.

mod support;

use hgvs_weaver::data::Strand;
use hgvs_weaver::equivalence::VariantEquivalence;
use hgvs_weaver::error::HgvsError;
use hgvs_weaver::mapper::VariantMapper;
use hgvs_weaver::{parse_hgvs_variant, CVariant, GVariant, PVariant, RVariant, SequenceVariant};
use support::{exon, transcript, Provider};

/// The transcript record: 5'UTR, exon 1's coding part (c.1_9, `ATG AAA CTG`),
/// exon 2 from c.10 (`GCC TAT CGC TAA`), then the 3'UTR. c.12 is a C.
const RECORD: &str = "GGGGGATGAAACTGGCCTATCGCTAACCGTATAAGTAAGG";
const INTRON: &str = "GTCCATTCAG";
const PROTEIN: &str = "MKLAYR";

/// The plus-strand genome: exon 1 at [0, 14), the intron at [14, 24), exon 2
/// at [24, 50). The genome has T where the record has C at c.12 (index 26),
/// so c.9+1 is g.15 (a G) and c.12 is g.27.
fn plus_genome() -> String {
    let mut g: Vec<char> = format!("{}{INTRON}{}", &RECORD[..14], &RECORD[14..])
        .chars()
        .collect();
    assert_eq!(g[26], 'C');
    g[26] = 'T';
    g.into_iter().collect()
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

/// NM_S.1 on the plus strand of NC_S.1 and NM_R.1, the same record, on the
/// minus strand of NC_R.1, the reverse complement of NC_S.1. Both encode
/// NP_S.1 / NP_R.1, `MKLAYR`; NP_T.1 is the same protein served with its stop.
/// NR_S.1 is the record as a non-coding transcript on NC_S.1, and NC_012920.1
/// the plus genome again, as a mitochondrion.
fn provider() -> Provider {
    let cds = Some((5, 25));
    let plus_exons = || {
        vec![
            exon((0, 14), (0, 13), Strand::Plus),
            exon((14, 40), (24, 49), Strand::Plus),
        ]
    };
    Provider::new()
        .sequence("NC_S.1", &plus_genome())
        .sequence("NC_R.1", &revcomp(&plus_genome()))
        .sequence("NC_012920.1", &plus_genome())
        .sequence("NM_S.1", RECORD)
        .sequence("NM_R.1", RECORD)
        .sequence("NR_S.1", RECORD)
        .sequence("NP_S.1", PROTEIN)
        .sequence("NP_R.1", PROTEIN)
        .sequence("NP_T.1", &format!("{PROTEIN}*"))
        .transcript(transcript(
            "NM_S.1",
            "NC_S.1",
            Strand::Plus,
            cds,
            plus_exons(),
        ))
        .transcript(transcript(
            "NR_S.1",
            "NC_S.1",
            Strand::Plus,
            None,
            plus_exons(),
        ))
        .transcript(transcript(
            "NM_R.1",
            "NC_R.1",
            Strand::Minus,
            cds,
            vec![
                exon((0, 14), (36, 49), Strand::Minus),
                exon((14, 40), (0, 25), Strand::Minus),
            ],
        ))
        .protein_for("NM_S.1", "NP_S.1")
        .protein_for("NM_R.1", "NP_R.1")
}

fn parse(s: &str) -> SequenceVariant {
    parse_hgvs_variant(s).unwrap()
}

fn coding(s: &str) -> CVariant {
    match parse(s) {
        SequenceVariant::Coding(c) => c,
        other => panic!("{other} is not c."),
    }
}

fn genomic(s: &str) -> GVariant {
    match parse(s) {
        SequenceVariant::Genomic(g) => g,
        other => panic!("{other} is not g."),
    }
}

fn protein(s: &str) -> PVariant {
    match parse(s) {
        SequenceVariant::Protein(p) => p,
        other => panic!("{other} is not p."),
    }
}

fn rna(s: &str) -> RVariant {
    match parse(s) {
        SequenceVariant::Rna(r) => r,
        other => panic!("{other} is not r."),
    }
}

/// `result` is the mismatch error naming the stated and the reference's bases
/// on `accession`.
#[track_caller]
fn assert_mismatch<T: std::fmt::Debug>(
    result: Result<T, HgvsError>,
    accession: &str,
    stated: &str,
    held: &str,
) {
    match result {
        Err(HgvsError::ReferenceMismatch {
            accession: ac,
            expected,
            found,
            ..
        }) => {
            assert_eq!(
                (ac.as_str(), expected.as_str(), found.as_str()),
                (accession, stated, held)
            );
        }
        other => {
            panic!("expected a mismatch of {stated} against {held} on {accession}, got {other:?}")
        }
    }
}

// --- validate ---

#[test]
fn validate_accepts_an_intronic_base_the_genome_carries() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert!(mapper.validate(&parse("NM_S.1:c.9+1G>A")).unwrap());
    assert!(mapper.validate(&parse("NM_R.1:c.9+1G>A")).unwrap());
}

#[test]
fn validate_refuses_an_intronic_base_the_genome_does_not_carry() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert!(!mapper.validate(&parse("NM_S.1:c.9+1C>A")).unwrap());
    assert!(!mapper.validate(&parse("NM_R.1:c.9+1C>A")).unwrap());
    assert!(!mapper.validate(&parse("NM_S.1:c.10-2_10-1delGG")).unwrap());
    assert!(!mapper.validate(&parse("NR_S.1:n.14+1C>A")).unwrap());
}

#[test]
fn validate_accepts_the_records_base_where_the_genome_differs() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert!(mapper.validate(&parse("NM_S.1:c.12C>A")).unwrap());
}

#[test]
fn validate_refuses_the_genomes_base_where_the_record_differs() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert!(!mapper.validate(&parse("NM_S.1:c.12T>A")).unwrap());
}

#[test]
fn validate_reads_a_range_across_an_exon_boundary_on_the_genome() {
    // The genome is the only sequence holding every base of c.10-1_12, and it
    // has T where the record has C at c.12.
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert!(mapper.validate(&parse("NM_S.1:c.10-1_12delGGCT")).unwrap());
}

#[test]
fn validate_refuses_a_cis_allele_with_a_member_whose_base_is_not_there() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert!(!mapper.validate(&parse("NM_S.1:c.[4A>T;5C>T]")).unwrap());
}

#[test]
fn validate_accepts_a_residue_the_protein_carries() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert!(mapper.validate(&parse("NP_S.1:p.Lys2Arg")).unwrap());
}

#[test]
fn validate_refuses_a_residue_the_protein_does_not_carry() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert!(!mapper.validate(&parse("NP_S.1:p.Arg2Leu")).unwrap());
}

#[test]
fn validate_reads_ter_one_past_a_protein_served_without_its_stop() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert!(mapper.validate(&parse("NP_S.1:p.Ter7Gln")).unwrap());
}

#[test]
fn validate_refuses_ter_past_a_protein_served_with_its_stop() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert!(!mapper.validate(&parse("NP_T.1:p.Ter8Gln")).unwrap());
}

#[test]
fn validate_accepts_a_repeat_unit_the_record_carries() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert!(mapper.validate(&parse("NM_S.1:c.4AAA[2]")).unwrap());
}

#[test]
fn validate_refuses_a_repeat_unit_the_record_does_not_carry() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert!(!mapper.validate(&parse("NM_S.1:c.4CTG[5]")).unwrap());
}

#[test]
fn validate_accepts_an_intronic_repeat_unit_read_on_the_minus_strand() {
    // c.9+1_9+2 is GT on the transcript's strand, AC on NC_R.1.
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert!(mapper.validate(&parse("NM_R.1:c.9+1GT[3]")).unwrap());
}

#[test]
fn validate_refuses_an_intronic_repeat_unit_the_genome_does_not_carry() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert!(!mapper.validate(&parse("NM_R.1:c.9+1AA[3]")).unwrap());
}

#[test]
fn validate_refuses_a_protein_repeat_unit_the_protein_does_not_carry() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert!(!mapper.validate(&parse("NP_S.1:p.2GG[3]")).unwrap());
}

// --- projection ---

#[test]
fn c_to_g_refuses_an_intronic_base_the_genome_does_not_carry() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_mismatch(
        mapper.c_to_g(&coding("NM_S.1:c.9+1C>A"), None),
        "NC_S.1",
        "C",
        "G",
    );
}

#[test]
fn c_to_g_carries_an_intronic_base_on_the_transcripts_strand() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_eq!(
        mapper
            .c_to_g(&coding("NM_R.1:c.9+1G>A"), None)
            .unwrap()
            .to_string(),
        "NC_R.1:g.36C>T"
    );
}

#[test]
fn c_to_g_refuses_an_intronic_base_read_on_the_minus_strand() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    // NC_R.1 holds C at g.36; read on the transcript's strand, a G.
    assert_mismatch(
        mapper.c_to_g(&coding("NM_R.1:c.9+1C>A"), None),
        "NC_R.1",
        "C",
        "G",
    );
}

#[test]
fn c_to_g_places_an_intronic_repeat_unit_on_the_minus_strand_where_the_genome_has_it() {
    // c.9+1_9+2 is GT on the transcript's strand: AC at g.35_36 on NC_R.1,
    // the whole run, written over both its bases.
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let g = mapper.c_to_g(&coding("NM_R.1:c.9+1GT[3]"), None).unwrap();
    assert_eq!(g.to_string(), "NC_R.1:g.35_36AC[3]");
}

#[test]
fn c_to_g_refuses_an_exonic_base_the_record_does_not_carry() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    // The genome's T at c.12 is not the record's base.
    assert_mismatch(
        mapper.c_to_g(&coding("NM_S.1:c.12T>A"), None),
        "NM_S.1",
        "T",
        "C",
    );
}

#[test]
fn c_to_g_refuses_bases_stated_after_del() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_mismatch(
        mapper.c_to_g(&coding("NM_S.1:c.4_5delCC"), None),
        "NM_S.1",
        "CC",
        "AA",
    );
}

#[test]
fn c_to_g_refuses_bases_stated_after_inv() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_mismatch(
        mapper.c_to_g(&coding("NM_S.1:c.4_5invTT"), None),
        "NM_S.1",
        "TT",
        "AA",
    );
}

#[test]
fn c_to_g_refuses_an_identity_whose_base_is_not_there() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_mismatch(
        mapper.c_to_g(&coding("NM_S.1:c.4C="), None),
        "NM_S.1",
        "C",
        "A",
    );
}

#[test]
fn c_to_g_refuses_a_repeat_unit_the_record_does_not_carry() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_mismatch(
        mapper.c_to_g(&coding("NM_S.1:c.4_6CTG[5]"), None),
        "NM_S.1",
        "CTG",
        "AAA",
    );
}

#[test]
fn c_to_g_writes_the_genomes_base_where_the_record_differs() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_eq!(
        mapper
            .c_to_g(&coding("NM_S.1:c.12C>A"), None)
            .unwrap()
            .to_string(),
        "NC_S.1:g.27T>A"
    );
}

fn cis(s: &str) -> hgvs_weaver::CisPhasedVariant {
    match parse(s) {
        SequenceVariant::CisPhased(c) => c,
        other => panic!("{other} is not an allele in cis"),
    }
}

/// `result` is the error of member `index` (0-based), and that error is the
/// mismatch of `stated` against `held` on `accession`.
#[track_caller]
fn assert_member_mismatch<T: std::fmt::Debug>(
    result: Result<T, HgvsError>,
    index: usize,
    accession: &str,
    stated: &str,
    held: &str,
) {
    match result {
        Err(HgvsError::CisMember {
            index: i, source, ..
        }) => {
            assert_eq!(i, index);
            assert_mismatch::<()>(Err(*source), accession, stated, held);
        }
        other => panic!("expected member {index} to fail, got {other:?}"),
    }
}

#[test]
fn c_to_g_cis_refuses_a_member_whose_base_the_record_does_not_carry() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_member_mismatch(
        mapper.c_to_g_cis(&cis("NM_S.1:c.[4A>T;5C>T]"), None),
        1,
        "NM_S.1",
        "C",
        "A",
    );
}

#[test]
fn g_to_c_cis_refuses_a_member_whose_base_the_genome_does_not_carry() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_member_mismatch(
        mapper.g_to_c_cis(&cis("NC_S.1:g.[15C>A;9A>T]"), "NM_S.1"),
        0,
        "NC_S.1",
        "C",
        "G",
    );
}

#[test]
fn n_to_g_refuses_a_base_the_record_does_not_carry() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let SequenceVariant::NonCoding(n) = parse("NR_S.1:n.6C>T") else {
        panic!("not n.");
    };
    assert_mismatch(mapper.n_to_g(&n, None), "NR_S.1", "C", "A");
}

#[test]
fn g_to_c_carries_an_intronic_base_the_genome_carries() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_eq!(
        mapper
            .g_to_c(&genomic("NC_S.1:g.15G>A"), "NM_S.1")
            .unwrap()
            .to_string(),
        "NM_S.1:c.9+1G>A"
    );
}

#[test]
fn g_to_c_refuses_an_intronic_base_the_genome_does_not_carry() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_mismatch(
        mapper.g_to_c(&genomic("NC_S.1:g.15C>A"), "NM_S.1"),
        "NC_S.1",
        "C",
        "G",
    );
}

#[test]
fn g_to_c_refuses_an_exonic_base_the_genome_does_not_carry() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    // The record's C at c.12 is not the genome's base.
    assert_mismatch(
        mapper.g_to_c(&genomic("NC_S.1:g.27C>A"), "NM_S.1"),
        "NC_S.1",
        "C",
        "T",
    );
}

#[test]
fn g_to_c_all_reports_the_mismatch_itself() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_mismatch(
        mapper.g_to_c_all(&genomic("NC_S.1:g.27C>A"), &hdp),
        "NC_S.1",
        "C",
        "T",
    );
}

// --- normalisation ---

#[test]
fn normalize_refuses_a_substitution_whose_base_is_not_there() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_mismatch(
        mapper.normalize_variant(parse("NM_S.1:c.4C>T")),
        "NM_S.1",
        "C",
        "A",
    );
}

#[test]
fn normalize_refuses_a_duplication_whose_base_is_not_there() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_mismatch(
        mapper.normalize_variant(parse("NC_S.1:g.6dupC")),
        "NC_S.1",
        "C",
        "A",
    );
}

#[test]
fn normalize_refuses_a_mitochondrial_base_that_is_not_there() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_mismatch(
        mapper.normalize_variant(parse("NC_012920.1:m.6C>T")),
        "NC_012920.1",
        "C",
        "A",
    );
}

#[test]
fn normalize_refuses_an_intronic_base_the_genome_does_not_carry() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_mismatch(
        mapper.normalize_variant(parse("NM_S.1:c.9+1C>A")),
        "NC_S.1",
        "C",
        "G",
    );
}

#[test]
fn normalize_refuses_a_residue_the_protein_does_not_carry() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_mismatch(
        mapper.normalize_variant(parse("NP_S.1:p.Arg2Leu")),
        "NP_S.1",
        "R",
        "K",
    );
}

// --- alleles, SPDI and VRS ---

#[test]
fn canonical_allele_refuses_a_base_the_genome_does_not_carry() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_mismatch(
        mapper.canonical_allele(&parse("NC_S.1:g.6C>T")),
        "NC_S.1",
        "C",
        "A",
    );
}

#[test]
fn canonical_allele_of_a_coding_variant_reads_the_record() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_mismatch(
        mapper.canonical_allele(&parse("NM_S.1:c.12T>A")),
        "NM_S.1",
        "T",
        "C",
    );
}

#[test]
fn spdi_refuses_a_base_the_genome_does_not_carry() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_mismatch(
        mapper.to_spdi(&parse("NC_S.1:g.6C>T"), false),
        "NC_S.1",
        "C",
        "A",
    );
}

#[test]
fn unambiguous_spdi_refuses_a_base_the_genome_does_not_carry() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_mismatch(
        mapper.to_spdi(&parse("NC_S.1:g.6C>T"), true),
        "NC_S.1",
        "C",
        "A",
    );
}

#[test]
fn from_spdi_refuses_a_deletion_the_sequence_does_not_carry() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_mismatch(mapper.from_spdi("NC_S.1:5:C:T"), "NC_S.1", "C", "A");
}

#[test]
fn from_spdi_reads_a_deletion_in_either_case() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_eq!(
        mapper.from_spdi("NC_S.1:5:a:T").unwrap().to_string(),
        "NC_S.1:g.6A>T"
    );
}

#[test]
fn to_vrs_refuses_a_base_the_genome_does_not_carry() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_mismatch(mapper.to_vrs(&parse("NC_S.1:g.6C>T")), "NC_S.1", "C", "A");
}

#[test]
fn to_vrs_refuses_a_duplication_the_record_does_not_carry() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_mismatch(mapper.to_vrs(&parse("NM_S.1:c.4dupC")), "NM_S.1", "C", "A");
}

#[test]
fn to_vrs_refuses_deleted_bases_of_a_length_insertion_the_genome_does_not_carry() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_mismatch(
        mapper.to_vrs(&parse("NC_S.1:g.6_7delCCinsN[3]")),
        "NC_S.1",
        "CC",
        "AT",
    );
}

#[test]
fn to_vrs_refuses_deleted_bases_of_a_coding_length_insertion_the_record_does_not_carry() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_mismatch(
        mapper.to_vrs(&parse("NM_S.1:c.4_5delCCinsN[3]")),
        "NM_S.1",
        "CC",
        "AA",
    );
}

#[test]
fn to_vrs_copy_number_change_refuses_deleted_bases_the_genome_does_not_carry() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_mismatch(
        mapper.to_vrs_copy_number_change(&parse("NC_S.1:g.6_7delCC")),
        "NC_S.1",
        "CC",
        "AT",
    );
}

#[test]
fn to_vrs_variation_refuses_a_cis_member_whose_base_is_not_there() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_mismatch(
        mapper.to_vrs_variation(&parse("NM_S.1:c.[4A>T;5C>T]")),
        "NM_S.1",
        "C",
        "A",
    );
}

#[test]
fn to_vrs_refuses_a_residue_the_protein_does_not_carry() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_mismatch(
        mapper.to_vrs(&parse("NP_S.1:p.Arg2Leu")),
        "NP_S.1",
        "R",
        "K",
    );
}

#[test]
fn to_vrs_refuses_a_protein_repeat_unit_the_protein_does_not_carry() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_mismatch(
        mapper.to_vrs(&parse("NP_S.1:p.2GG[3]")),
        "NP_S.1",
        "GG",
        "KL",
    );
}

// --- protein consequences ---

#[test]
fn c_to_p_refuses_an_exonic_base_the_record_does_not_carry() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_mismatch(
        mapper.c_to_p(&coding("NM_S.1:c.4C>T"), None),
        "NM_S.1",
        "C",
        "A",
    );
}

#[test]
fn c_to_p_refuses_a_three_prime_utr_base_it_predicts_no_change_from() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_mismatch(
        mapper.c_to_p(&coding("NM_S.1:c.*3A>T"), None),
        "NM_S.1",
        "A",
        "G",
    );
}

#[test]
fn c_to_p_refuses_an_intronic_base_it_predicts_an_unknown_effect_from() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_mismatch(
        mapper.c_to_p(&coding("NM_S.1:c.9+1C>A"), None),
        "NC_S.1",
        "C",
        "G",
    );
}

#[test]
fn c_to_p_of_an_edit_before_the_record_is_refused() {
    // c.-8 lies before the record's first base: there is no sequence to place
    // the change on.
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let result = mapper.c_to_p(&coding("NM_S.1:c.-8del"), None);
    assert!(
        matches!(&result, Err(HgvsError::ValidationError(m)) if m.contains("before transcript start")),
        "{result:?}"
    );
}

#[test]
fn r_to_p_refuses_a_base_the_record_does_not_carry() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_mismatch(
        mapper.r_to_p(&rna("NM_S.1:r.4c>u"), None),
        "NM_S.1",
        "C",
        "A",
    );
}

#[test]
fn protein_allele_refuses_a_base_the_record_does_not_carry() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_mismatch(
        mapper.protein_allele(&coding("NM_S.1:c.4C>T"), None),
        "NM_S.1",
        "C",
        "A",
    );
}

#[test]
fn p_to_c_refuses_a_residue_the_protein_does_not_carry() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_mismatch(
        mapper.p_to_c(&protein("NP_S.1:p.Arg2Leu"), None),
        "NP_S.1",
        "R",
        "K",
    );
}

// --- equivalence ---

#[test]
fn equivalence_refuses_a_base_the_genome_does_not_carry() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let eq = VariantEquivalence::new(&mapper, &hdp);
    assert_mismatch(
        eq.equivalent(&parse("NC_S.1:g.6C>T"), &parse("NM_S.1:c.1A>T")),
        "NC_S.1",
        "C",
        "A",
    );
}

#[test]
fn equivalence_refuses_a_residue_the_protein_does_not_carry() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let eq = VariantEquivalence::new(&mapper, &hdp);
    assert_mismatch(
        eq.equivalent(
            &parse("NP_S.1:p.Arg2Leu"),
            &parse("NP_S.1:p.Lys2_Leu3delinsLeu"),
        ),
        "NP_S.1",
        "R",
        "K",
    );
}
