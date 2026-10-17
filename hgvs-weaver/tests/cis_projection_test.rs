//! Projecting an allele in cis between the genome and its transcripts, and
//! between a transcript's `c.`, `n.` and `r.` spellings: member by member,
//! the members reversed across strands, all or nothing.

mod support;

use hgvs_weaver::data::Strand;
use hgvs_weaver::error::HgvsError;
use hgvs_weaver::mapper::VariantMapper;
use hgvs_weaver::utils::reverse_complement;
use hgvs_weaver::{parse_hgvs_variant, CisPhasedVariant, SequenceVariant};
use support::{exon, single_exon_transcript, transcript, Provider};

/// c.1 is index 5 of `NM_X.1`, which is genomic index 10..50 of `NC_X.1`, so
/// c.7 is g.22.
const UTR5: &str = "GGGGG";
/// M K L A Y R *
const CDS: &str = "ATGAAACTGGCCTATCGCTAA";
const UTR3: &str = "CCGTATAAGTAAGG";
/// `NM_S.1` is `NM_X.1` with exon 1 ending at c.5, before this intron.
const INTRON: &str = "GTAAGTCCAG";

fn mrna() -> String {
    format!("{UTR5}{CDS}{UTR3}")
}

/// `NC_X.1`: `NM_X.1` on the plus strand at index 10, and `NR_XR.1`, the
/// same forty bases read on the minus strand, so n.23 is g.28 and n.29 is
/// g.22; `NM_XR.1` is `NR_XR.1` with a CDS from n.4 to n.30.
fn provider() -> Provider {
    let mrna = mrna();
    let genome = format!("TTTTTTTTTT{mrna}CCCCCCCCCC");
    let split = format!("TTTTTTTTTT{}{INTRON}{}CCCCCCCCCC", &mrna[..10], &mrna[10..]);
    let len = mrna.len() as i32;
    let minus = |ac: &str, cds| {
        transcript(
            ac,
            "NC_X.1",
            Strand::Minus,
            cds,
            vec![exon((0, len), (10, 10 + len - 1), Strand::Minus)],
        )
    };
    Provider::new()
        .sequence("NC_X.1", &genome)
        .sequence("NC_S.1", &split)
        .sequence("NM_X.1", &mrna)
        .sequence("NM_S.1", &mrna)
        .sequence("NR_XR.1", &reverse_complement(&mrna))
        .sequence("NM_XR.1", &reverse_complement(&mrna))
        .transcript(single_exon_transcript(
            "NM_X.1",
            "NC_X.1",
            10,
            Strand::Plus,
            UTR5.len() as i32,
            (UTR5.len() + CDS.len()) as i32 - 1,
            len,
        ))
        .transcript(transcript(
            "NM_S.1",
            "NC_S.1",
            Strand::Plus,
            Some((UTR5.len() as i32, (UTR5.len() + CDS.len()) as i32 - 1)),
            vec![
                exon((0, 10), (10, 19), Strand::Plus),
                exon((10, len), (30, 30 + len - 10 - 1), Strand::Plus),
            ],
        ))
        .transcript(minus("NR_XR.1", None))
        // the first ten bases of NC_X.1 on the minus strand, so n.10 is g.1
        .sequence("NR_EDGE.1", "AAAAAAAAAA")
        .transcript(transcript(
            "NR_EDGE.1",
            "NC_X.1",
            Strand::Minus,
            None,
            vec![exon((0, 10), (0, 9), Strand::Minus)],
        ))
        .transcript(minus("NM_XR.1", Some((3, 29))))
}

fn cis(s: &str) -> CisPhasedVariant {
    match parse_hgvs_variant(s).unwrap() {
        SequenceVariant::CisPhased(c) => c,
        other => panic!("{other} is not a cis allele"),
    }
}

#[test]
fn members_project_to_the_genome_one_by_one() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let g = mapper
        .c_to_g_cis(&cis("NM_X.1:c.[7C>T;13T>G]"), None)
        .unwrap();
    assert_eq!(g.to_string(), "NC_X.1:g.[22C>T;28T>G]");
}

#[test]
fn members_project_from_the_genome_one_by_one() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let c = mapper
        .g_to_c_cis(&cis("NC_X.1:g.[22C>T;28T>G]"), "NM_X.1")
        .unwrap();
    assert_eq!(c.to_string(), "NM_X.1:c.[7C>T;13T>G]");
}

#[test]
fn the_minus_strand_reverses_the_members_on_the_genome() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let g = mapper
        .n_to_g_cis(&cis("NR_XR.1:n.[23A>C;29G>A]"), None)
        .unwrap();
    assert_eq!(g.to_string(), "NC_X.1:g.[22C>T;28T>G]");
}

#[test]
fn the_minus_strand_reverses_the_members_on_the_transcript() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    // NM_XR.1's CDS starts at n.4, so n.23 is c.20 and n.29 is c.26
    let c = mapper
        .g_to_c_cis(&cis("NC_X.1:g.[22C>T;28T>G]"), "NM_XR.1")
        .unwrap();
    assert_eq!(c.to_string(), "NM_XR.1:c.[20A>C;26G>A]");
}

#[test]
fn members_in_different_exons_project_to_each() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let g = mapper
        .c_to_g_cis(&cis("NM_S.1:c.[4A>T;6A>T]"), None)
        .unwrap();
    assert_eq!(g.to_string(), "NC_S.1:g.[19A>T;31A>T]");
}

#[test]
fn a_member_in_an_intron_projects_as_given() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    // c.5+2 is the intron's second base, a T, at g.22
    let allele = cis("NM_S.1:c.[5+2T>A;7C>T]");
    let g = mapper.c_to_g_cis(&allele, None).unwrap();
    assert_eq!(g.to_string(), "NC_S.1:g.[22T>A;32C>T]");
    let back = mapper.g_to_c_cis(&g, "NM_S.1").unwrap();
    assert_eq!(back, allele);
}

#[test]
fn a_member_that_fails_to_project_fails_the_allele_and_is_named() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    // r.5_6del spans exon 1's junction, so it has no genomic form
    let err = mapper
        .r_to_g_cis(&cis("NM_S.1:r.[4a>u;5_6del]"), None)
        .unwrap_err();
    let HgvsError::CisMember {
        index,
        member,
        source,
        ..
    } = &err
    else {
        panic!("not a member's error: {err}");
    };
    assert_eq!((*index, member.as_str()), (1, "NM_S.1:r.5_6del"));
    assert!(
        matches!(source.as_ref(), HgvsError::UnsupportedOperation(m) if m.contains("splice junction")),
        "{source}"
    );
}

#[test]
fn a_member_outside_every_exon_fails_the_allele() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let err = mapper
        .c_to_g_cis(&cis("NM_X.1:c.[7C>T;*20del]"), None)
        .unwrap_err();
    assert!(
        matches!(&err, HgvsError::CisMember { index: 1, source, .. }
            if matches!(source.as_ref(), HgvsError::ValidationError(_))),
        "{err}"
    );
}

#[test]
fn a_one_member_allele_stays_an_allele() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let g = mapper.c_to_g_cis(&cis("NM_X.1:c.[7C>T]"), None).unwrap();
    assert_eq!(g.to_string(), "NC_X.1:g.[22C>T]");
}

#[test]
fn a_minus_strand_duplication_alone_stays_a_duplication() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let g = mapper
        .n_to_g_cis(&cis("NR_XR.1:n.[23_25dup;29G>A]"), None)
        .unwrap();
    assert_eq!(g.to_string(), "NC_X.1:g.[22C>T;26_28dup]");
}

#[test]
fn a_minus_strand_duplication_beside_an_insertion_is_written_as_the_insertion() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    // n.23_25 is g.26_28, CCT; the A goes before what is copied
    let g = mapper
        .n_to_g_cis(&cis("NR_XR.1:n.[22_23insA;23_25dup]"), None)
        .unwrap();
    assert_eq!(g.to_string(), "NC_X.1:g.[25_26insCCT;28_29insT]");
}

#[test]
fn a_genomic_duplication_beside_an_insertion_is_written_as_the_insertion_on_a_minus_strand_transcript(
) {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    // g.26_28 is c.20_22 of NM_XR.1, AGG; the A goes before what is copied
    let c = mapper
        .g_to_c_cis(&cis("NC_X.1:g.[25_26insA;26_28dup]"), "NM_XR.1")
        .unwrap();
    assert_eq!(c.to_string(), "NM_XR.1:c.[19_20insAGG;22_23insT]");
}

#[test]
fn a_coding_allele_is_respelled_in_rna_and_back() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let allele = cis("NM_X.1:c.[7C>T;13T>G]");
    let r = mapper.c_to_r_cis(&allele).unwrap();
    assert_eq!(r.to_string(), "NM_X.1:r.[7c>u;13u>g]");
    assert_eq!(mapper.r_to_c_cis(&r).unwrap(), allele);
}

#[test]
fn a_non_coding_allele_is_respelled_in_rna_and_back() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let allele = cis("NR_XR.1:n.[23A>C;29G>A]");
    let r = mapper.n_to_r_cis(&allele).unwrap();
    assert_eq!(r.to_string(), "NR_XR.1:r.[23a>c;29g>a]");
    assert_eq!(mapper.r_to_n_cis(&r).unwrap(), allele);
}

#[test]
fn an_rna_allele_projects_to_the_genome() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let g = mapper
        .r_to_g_cis(&cis("NR_XR.1:r.[23a>c;29g>a]"), None)
        .unwrap();
    assert_eq!(g.to_string(), "NC_X.1:g.[22C>T;28T>G]");
}

#[test]
fn every_transcript_over_the_members_gets_the_allele() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    // NM_S.1 is the one transcript on NC_S.1
    let found: Vec<String> = mapper
        .g_to_c_all_cis(&cis("NC_S.1:g.[19A>T;31A>T]"), &hdp)
        .unwrap()
        .iter()
        .map(ToString::to_string)
        .collect();
    assert_eq!(found, ["NM_S.1:c.[4A>T;6A>T]"]);
}

#[test]
fn an_allele_in_another_system_is_refused() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let err = mapper
        .c_to_g_cis(&cis("NC_X.1:g.[22C>T;28T>G]"), None)
        .unwrap_err();
    assert!(
        matches!(&err, HgvsError::UnsupportedOperation(m) if m.contains("c. members")),
        "{err}"
    );
}

#[test]
fn a_coding_minus_strand_duplication_beside_an_insertion_is_written_as_the_insertion() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    // c.19 is n.22 and c.20_22 is n.23_25 on NM_XR.1
    let g = mapper
        .c_to_g_cis(&cis("NM_XR.1:c.[19_20insA;20_22dup]"), None)
        .unwrap();
    assert_eq!(g.to_string(), "NC_X.1:g.[25_26insCCT;28_29insT]");
}

#[test]
fn an_rna_minus_strand_duplication_beside_an_insertion_is_written_as_the_insertion() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let g = mapper
        .r_to_g_cis(&cis("NR_XR.1:r.[22_23insa;23_25dup]"), None)
        .unwrap();
    assert_eq!(g.to_string(), "NC_X.1:g.[25_26insCCT;28_29insT]");
}

#[test]
fn an_insertion_of_a_length_beside_a_duplication_is_met_like_one_of_bases() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let g = mapper
        .n_to_g_cis(&cis("NR_XR.1:n.[22_23insN[5];23_25dup]"), None)
        .unwrap();
    assert_eq!(g.to_string(), "NC_X.1:g.[25_26insCCT;28_29insN[5]]");
}

#[test]
fn a_rewritten_duplication_can_meet_the_next_whatever_the_order_written() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    // n.26_28 is g.23_25, TGG: its copy would go where n.23_25's now goes
    let g = mapper
        .n_to_g_cis(&cis("NR_XR.1:n.[26_28dup;22_23insA;23_25dup]"), None)
        .unwrap();
    assert_eq!(
        g.to_string(),
        "NC_X.1:g.[25_26insCCT;28_29insT;22_23insTGG]"
    );
}

#[test]
fn a_rewritten_genomic_duplication_can_meet_the_next() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    // g.29_31 is c.17_19 of NM_XR.1, GAT
    let c = mapper
        .g_to_c_cis(&cis("NC_X.1:g.[25_26insA;26_28dup;29_31dup]"), "NM_XR.1")
        .unwrap();
    assert_eq!(
        c.to_string(),
        "NM_XR.1:c.[16_17insGAT;19_20insAGG;22_23insT]"
    );
}

#[test]
fn a_rewritten_duplication_written_in_another_order_meets_the_next() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let g = mapper
        .n_to_g_cis(&cis("NR_XR.1:n.[23_25dup;26_28dup;22_23insA]"), None)
        .unwrap();
    assert_eq!(
        g.to_string(),
        "NC_X.1:g.[28_29insT;22_23insTGG;25_26insCCT]"
    );
}

#[test]
fn an_uncertain_duplication_beside_an_insertion_across_strands_is_refused() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let err = mapper
        .g_to_c_cis(&cis("NC_X.1:g.[25_26insA;(26_28)dup]"), "NM_XR.1")
        .unwrap_err();
    assert!(
        matches!(&err, HgvsError::CisMember { index: 1, source, .. }
            if matches!(source.as_ref(), HgvsError::UnsupportedOperation(m) if m.contains("uncertain"))),
        "{err}"
    );
}

#[test]
fn an_uncertain_transcript_duplication_beside_an_insertion_is_refused() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let err = mapper
        .n_to_g_cis(&cis("NR_XR.1:n.[22_23insA;(23_25)dup]"), None)
        .unwrap_err();
    assert!(
        matches!(&err, HgvsError::CisMember { index: 1, source, .. }
            if matches!(source.as_ref(), HgvsError::UnsupportedOperation(m) if m.contains("uncertain"))),
        "{err}"
    );
}

#[test]
fn a_copy_of_the_genomes_first_base_has_no_junction_before_it() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    // n.9_10 is g.1_2; the insertion is between g.2 and g.3
    let err = mapper
        .n_to_g_cis(&cis("NR_EDGE.1:n.[8_9insG;9_10dup]"), None)
        .unwrap_err();
    assert!(
        matches!(&err, HgvsError::CisMember { index: 1, source, .. }
            if matches!(source.as_ref(), HgvsError::ValidationError(m) if m.contains("first base"))),
        "{err}"
    );
}
