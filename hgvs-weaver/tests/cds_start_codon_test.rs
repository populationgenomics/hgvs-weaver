//! A CDS numbered from a first codon at which translation does not initiate, as MAVE-HGVS
//! numbers a target that is part of a protein: codon 1 is read by the standard table, and no
//! rule about the start codon applies.

mod support;

use hgvs_weaver::data::{CdsStart, Strand, TranscriptData};
use hgvs_weaver::equivalence::EquivalenceLevel;
use hgvs_weaver::mapper::VariantMapper;
use hgvs_weaver::{parse_hgvs_variant, SequenceVariant, VariantEquivalence};
use support::{exon, transcript, Provider};

const GENOME_AC: &str = "NC_000099.1";

/// Nine codons and a stop, GAT GCT GAA CAA CCA CTT TCT GAA GCT TAA: Asp Ala Glu Gln Pro Leu
/// Ser Glu Ala.
const CDS: &str = "GATGCTGAACAACCACTTTCTGAAGCTTAA";
const UTR3: &str = "GGCCGGCCGGC";

/// - `NM_FRAG.1` is the CDS and its 3'UTR, numbered from GAT, at which translation does not
///   initiate; `NP_FRAG.1` is its translation.
/// - `NM_INIT.1` is the same record with GAT as its start codon; `NP_INIT.1` reads it as Met.
/// - `NM_FRAGUTR.1` is `NM_FRAG.1` with two bases before c.1.
fn provider() -> Provider {
    let record = format!("{CDS}{UTR3}");
    let genome = format!("GG{record}ACGTACGTACGTACGTACGT");
    let mut provider = Provider::new()
        .sequence(GENOME_AC, &genome)
        .identifier_type(
            GENOME_AC,
            hgvs_weaver::data::IdentifierType::GenomicAccession,
        );
    let len = record.len() as i32;
    for (ac, start, utr5, protein) in [
        ("NM_FRAG.1", CdsStart::Internal, "", "DAEQPLSEA"),
        ("NM_INIT.1", CdsStart::Initiator, "", "MAEQPLSEA"),
        ("NM_FRAGUTR.1", CdsStart::Internal, "GG", "DAEQPLSEA"),
    ] {
        let first = utr5.len() as i32;
        let mut model: TranscriptData = transcript(
            ac,
            GENOME_AC,
            Strand::Plus,
            Some((first, first + CDS.len() as i32 - 1)),
            vec![exon((0, first + len), (2 - first, 1 + len), Strand::Plus)],
        );
        model.cds_start = start;
        let np = ac.replace("NM_", "NP_");
        provider = provider
            .sequence(ac, &format!("{utr5}{record}"))
            .sequence(&np, protein)
            .transcript(model)
            .protein_for(ac, &np);
    }
    provider
}

fn c_to_p(hdp: &Provider, s: &str) -> String {
    let mapper = VariantMapper::new(hdp);
    match parse_hgvs_variant(s).unwrap() {
        SequenceVariant::Coding(c) => mapper.c_to_p(&c, None),
        other => panic!("{other} is not c."),
    }
    .map(|p| p.to_string())
    .unwrap_or_else(|e| panic!("{s}: {e}"))
}

#[test]
fn a_substitution_in_codon_1_is_a_substitution_of_its_own_residue() {
    assert_eq!(
        c_to_p(&provider(), "NM_FRAG.1:c.1G>A"),
        "NP_FRAG.1:p.(Asp1Asn)"
    );
}

#[test]
fn a_substitution_in_an_initiator_of_the_same_bases_is_an_unknown_change_from_met1() {
    assert_eq!(
        c_to_p(&provider(), "NM_INIT.1:c.1G>A"),
        "NP_INIT.1:p.(Met1?)"
    );
}

#[test]
fn a_delins_over_codon_1_replaces_its_residue() {
    assert_eq!(
        c_to_p(&provider(), "NM_FRAG.1:c.1_3delinsTGG"),
        "NP_FRAG.1:p.(Asp1Trp)"
    );
}

#[test]
fn a_deletion_of_codon_1_deletes_its_residue() {
    assert_eq!(
        c_to_p(&provider(), "NM_FRAG.1:c.1_3del"),
        "NP_FRAG.1:p.(Asp1del)"
    );
}

#[test]
fn an_insertion_inside_codon_1_before_its_residue_replaces_residue_1_with_both() {
    // GAT made GGTGAT, read on its 3'-most spelling c.1_2insGTG: Gly before Asp.
    assert_eq!(
        c_to_p(&provider(), "NM_FRAG.1:c.1_3delinsGGTGAT"),
        "NP_FRAG.1:p.(Asp1delinsGlyAsp)"
    );
}

#[test]
fn an_insertion_just_before_codon_1_replaces_residue_1_with_both() {
    // GAT made TTTGAT is TTT inserted before c.1, which translation reaches from upstream.
    assert_eq!(
        c_to_p(&provider(), "NM_FRAG.1:c.1_3delinsTTTGAT"),
        "NP_FRAG.1:p.(Asp1delinsPheAsp)"
    );
}

#[test]
fn an_insertion_between_c_minus_1_and_c1_replaces_residue_1_with_both() {
    assert_eq!(
        c_to_p(&provider(), "NM_FRAGUTR.1:c.-1_1insTTT"),
        "NP_FRAGUTR.1:p.(Asp1delinsPheAsp)"
    );
}

#[test]
fn an_insertion_of_one_base_just_before_codon_1_shifts_the_frame_from_residue_1() {
    // T before GAT GCT GAA reads TGA: the new frame stops at once.
    assert_eq!(
        c_to_p(&provider(), "NM_FRAGUTR.1:c.-1_1insT"),
        "NP_FRAGUTR.1:p.(Asp1Ter)"
    );
}

#[test]
fn a_cis_allele_with_an_insertion_just_before_codon_1_reads_it_with_its_other_members() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let SequenceVariant::CisPhased(cis) =
        parse_hgvs_variant("NM_FRAGUTR.1:c.[-1_1insTTT;7G>C]").unwrap()
    else {
        panic!("not a cis allele");
    };
    assert_eq!(
        mapper.c_to_p_cis(&cis, None).unwrap().to_string(),
        "NP_FRAGUTR.1:p.[(Asp1delinsPheAsp;Glu3Gln)]"
    );
}

#[test]
fn a_deletion_of_one_base_in_codon_1_is_a_frameshift_from_residue_1() {
    // ATG CTG AAC: Met then Leu, so the first residue changed is residue 1.
    assert_eq!(
        c_to_p(&provider(), "NM_FRAG.1:c.1del"),
        "NP_FRAG.1:p.(Asp1MetfsTer?)"
    );
}

#[test]
fn a_duplication_of_codon_1_duplicates_its_residue() {
    assert_eq!(
        c_to_p(&provider(), "NM_FRAG.1:c.1_3dup"),
        "NP_FRAG.1:p.(Asp1dup)"
    );
}

#[test]
fn an_identity_over_codon_1_names_its_own_residue_unchanged() {
    assert_eq!(
        c_to_p(&provider(), "NM_FRAG.1:c.1_3="),
        "NP_FRAG.1:p.(Asp1=)"
    );
}

#[test]
fn a_cis_allele_with_a_member_in_codon_1_predicts_each_residue_changed() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let SequenceVariant::CisPhased(cis) = parse_hgvs_variant("NM_FRAG.1:c.[1G>A;7G>C]").unwrap()
    else {
        panic!("not a cis allele");
    };
    assert_eq!(
        mapper.c_to_p_cis(&cis, None).unwrap().to_string(),
        "NP_FRAG.1:p.[(Asp1Asn;Glu3Gln)]"
    );
}

#[test]
fn a_change_in_codon_1_is_identical_to_the_substitution_it_makes() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let parse = |s: &str| parse_hgvs_variant(s).unwrap();
    assert_eq!(
        VariantEquivalence::new(&mapper, &hdp)
            .equivalent_level(&parse("NM_FRAG.1:c.1G>A"), &parse("NP_FRAG.1:p.Asp1Asn"))
            .unwrap(),
        EquivalenceLevel::Identity
    );
}

#[test]
fn a_change_in_codon_1_leaves_the_protein_with_its_residue_replaced() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let SequenceVariant::Coding(c) = parse_hgvs_variant("NM_FRAG.1:c.1G>A").unwrap() else {
        panic!("not c.");
    };
    assert_eq!(
        mapper
            .predicted_protein(&c, None)
            .unwrap()
            .unwrap()
            .residues,
        "NAEQPLSEA"
    );
}

#[test]
fn codon_1_is_numbered_c1_on_the_genome() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let SequenceVariant::Coding(c) = parse_hgvs_variant("NM_FRAGUTR.1:c.1G>A").unwrap() else {
        panic!("not c.");
    };
    assert_eq!(
        mapper.c_to_g(&c, None).unwrap().to_string(),
        "NC_000099.1:g.3G>A"
    );
}

#[test]
fn a_change_wholly_before_codon_1_is_an_unknown_effect() {
    assert_eq!(
        c_to_p(&provider(), "NM_FRAGUTR.1:c.-1G>T"),
        "NP_FRAGUTR.1:p.?"
    );
}

#[test]
fn a_change_from_before_codon_1_into_it_is_an_unknown_change_from_its_residue() {
    assert_eq!(
        c_to_p(&provider(), "NM_FRAGUTR.1:c.-1_2del"),
        "NP_FRAGUTR.1:p.(Asp1?)"
    );
}

#[test]
fn a_cis_allele_with_a_member_from_before_codon_1_into_it_is_an_unknown_change_from_its_residue() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let SequenceVariant::CisPhased(cis) =
        parse_hgvs_variant("NM_FRAGUTR.1:c.[-1_2del;7G>C]").unwrap()
    else {
        panic!("not a cis allele");
    };
    assert_eq!(
        mapper.c_to_p_cis(&cis, None).unwrap().to_string(),
        "NP_FRAGUTR.1:p.(Asp1?)"
    );
}

#[test]
fn unspelled_bases_from_before_codon_1_into_it_are_an_unknown_change_from_its_residue() {
    assert_eq!(
        c_to_p(&provider(), "NM_FRAGUTR.1:c.-1_4delinsN[10]"),
        "NP_FRAGUTR.1:p.(Asp1?)"
    );
}

#[test]
fn an_insertion_of_unspelled_bases_just_before_codon_1_is_an_unknown_change_from_its_residue() {
    assert_eq!(
        c_to_p(&provider(), "NM_FRAGUTR.1:c.-1_1ins(10)"),
        "NP_FRAGUTR.1:p.(Asp1?)"
    );
}

#[test]
fn a_provider_says_translation_does_not_initiate_at_codon_1_with_cds_start_initiates() {
    let model: TranscriptData = serde_json::from_str(
        r#"{"ac":"NM_X.1","gene":"X","cds_start_index":0,"cds_end_index":29,"cds_start_initiates":false,"strand":1,"reference_accession":"NC_000099.1","exons":[]}"#,
    )
    .unwrap();
    assert_eq!(model.cds_start, CdsStart::Internal);
}

#[test]
fn an_open_cds_start_that_initiates_is_refused() {
    let err = serde_json::from_str::<TranscriptData>(
        r#"{"ac":"NM_X.1","gene":"X","cds_start_index":0,"cds_end_index":29,"cds_start_open":true,"cds_start_initiates":true,"strand":1,"reference_accession":"NC_000099.1","exons":[]}"#,
    )
    .unwrap_err();
    assert!(err.to_string().contains("cds_start_initiates"), "{err}");
}

#[test]
fn a_null_cds_start_initiates_is_refused() {
    let err = serde_json::from_str::<TranscriptData>(
        r#"{"ac":"NM_X.1","gene":"X","cds_start_index":0,"cds_end_index":29,"cds_start_initiates":null,"strand":1,"reference_accession":"NC_000099.1","exons":[]}"#,
    )
    .unwrap_err();
    assert!(err.to_string().contains("boolean"), "{err}");
}

#[test]
fn a_model_whose_first_codon_does_not_initiate_reads_back_as_written() {
    let model: TranscriptData = serde_json::from_str(
        r#"{"ac":"NM_X.1","gene":"X","cds_start_index":0,"cds_end_index":29,"cds_start_initiates":false,"strand":1,"reference_accession":"NC_000099.1","exons":[]}"#,
    )
    .unwrap();
    let again: TranscriptData =
        serde_json::from_str(&serde_json::to_string(&model).unwrap()).unwrap();
    assert_eq!(again.cds_start, CdsStart::Internal);
}
