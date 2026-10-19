//! A CDS that runs off an end of the record has no start or stop codon there. A position
//! numbered from that end has no meaning on the record, in either direction, and a protein
//! consequence that depends on it cannot be predicted. Issue #47.

mod support;

use hgvs_weaver::data::{CdsStart, Strand, TranscriptData};
use hgvs_weaver::equivalence::{EquivalenceLevel, VariantEquivalence};
use hgvs_weaver::mapper::VariantMapper;
use hgvs_weaver::parse_hgvs_variant;
use hgvs_weaver::SequenceVariant;
use support::{exon, transcript, Provider};

const GENOME_AC: &str = "NC_000099.1";

/// One exon of 43 bases at genome 0..=42. The record reads GG, then ATG and eight codons, then
/// TAA, then eleven bases of 3'UTR: the CDS is indices 2..=31, ten codons.
///
/// - `NM_WHOLE.1` states the CDS whole.
/// - `NM_OPEN5.1` says the CDS runs off the 5' end: index 2 is its first coding base, not a start.
/// - `NM_OPEN3.1` says the CDS runs off the 3' end: index 31 is its last coding base, not a stop.
fn provider() -> Provider {
    let record = format!(
        "GG{}{}{}{}",
        "ATG", "GCTGAACAACCACTTTCTGAAGCT", "TAA", "GGCCGGCCGGC"
    );
    assert_eq!(record.len(), 43);
    let genome = format!("{record}ACGTACGTACGTACGTACGT");
    let mut provider = Provider::new()
        .sequence(GENOME_AC, &genome)
        .identifier_type(
            GENOME_AC,
            hgvs_weaver::data::IdentifierType::GenomicAccession,
        );
    for (ac, start_open, end_open) in [
        ("NM_WHOLE.1", false, false),
        ("NM_OPEN5.1", true, false),
        ("NM_OPEN3.1", false, true),
    ] {
        let mut model: TranscriptData = transcript(
            ac,
            GENOME_AC,
            Strand::Plus,
            Some((2, 31)),
            vec![exon((0, 43), (0, 42), Strand::Plus)],
        );
        if start_open {
            model.cds_start = CdsStart::Open;
        }
        model.cds_end_open = end_open;
        provider = provider
            .sequence(ac, &record)
            .sequence(&ac.replace("NM_", "NP_"), "MAEQPLSEA")
            .transcript(model)
            .protein_for(ac, &ac.replace("NM_", "NP_"));
    }
    provider
}

fn c_to_g(mapper: &VariantMapper, s: &str) -> Result<String, String> {
    match parse_hgvs_variant(s).unwrap() {
        SequenceVariant::Coding(c) => mapper
            .c_to_g(&c, None)
            .map(|g| g.to_string())
            .map_err(|e| e.to_string()),
        SequenceVariant::NonCoding(n) => mapper
            .n_to_g(&n, None)
            .map(|g| g.to_string())
            .map_err(|e| e.to_string()),
        other => panic!("{other} is not c. or n."),
    }
}

fn g_to_c(mapper: &VariantMapper, s: &str, tx: &str) -> Result<String, String> {
    match parse_hgvs_variant(s).unwrap() {
        SequenceVariant::Genomic(g) => mapper
            .g_to_c(&g, tx)
            .map(|c| c.to_string())
            .map_err(|e| e.to_string()),
        other => panic!("{other} is not g."),
    }
}

fn c_to_p(mapper: &VariantMapper, s: &str) -> Result<String, String> {
    match parse_hgvs_variant(s).unwrap() {
        SequenceVariant::Coding(c) => mapper
            .c_to_p(&c, None)
            .map(|p| p.to_string())
            .map_err(|e| e.to_string()),
        other => panic!("{other} is not c."),
    }
}

#[test]
fn a_whole_cds_numbers_from_both_codons() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_eq!(
        c_to_g(&mapper, "NM_WHOLE.1:c.1A>T").unwrap(),
        "NC_000099.1:g.3A>T"
    );
    assert_eq!(
        c_to_g(&mapper, "NM_WHOLE.1:c.-1G>T").unwrap(),
        "NC_000099.1:g.2G>T"
    );
    assert_eq!(
        c_to_g(&mapper, "NM_WHOLE.1:c.*1G>T").unwrap(),
        "NC_000099.1:g.33G>T"
    );
}

#[test]
fn a_5_open_cds_refuses_every_position_numbered_from_the_start_codon() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    for name in [
        "NM_OPEN5.1:c.1A>T",
        "NM_OPEN5.1:c.-1G>T",
        "NM_OPEN5.1:c.10C>T",
        "NM_OPEN5.1:c.4+2A>T",
    ] {
        let err = c_to_g(&mapper, name).unwrap_err();
        assert!(
            err.contains("open at the 5' end") && err.contains("start codon"),
            "{name}: {err}"
        );
    }
}

#[test]
fn a_5_open_cds_still_numbers_from_the_stop_codon_and_as_n() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_eq!(
        c_to_g(&mapper, "NM_OPEN5.1:c.*1G>T").unwrap(),
        "NC_000099.1:g.33G>T"
    );
    assert_eq!(
        c_to_g(&mapper, "NM_OPEN5.1:n.3A>T").unwrap(),
        "NC_000099.1:g.3A>T"
    );
}

#[test]
fn a_3_open_cds_refuses_every_position_numbered_from_the_stop_codon() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let err = c_to_g(&mapper, "NM_OPEN3.1:c.*1G>T").unwrap_err();
    assert!(
        err.contains("open at the 3' end") && err.contains("stop codon"),
        "{err}"
    );
    assert_eq!(
        c_to_g(&mapper, "NM_OPEN3.1:c.1A>T").unwrap(),
        "NC_000099.1:g.3A>T"
    );
    assert_eq!(
        c_to_g(&mapper, "NM_OPEN3.1:c.-1G>T").unwrap(),
        "NC_000099.1:g.2G>T"
    );
}

#[test]
fn a_genome_position_is_refused_when_its_transcript_position_would_be_numbered_from_an_open_end() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    // genome 33 is the first base after the CDS: c.*1 on a whole CDS, undefined on a 3'-open one
    assert_eq!(
        g_to_c(&mapper, "NC_000099.1:g.33G>T", "NM_WHOLE.1").unwrap(),
        "NM_WHOLE.1:c.*1G>T"
    );
    assert!(g_to_c(&mapper, "NC_000099.1:g.33G>T", "NM_OPEN3.1")
        .unwrap_err()
        .contains("open at the 3' end"));
    assert_eq!(
        g_to_c(&mapper, "NC_000099.1:g.33G>T", "NM_OPEN5.1").unwrap(),
        "NM_OPEN5.1:c.*1G>T"
    );
    // genome 3 is c.1: defined on a 3'-open CDS, undefined on a 5'-open one
    assert_eq!(
        g_to_c(&mapper, "NC_000099.1:g.3A>T", "NM_OPEN3.1").unwrap(),
        "NM_OPEN3.1:c.1A>T"
    );
    assert!(g_to_c(&mapper, "NC_000099.1:g.3A>T", "NM_OPEN5.1")
        .unwrap_err()
        .contains("open at the 5' end"));
}

#[test]
fn a_5_open_cds_predicts_no_protein_consequence() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    for name in ["NM_OPEN5.1:c.4G>A", "NM_OPEN5.1:c.*1G>T"] {
        let err = c_to_p(&mapper, name).unwrap_err();
        assert!(
            err.contains("open at the 5' end") && err.contains("no protein consequence"),
            "{name}: {err}"
        );
    }
}

#[test]
fn a_3_open_cds_predicts_an_in_frame_change_inside_the_coding_bases_it_carries() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    // codon 2, GCT (Ala), becomes GTT (Val); the same on the whole CDS
    assert_eq!(
        c_to_p(&mapper, "NM_OPEN3.1:c.5C>T").unwrap(),
        "NP_OPEN3.1:p.(Ala2Val)"
    );
    assert_eq!(
        c_to_p(&mapper, "NM_WHOLE.1:c.5C>T").unwrap(),
        "NP_WHOLE.1:p.(Ala2Val)"
    );
    // codon 3, GAA (Glu), becomes TAA: a stop the record carries, needing no stop codon of its own
    assert_eq!(
        c_to_p(&mapper, "NM_OPEN3.1:c.7G>T").unwrap(),
        "NP_OPEN3.1:p.(Glu3Ter)"
    );
    // an in-frame deletion of codon 4
    assert_eq!(
        c_to_p(&mapper, "NM_OPEN3.1:c.10_12del").unwrap(),
        "NP_OPEN3.1:p.(Gln4del)"
    );
}

/// A protein fragment, MAEQPLSEAG, as MaveDB serves one: `NM_FRAG.1` is its
/// ten codons and nothing more, its CDS open at the 3' end. `NM_SHUT.1` is the
/// same codons followed by TAA and a 3'UTR, its CDS closed. The last codon,
/// c.28_30, is GGA (Gly), so one base makes it Glu, Gly again, or a stop.
fn fragment() -> Provider {
    let coding = "ATGGCTGAACAACCACTTTCTGAAGCTGGA";
    let shut = format!("{coding}TAAGGCCGG");
    let mut open: TranscriptData = transcript(
        "NM_FRAG.1",
        GENOME_AC,
        Strand::Plus,
        Some((0, 29)),
        vec![exon((0, 30), (0, 29), Strand::Plus)],
    );
    open.cds_end_open = true;
    let closed: TranscriptData = transcript(
        "NM_SHUT.1",
        GENOME_AC,
        Strand::Plus,
        Some((0, 32)),
        vec![exon((0, 39), (0, 38), Strand::Plus)],
    );
    Provider::new()
        .sequence(GENOME_AC, &shut)
        .identifier_type(
            GENOME_AC,
            hgvs_weaver::data::IdentifierType::GenomicAccession,
        )
        .sequence("NM_FRAG.1", coding)
        .sequence("NM_SHUT.1", &shut)
        .sequence("NP_FRAG.1", "MAEQPLSEAG")
        .sequence("NP_SHUT.1", "MAEQPLSEAG")
        .transcript(open)
        .transcript(closed)
        .protein_for("NM_FRAG.1", "NP_FRAG.1")
        .protein_for("NM_SHUT.1", "NP_SHUT.1")
}

fn cis(s: &str) -> hgvs_weaver::CisPhasedVariant {
    match parse_hgvs_variant(s).unwrap() {
        SequenceVariant::CisPhased(c) => c,
        other => panic!("{other} is not a cis allele"),
    }
}

fn c_to_p_cis(mapper: &VariantMapper, s: &str) -> Result<String, String> {
    mapper
        .c_to_p_cis(&cis(s), None)
        .map(|p| p.to_string())
        .map_err(|e| e.to_string())
}

fn refused_as_open(result: Result<String, String>) -> bool {
    result.is_err_and(|e| e.contains("open at the 3' end"))
}

#[test]
fn a_3_open_cds_predicts_a_missense_in_its_last_codon() {
    let hdp = fragment();
    let mapper = VariantMapper::new(&hdp);
    assert_eq!(
        c_to_p(&mapper, "NM_FRAG.1:c.29G>A").unwrap(),
        "NP_FRAG.1:p.(Gly10Glu)"
    );
}

#[test]
fn a_3_open_cds_predicts_a_silent_change_in_its_last_codon() {
    let hdp = fragment();
    let mapper = VariantMapper::new(&hdp);
    assert_eq!(
        c_to_p(&mapper, "NM_FRAG.1:c.30A>G").unwrap(),
        "NP_FRAG.1:p.(Gly10=)"
    );
}

#[test]
fn a_3_open_cds_predicts_a_nonsense_change_in_its_last_codon() {
    let hdp = fragment();
    let mapper = VariantMapper::new(&hdp);
    assert_eq!(
        c_to_p(&mapper, "NM_FRAG.1:c.28G>T").unwrap(),
        "NP_FRAG.1:p.(Gly10Ter)"
    );
}

#[test]
fn a_3_open_cds_predicts_a_delins_inside_its_last_codon() {
    let hdp = fragment();
    let mapper = VariantMapper::new(&hdp);
    // GGA becomes CCA, Pro
    assert_eq!(
        c_to_p(&mapper, "NM_FRAG.1:c.28_29delinsCC").unwrap(),
        "NP_FRAG.1:p.(Gly10Pro)"
    );
}

#[test]
fn a_3_open_cds_refuses_deleting_its_last_residue() {
    let hdp = fragment();
    let mapper = VariantMapper::new(&hdp);
    // the residue after Gly10, which the record does not carry, could be Gly
    assert!(refused_as_open(c_to_p(&mapper, "NM_FRAG.1:c.28_30del")));
}

#[test]
fn a_3_open_cds_predicts_deleting_a_residue_whose_3_prime_spelling_stays_inside_it() {
    let hdp = fragment();
    let mapper = VariantMapper::new(&hdp);
    // Ala9 is followed by Gly10, which differs, so the deletion is Ala9's
    assert_eq!(
        c_to_p(&mapper, "NM_FRAG.1:c.25_27del").unwrap(),
        "NP_FRAG.1:p.(Ala9del)"
    );
}

#[test]
fn a_3_open_cds_refuses_duplicating_its_last_residue() {
    let hdp = fragment();
    let mapper = VariantMapper::new(&hdp);
    assert!(refused_as_open(c_to_p(&mapper, "NM_FRAG.1:c.28_30dup")));
}

#[test]
fn a_3_open_cds_refuses_inserting_before_its_last_residue_a_residue_that_repeats_it() {
    let hdp = fragment();
    let mapper = VariantMapper::new(&hdp);
    // GGG is Gly: shifted 3', the insertion is a duplication of the last residue
    assert!(refused_as_open(c_to_p(&mapper, "NM_FRAG.1:c.27_28insGGG")));
}

#[test]
fn a_3_open_cds_predicts_inserting_before_its_last_residue_a_residue_that_does_not_repeat_it() {
    let hdp = fragment();
    let mapper = VariantMapper::new(&hdp);
    assert_eq!(
        c_to_p(&mapper, "NM_FRAG.1:c.27_28insCCC").unwrap(),
        "NP_FRAG.1:p.(Ala9_Gly10insPro)"
    );
}

#[test]
fn a_3_open_cds_refuses_a_frameshift_that_reads_to_its_end() {
    let hdp = fragment();
    let mapper = VariantMapper::new(&hdp);
    // ATG GTG AAC AAC CAC TTT CTG AAG CTG GA: no stop in the new frame
    assert!(refused_as_open(c_to_p(&mapper, "NM_FRAG.1:c.5del")));
}

#[test]
fn a_3_open_cds_predicts_a_frameshift_that_stops_inside_it() {
    let hdp = fragment();
    let mapper = VariantMapper::new(&hdp);
    // ATG GCT GAA CAA CCC ACT TTC TGA: the new frame stops at codon 8
    assert_eq!(
        c_to_p(&mapper, "NM_FRAG.1:c.14dup").unwrap(),
        "NP_FRAG.1:p.(Leu6ThrfsTer3)"
    );
}

#[test]
fn a_3_open_cds_refuses_a_frameshift_that_leaves_its_last_residue_as_it_was() {
    let hdp = fragment();
    let mapper = VariantMapper::new(&hdp);
    // GGA becomes GGT and an A: Gly10 is unchanged, and the shifted frame
    // reads on past the record
    assert!(refused_as_open(c_to_p(&mapper, "NM_FRAG.1:c.29_30insT")));
}

#[test]
fn a_3_open_cds_refuses_losing_a_stop_its_record_carries_with_no_stop_after() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    // c.28_30 is TAA, the last codon the record carries; TCA reads on past it
    assert!(refused_as_open(c_to_p(&mapper, "NM_OPEN3.1:c.29A>C")));
}

#[test]
fn a_3_open_cds_predicts_a_change_to_its_start_codon_as_unknown_from_met1() {
    let hdp = fragment();
    let mapper = VariantMapper::new(&hdp);
    assert_eq!(
        c_to_p(&mapper, "NM_FRAG.1:c.1A>G").unwrap(),
        "NP_FRAG.1:p.(Met1?)"
    );
}

/// A selenoprotein fragment, MUAEQPLSEG: `NM_SECFRAG.1` is ten codons, the
/// second a TGA read as selenocysteine, its CDS open at the 3' end.
fn selenoprotein_fragment() -> Provider {
    let coding = "ATGTGAGCTGAACAACCACTTTCTGAAGGA";
    let mut model: TranscriptData = transcript(
        "NM_SECFRAG.1",
        GENOME_AC,
        Strand::Plus,
        Some((0, 29)),
        vec![exon((0, 30), (0, 29), Strand::Plus)],
    );
    model.cds_end_open = true;
    Provider::new()
        .sequence(GENOME_AC, coding)
        .sequence("NM_SECFRAG.1", coding)
        .sequence("NP_SECFRAG.1", "MUAEQPLSEG")
        .transcript(model)
        .protein_for("NM_SECFRAG.1", "NP_SECFRAG.1")
}

#[test]
fn a_3_open_cds_reads_an_in_frame_tga_as_selenocysteine_not_its_end() {
    let hdp = selenoprotein_fragment();
    let mapper = VariantMapper::new(&hdp);
    assert_eq!(
        c_to_p(&mapper, "NM_SECFRAG.1:c.8C>T").unwrap(),
        "NP_SECFRAG.1:p.(Ala3Val)"
    );
}

#[test]
fn a_3_open_cds_refuses_a_frameshift_after_a_selenocysteine() {
    let hdp = selenoprotein_fragment();
    let mapper = VariantMapper::new(&hdp);
    // ATG TGA GCT GAA CAA CCA CTT TCT AAG GA: no stop in the new frame
    assert!(refused_as_open(c_to_p(&mapper, "NM_SECFRAG.1:c.20del")));
}

#[test]
fn a_closed_cds_predicts_deleting_its_last_residue() {
    let hdp = fragment();
    let mapper = VariantMapper::new(&hdp);
    assert_eq!(
        c_to_p(&mapper, "NM_SHUT.1:c.28_30del").unwrap(),
        "NP_SHUT.1:p.(Gly10del)"
    );
}

#[test]
fn a_closed_cds_predicts_duplicating_its_last_residue() {
    let hdp = fragment();
    let mapper = VariantMapper::new(&hdp);
    assert_eq!(
        c_to_p(&mapper, "NM_SHUT.1:c.28_30dup").unwrap(),
        "NP_SHUT.1:p.(Gly10dup)"
    );
}

#[test]
fn a_closed_cds_predicts_a_frameshift_that_reads_past_its_codons() {
    let hdp = fragment();
    let mapper = VariantMapper::new(&hdp);
    assert_eq!(
        c_to_p(&mapper, "NM_SHUT.1:c.5del").unwrap(),
        "NP_SHUT.1:p.(Ala2ValfsTer?)"
    );
}

#[test]
fn a_closed_cds_predicts_a_missense_in_its_last_codon() {
    let hdp = fragment();
    let mapper = VariantMapper::new(&hdp);
    assert_eq!(
        c_to_p(&mapper, "NM_SHUT.1:c.29G>A").unwrap(),
        "NP_SHUT.1:p.(Gly10Glu)"
    );
}

#[test]
fn a_closed_cds_predicts_a_nonsense_change_in_its_last_codon() {
    let hdp = fragment();
    let mapper = VariantMapper::new(&hdp);
    assert_eq!(
        c_to_p(&mapper, "NM_SHUT.1:c.28G>T").unwrap(),
        "NP_SHUT.1:p.(Gly10Ter)"
    );
}

#[test]
fn a_closed_cds_predicts_deleting_a_residue_before_its_last() {
    let hdp = fragment();
    let mapper = VariantMapper::new(&hdp);
    assert_eq!(
        c_to_p(&mapper, "NM_SHUT.1:c.25_27del").unwrap(),
        "NP_SHUT.1:p.(Ala9del)"
    );
}

#[test]
fn a_3_open_cds_predicts_an_allele_of_two_in_frame_changes_inside_it() {
    let hdp = fragment();
    let mapper = VariantMapper::new(&hdp);
    assert_eq!(
        c_to_p_cis(&mapper, "NM_FRAG.1:c.[5C>T;13C>T]").unwrap(),
        "NP_FRAG.1:p.[(Ala2Val;Pro5Ser)]"
    );
}

#[test]
fn a_3_open_cds_predicts_an_allele_with_a_member_in_its_last_codon() {
    let hdp = fragment();
    let mapper = VariantMapper::new(&hdp);
    assert_eq!(
        c_to_p_cis(&mapper, "NM_FRAG.1:c.[5C>T;29G>A]").unwrap(),
        "NP_FRAG.1:p.[(Ala2Val;Gly10Glu)]"
    );
}

#[test]
fn a_3_open_cds_predicts_an_allele_with_a_silent_member_in_its_last_codon() {
    let hdp = fragment();
    let mapper = VariantMapper::new(&hdp);
    // the silent member says nothing beside the change
    assert_eq!(
        c_to_p_cis(&mapper, "NM_FRAG.1:c.[5C>T;30A>G]").unwrap(),
        "NP_FRAG.1:p.(Ala2Val)"
    );
}

#[test]
fn a_3_open_cds_refuses_an_allele_that_deletes_its_last_residue() {
    let hdp = fragment();
    let mapper = VariantMapper::new(&hdp);
    assert!(refused_as_open(c_to_p_cis(
        &mapper,
        "NM_FRAG.1:c.[5C>T;28_30del]"
    )));
}

#[test]
fn a_3_open_cds_refuses_an_allele_with_a_frameshift_member_as_it_does_the_member() {
    let hdp = fragment();
    let mapper = VariantMapper::new(&hdp);
    assert!(refused_as_open(c_to_p_cis(
        &mapper,
        "NM_FRAG.1:c.[5del;13C>T]"
    )));
}

#[test]
fn a_3_open_cds_predicts_an_allele_whose_frameshifts_cancel() {
    let hdp = fragment();
    let mapper = VariantMapper::new(&hdp);
    // ATG G|TG AAC A|A|A CCA: the frame is restored before codon 5
    assert_eq!(
        c_to_p_cis(&mapper, "NM_FRAG.1:c.[5del;11_12insA]").unwrap(),
        "NP_FRAG.1:p.(Ala2_Gln4delinsValAsnLys)"
    );
}

#[test]
fn a_closed_cds_predicts_an_allele_with_a_member_in_its_last_codon() {
    let hdp = fragment();
    let mapper = VariantMapper::new(&hdp);
    assert_eq!(
        c_to_p_cis(&mapper, "NM_SHUT.1:c.[5C>T;29G>A]").unwrap(),
        "NP_SHUT.1:p.[(Ala2Val;Gly10Glu)]"
    );
}

#[test]
fn a_closed_cds_predicts_an_allele_that_deletes_its_last_residue() {
    let hdp = fragment();
    let mapper = VariantMapper::new(&hdp);
    assert_eq!(
        c_to_p_cis(&mapper, "NM_SHUT.1:c.[5C>T;28_30del]").unwrap(),
        "NP_SHUT.1:p.[(Ala2Val;Gly10del)]"
    );
}

#[test]
fn a_3_open_cds_predicts_an_allele_whose_changes_meet_once_described() {
    // MAEQPLGGGG: the deletion of codon 7, alone, shifts 3' to the last Gly
    // and would need the residue after it; with Gly10 made Glu the two meet
    // as one delins that does not.
    let record = "ATGGCTGAACAACCACTTGGTGGCGGAGGG";
    let mut model: TranscriptData = transcript(
        "NM_GLYFRAG.1",
        GENOME_AC,
        Strand::Plus,
        Some((0, 29)),
        vec![exon((0, 30), (0, 29), Strand::Plus)],
    );
    model.cds_end_open = true;
    let hdp = Provider::new()
        .sequence(GENOME_AC, record)
        .sequence("NM_GLYFRAG.1", record)
        .sequence("NP_GLYFRAG.1", "MAEQPLGGGG")
        .transcript(model)
        .protein_for("NM_GLYFRAG.1", "NP_GLYFRAG.1");
    let mapper = VariantMapper::new(&hdp);
    assert_eq!(
        c_to_p_cis(&mapper, "NM_GLYFRAG.1:c.[19_21del;29G>A]").unwrap(),
        "NP_GLYFRAG.1:p.(Gly9_Gly10delinsGlu)"
    );
}

#[test]
fn a_3_open_cds_leaves_the_residues_of_an_allele_inside_it() {
    let hdp = fragment();
    let left = VariantMapper::new(&hdp)
        .predicted_protein_cis(&cis("NM_FRAG.1:c.[5C>T;29G>A]"), None)
        .unwrap()
        .expect("an in-frame allele leaves residues");
    assert_eq!((left.residues.as_str(), left.stops), ("MVEQPLSEAE", false));
}

#[test]
fn a_3_open_cds_leaves_no_residues_it_cannot_read_for_an_allele() {
    let hdp = fragment();
    let result =
        VariantMapper::new(&hdp).predicted_protein_cis(&cis("NM_FRAG.1:c.[5C>T;28_30del]"), None);
    assert!(
        matches!(&result, Err(e) if e.to_string().contains("open at the 3' end")),
        "{result:?}"
    );
}

#[test]
fn a_3_open_cds_allele_is_equivalent_to_the_protein_allele_it_predicts() {
    let hdp = fragment();
    let mapper = VariantMapper::new(&hdp);
    let level = VariantEquivalence::new(&mapper, &hdp)
        .equivalent_level(
            &parse_hgvs_variant("NM_FRAG.1:c.[5C>T;29G>A]").unwrap(),
            &parse_hgvs_variant("NP_FRAG.1:p.[Ala2Val;Gly10Glu]").unwrap(),
        )
        .unwrap();
    assert_eq!(level, EquivalenceLevel::Identity);
}

#[test]
fn a_3_open_cds_allele_that_needs_a_base_past_it_is_judged_against_no_protein() {
    let hdp = fragment();
    let mapper = VariantMapper::new(&hdp);
    let result = VariantEquivalence::new(&mapper, &hdp).equivalent_level(
        &parse_hgvs_variant("NM_FRAG.1:c.[5C>T;28_30del]").unwrap(),
        &parse_hgvs_variant("NP_FRAG.1:p.[Ala2Val;Gly10del]").unwrap(),
    );
    assert!(
        matches!(&result, Err(e) if e.to_string().contains("open at the 3' end")),
        "{result:?}"
    );
}

#[test]
fn a_3_open_cds_deleting_a_residue_whose_3_prime_spelling_reaches_its_last_is_refused() {
    // MKLAA, the record ending at the last coding base: deleting codon 4
    // leaves MKLA, and shifted 3' the deletion is of Ala5, the last residue,
    // which the unseen next residue could repeat.
    let record = "ATGAAACTGGCCGCC";
    let mut model: TranscriptData = transcript(
        "NM_OPENEND.1",
        GENOME_AC,
        Strand::Plus,
        Some((0, 14)),
        vec![exon((0, 15), (0, 14), Strand::Plus)],
    );
    model.cds_end_open = true;
    let hdp = Provider::new()
        .sequence(GENOME_AC, record)
        .sequence("NM_OPENEND.1", record)
        .sequence("NP_OPENEND.1", "MKLAA")
        .transcript(model)
        .protein_for("NM_OPENEND.1", "NP_OPENEND.1");
    let mapper = VariantMapper::new(&hdp);
    assert!(refused_as_open(c_to_p(&mapper, "NM_OPENEND.1:c.10_12del")));
}

#[test]
fn a_cds_with_no_stop_shifted_to_the_end_of_its_record_says_nothing() {
    // MKLAA and two bases, no stop, the end not declared open: each deletion
    // shifts the frame and runs to the end of the record, so no codon is left
    // where the protein first changes, and none of them makes a stop.
    let record = "ATGAAACTGGCCGCCGG";
    let model: TranscriptData = transcript(
        "NM_NOSTOP.1",
        GENOME_AC,
        Strand::Plus,
        Some((0, 14)),
        vec![exon((0, 17), (0, 16), Strand::Plus)],
    );
    let hdp = Provider::new()
        .sequence(GENOME_AC, record)
        .sequence("NM_NOSTOP.1", record)
        .sequence("NP_NOSTOP.1", "MKLAA")
        .transcript(model)
        .protein_for("NM_NOSTOP.1", "NP_NOSTOP.1");
    let mapper = VariantMapper::new(&hdp);
    assert_eq!(
        c_to_p(&mapper, "NM_NOSTOP.1:c.10_*2del").unwrap(),
        "NP_NOSTOP.1:p.?"
    );
    assert_eq!(
        c_to_p(&mapper, "NM_NOSTOP.1:c.13_*2del").unwrap(),
        "NP_NOSTOP.1:p.?"
    );
}

#[test]
fn a_provider_that_says_nothing_about_open_ends_means_both_are_whole() {
    let model: TranscriptData = serde_json::from_str(
        r#"{"ac":"NM_X.1","gene":"X","cds_start_index":2,"cds_end_index":31,"strand":1,"reference_accession":"NC_000099.1","exons":[]}"#,
    )
    .unwrap();
    assert!(model.cds_start == CdsStart::Initiator && !model.cds_end_open);
}
