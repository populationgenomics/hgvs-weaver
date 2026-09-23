//! Regression tests for c./n. -> genomic coordinate resolution on the paths
//! that used to route through `DataProvider::c_to_g`.
//!
//! Those paths (SPDI intervals and c.-vs-c. equivalence) handed the provider a
//! bare index with no anchor, so `c.*1` was indistinguishable from `c.1` and no
//! adapter could get the minus strand or a non-zero CDS start right. The
//! fixtures here deliberately use a CDS that does not start at index 0 and a
//! minus-strand transcript.

mod support;

use hgvs_weaver::data::{Strand, TranscriptData};
use hgvs_weaver::equivalence::{EquivalenceLevel, VariantEquivalence};
use hgvs_weaver::error::HgvsError;
use hgvs_weaver::mapper::VariantMapper;
use hgvs_weaver::parse_hgvs_variant;
use hgvs_weaver::structs::IntervalSpdi;
use hgvs_weaver::SequenceVariant;
use support::{single_exon_transcript, Provider};

const REF_AC: &str = "NC_TEST.1";
fn genome() -> String {
    "ACGT".repeat(500)
}

/// A second reference: 5 T, three copies of GCCATT, 5 A, then C to 100 bases.
/// GCCATT is not its own reverse complement (that is AATGGC), so a repeat on
/// the minus strand reads differently from the genome.
const REP_AC: &str = "NC_REP.1";
fn rep_genome() -> String {
    format!("TTTTT{}AAAAA{}", "GCCATT".repeat(3), "C".repeat(72))
}

/// One exon: transcript indices 0..=99 <-> genomic indices 1000..=1099.
fn transcript(ac: &str, strand: Strand, cds_start: i32, cds_end: i32) -> TranscriptData {
    transcript_on(ac, REF_AC, 1000, strand, cds_start, cds_end)
}

/// One 100-base exon on `reference` starting at genomic index `g0`.
fn transcript_on(
    ac: &str,
    reference: &str,
    g0: i32,
    strand: Strand,
    cds_start: i32,
    cds_end: i32,
) -> TranscriptData {
    single_exon_transcript(ac, reference, g0, strand, cds_start, cds_end, 100)
}

/// The bases of `transcript`'s one exon, read from its genome on its strand.
fn transcript_sequence(transcript: &TranscriptData) -> String {
    let genome = if transcript.reference_accession == REP_AC {
        rep_genome()
    } else {
        genome()
    };
    let exon = &transcript.exons[0];
    let exonic = &genome[exon.reference_start.0 as usize..=exon.reference_end.0 as usize];
    match transcript.strand {
        Strand::Plus => exonic.to_string(),
        Strand::Minus => exonic
            .chars()
            .rev()
            .map(|c| match c {
                'A' => 'T',
                'C' => 'G',
                'G' => 'C',
                'T' => 'A',
                other => other,
            })
            .collect(),
    }
}

fn provider() -> Provider {
    let transcripts = [
        // CDS at transcript indices 10..=39 on the plus strand.
        transcript("NM_PLUS10.1", Strand::Plus, 10, 39),
        // CDS spanning the whole transcript on the plus strand.
        transcript("NM_PLUS0.1", Strand::Plus, 0, 99),
        // CDS at transcript indices 10..=39 on the minus strand.
        transcript("NM_MINUS10.1", Strand::Minus, 10, 39),
        // The whole of NC_REP.1, read on the minus strand, CDS the whole transcript.
        transcript_on("NM_REP_MINUS.1", REP_AC, 0, Strand::Minus, 0, 99),
    ];
    let mut provider = Provider::new()
        .sequence(REF_AC, &genome())
        .sequence(REP_AC, &rep_genome());
    for transcript in transcripts {
        provider = provider
            .sequence(&transcript.ac, &transcript_sequence(&transcript))
            .transcript(transcript);
    }
    provider
}

fn coding_interval(hgvs: &str) -> hgvs_weaver::structs::BaseOffsetInterval {
    match parse_hgvs_variant(hgvs).unwrap() {
        SequenceVariant::Coding(v) => v.posedit.pos.unwrap(),
        other => panic!("expected a c. variant, got {:?}", other),
    }
}

#[test]
fn spdi_interval_honours_cds_end_anchor() {
    // c.*1 is the base after the stop codon: transcript index 40, genomic 1040.
    let iv = coding_interval("NM_PLUS10.1:c.*1A>G");
    let got = iv.spdi_interval("NM_PLUS10.1", &provider()).unwrap();
    assert_eq!(got, (1040, 1041, REF_AC.to_string()));
}

#[test]
fn spdi_interval_honours_non_zero_cds_start() {
    // c.1 is transcript index 10, genomic 1010.
    let iv = coding_interval("NM_PLUS10.1:c.1A>G");
    let got = iv.spdi_interval("NM_PLUS10.1", &provider()).unwrap();
    assert_eq!(got, (1010, 1011, REF_AC.to_string()));
}

#[test]
fn spdi_interval_on_minus_strand() {
    // Minus strand: transcript index 0 is genomic 1099, so c.1 (index 10) is 1089.
    let iv = coding_interval("NM_MINUS10.1:c.1A>G");
    let got = iv.spdi_interval("NM_MINUS10.1", &provider()).unwrap();
    assert_eq!(got, (1089, 1090, REF_AC.to_string()));

    // A multi-base interval is reported low-to-high on the genome.
    let iv = coding_interval("NM_MINUS10.1:c.1_3del");
    let got = iv.spdi_interval("NM_MINUS10.1", &provider()).unwrap();
    assert_eq!(got, (1087, 1090, REF_AC.to_string()));
}

#[test]
fn spdi_interval_applies_intronic_offset_by_strand() {
    // On the minus strand a +5 intronic offset moves towards lower genomic indices.
    let iv = coding_interval("NM_MINUS10.1:c.1+5A>G");
    let got = iv.spdi_interval("NM_MINUS10.1", &provider()).unwrap();
    assert_eq!(got, (1084, 1085, REF_AC.to_string()));
}

#[test]
fn coding_variants_on_transcripts_with_different_cds_starts_are_equivalent() {
    // Both name genomic index 1040 (reference base A): c.*1 on the 10..=39 CDS
    // and c.41 on the whole-transcript CDS.
    let v1 = parse_hgvs_variant("NM_PLUS10.1:c.*1A>G").unwrap();
    let v2 = parse_hgvs_variant("NM_PLUS0.1:c.41A>G").unwrap();
    let hdp = provider();
    let eq_mapper = VariantMapper::new(&hdp);
    let eq = VariantEquivalence::new(&eq_mapper, &hdp);
    assert_eq!(
        eq.equivalent_level(&v1, &v2).unwrap(),
        EquivalenceLevel::Analogous
    );

    // And a different base is not.
    let v3 = parse_hgvs_variant("NM_PLUS0.1:c.42C>G").unwrap();
    assert_eq!(
        eq.equivalent_level(&v1, &v3).unwrap(),
        EquivalenceLevel::Different
    );
}

#[test]
fn validate_checks_stated_reference_through_transcript_coordinates() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let ok = |hgvs: &str| mapper.validate(&parse_hgvs_variant(hgvs).unwrap()).unwrap();

    // Plus strand, CDS at index 10: c.1 is transcript index 10, genome[1010] = G.
    assert!(ok("NM_PLUS10.1:c.1G>A"));
    assert!(!ok("NM_PLUS10.1:c.1A>G"));
    // c.*1 is transcript index 40, genome[1040] = A.
    assert!(ok("NM_PLUS10.1:c.*1A>G"));
    // Minus strand: transcript index 10 is the complement of genome[1089] = C.
    assert!(ok("NM_MINUS10.1:c.1G>A"));
    assert!(!ok("NM_MINUS10.1:c.1C>A"));
    // Intronic positions are accepted unchecked; genomic goes straight to the reference.
    assert!(ok("NM_MINUS10.1:c.1+5T>A"));
    assert!(ok("NC_TEST.1:g.1011G>A"));
    assert!(!ok("NC_TEST.1:g.1011A>G"));
}

#[test]
fn normalize_leaves_intronic_coding_variants_as_written() {
    // An insertion straddling an exon boundary has no transcript-space
    // normalisation. It must come back unchanged, not as an error, so that
    // c_to_p (p.?) and to_spdi (via g.) still run on it.
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    for hgvs in [
        "NM_PLUS10.1:c.30_30+1insA",
        "NM_MINUS10.1:c.1+5del",
        "NM_PLUS10.1:n.5-2_5del",
    ] {
        let v = parse_hgvs_variant(hgvs).unwrap();
        assert_eq!(mapper.normalize_variant(v).unwrap().to_string(), hgvs);
    }
}

#[test]
fn normalize_converts_genomic_and_noncoding_insertions_to_duplications() {
    // Reference is ACGT repeated: genomic index 1010 is G, 1011 is T.
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let norm = |hgvs: &str| {
        mapper
            .normalize_variant(parse_hgvs_variant(hgvs).unwrap())
            .unwrap()
            .to_string()
    };
    // Inserting T after g.1012 (index 1011, a T): the run of one T ends there, so it is a dup.
    assert_eq!(norm("NC_TEST.1:g.1012_1013insT"), "NC_TEST.1:g.1012dup");
    // Inserting ACGT into the repeat shifts to the end of the reference and duplicates.
    assert_eq!(
        norm("NC_TEST.1:g.1012_1013insACGT"),
        "NC_TEST.1:g.1997_2000dup"
    );
    // Non-coding on the plus strand: transcript index 11 is genomic 1011 (T).
    assert_eq!(norm("NM_PLUS0.1:n.12_13insT"), "NM_PLUS0.1:n.12dup");
    // An insertion that repeats nothing stays an insertion; here it slides one
    // base 3' because inserting AA before an A equals inserting it after.
    assert_eq!(
        norm("NC_TEST.1:g.1012_1013insAA"),
        "NC_TEST.1:g.1013_1014insAA"
    );
    assert_eq!(
        norm("NC_TEST.1:g.1012_1013insCC"),
        "NC_TEST.1:g.1012_1013insCC"
    );
}

#[test]
fn plain_spdi_places_an_insertion_at_its_second_flank() {
    // SPDI counts the bases before the change, so an insertion between
    // g.1012 and g.1013 sits at 0-based position 1012, like a duplication of
    // g.1012 does. Reference here: index 1011 is T, 1012 is A.
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let spdi = |hgvs: &str| {
        mapper
            .to_spdi(&parse_hgvs_variant(hgvs).unwrap(), false)
            .unwrap()
    };
    assert_eq!(spdi("NC_TEST.1:g.1012_1013insCC"), "NC_TEST.1:1012::CC");
    assert_eq!(spdi("NC_TEST.1:g.1012dupT"), "NC_TEST.1:1012::T");
    assert_eq!(spdi("NC_TEST.1:g.1011_1012del"), "NC_TEST.1:1010:GT:");
    assert_eq!(spdi("NC_TEST.1:g.1013A>G"), "NC_TEST.1:1012:A:G");
}

#[test]
fn mitochondrial_variants_share_the_genomic_implementation() {
    // m. is g. on the mitochondrial reference: same normalisation, SPDI,
    // validation and equivalence, written back with the m. letter.
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let norm = |hgvs: &str| {
        mapper
            .normalize_variant(parse_hgvs_variant(hgvs).unwrap())
            .unwrap()
            .to_string()
    };
    assert_eq!(norm("NC_TEST.1:m.1012_1013insT"), "NC_TEST.1:m.1012dup");
    assert_eq!(norm("NC_TEST.1:m.1008_1010del"), "NC_TEST.1:m.1008_1010del");
    // Stated bases stay when the edit does not move.
    assert_eq!(
        norm("NC_TEST.1:m.1008_1010delTAC"),
        "NC_TEST.1:m.1008_1010delTAC"
    );

    let spdi = |hgvs: &str| {
        mapper
            .to_spdi(&parse_hgvs_variant(hgvs).unwrap(), true)
            .unwrap()
    };
    assert_eq!(spdi("NC_TEST.1:m.1013A>G"), spdi("NC_TEST.1:g.1013A>G"));

    let valid = |hgvs: &str| mapper.validate(&parse_hgvs_variant(hgvs).unwrap()).unwrap();
    assert!(valid("NC_TEST.1:m.1013A>G"));
    assert!(!valid("NC_TEST.1:m.1013C>G"));

    let eq_mapper = VariantMapper::new(&hdp);

    let eq = VariantEquivalence::new(&eq_mapper, &hdp);
    let m = parse_hgvs_variant("NC_TEST.1:m.1013A>G").unwrap();
    let g = parse_hgvs_variant("NC_TEST.1:g.1013A>G").unwrap();
    assert_eq!(
        eq.equivalent_level(&m, &g).unwrap(),
        EquivalenceLevel::Identity
    );
}

#[test]
fn canonical_alleles_make_spdi_vrs_and_equivalence_one_value() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let parse = |h: &str| parse_hgvs_variant(h).unwrap();
    // Two spellings of one change in the ACGT repeat: inserting ACGT anywhere
    // in the run, or duplicating any copy, is the same allele.
    let a = mapper
        .canonical_allele(&parse("NC_TEST.1:g.1012_1013insACGT"))
        .unwrap();
    let b = mapper
        .canonical_allele(&parse("NC_TEST.1:g.1005_1008dup"))
        .unwrap();
    let c = mapper
        .canonical_allele(&parse("NM_PLUS0.1:c.12_13insACGT"))
        .unwrap();
    assert_eq!(a, b);
    assert_eq!(a, c);
    // The run is the whole 2000-base repeat, so the allele spans all of it.
    assert_eq!((a.start, a.end), (0, 2000));
    assert_eq!(a.repeat_subunit, Some(4));
    assert_eq!(
        mapper
            .to_spdi_unambiguous(&parse("NC_TEST.1:g.1005_1008dup"))
            .unwrap(),
        a.spdi()
    );

    let vrs_a = mapper
        .to_vrs(&parse("NC_TEST.1:g.1012_1013insACGT"))
        .unwrap();
    let vrs_c = mapper.to_vrs(&parse("NM_PLUS0.1:c.12_13insACGT")).unwrap();
    assert_eq!(vrs_a.id, vrs_c.id);
    assert!(vrs_a.id.starts_with("ga4gh:VA."));
    // No provider hook here, so the refget accession is computed from the sequence.
    assert_eq!(
        vrs_a.location.sequence_reference.refget_accession,
        hgvs_weaver::vrs::refget_accession(&genome())
    );
    assert_eq!(vrs_a.expressions[0].syntax, "hgvs.g");
    assert_eq!(vrs_c.expressions[0].value, "NM_PLUS0.1:c.12_13insACGT");

    // Equivalence rides on the same value.
    let eq_mapper = VariantMapper::new(&hdp);
    let eq = VariantEquivalence::new(&eq_mapper, &hdp);
    assert!(eq
        .equivalent(
            &parse("NC_TEST.1:g.1012_1013insACGT"),
            &parse("NC_TEST.1:g.1005_1008dup")
        )
        .unwrap());
    // A substitution is trimmed to the base that changes and is not a repeat.
    let s = mapper
        .canonical_allele(&parse("NC_TEST.1:g.1011_1013delGTAinsGTC"))
        .unwrap();
    assert_eq!(s.spdi(), "NC_TEST.1:1012:A:C");
    assert_eq!(s.repeat_subunit, None);
}

#[test]
fn repeat_on_the_minus_strand_projects_to_its_whole_run() {
    // On the transcript (reverse complement of NC_REP.1), the GCCATT run reads
    // as AATGGC x3 starting at transcript index 77, i.e. c.78.
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let parse = |h: &str| parse_hgvs_variant(h).unwrap();
    let SequenceVariant::Coding(c) = parse("NM_REP_MINUS.1:c.78AATGGC[4]") else {
        panic!()
    };
    // Projected to the genome the repeat covers the whole run, indices 5..=22,
    // written in plus-strand orientation.
    let g = mapper.c_to_g(&c, None).unwrap();
    assert_eq!(g.to_string(), "NC_REP.1:g.6_23GCCATT[4]");

    // One more copy, however it is written, is one allele.
    let expand = mapper
        .canonical_allele(&parse("NM_REP_MINUS.1:c.78AATGGC[4]"))
        .unwrap();
    let genomic = mapper
        .canonical_allele(&parse("NC_REP.1:g.6GCCATT[4]"))
        .unwrap();
    let dup = mapper
        .canonical_allele(&parse("NC_REP.1:g.18_23dup"))
        .unwrap();
    let ins = mapper
        .canonical_allele(&parse("NC_REP.1:g.23_24insGCCATT"))
        .unwrap();
    assert_eq!(expand, genomic);
    assert_eq!(expand, dup);
    assert_eq!(expand, ins);
    assert_eq!(expand.repeat_subunit, Some(6));
    assert_eq!(expand.alternate.len() - expand.reference.len(), 6);

    // And the reverse projection gives back the run on the transcript.
    let SequenceVariant::Genomic(gv) = parse("NC_REP.1:g.6GCCATT[4]") else {
        panic!()
    };
    let back = mapper.g_to_c(&gv, "NM_REP_MINUS.1").unwrap();
    assert_eq!(back.to_string(), "NM_REP_MINUS.1:c.78_95AATGGC[4]");
}

#[test]
fn a_position_in_no_exon_is_an_error_not_an_extrapolation() {
    // NM_PLUS10.1 is one exon of 100 bases and c.*60 is its last base, so
    // c.*61 lies in no exon. It is refused wherever it is resolved, never
    // projected by extending the exon past its end.
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let var = parse_hgvs_variant("NM_PLUS10.1:c.*61A>G").unwrap();
    let SequenceVariant::Coding(c) = &var else {
        panic!()
    };
    for err in [
        mapper.c_to_g(c, None).unwrap_err(),
        mapper.validate(&var).unwrap_err(),
        mapper.to_spdi_unambiguous(&var).unwrap_err(),
    ] {
        assert!(matches!(err, HgvsError::ValidationError(_)), "{err}");
    }
    // The last base is fine.
    assert!(mapper
        .c_to_g(&coding_variant("NM_PLUS10.1:c.*60T>A"), None)
        .is_ok());
}

fn coding_variant(hgvs: &str) -> hgvs_weaver::structs::CVariant {
    match parse_hgvs_variant(hgvs).unwrap() {
        SequenceVariant::Coding(c) => c,
        other => panic!("{other} is not c."),
    }
}
