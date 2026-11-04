//! A transcript record its provider places on no genome, as a GenBank mRNA
//! is served: read on its own sequence, and refused wherever an operation
//! needs the genome.

mod support;

use hgvs_weaver::coords::TranscriptPos;
use hgvs_weaver::data::{IdentifierType, Strand, TranscriptData};
use hgvs_weaver::equivalence::VariantEquivalence;
use hgvs_weaver::error::HgvsError;
use hgvs_weaver::mapper::VariantMapper;
use hgvs_weaver::{parse_hgvs_variant, SequenceVariant};
use support::{exon, transcript, Provider};

/// `NM_U.1`, `ATGAA|AGGCCTGTAACC` with the CDS at transcript indices 0..=14,
/// encoding `NP_U.1`, Met-Lys-Gly-Leu. Placed, its junction would fall
/// between c.5 and c.6, so c.4, c.5 and c.6, all A, are a run across it.
fn unplaced() -> Provider {
    Provider::new()
        .sequence("NM_U.1", "ATGAAAGGCCTGTAACC")
        .transcript(TranscriptData::unplaced(
            "NM_U.1",
            "UNPLACED",
            Some((TranscriptPos(0), TranscriptPos(14))),
        ))
        .sequence("NP_U.1", "MKGL")
        .protein_for("NM_U.1", "NP_U.1")
}

/// `NM_U.1` beside `NM_J.1`, the same sequence placed on `NC_J.1` as two
/// exons, c.1 to c.5 at g.16 to g.20 and c.6 to c.*2 at g.31 to g.42, the
/// intron between them `GTCCCCCCAG`; both are transcripts of the gene `TEST`.
fn beside_a_placed_transcript() -> Provider {
    let record = "ATGAAAGGCCTGTAACC";
    let flank = "CCCCCCCCCCCCCCC";
    Provider::new()
        .sequence(
            "NC_J.1",
            &format!("{flank}ATGAAGTCCCCCCAGAGGCCTGTAACC{flank}"),
        )
        .sequence("NM_J.1", record)
        .transcript(transcript(
            "NM_J.1",
            "NC_J.1",
            Strand::Plus,
            Some((0, 14)),
            vec![
                exon((0, 5), (15, 19), Strand::Plus),
                exon((5, 17), (30, 41), Strand::Plus),
            ],
        ))
        .sequence("NM_U.1", record)
        .transcript(TranscriptData::unplaced(
            "NM_U.1",
            "TEST",
            Some((TranscriptPos(0), TranscriptPos(14))),
        ))
        .identifier_type("TEST", IdentifierType::GeneSymbol)
}

fn parse(s: &str) -> SequenceVariant {
    parse_hgvs_variant(s).unwrap_or_else(|e| panic!("{s}: {e}"))
}

fn coding(s: &str) -> hgvs_weaver::structs::CVariant {
    match parse(s) {
        SequenceVariant::Coding(c) => c,
        other => panic!("{other} is not c."),
    }
}

/// Whether `result` is the refusal of `NM_U.1` for having no placement.
fn unplaced_refusal<T: std::fmt::Debug>(result: &Result<T, HgvsError>) -> bool {
    matches!(result, Err(HgvsError::Unplaced { accession }) if accession == "NM_U.1")
}

#[test]
fn an_exonic_change_is_validated_on_the_records_bases() {
    let hdp = unplaced();

    assert!(VariantMapper::new(&hdp)
        .validate(&parse("NM_U.1:c.7G>T"))
        .unwrap());
}

#[test]
fn a_base_the_record_does_not_hold_is_invalid() {
    let hdp = unplaced();

    assert!(!VariantMapper::new(&hdp)
        .validate(&parse("NM_U.1:c.7C>T"))
        .unwrap());
}

#[test]
fn two_spellings_on_the_record_are_not_compared_without_the_genome() {
    let hdp = unplaced();
    let mapper = VariantMapper::new(&hdp);
    let equivalence = VariantEquivalence::new(&mapper, &hdp);

    let result = equivalence.equivalent_level(&parse("NM_U.1:c.4del"), &parse("NM_U.1:c.6del"));

    assert!(unplaced_refusal(&result), "{result:?}");
}

#[test]
fn a_deletion_shifts_across_where_a_junction_would_be() {
    let hdp = unplaced();

    let normalized = VariantMapper::new(&hdp)
        .normalize_variant(parse("NM_U.1:c.4del"))
        .unwrap();

    assert_eq!(normalized.to_string(), "NM_U.1:c.6del");
}

#[test]
fn an_rna_deletion_shifts_across_where_a_junction_would_be() {
    let hdp = unplaced();

    let normalized = VariantMapper::new(&hdp)
        .normalize_variant(parse("NM_U.1:r.4del"))
        .unwrap();

    assert_eq!(normalized.to_string(), "NM_U.1:r.6del");
}

#[test]
fn a_repeat_at_the_records_count_normalises_on_the_record() {
    let hdp = unplaced();

    let normalized = VariantMapper::new(&hdp)
        .normalize_variant(parse("NM_U.1:c.4A[3]"))
        .unwrap();

    assert_eq!(normalized.to_string(), "NM_U.1:c.4_6=");
}

#[test]
fn a_position_past_the_records_end_is_outside_it_not_unplaced() {
    let hdp = unplaced();

    let result = VariantMapper::new(&hdp).normalize_variant(parse("NM_U.1:c.*5del"));

    assert!(
        matches!(result, Err(HgvsError::PositionOutsideSequence { .. })),
        "{result:?}"
    );
}

#[test]
fn a_coding_change_has_its_protein_consequence() {
    let hdp = unplaced();

    let p = VariantMapper::new(&hdp)
        .c_to_p(&coding("NM_U.1:c.5A>G"), None)
        .unwrap();

    assert_eq!(p.to_string(), "NP_U.1:p.(Lys2Arg)");
}

#[test]
fn a_coding_change_has_its_rna_form() {
    let hdp = unplaced();

    let r = VariantMapper::new(&hdp)
        .c_to_r(&coding("NM_U.1:c.5A>G"))
        .unwrap();

    assert_eq!(r.to_string(), "NM_U.1:r.5a>g");
}

#[test]
fn an_intronic_position_is_refused_for_want_of_a_placement() {
    let hdp = unplaced();

    let result = VariantMapper::new(&hdp).validate(&parse("NM_U.1:c.5+1G>A"));

    assert!(unplaced_refusal(&result), "{result:?}");
}

#[test]
fn an_intronic_position_is_not_normalised() {
    let hdp = unplaced();

    let result = VariantMapper::new(&hdp).normalize_variant(parse("NM_U.1:c.5+1del"));

    assert!(unplaced_refusal(&result), "{result:?}");
}

#[test]
fn an_intronic_position_has_no_protein_consequence() {
    let hdp = unplaced();

    let result = VariantMapper::new(&hdp).c_to_p(&coding("NM_U.1:c.5+1G>A"), None);

    assert!(unplaced_refusal(&result), "{result:?}");
}

#[test]
fn a_coding_change_has_no_genomic_form() {
    let hdp = unplaced();

    let result = VariantMapper::new(&hdp).c_to_g(&coding("NM_U.1:c.5A>G"), None);

    assert!(unplaced_refusal(&result), "{result:?}");
}

#[test]
fn an_rna_change_has_no_genomic_form() {
    let hdp = unplaced();
    let r = match parse("NM_U.1:r.5a>g") {
        SequenceVariant::Rna(r) => r,
        other => panic!("{other} is not r."),
    };

    let result = VariantMapper::new(&hdp).r_to_g(&r, None);

    assert!(unplaced_refusal(&result), "{result:?}");
}

#[test]
fn a_coding_change_has_no_spdi() {
    let hdp = unplaced();

    let result = VariantMapper::new(&hdp).to_spdi(&parse("NM_U.1:c.5A>G"), false);

    assert!(unplaced_refusal(&result), "{result:?}");
}

#[test]
fn a_gene_symbols_intronic_change_is_read_on_its_placed_transcripts_alone() {
    let hdp = beside_a_placed_transcript();

    let consistent = VariantMapper::new(&hdp)
        .consistent_variants(&parse("TEST:c.5+1G>A"), None)
        .unwrap();

    let found: Vec<String> = consistent.iter().map(|c| c.variant.to_string()).collect();
    assert_eq!(found, ["NM_J.1:c.5+1G>A"]);
}

#[test]
fn a_gene_symbols_change_no_placed_transcript_holds_is_unplaced_rather_than_inconsistent() {
    let hdp = beside_a_placed_transcript();

    let result = VariantMapper::new(&hdp).consistent_variants(&parse("TEST:c.5+1C>A"), None);

    assert!(unplaced_refusal(&result), "{result:?}");
}

#[test]
fn a_gene_symbol_is_equivalent_to_the_genome_through_its_placed_transcript() {
    let hdp = beside_a_placed_transcript();
    let mapper = VariantMapper::new(&hdp);
    let equivalence = VariantEquivalence::new(&mapper, &hdp);

    let equivalent = equivalence
        .equivalent(&parse("TEST:c.4del"), &parse("NC_J.1:g.20del"))
        .unwrap();

    assert!(equivalent);
}

#[test]
fn a_gene_symbols_intronic_change_is_compared_on_its_placed_transcripts() {
    let hdp = beside_a_placed_transcript();
    let mapper = VariantMapper::new(&hdp);
    let equivalence = VariantEquivalence::new(&mapper, &hdp);

    let equivalent = equivalence
        .equivalent(&parse("TEST:c.5+1G>A"), &parse("NC_J.1:g.21G>A"))
        .unwrap();

    assert!(equivalent);
}

#[test]
fn a_region_search_that_names_an_unplaced_record_is_refused() {
    let hdp = unplaced();
    let g = match parse("NC_J.1:g.20del") {
        SequenceVariant::Genomic(g) => g,
        other => panic!("{other} is not g."),
    };
    let searcher = OnlyUnplaced;

    let result = VariantMapper::new(&hdp).g_to_c_all(&g, &searcher);

    assert!(unplaced_refusal(&result), "{result:?}");
}

/// A region search that, against its contract, names the unplaced `NM_U.1`.
struct OnlyUnplaced;

impl hgvs_weaver::data::TranscriptSearch for OnlyUnplaced {
    fn get_transcripts_for_region(
        &self,
        _: &str,
        _: i32,
        _: i32,
    ) -> Result<Vec<String>, HgvsError> {
        Ok(vec!["NM_U.1".to_string()])
    }
}

#[test]
fn a_record_without_placement_keys_is_unplaced() {
    let record: TranscriptData = serde_json::from_str(
        r#"{"ac":"NM_U.1","gene":"UNPLACED","cds_start_index":0,"cds_end_index":14}"#,
    )
    .unwrap();

    assert!(record.placement.is_none());
}

#[test]
fn a_record_with_part_of_a_placement_is_refused() {
    let result: Result<TranscriptData, _> = serde_json::from_str(
        r#"{"ac":"NM_U.1","gene":"UNPLACED","cds_start_index":0,"cds_end_index":14,"strand":1}"#,
    );

    let message = result.unwrap_err().to_string();
    assert!(
        message.contains("reference_accession") && message.contains("exons"),
        "{message}"
    );
}

#[test]
fn a_null_placement_key_is_refused_rather_than_read_as_absent() {
    let result: Result<TranscriptData, _> = serde_json::from_str(
        r#"{"ac":"NM_U.1","gene":"UNPLACED","cds_start_index":0,"cds_end_index":14,"strand":null,"reference_accession":null,"exons":null}"#,
    );

    let message = result.unwrap_err().to_string();
    assert!(message.contains("null"), "{message}");
}

#[test]
fn an_unplaced_record_round_trips_through_its_serialised_form() {
    let record = TranscriptData::unplaced("NM_U.1", "UNPLACED", None);

    let json = serde_json::to_value(&record).unwrap();
    let read: TranscriptData = serde_json::from_value(json.clone()).unwrap();

    assert!(read.placement.is_none(), "{json}");
}
