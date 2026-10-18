//! `c.=`: the whole coding transcript analysed and found unchanged. It names
//! no position, so it has a protein consequence (`p.(=)`, on the transcript's
//! own protein) but no genomic form, no allele and no `r.` spelling: `r.=`
//! says the RNA was analysed, a different claim.

mod support;

use hgvs_weaver::data::Strand;
use hgvs_weaver::equivalence::{EquivalenceLevel, VariantEquivalence};
use hgvs_weaver::error::HgvsError;
use hgvs_weaver::mapper::VariantMapper;
use hgvs_weaver::{parse_hgvs_variant, SequenceVariant};
use support::{exon, transcript, Provider};

/// A 200-base genome: a CGT repeat with ATG at 15 and TAA at 142.
fn genome() -> String {
    let mut g: Vec<u8> = b"CGT".iter().copied().cycle().take(200).collect();
    g[15..18].copy_from_slice(b"ATG");
    g[30..34].copy_from_slice(b"AAAA");
    g[142..145].copy_from_slice(b"TAA");
    String::from_utf8(g).unwrap()
}

/// The coding NM_W.1 (CDS at transcript indices 5..=94, encoding NP_W.1), the
/// same with its CDS open at the 5' end (NM_OPEN.1), and the non-coding
/// NR_W.1, all spliced from genome[10, 60) and genome[100, 160). NP_OTHER.1
/// is a protein no transcript here encodes.
fn provider() -> Provider {
    let g = genome();
    let spliced = format!("{}{}", &g[10..60], &g[100..160]);
    let exons = || {
        vec![
            exon((0, 50), (10, 59), Strand::Plus),
            exon((50, 110), (100, 159), Strand::Plus),
        ]
    };
    Provider::new()
        .sequence("NC_W.1", &g)
        .sequence("NM_W.1", &spliced)
        .sequence("NR_W.1", &spliced)
        .sequence("NM_OPEN.1", &spliced)
        .sequence("NP_W.1", "MRRRRKSRRRRRRRRVVVVVVVVVVVVVV")
        .sequence("NP_OTHER.1", "MRRRRKSRRRRRRRRVVVVVVVVVVVVVV")
        .transcript(transcript(
            "NM_W.1",
            "NC_W.1",
            Strand::Plus,
            Some((5, 94)),
            exons(),
        ))
        .transcript(transcript("NR_W.1", "NC_W.1", Strand::Plus, None, exons()))
        .transcript(hgvs_weaver::data::TranscriptData {
            cds_start_open: true,
            ..transcript("NM_OPEN.1", "NC_W.1", Strand::Plus, Some((5, 94)), exons())
        })
        .protein_for("NM_W.1", "NP_W.1")
        .protein_for("NM_OPEN.1", "NP_OPEN.1")
}

fn parse(s: &str) -> SequenceVariant {
    parse_hgvs_variant(s).unwrap_or_else(|e| panic!("{s}: {e}"))
}

fn coding(s: &str) -> hgvs_weaver::CVariant {
    match parse(s) {
        SequenceVariant::Coding(c) => c,
        other => panic!("{other} is not c."),
    }
}

fn rna(s: &str) -> hgvs_weaver::RVariant {
    match parse(s) {
        SequenceVariant::Rna(r) => r,
        other => panic!("{other} is not r."),
    }
}

fn level(a: &str, b: &str) -> Result<EquivalenceLevel, HgvsError> {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    VariantEquivalence::new(&mapper, &hdp).equivalent_level(&parse(a), &parse(b))
}

fn is_parse_error<T>(r: Result<T, HgvsError>) -> bool {
    matches!(r, Err(HgvsError::PestError(_)))
}

fn is_unsupported<T>(r: Result<T, HgvsError>) -> bool {
    matches!(r, Err(HgvsError::UnsupportedOperation(_)))
}

fn from_json(json: &str) -> Result<SequenceVariant, serde_json::Error> {
    serde_json::from_str(json)
}

/// Whether `json` was refused because `system` cannot write what it describes.
fn refused_as(json: &str, system: &str) -> bool {
    from_json(json).is_err_and(|e| e.to_string().contains(&format!("{system} cannot write")))
}

// --- Grammar ---

#[test]
fn c_identity_prints_back_as_parsed() {
    assert_eq!(parse("NM_W.1(GENE):c.=").to_string(), "NM_W.1(GENE):c.=");
}

#[test]
fn c_identity_alone_in_brackets_prints_back_as_parsed() {
    assert_eq!(parse("NM_W.1:c.[=]").to_string(), "NM_W.1:c.[=]");
}

#[test]
fn n_identity_does_not_parse() {
    assert!(is_parse_error(parse_hgvs_variant("NR_W.1:n.=")));
}

#[test]
fn g_identity_does_not_parse() {
    assert!(is_parse_error(parse_hgvs_variant("NC_W.1:g.=")));
}

#[test]
fn m_identity_does_not_parse() {
    assert!(is_parse_error(parse_hgvs_variant("NC_012920.1:m.=")));
}

#[test]
fn a_predicted_c_identity_does_not_parse() {
    assert!(is_parse_error(parse_hgvs_variant("NM_W.1:c.(=)")));
}

// --- Alleles in cis ---

#[test]
fn c_identity_beside_a_change_is_refused() {
    assert!(matches!(
        parse_hgvs_variant("NM_W.1:c.[=;10C>G]"),
        Err(HgvsError::ValidationError(m)) if m.contains("NM_W.1:c.= states the whole sequence unchanged")
    ));
}

#[test]
fn r_identity_beside_a_change_is_refused() {
    assert!(matches!(
        parse_hgvs_variant("NM_W.1:r.[10c>g;=]"),
        Err(HgvsError::ValidationError(m)) if m.contains("NM_W.1:r.= states the whole sequence unchanged")
    ));
}

#[test]
fn p_identity_beside_a_change_is_refused() {
    assert!(matches!(
        parse_hgvs_variant("NP_W.1:p.[=;Arg2Gly]"),
        Err(HgvsError::ValidationError(m)) if m.contains("NP_W.1:p.= states the whole sequence unchanged")
    ));
}

/// A statement about the transcript other than the whole-sequence identity,
/// and an identity at a position, keep their place in a cis allele.
macro_rules! statement_in_cis_parses {
    ($($name:ident: $hgvs:literal,)*) => {$(
        #[test]
        fn $name() {
            assert_eq!(parse($hgvs).to_string(), $hgvs);
        }
    )*};
}

statement_in_cis_parses! {
    r_splicing_beside_a_change_parses: "NM_W.1:r.[spl;10c>g]",
    r_unknown_beside_a_change_parses: "NM_W.1:r.[?;10c>g]",
    r_no_rna_beside_a_change_parses: "NM_W.1:r.[0;10c>g]",
    r_splicing_twice_parses: "NM_W.1:r.[spl;spl]",
    p_unknown_beside_a_change_parses: "NP_W.1:p.[?;Arg2Gly]",
    p_no_protein_beside_a_change_parses: "NP_W.1:p.[0;Arg2Gly]",
    p_one_residue_unchanged_beside_a_change_parses: "NP_W.1:p.[Arg2=;Arg3Gly]",
    c_one_base_unchanged_beside_a_change_parses: "NM_W.1:c.[10=;12C>G]",
}

// --- Deserialisation admits only what parses ---

#[test]
fn c_identity_round_trips_through_json() {
    let v = parse("NM_W.1:c.=");
    let back = from_json(&serde_json::to_string(&v).unwrap()).unwrap();
    assert_eq!(back.to_string(), "NM_W.1:c.=");
}

#[test]
fn predicted_r_identity_round_trips_through_json() {
    let v = parse("NM_W.1:r.(=)");
    let back = from_json(&serde_json::to_string(&v).unwrap()).unwrap();
    assert_eq!(back.to_string(), "NM_W.1:r.(=)");
}

#[test]
fn a_deserialised_c_splicing_statement_is_refused() {
    let json = r#"{"variant_type": "Coding", "ac": "NM_W.1", "gene": null, "posedit": {"pos": null,
        "edit": {"type": "Special", "value": "spl", "uncertain": false}, "uncertain": false, "predicted": false}}"#;
    assert!(refused_as(json, "c."));
}

#[test]
fn a_deserialised_predicted_c_identity_is_refused() {
    let json = r#"{"variant_type": "Coding", "ac": "NM_W.1", "gene": null, "posedit": {"pos": null,
        "edit": {"type": "WholeIdentity"}, "uncertain": false, "predicted": true}}"#;
    assert!(refused_as(json, "c."));
}

#[test]
fn a_deserialised_g_identity_is_refused() {
    let json = r#"{"variant_type": "Genomic", "ac": "NC_W.1", "gene": null, "posedit": {"pos": null,
        "edit": {"type": "WholeIdentity"}, "uncertain": false, "predicted": false}}"#;
    assert!(refused_as(json, "g."));
}

#[test]
fn a_deserialised_n_identity_is_refused() {
    let json = r#"{"variant_type": "NonCoding", "ac": "NR_W.1", "gene": null, "posedit": {"pos": null,
        "edit": {"type": "WholeIdentity"}, "uncertain": false, "predicted": false}}"#;
    assert!(refused_as(json, "n."));
}

#[test]
fn a_deserialised_g_change_without_a_position_is_refused() {
    // It would print as g.=.
    let json = r#"{"variant_type": "Genomic", "ac": "NC_W.1", "gene": null, "posedit": {"pos": null,
        "edit": {"type": "None"}, "uncertain": false, "predicted": false}}"#;
    assert!(refused_as(json, "g."));
}

#[test]
fn an_r_identity_serialised_as_a_statement_string_still_loads() {
    let json = r#"{"variant_type": "Rna", "ac": "NM_W.1", "gene": null, "posedit": {"pos": null,
        "edit": {"type": "Special", "value": "=", "uncertain": false}, "uncertain": false, "predicted": true}}"#;
    let v = from_json(json).unwrap();
    assert_eq!(v.to_string(), "NM_W.1:r.(=)");
}

#[test]
fn an_r_statement_with_inner_whitespace_round_trips_through_json() {
    let v = parse("NM_W.1:r.spl ?");
    let back = from_json(&serde_json::to_string(&v).unwrap()).unwrap();
    assert_eq!(back, v);
}

#[test]
fn a_deserialised_uncertain_g_change_is_refused() {
    let json = r#"{"variant_type": "Genomic", "ac": "NC_W.1", "gene": null, "posedit": {
        "pos": {"start": {"base": 10, "end": null, "uncertain": false}, "end": null, "uncertain": false},
        "edit": {"type": "RefAlt", "ref_": "C", "alt": "G", "uncertain": false},
        "uncertain": true, "predicted": false}}"#;
    assert!(refused_as(json, "g."));
}

#[test]
fn a_deserialised_predicted_c_change_is_refused() {
    let json = r#"{"variant_type": "Coding", "ac": "NM_W.1", "gene": null, "posedit": {
        "pos": {"start": {"base": 10, "offset": null, "anchor": "CdsStart", "uncertain": false},
                "end": null, "uncertain": false},
        "edit": {"type": "RefAlt", "ref_": "C", "alt": "G", "uncertain": false},
        "uncertain": false, "predicted": true}}"#;
    assert!(refused_as(json, "c."));
}

#[test]
fn a_deserialised_cis_allele_with_c_identity_beside_a_change_is_refused() {
    let json = r#"{"variant_type": "CisPhased", "ac": "NM_W.1", "gene": null, "members": [
        {"variant_type": "Coding", "ac": "NM_W.1", "gene": null, "posedit": {"pos": null,
            "edit": {"type": "WholeIdentity"}, "uncertain": false, "predicted": false}},
        {"variant_type": "Coding", "ac": "NM_W.1", "gene": null, "posedit": {
            "pos": {"start": {"base": 10, "offset": null, "anchor": "CdsStart", "uncertain": false},
                    "end": null, "uncertain": false},
            "edit": {"type": "RefAlt", "ref_": "C", "alt": "G", "uncertain": false},
            "uncertain": false, "predicted": false}}]}"#;
    assert!(from_json(json).is_err_and(|e| e
        .to_string()
        .contains("states the whole sequence unchanged")));
}

// --- Protein ---

#[test]
fn c_identity_predicts_its_transcripts_protein_unchanged() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let p = mapper.c_to_p(&coding("NM_W.1:c.="), None).unwrap();
    assert_eq!(p.to_string(), "NP_W.1:p.(=)");
}

#[test]
fn c_identity_predicts_the_protein_unchanged_when_that_protein_is_named() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let p = mapper
        .c_to_p(&coding("NM_W.1:c.="), Some("NP_W.1"))
        .unwrap();
    assert_eq!(p.to_string(), "NP_W.1:p.(=)");
}

#[test]
fn c_identity_says_nothing_about_a_protein_its_transcript_does_not_encode() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert!(matches!(
        mapper.c_to_p(&coding("NM_W.1:c.="), Some("NP_OTHER.1")),
        Err(HgvsError::ValidationError(m)) if m.contains("NM_W.1 encodes NP_W.1, not NP_OTHER.1")
    ));
}

#[test]
fn c_identity_on_a_transcript_without_a_cds_has_no_protein_consequence() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert!(matches!(
        mapper.c_to_p(&coding("NR_W.1:c.="), Some("NP_W.1")),
        Err(HgvsError::ValidationError(m)) if m.contains("NR_W.1 has no CDS")
    ));
}

#[test]
fn c_identity_predicts_the_protein_unchanged_on_a_cds_open_at_its_start() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let p = mapper.c_to_p(&coding("NM_OPEN.1:c.="), None).unwrap();
    assert_eq!(p.to_string(), "NP_OPEN.1:p.(=)");
}

#[test]
fn c_identity_has_no_protein_allele() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert!(is_unsupported(
        mapper.protein_allele(&coding("NM_W.1:c.="), None)
    ));
}

#[test]
fn r_identity_on_a_transcript_without_a_cds_has_no_protein_consequence() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert!(matches!(
        mapper.r_to_p(&rna("NR_W.1:r.="), Some("NP_W.1")),
        Err(HgvsError::ValidationError(m)) if m.contains("NR_W.1 has no CDS")
    ));
}

#[test]
fn r_splicing_on_a_transcript_without_a_cds_has_no_protein_consequence() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert!(matches!(
        mapper.r_to_p(&rna("NR_W.1:r.spl"), Some("NP_W.1")),
        Err(HgvsError::ValidationError(m)) if m.contains("NR_W.1 has no CDS")
    ));
}

// --- No other spelling, no genomic form, no allele ---

#[test]
fn c_identity_has_no_r_spelling() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert!(is_unsupported(mapper.c_to_r(&coding("NM_W.1:c.="))));
}

#[test]
fn r_identity_has_no_c_spelling() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert!(is_unsupported(mapper.r_to_c(&rna("NM_W.1:r.="))));
}

#[test]
fn c_identity_has_no_genomic_form() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert!(matches!(
        mapper.c_to_g(&coding("NM_W.1:c.="), None),
        Err(HgvsError::UnsupportedOperation(m)) if m.contains("NM_W.1:c.= describes the whole transcript")
    ));
}

#[test]
fn c_identity_has_no_canonical_allele() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert!(is_unsupported(
        mapper.canonical_allele(&parse("NM_W.1:c.="))
    ));
}

#[test]
fn c_identity_has_no_spdi() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert!(is_unsupported(mapper.to_spdi(&parse("NM_W.1:c.="), false)));
}

#[test]
fn c_identity_has_no_unambiguous_spdi() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert!(is_unsupported(mapper.to_spdi(&parse("NM_W.1:c.="), true)));
}

#[test]
fn c_identity_has_no_vrs_allele() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert!(is_unsupported(mapper.to_vrs(&parse("NM_W.1:c.="))));
}

// --- Validation and normalisation ---

#[test]
fn c_identity_on_a_coding_transcript_validates() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert!(mapper.validate(&parse("NM_W.1:c.=")).unwrap());
}

#[test]
fn c_identity_on_a_transcript_without_a_cds_does_not_validate() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert!(matches!(
        mapper.validate(&parse("NR_W.1:c.=")),
        Err(HgvsError::ValidationError(m)) if m.contains("NR_W.1 has no CDS")
    ));
}

#[test]
fn c_identity_on_a_transcript_the_provider_lacks_does_not_validate() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert!(matches!(
        mapper.validate(&parse("NM_ABSENT.1:c.=")),
        Err(HgvsError::DataProviderError(m)) if m.contains("NM_ABSENT.1")
    ));
}

#[test]
fn normalising_c_identity_leaves_it_as_written() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let normalised = mapper.normalize_variant(parse("NM_W.1:c.=")).unwrap();
    assert_eq!(normalised.to_string(), "NM_W.1:c.=");
}

// --- Equivalence ---

#[test]
fn c_identity_is_its_one_member_cis_allele() {
    assert_eq!(
        level("NM_W.1:c.=", "NM_W.1:c.[=]").unwrap(),
        EquivalenceLevel::Identity
    );
}

#[test]
fn r_identity_is_its_one_member_cis_allele() {
    assert_eq!(
        level("NM_W.1:r.=", "NM_W.1:r.[=]").unwrap(),
        EquivalenceLevel::Identity
    );
}

#[test]
fn the_rna_predicted_unchanged_is_the_rna_found_unchanged() {
    // A prediction is taken as its confirmation, as p.(=) is p.=.
    assert_eq!(
        level("NM_W.1:r.(=)", "NM_W.1:r.=").unwrap(),
        EquivalenceLevel::Identity
    );
}

#[test]
fn c_identity_is_not_r_identity() {
    // One says the DNA was analysed, the other the RNA.
    assert_eq!(
        level("NM_W.1:c.=", "NM_W.1:r.=").unwrap(),
        EquivalenceLevel::Different
    );
}

#[test]
fn c_identity_is_the_predicted_protein_unchanged_on_its_own_protein() {
    assert_eq!(
        level("NM_W.1:c.=", "NP_W.1:p.(=)").unwrap(),
        EquivalenceLevel::Identity
    );
}

#[test]
fn c_identity_is_analogous_to_its_own_protein_unchanged() {
    assert_eq!(
        level("NM_W.1:c.=", "NP_W.1:p.=").unwrap(),
        EquivalenceLevel::Analogous
    );
}

#[test]
fn c_identity_is_not_another_protein_unchanged() {
    assert_eq!(
        level("NM_W.1:c.=", "NP_OTHER.1:p.(=)").unwrap(),
        EquivalenceLevel::Different
    );
}

#[test]
fn c_identity_is_not_one_base_unchanged() {
    assert_eq!(
        level("NM_W.1:c.=", "NM_W.1:c.10=").unwrap(),
        EquivalenceLevel::Different
    );
}

#[test]
fn c_identity_is_not_a_change() {
    assert_eq!(
        level("NM_W.1:c.=", "NM_W.1:c.10C>G").unwrap(),
        EquivalenceLevel::Different
    );
}
