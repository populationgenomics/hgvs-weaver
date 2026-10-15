//! `c.=`: the whole transcript analysed and found unchanged. It is a statement
//! about the molecule, like `r.=` and `p.=`, so it has another spelling (r.=,
//! and n.= on a non-coding transcript) and a protein consequence (p.(=)), but
//! no genomic form and no allele.

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

/// The coding NM_W.1 (CDS at transcript indices 5..=94) and the non-coding
/// NR_W.1, both spliced from genome[10, 60) and genome[100, 160).
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
        .sequence("NP_W.1", "MRRRRKSRRRRRRRRVVVVVVVVVVVVVV")
        .transcript(transcript(
            "NM_W.1",
            "NC_W.1",
            Strand::Plus,
            Some((5, 94)),
            exons(),
        ))
        .transcript(transcript("NR_W.1", "NC_W.1", Strand::Plus, None, exons()))
        .protein_for("NM_W.1", "NP_W.1")
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

fn level(a: &str, b: &str) -> EquivalenceLevel {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    VariantEquivalence::new(&mapper, &hdp)
        .equivalent_level(&parse(a), &parse(b))
        .unwrap_or_else(|e| panic!("{a} vs {b}: {e}"))
}

#[test]
fn whole_sequence_identity_parses_and_prints_back() {
    for s in [
        "NM_W.1:c.=",
        "NM_W.1(GENE):c.=",
        "NR_W.1:n.=",
        "NM_W.1:c.[=]",
    ] {
        assert_eq!(parse(s).to_string(), s);
    }
}

#[test]
fn identity_without_a_position_is_refused_where_hgvs_does_not_define_it() {
    // A genome has no whole-sequence statement, and c. has no predicted form.
    for s in ["NC_W.1:g.=", "NC_012920.1:m.=", "NM_W.1:c.(=)"] {
        assert!(
            matches!(parse_hgvs_variant(s), Err(HgvsError::PestError(_))),
            "{s} parsed"
        );
    }
}

#[test]
fn a_whole_molecule_statement_shares_no_cis_allele_with_a_change() {
    for s in [
        "NM_W.1:c.[=;10C>G]",
        "NM_W.1:c.[10C>G;=]",
        "NM_W.1:r.[10c>g;=]",
        "NM_W.1:r.[spl;10c>g]",
    ] {
        assert!(
            matches!(parse_hgvs_variant(s), Err(HgvsError::ValidationError(m)) if m.contains("whole molecule")),
            "{s} was accepted"
        );
    }
}

#[test]
fn c_identity_is_r_identity_in_other_letters() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_eq!(
        mapper.c_to_r(&coding("NM_W.1:c.=")).unwrap().to_string(),
        "NM_W.1:r.="
    );
    assert_eq!(
        mapper.r_to_c(&rna("NM_W.1:r.=")).unwrap().to_string(),
        "NM_W.1:c.="
    );
}

#[test]
fn r_identity_on_a_non_coding_transcript_is_n_identity() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_eq!(
        mapper.r_to_n(&rna("NR_W.1:r.=")).unwrap().to_string(),
        "NR_W.1:n.="
    );
}

#[test]
fn a_predicted_r_identity_loses_its_parentheses_in_dna_letters() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_eq!(
        mapper.r_to_c(&rna("NM_W.1:r.(=)")).unwrap().to_string(),
        "NM_W.1:c.="
    );
}

#[test]
fn the_whole_transcript_unchanged_predicts_the_protein_unchanged() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_eq!(
        mapper
            .c_to_p(&coding("NM_W.1:c.="), None)
            .unwrap()
            .to_string(),
        "NP_W.1:p.(=)"
    );
}

#[test]
fn a_whole_molecule_statement_validates_with_nothing_to_check() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    for s in [
        "NM_W.1:c.=",
        "NR_W.1:n.=",
        "NM_W.1:r.=",
        "NR_W.1:r.=",
        "NM_W.1:r.spl",
        "NP_W.1:p.=",
        "NP_W.1:p.0",
        "NM_W.1:c.[=]",
    ] {
        assert!(mapper.validate(&parse(s)).unwrap(), "{s}");
    }
}

#[test]
fn validating_a_whole_molecule_statement_needs_its_reference() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    for s in ["NM_ABSENT.1:c.=", "NP_ABSENT.1:p.="] {
        assert!(
            matches!(mapper.validate(&parse(s)), Err(HgvsError::DataProviderError(m)) if m.contains("ABSENT")),
            "{s} validated"
        );
    }
}

#[test]
fn validating_c_identity_needs_a_coding_transcript() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert!(matches!(
        mapper.validate(&parse("NR_W.1:c.=")),
        Err(HgvsError::ValidationError(m)) if m.contains("no CDS")
    ));
}

#[test]
fn validating_n_identity_needs_a_non_coding_transcript() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert!(matches!(
        mapper.validate(&parse("NM_W.1:n.=")),
        Err(HgvsError::ValidationError(m)) if m.contains("has a CDS")
    ));
}

#[test]
fn normalising_whole_sequence_identity_leaves_it_as_written() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    for s in ["NM_W.1:c.=", "NR_W.1:n.=", "NM_W.1:c.[=]"] {
        assert_eq!(mapper.normalize_variant(parse(s)).unwrap().to_string(), s);
    }
}

#[test]
fn whole_sequence_identity_has_no_genomic_form() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert!(matches!(
        mapper.c_to_g(&coding("NM_W.1:c.="), None),
        Err(HgvsError::UnsupportedOperation(m)) if m.contains("NM_W.1:c.=")
    ));
    assert!(matches!(
        mapper.r_to_g(&rna("NR_W.1:r.="), None),
        Err(HgvsError::UnsupportedOperation(m)) if m.contains("NR_W.1:r.= describes the transcript as a whole")
    ));
}

#[test]
fn whole_sequence_identity_has_no_allele() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    for s in ["NM_W.1:c.=", "NR_W.1:n.=", "NM_W.1:r.=", "NP_W.1:p.="] {
        let v = parse(s);
        let unsupported =
            |r: Result<String, HgvsError>| matches!(r, Err(HgvsError::UnsupportedOperation(_)));
        assert!(unsupported(mapper.to_spdi(&v, false)), "{s} to_spdi");
        assert!(unsupported(mapper.to_spdi(&v, true)), "{s} unambiguous");
        assert!(unsupported(mapper.to_vrs(&v).map(|a| a.id)), "{s} to_vrs");
    }
}

#[test]
fn c_identity_and_r_identity_on_one_transcript_are_equivalent() {
    assert_eq!(
        level("NM_W.1:c.=", "NM_W.1:r.="),
        EquivalenceLevel::Identity
    );
}

#[test]
fn whole_sequence_identity_is_its_one_member_cis_allele() {
    assert_eq!(
        level("NM_W.1:c.=", "NM_W.1:c.[=]"),
        EquivalenceLevel::Identity
    );
}

#[test]
fn whole_sequence_identity_implies_the_protein_unchanged() {
    assert_eq!(
        level("NM_W.1:c.=", "NP_W.1:p.(=)"),
        EquivalenceLevel::Identity
    );
    assert_eq!(
        level("NM_W.1:c.=", "NP_W.1:p.="),
        EquivalenceLevel::Analogous
    );
}

#[test]
fn whole_sequence_identity_differs_from_a_change_or_one_unchanged_base() {
    for other in [
        "NM_W.1:c.10C>G",
        "NM_W.1:c.10=",
        "NC_W.1:g.30=",
        "NP_W.1:p.Arg2Gly",
    ] {
        assert_eq!(
            level("NM_W.1:c.=", other),
            EquivalenceLevel::Different,
            "{other}"
        );
    }
}

#[test]
fn non_coding_identity_says_nothing_about_a_protein() {
    for (a, b) in [("NR_W.1:n.=", "NP_W.1:p.="), ("NR_W.1:r.=", "NP_W.1:p.(=)")] {
        assert_eq!(level(a, b), EquivalenceLevel::Different, "{a} vs {b}");
    }
}
