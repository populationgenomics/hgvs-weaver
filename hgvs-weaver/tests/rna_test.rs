//! r. variants: the transcript in RNA letters. Every operation is a
//! conversion to the c. or n. spelling and back, except that a change
//! spanning a splice junction describes the spliced RNA and cannot be
//! projected to the genome.

mod support;

use hgvs_weaver::data::Strand;
use hgvs_weaver::equivalence::VariantEquivalence;
use hgvs_weaver::error::HgvsError;
use hgvs_weaver::mapper::VariantMapper;
use hgvs_weaver::{parse_hgvs_variant, SequenceVariant};
use support::{exon, transcript, Provider};

/// A 200-base genome: a CGT repeat with ATG at 15, AAAA at 30 and TAA at 142.
fn genome() -> String {
    let mut g: Vec<u8> = b"CGT".iter().copied().cycle().take(200).collect();
    g[15..18].copy_from_slice(b"ATG");
    g[30..34].copy_from_slice(b"AAAA");
    g[142..145].copy_from_slice(b"TAA");
    String::from_utf8(g).unwrap()
}

/// Exon 1 is genome[10, 60), exon 2 genome[100, 160): a 110-base transcript.
/// The coding transcript's CDS is transcript indices 5..=94 (c.1 = index 5,
/// so c.N is index N + 4; the junction falls between c.45 and c.46).
fn spliced() -> String {
    let g = genome();
    format!("{}{}", &g[10..60], &g[100..160])
}

/// The coding NM_R.1 and the non-coding NR_R.1, both spliced from the genome.
fn provider() -> Provider {
    let exons = || {
        vec![
            exon((0, 50), (10, 59), Strand::Plus),
            exon((50, 110), (100, 159), Strand::Plus),
        ]
    };
    Provider::new()
        .sequence("NC_R.1", &genome())
        .sequence("NM_R.1", &spliced())
        .sequence("NR_R.1", &spliced())
        .transcript(transcript(
            "NM_R.1",
            "NC_R.1",
            Strand::Plus,
            Some((5, 94)),
            exons(),
        ))
        .transcript(transcript("NR_R.1", "NC_R.1", Strand::Plus, None, exons()))
        .protein_for("NM_R.1", "NP_R.1")
}

fn parse(s: &str) -> SequenceVariant {
    parse_hgvs_variant(s).unwrap_or_else(|e| panic!("{s}: {e}"))
}

fn rna(s: &str) -> hgvs_weaver::RVariant {
    match parse(s) {
        SequenceVariant::Rna(r) => r,
        other => panic!("{other} is not r."),
    }
}

fn coding(s: &str) -> hgvs_weaver::CVariant {
    match parse(s) {
        SequenceVariant::Coding(c) => c,
        other => panic!("{other} is not c."),
    }
}

#[test]
fn statements_about_the_transcript_round_trip() {
    for s in [
        "NM_R.1:r.0",
        "NM_R.1:r.0?",
        "NM_R.1:r.?",
        "NM_R.1:r.spl",
        "NM_R.1:r.spl?",
        "NM_R.1:r.(=)",
        "NM_R.1:r.=",
    ] {
        assert_eq!(parse(s).to_string(), s);
    }
}

#[test]
fn predicted_changes_keep_their_parentheses() {
    for s in ["NM_R.1:r.(10c>g)", "NM_R.1:r.(16_20del)"] {
        assert_eq!(parse(s).to_string(), s);
    }
    // The parentheses set only the predicted flag.
    let predicted = rna("NM_R.1:r.(10c>g)").posedit;
    let observed = rna("NM_R.1:r.10c>g").posedit;
    assert!(predicted.predicted);
    assert!(!observed.predicted);
    assert_eq!(
        hgvs_weaver::structs::PosEdit {
            predicted: false,
            ..predicted
        },
        observed
    );
    // An uncertain interval is not a predicted change.
    assert!(!rna("NM_R.1:r.(10_20)del").posedit.predicted);
    // c. has no parenthesised form in this grammar, so the flag is dropped on
    // the way to c. and the result still parses.
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let c = mapper.r_to_c(&rna("NM_R.1:r.(10c>g)")).unwrap();
    assert_eq!(c.to_string(), "NM_R.1:c.10C>G");
    assert_eq!(parse(&c.to_string()).to_string(), c.to_string());
}

#[test]
fn r_is_c_in_rna_letters() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    for (r, c) in [
        ("NM_R.1:r.10c>g", "NM_R.1:c.10C>G"),
        ("NM_R.1:r.-3_5delinsuu", "NM_R.1:c.-3_5delinsTT"),
        ("NM_R.1:r.*5u>a", "NM_R.1:c.*5T>A"),
        ("NM_R.1:r.45+2a>g", "NM_R.1:c.45+2A>G"),
        ("NM_R.1:r.16_17insaug", "NM_R.1:c.16_17insATG"),
        ("NM_R.1:r.16dup", "NM_R.1:c.16dup"),
        ("NM_R.1:r.10_12inv", "NM_R.1:c.10_12inv"),
        ("NM_R.1:r.16_20del", "NM_R.1:c.16_20del"),
        ("NM_R.1:r.16a[4]", "NM_R.1:c.16A[4]"),
        ("NM_R.1:r.16=", "NM_R.1:c.16="),
    ] {
        assert_eq!(mapper.r_to_c(&rna(r)).unwrap().to_string(), c, "{r}");
        assert_eq!(mapper.c_to_r(&coding(c)).unwrap().to_string(), r, "{c}");
    }
    // r.0 and friends have no c. spelling.
    let err = mapper.r_to_c(&rna("NM_R.1:r.spl")).unwrap_err();
    assert!(matches!(err, HgvsError::UnsupportedOperation(_)), "{err}");
}

#[test]
fn r_on_a_non_coding_transcript_is_n() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let n = mapper.r_to_n(&rna("NR_R.1:r.10c>g")).unwrap();
    assert_eq!(n.to_string(), "NR_R.1:n.10C>G");
    assert_eq!(mapper.n_to_r(&n).unwrap().to_string(), "NR_R.1:r.10c>g");
    assert_eq!(
        mapper
            .r_as_transcript(&rna("NR_R.1:r.10c>g"))
            .unwrap()
            .to_string(),
        "NR_R.1:n.10C>G"
    );
    assert_eq!(
        mapper
            .r_as_transcript(&rna("NM_R.1:r.10c>g"))
            .unwrap()
            .to_string(),
        "NM_R.1:c.10C>G"
    );
    // The wrong numbering for the transcript is refused, as are CDS-relative
    // positions where there is no CDS.
    assert!(mapper.r_to_c(&rna("NR_R.1:r.10c>g")).is_err());
    assert!(mapper.r_to_n(&rna("NM_R.1:r.10c>g")).is_err());
    assert!(mapper.r_to_n(&rna("NR_R.1:r.*5u>a")).is_err());
    assert!(mapper.r_to_n(&rna("NR_R.1:r.-3c>g")).is_err());
}

#[test]
fn r_projects_to_the_genome_within_one_exon_only() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let g = |s: &str| mapper.as_genomic(&parse(s)).unwrap().map(|v| v.to_string());
    // c.10 is transcript index 14, genome index 24: g.25.
    assert_eq!(g("NM_R.1:r.10c>g").unwrap(), "NC_R.1:g.25C>G");
    assert_eq!(g("NM_R.1:r.10c>g").unwrap(), g("NM_R.1:c.10C>G").unwrap());
    // The last base of exon 1 and the first of exon 2 project; a change
    // covering both does not, though its c. spelling does.
    assert!(g("NM_R.1:r.45del").is_ok());
    assert!(g("NM_R.1:r.46del").is_ok());
    for s in [
        "NM_R.1:r.44_47del",
        "NM_R.1:r.45_46insa",
        "NM_R.1:r.45_46del",
        "NR_R.1:r.49_52del",
    ] {
        let err = g(s).unwrap_err();
        assert!(
            matches!(&err, HgvsError::UnsupportedOperation(m) if m.contains("splice junction")),
            "{s}: {err}"
        );
    }
    assert!(g("NM_R.1:c.44_47del").is_ok());
    // An intronic position names the genome directly.
    assert_eq!(
        g("NM_R.1:r.45+2u>a").unwrap(),
        g("NM_R.1:c.45+2T>A").unwrap()
    );
    assert_eq!(g("NR_R.1:r.10c>g").unwrap(), "NC_R.1:g.20C>G");
}

#[test]
fn r_predicts_the_protein_like_c() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let p_from_r = |s: &str| mapper.r_to_p(&rna(s), None).map(|p| p.to_string());
    let p_from_c = |s: &str| mapper.c_to_p(&coding(s), None).unwrap().to_string();
    assert_eq!(
        p_from_r("NM_R.1:r.10c>g").unwrap(),
        p_from_c("NM_R.1:c.10C>G")
    );
    // Across the junction the RNA has no genomic form but a protein all the same.
    assert_eq!(
        p_from_r("NM_R.1:r.44_47del").unwrap(),
        p_from_c("NM_R.1:c.44_47del")
    );
    assert!(p_from_c("NM_R.1:c.44_47del").contains("fs"));
    assert_eq!(p_from_r("NM_R.1:r.0").unwrap(), "NP_R.1:p.0");
    assert_eq!(p_from_r("NM_R.1:r.0?").unwrap(), "NP_R.1:p.0?");
    assert_eq!(p_from_r("NM_R.1:r.spl").unwrap(), "NP_R.1:p.?");
    assert_eq!(p_from_r("NM_R.1:r.?").unwrap(), "NP_R.1:p.?");
    assert_eq!(p_from_r("NM_R.1:r.(=)").unwrap(), "NP_R.1:p.(=)");
    assert_eq!(
        mapper
            .r_to_p(&rna("NM_R.1:r.spl"), Some("NP_X.1"))
            .unwrap()
            .to_string(),
        "NP_X.1:p.?"
    );
    assert!(
        p_from_r("NR_R.1:r.10c>g").is_err(),
        "a non-coding transcript makes no protein"
    );
}

#[test]
fn r_normalises_validates_and_has_alleles_like_c() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    // c.16..c.20 are AAAA A? No: transcript indices 20..24 are A, i.e. c.16..c.19,
    // so a deletion of one A rolls to c.19.
    let n = |s: &str| mapper.normalize_variant(parse(s)).unwrap().to_string();
    // The r. result is the c. (or n.) result in RNA letters.
    let as_r = |s: &str| match parse(s) {
        SequenceVariant::Coding(c) => mapper.c_to_r(&c).unwrap().to_string(),
        SequenceVariant::NonCoding(v) => mapper.n_to_r(&v).unwrap().to_string(),
        other => panic!("{other}"),
    };
    assert!(
        n("NM_R.1:r.16del").starts_with("NM_R.1:r.19del"),
        "{}",
        n("NM_R.1:r.16del")
    );
    assert_eq!(n("NM_R.1:r.16del"), as_r(&n("NM_R.1:c.16del")));
    assert_eq!(n("NM_R.1:r.15_16insa"), as_r(&n("NM_R.1:c.15_16insA")));
    assert!(
        n("NM_R.1:r.15_16insa").contains("r.19dup"),
        "{}",
        n("NM_R.1:r.15_16insa")
    );
    assert_eq!(n("NM_R.1:r.spl"), "NM_R.1:r.spl");
    // Without a CDS the A run is n.21..n.24.
    assert_eq!(n("NR_R.1:r.21del"), as_r(&n("NR_R.1:n.21del")));
    assert!(
        n("NR_R.1:r.21del").starts_with("NR_R.1:r.24del"),
        "{}",
        n("NR_R.1:r.21del")
    );

    assert!(mapper.validate(&parse("NM_R.1:r.16a>g")).unwrap());
    assert!(!mapper.validate(&parse("NM_R.1:r.16c>g")).unwrap());
    assert!(mapper.validate(&parse("NM_R.1:r.16_19delaaaa")).unwrap());

    let r = parse("NM_R.1:r.16del");
    let c = parse("NM_R.1:c.16del");
    assert_eq!(
        mapper.canonical_allele(&r).unwrap(),
        mapper.canonical_allele(&c).unwrap()
    );
    assert_eq!(
        mapper.to_spdi_unambiguous(&r).unwrap(),
        mapper.to_spdi_unambiguous(&c).unwrap()
    );
    let vrs = mapper.to_vrs(&r).unwrap();
    assert_eq!(vrs.id, mapper.to_vrs(&c).unwrap().id);
    assert_eq!(vrs.expressions[0].syntax, "hgvs.r");

    let eq = VariantEquivalence::new(&mapper, &hdp);
    assert!(eq.equivalent_level(&r, &c).unwrap().is_equivalent());
    assert!(eq
        .equivalent_level(&r, &parse("NM_R.1:c.19del"))
        .unwrap()
        .is_equivalent());
    assert!(!eq
        .equivalent_level(&r, &parse("NM_R.1:c.10del"))
        .unwrap()
        .is_equivalent());
}
