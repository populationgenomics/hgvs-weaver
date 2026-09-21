//! Alleles in cis, `c.[145C>T;147C>G]`: parsed, formatted, normalised and
//! validated member by member, rendered as VRS `CisPhasedBlock` objects and
//! read back.

mod support;

use hgvs_weaver::data::Strand;
use hgvs_weaver::error::HgvsError;
use hgvs_weaver::mapper::VariantMapper;
use hgvs_weaver::{parse_hgvs_variant, SequenceVariant, Variant};
use support::{single_exon_transcript, Provider};

/// The transcript sits at genomic indices 10..50; c.1 is index 5 of it, so
/// g. = c. + 15 within the CDS.
const UTR5: &str = "GGGGG";
/// M K L A Y R *
const CDS: &str = "ATGAAACTGGCCTATCGCTAA";
const UTR3: &str = "CCGTATAAGTAAGG";
const PROTEIN: &str = "MKLAYR";

fn transcript() -> String {
    format!("{UTR5}{CDS}{UTR3}")
}

fn genome() -> String {
    format!("TTTTTTTTTT{}CCCCCCCCCC", transcript())
}

fn provider() -> Provider {
    let transcript = transcript();
    Provider::new()
        .sequence("NC_X.1", &genome())
        .sequence("NM_X.1", &transcript)
        .sequence("NP_X.1", PROTEIN)
        .transcript(single_exon_transcript(
            "NM_X.1",
            "NC_X.1",
            10,
            Strand::Plus,
            UTR5.len() as i32,
            (UTR5.len() + CDS.len()) as i32 - 1,
            transcript.len() as i32,
        ))
        .protein_for("NM_X.1", "NP_X.1")
}

fn cis(s: &str) -> hgvs_weaver::CisPhasedVariant {
    match parse_hgvs_variant(s).unwrap() {
        SequenceVariant::CisPhased(c) => c,
        other => panic!("{other} is not a cis allele"),
    }
}

#[test]
fn cis_alleles_parse_and_print_back_in_every_coordinate_system() {
    for hgvs in [
        "NC_X.1:g.[22C>T;28T>G]",
        "NC_012920.1:m.[8993T>G;9000del]",
        "NM_X.1:c.[7C>T;13T>G]",
        "NM_X.1:c.[7C>T;13T>G;20_21insA]",
        "NM_X.1(TEST):c.[122-6T>A;153C>T]",
        "NR_X.1:n.[7C>T;13T>G]",
        "NM_X.1:r.[7c>u;13u>g]",
        "NP_X.1:p.[Lys2Leu;Ala4del]",
        "NP_X.1:p.[Glu27Trp;Lys212fs]",
        "NM_X.1:c.[7C>T]",
    ] {
        let var = parse_hgvs_variant(hgvs).unwrap_or_else(|e| panic!("{hgvs}: {e}"));
        assert_eq!(var.to_string(), hgvs);
        let SequenceVariant::CisPhased(c) = &var else {
            panic!("{hgvs} did not parse as a cis allele");
        };
        assert_eq!(c.members.len(), hgvs.matches(';').count() + 1);
        for m in &c.members {
            assert_eq!(m.ac(), c.ac);
            assert_eq!(m.gene(), c.gene.as_deref());
            assert_eq!(m.coordinate_type(), var.coordinate_type());
        }
    }
    let c = cis("NM_X.1(TEST):c.[7C>T;13T>G]");
    assert_eq!(
        (c.ac(), c.gene(), c.coordinate_type()),
        ("NM_X.1", Some("TEST"), "c")
    );
    assert_eq!(c.members[0].to_string(), "NM_X.1(TEST):c.7C>T");
    assert_eq!(c.members[1].to_string(), "NM_X.1(TEST):c.13T>G");
}

#[test]
fn set_ac_reaches_the_members() {
    let mut var = parse_hgvs_variant("NM_X.1:c.[7C>T;13T>G]").unwrap();
    var.set_ac("NM_Y.2".into());
    assert_eq!(var.to_string(), "NM_Y.2:c.[7C>T;13T>G]");
    let SequenceVariant::CisPhased(c) = var else {
        unreachable!()
    };
    assert!(c.members.iter().all(|m| m.ac() == "NM_Y.2"));
}

#[test]
fn alleles_in_trans_are_two_molecules_and_are_refused() {
    for hgvs in [
        "NM_X.1:c.[7C>T];[13T>G]",
        "NM_X.1:c.[7C>T];[13T>G];[20del]",
        "NP_X.1:p.[Lys2Leu];[Ala4del]",
    ] {
        let err = parse_hgvs_variant(hgvs).unwrap_err();
        assert!(
            matches!(&err, HgvsError::UnsupportedOperation(m) if m.contains("two molecules")),
            "{hgvs}: {err}"
        );
    }
    // Malformed brackets are plain parse errors.
    for hgvs in [
        "NM_X.1:c.[]",
        "NM_X.1:c.[7C>T;]",
        "NM_X.1:c.[7C>T",
        "NM_X.1:c.7C>T]",
    ] {
        assert!(
            matches!(parse_hgvs_variant(hgvs), Err(HgvsError::PestError(_))),
            "{hgvs}"
        );
    }
}

#[test]
fn members_are_one_system_on_one_accession() {
    let a = parse_hgvs_variant("NM_X.1:c.7C>T").unwrap();
    let b = parse_hgvs_variant("NC_X.1:g.22C>T").unwrap();
    let p = parse_hgvs_variant("NM_X.1:p.Lys2Leu").unwrap();
    assert!(hgvs_weaver::CisPhasedVariant::new("NM_X.1".into(), None, vec![a.clone(), b]).is_err());
    assert!(hgvs_weaver::CisPhasedVariant::new("NM_X.1".into(), None, vec![a.clone(), p]).is_err());
    assert!(hgvs_weaver::CisPhasedVariant::new("NM_X.1".into(), None, vec![]).is_err());
    assert!(hgvs_weaver::CisPhasedVariant::new("NM_X.1".into(), None, vec![a]).is_ok());
}

#[test]
fn json_round_trips_with_the_members_nested() {
    let var = parse_hgvs_variant("NM_X.1:c.[7C>T;13T>G]").unwrap();
    let json = serde_json::to_string(&var).unwrap();
    assert!(
        json.starts_with(r#"{"variant_type":"CisPhased","ac":"NM_X.1""#),
        "{json}"
    );
    assert_eq!(json.matches(r#""variant_type":"Coding""#).count(), 2);
    let back: SequenceVariant = serde_json::from_str(&json).unwrap();
    assert_eq!(back, var);
}

#[test]
fn members_normalise_and_validate_one_by_one() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    // c.4 is the first A of AAA (c.4_6): the deletion shifts to c.6.
    let var = parse_hgvs_variant("NM_X.1:c.[4del;13T>G]").unwrap();
    assert_eq!(
        mapper.normalize_variant(var).unwrap().to_string(),
        "NM_X.1:c.[6del;13T>G]"
    );
    let var = parse_hgvs_variant("NC_X.1:g.[19del;28T>G]").unwrap();
    assert_eq!(
        mapper.normalize_variant(var).unwrap().to_string(),
        "NC_X.1:g.[21del;28T>G]"
    );

    let ok = parse_hgvs_variant("NM_X.1:c.[7C>T;13T>G]").unwrap();
    assert!(mapper.validate(&ok).unwrap());
    let wrong_second = parse_hgvs_variant("NM_X.1:c.[7C>T;13A>G]").unwrap();
    assert!(!mapper.validate(&wrong_second).unwrap());
    let protein = parse_hgvs_variant("NP_X.1:p.[Lys2Leu;Ala4del]").unwrap();
    assert!(mapper.validate(&protein).unwrap());
    let wrong_residue = parse_hgvs_variant("NP_X.1:p.[Lys2Leu;Gly4del]").unwrap();
    assert!(!mapper.validate(&wrong_residue).unwrap());
}

#[test]
fn a_cis_allele_has_no_single_canonical_allele_or_spdi() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let var = parse_hgvs_variant("NM_X.1:c.[7C>T;13T>G]").unwrap();
    let err = mapper.canonical_allele(&var).unwrap_err();
    assert!(
        matches!(&err, HgvsError::UnsupportedOperation(m) if m.contains("per member")),
        "{err}"
    );
    assert!(matches!(
        mapper.to_spdi(&var, false),
        Err(HgvsError::UnsupportedOperation(_))
    ));
    assert!(matches!(
        mapper.to_spdi_unambiguous(&var),
        Err(HgvsError::UnsupportedOperation(_))
    ));
    assert!(mapper.as_genomic(&var).is_none());
    // The members themselves are ordinary variants.
    let SequenceVariant::CisPhased(c) = &var else {
        unreachable!()
    };
    assert_eq!(
        mapper.to_spdi_unambiguous(&c.members[0]).unwrap(),
        "NC_X.1:21:C:T"
    );
}
