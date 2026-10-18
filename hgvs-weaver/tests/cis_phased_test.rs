//! Alleles in cis, `c.[145C>T;147C>G]`: parsed, formatted, normalised and
//! validated member by member, rendered as VRS `CisPhasedBlock` objects and
//! read back; and read into the one sequence their members make.

mod support;

use std::collections::BTreeSet;

use hgvs_weaver::data::Strand;
use hgvs_weaver::equivalence::{EquivalenceLevel, VariantEquivalence};
use hgvs_weaver::error::HgvsError;
use hgvs_weaver::mapper::VariantMapper;
use hgvs_weaver::reading::CisReading;
use hgvs_weaver::vrs::VrsBound::Exact;
use hgvs_weaver::vrs::{refget_accession, VrsExpression, VrsVariation};
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
        .sequence("NC_OTHER.1", "GATTACAGATTACAGATTACAGATTACAGATTACAGATTACA")
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
fn members_that_meet_no_other_normalise_alone_and_validate_one_by_one() {
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
fn a_single_variant_that_makes_the_reference_is_its_statement() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let normalised = |s: &str| {
        mapper
            .normalize_variant(parse_hgvs_variant(s).unwrap())
            .unwrap_or_else(|e| panic!("{s}: {e}"))
            .to_string()
    };
    // c.9_12 is GGCC, its own reverse complement.
    assert_eq!(normalised("NM_X.1:c.9_12inv"), "NM_X.1:c.9_12=");
    assert_eq!(normalised("NM_X.1:c.7_9delinsCTG"), "NM_X.1:c.7_9=");
    // Three copies of the A that c.4_6 already holds: the run it resolves to.
    assert_eq!(normalised("NM_X.1:c.4A[3]"), "NM_X.1:c.4_6=");
}

#[test]
fn statements_about_different_bases_differ() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let equivalence = VariantEquivalence::new(&mapper, &hdp);
    let level = |a: &str, b: &str| {
        equivalence
            .equivalent_level(
                &parse_hgvs_variant(a).unwrap(),
                &parse_hgvs_variant(b).unwrap(),
            )
            .unwrap_or_else(|e| panic!("{a} vs {b}: {e}"))
    };
    // Both make the reference; they state different bases.
    assert_eq!(
        level("NM_X.1:c.10=", "NM_X.1:c.20="),
        EquivalenceLevel::Different
    );
    // c.10 is g.25: the projection states the same base.
    assert_eq!(
        level("NM_X.1:c.10=", "NC_X.1:g.25="),
        EquivalenceLevel::Identity
    );
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

#[test]
fn a_cis_allele_is_a_cis_phased_block_of_its_members_alleles() {
    let hdp = provider();
    let mapper = VariantMapper::with_refget(&hdp, &hdp);
    let var = parse_hgvs_variant("NM_X.1:c.[7C>T;13T>G]").unwrap();
    let block = mapper.to_vrs_cis_phased(&var).unwrap();
    assert_eq!(block.type_, "CisPhasedBlock");
    assert!(block.id.starts_with("ga4gh:CPB."), "{}", block.id);
    assert_eq!(block.id, format!("ga4gh:CPB.{}", block.digest));

    // Each member is the Allele the member variant renders as on its own,
    // written as its projection to the genome.
    let c = cis("NM_X.1:c.[7C>T;13T>G]");
    assert_eq!(block.members.len(), 2);
    for (allele, member) in block.members.iter().zip(&c.members) {
        assert_eq!(allele.id, mapper.to_vrs(member).unwrap().id);
    }
    assert_eq!(block.members[0].expressions[0].value, "NC_X.1:g.22C>T");
    // c.7 is g.22, interbase [21, 22); c.13 is g.28.
    let location = &block.members[0].location;
    assert_eq!((location.start, location.end), (Exact(21), Exact(22)));
    let location = &block.members[1].location;
    assert_eq!((location.start, location.end), (Exact(27), Exact(28)));

    // All on the genome, which the block names too.
    let refget = refget_accession(&genome());
    let reference = block.sequence_reference.as_ref().unwrap();
    assert_eq!(reference.refget_accession, refget);
    assert_eq!(reference.residue_alphabet, "na");
    assert_eq!(*reference, block.members[0].location.sequence_reference);
    assert_eq!(
        block.expressions,
        vec![VrsExpression {
            syntax: "hgvs.c".into(),
            value: "NM_X.1:c.[7C>T;13T>G]".into(),
        }]
    );

    let variation = mapper.to_vrs_variation(&var).unwrap();
    assert_eq!(variation.id(), block.id);
    match variation {
        VrsVariation::CisPhasedBlock(b) => assert_eq!(b, block),
        other => panic!("expected a CisPhasedBlock, got {other:?}"),
    }

    let json = block.to_json();
    for field in [
        r#""type":"CisPhasedBlock""#,
        r#""members":[{"id":"ga4gh:VA."#,
        r#""sequenceReference":{"type":"SequenceReference","refgetAccession":"SQ."#,
        r#""expressions":[{"syntax":"hgvs.c","value":"NM_X.1:c.[7C>T;13T>G]"}]"#,
    ] {
        assert!(json.contains(field), "{field} missing from {json}");
    }
}

#[test]
fn the_identifier_is_the_same_whichever_way_the_members_are_written() {
    let hdp = provider();
    let mapper = VariantMapper::with_refget(&hdp, &hdp);
    let block = |hgvs: &str| {
        mapper
            .to_vrs_cis_phased(&parse_hgvs_variant(hgvs).unwrap())
            .unwrap()
    };
    let forward = block("NM_X.1:c.[7C>T;13T>G]");
    let reversed = block("NM_X.1:c.[13T>G;7C>T]");
    assert_eq!(forward.id, reversed.id);
    // The same changes spelled on the genome, or written unnormalised.
    assert_eq!(block("NC_X.1:g.[22C>T;28T>G]").id, forward.id);
    assert_eq!(
        block("NM_X.1:c.[4del;13T>G]").id,
        block("NC_X.1:g.[21del;28T>G]").id
    );
    // Different members, different block.
    assert_ne!(block("NM_X.1:c.[7C>T;13T>A]").id, forward.id);
}

#[test]
fn to_vrs_is_for_alleles_and_to_vrs_cis_phased_for_cis_alleles() {
    let hdp = provider();
    let mapper = VariantMapper::with_refget(&hdp, &hdp);
    let cis = parse_hgvs_variant("NM_X.1:c.[7C>T;13T>G]").unwrap();
    let err = mapper.to_vrs(&cis).unwrap_err();
    assert!(
        matches!(&err, HgvsError::UnsupportedOperation(m) if m.contains("to_vrs_variation")),
        "{err}"
    );
    let plain = parse_hgvs_variant("NM_X.1:c.7C>T").unwrap();
    assert!(matches!(
        mapper.to_vrs_cis_phased(&plain),
        Err(HgvsError::UnsupportedOperation(_))
    ));
    // A member without an allele fails the block.
    let with_copy = parse_hgvs_variant("NC_X.1:g.[22C>T;30_35copy3]").unwrap();
    assert!(matches!(
        mapper.to_vrs_cis_phased(&with_copy),
        Err(HgvsError::UnsupportedOperation(_))
    ));
}

#[test]
fn cis_phased_blocks_read_back_as_cis_alleles_on_their_own_sequence() {
    let hdp = provider();
    let mapper = VariantMapper::with_refget(&hdp, &hdp);
    for (hgvs, expected) in [
        ("NM_X.1:c.[7C>T;13T>G]", "NC_X.1:g.[22C>T;28T>G]"),
        // Each member comes back 3'-normalised.
        ("NM_X.1:c.[4del;13T>G]", "NC_X.1:g.[21del;28T>G]"),
        ("NC_X.1:g.[19del;28T>G]", "NC_X.1:g.[21del;28T>G]"),
        ("NC_X.1:m.[22C>T;28T>G]", "NC_X.1:g.[22C>T;28T>G]"),
        ("NP_X.1:p.[Lys2Leu;Ala4del]", "NP_X.1:p.[Lys2Leu;Ala4del]"),
    ] {
        let var = parse_hgvs_variant(hgvs).unwrap();
        let block = mapper.to_vrs_cis_phased(&var).unwrap();
        let back = mapper.from_vrs(&block.to_json(), None).unwrap();
        assert_eq!(back.to_string(), expected, "{hgvs}");
        assert_eq!(
            mapper.to_vrs_cis_phased(&back).unwrap().id,
            block.id,
            "{hgvs}"
        );
    }
    // Without a refget lookup the accession is passed, and checked.
    let mapper = VariantMapper::new(&hdp);
    let var = parse_hgvs_variant("NM_X.1:c.[7C>T;13T>G]").unwrap();
    let json = mapper.to_vrs_cis_phased(&var).unwrap().to_json();
    assert!(matches!(
        mapper.from_vrs(&json, None),
        Err(HgvsError::DataProviderError(_))
    ));
    assert_eq!(
        mapper.from_vrs(&json, Some("NC_X.1")).unwrap().to_string(),
        "NC_X.1:g.[22C>T;28T>G]"
    );
    assert!(matches!(
        mapper.from_vrs(&json, Some("NC_OTHER.1")),
        Err(HgvsError::ValidationError(_))
    ));
}

#[test]
fn blocks_from_other_producers_parse() {
    let hdp = provider();
    let mapper = VariantMapper::with_refget(&hdp, &hdp);
    let refget = refget_accession(&genome());
    let allele = |start: usize, alt: &str, reference: &str| {
        format!(
            r#"{{"type":"Allele","location":{{"type":"SequenceLocation"{reference},
                "start":{start},"end":{}}},
                "state":{{"type":"LiteralSequenceExpression","sequence":"{alt}"}}}}"#,
            start + 1
        )
    };
    // The spec's shape: the sequence stated once on the block, no ids.
    let json = format!(
        r#"{{"type":"CisPhasedBlock","members":[{},{}],
            "sequenceReference":{{"type":"SequenceReference","refgetAccession":"{refget}"}}}}"#,
        allele(21, "T", ""),
        allele(27, "G", ""),
    );
    assert_eq!(
        mapper.from_vrs(&json, None).unwrap().to_string(),
        "NC_X.1:g.[22C>T;28T>G]"
    );
    // The same with the sequence on each member and none on the block.
    let on_member = format!(
        r#","sequenceReference":{{"type":"SequenceReference","refgetAccession":"{refget}"}}"#
    );
    let json = format!(
        r#"{{"type":"CisPhasedBlock","members":[{},{}]}}"#,
        allele(21, "T", &on_member),
        allele(27, "G", &on_member),
    );
    assert_eq!(
        mapper.from_vrs(&json, None).unwrap().to_string(),
        "NC_X.1:g.[22C>T;28T>G]"
    );
    // A member on another sequence is refused.
    let other = format!(
        r#","sequenceReference":{{"type":"SequenceReference","refgetAccession":"{}"}}"#,
        refget_accession("GATTACAGATTACAGATTACAGATTACAGATTACAGATTACA")
    );
    let json = format!(
        r#"{{"type":"CisPhasedBlock","members":[{},{}]}}"#,
        allele(21, "T", &on_member),
        allele(3, "C", &other),
    );
    assert!(matches!(
        mapper.from_vrs(&json, None),
        Err(HgvsError::ValidationError(_))
    ));
    // As is a block whose stated sequence is not its members'.
    let json = format!(
        r#"{{"type":"CisPhasedBlock","members":[{}],
            "sequenceReference":{{"type":"SequenceReference","refgetAccession":"SQ.x"}}}}"#,
        allele(21, "T", &on_member),
    );
    assert!(matches!(
        mapper.from_vrs(&json, None),
        Err(HgvsError::ValidationError(_))
    ));
}

#[test]
fn cis_alleles_compare_as_sets_of_their_members_alleles() {
    let hdp = provider();
    let mapper = VariantMapper::with_refget(&hdp, &hdp);
    let equivalence = VariantEquivalence::new(&mapper, &hdp);
    let level = |a: &str, b: &str| {
        equivalence
            .equivalent_level(
                &parse_hgvs_variant(a).unwrap(),
                &parse_hgvs_variant(b).unwrap(),
            )
            .unwrap_or_else(|e| panic!("{a} vs {b}: {e}"))
    };
    use EquivalenceLevel::{Analogous, Different, Identity};
    assert_eq!(
        level("NM_X.1:c.[7C>T;13T>G]", "NM_X.1:c.[7C>T;13T>G]"),
        Identity
    );
    // The same members in the other order are the same spelling; spelled on
    // the genome, or written unnormalised, they name the same molecule.
    assert_eq!(
        level("NM_X.1:c.[7C>T;13T>G]", "NM_X.1:c.[13T>G;7C>T]"),
        Identity
    );
    assert_eq!(
        level("NM_X.1:c.[7C>T;13T>G]", "NC_X.1:g.[22C>T;28T>G]"),
        Analogous
    );
    assert_eq!(
        level("NM_X.1:c.[4del;13T>G]", "NC_X.1:g.[21del;28T>G]"),
        Analogous
    );
    assert_eq!(
        level("NC_X.1:m.[22C>T;28T>G]", "NC_X.1:g.[28T>G;22C>T]"),
        Analogous
    );
    // A different member, a missing member, an extra member.
    assert_eq!(
        level("NM_X.1:c.[7C>T;13T>G]", "NM_X.1:c.[7C>T;13T>A]"),
        Different
    );
    assert_eq!(level("NM_X.1:c.[7C>T;13T>G]", "NM_X.1:c.[7C>T]"), Different);
    assert_eq!(
        level("NM_X.1:c.[7C>T;13T>G]", "NM_X.1:c.[7C>T;13T>G;20del]"),
        Different
    );
    // Against a plain variant a cis allele of several members is Different,
    // one of a single member is that member, to the letter when it is.
    assert_eq!(level("NM_X.1:c.[7C>T;13T>G]", "NM_X.1:c.7C>T"), Different);
    assert_eq!(level("NM_X.1:c.7C>T", "NM_X.1:c.[7C>T;13T>G]"), Different);
    assert_eq!(level("NM_X.1:c.[7C>T]", "NM_X.1:c.7C>T"), Identity);
    assert_eq!(level("NM_X.1:c.[4del]", "NC_X.1:g.21del"), Analogous);
    assert_eq!(level("NC_X.1:g.22C>T", "NM_X.1:c.[7C>T]"), Identity);
    assert_eq!(level("NM_X.1:c.[7C>T]", "NM_X.1:c.13T>G"), Different);
    // Protein cis alleles compare by their members' protein alleles; no
    // nucleotide cis allele is projected to protein.
    assert_eq!(
        level("NP_X.1:p.[Lys2Leu;Ala4del]", "NP_X.1:p.[Ala4del;Lys2Leu]"),
        Identity
    );
    assert_eq!(
        level("NP_X.1:p.[Lys2Leu;Ala4del]", "NP_X.1:p.[K2L;A4del]"),
        Identity
    );
    assert_eq!(
        level(
            "NP_X.1:p.[Lys2Leu;Ala4del]",
            "NP_X.1:p.[Ala4del;Lys2delinsLeu]"
        ),
        Analogous
    );
    assert_eq!(
        level("NP_X.1:p.[Lys2Leu;Ala4del]", "NP_X.1:p.[Lys2Leu;Tyr5del]"),
        Different
    );
    assert_eq!(
        level("NM_X.1:c.[7C>T;13T>G]", "NP_X.1:p.Leu3Phe"),
        Different
    );
    // With a member that has no canonical allele, the same members in
    // another order are the same spelling, and others are Different, not an
    // error.
    assert_eq!(
        level("NC_X.1:g.[22C>T;30_35copy3]", "NC_X.1:g.[30_35copy3;22C>T]"),
        Identity
    );
    assert_eq!(
        level("NC_X.1:g.[22C>T;30_35copy3]", "NC_X.1:g.[22C>T;30_35copy4]"),
        Different
    );
}

// Reading an allele into the one sequence its members make: the overlap rule
// of choices.md, "Overlapping members of a cis allele need one reading", and
// the ambiguity it can end in.

/// `NM_X.1` and `NC_X.1` as above; `NM_S.1` is `NM_X.1` with exon 1 ending at
/// c.5 and a ten-base intron after it in `NC_S.1`, so c.4 is g.19 and c.6 is
/// g.31; `NR_MS.1` is `NM_S.1`'s two exons read on the minus strand of
/// `NC_S.1`, so its n.30 is g.31 and n.30+2 the intron's g.29; `NR_XR.1`
/// reads `NM_X.1`'s forty bases of `NC_X.1` on the minus strand, so n.23 is
/// g.28; `NC_T.1` is ten Ts; `NC_P.1` is seventeen `GAAGCCC` and a `GAAACCC`;
/// `NC_Q.1` is eight `ACGTTGCAAGCTTCGA`; `NC_U.1` and `NC_V.1` are short AT
/// sequences with runs, as are `NC_W.1`, `NC_Y.1` and `NC_Z.1`.
fn reading_provider() -> Provider {
    let transcript = transcript();
    let len = transcript.len() as i32;
    let split = format!(
        "TTTTTTTTTT{}GTAAGTCCAG{}CCCCCCCCCC",
        &transcript[..10],
        &transcript[10..]
    );
    let minus_spliced = hgvs_weaver::utils::reverse_complement(&transcript);
    provider()
        .sequence("NC_S.1", &split)
        .sequence("NM_S.1", &transcript)
        .sequence("NR_MS.1", &minus_spliced)
        .sequence(
            "NR_XR.1",
            &hgvs_weaver::utils::reverse_complement(&transcript),
        )
        .sequence("NC_T.1", "TTTTTTTTTT")
        .sequence("NC_P.1", &format!("{}GAAACCC", "GAAGCCC".repeat(17)))
        .sequence("NC_Q.1", &"ACGTTGCAAGCTTCGA".repeat(8))
        .sequence("NC_U.1", "ATTTATTTTAATT")
        .sequence("NC_V.1", "ATAAATATTA")
        .sequence("NC_W.1", "TTTAAATTTAAATTTAAA")
        .sequence("NC_Y.1", "TTTTAAAT")
        .sequence("NC_Z.1", "TATTATTA")
        .transcript(support::transcript(
            "NM_S.1",
            "NC_S.1",
            Strand::Plus,
            Some((UTR5.len() as i32, (UTR5.len() + CDS.len()) as i32 - 1)),
            vec![
                support::exon((0, 10), (10, 19), Strand::Plus),
                support::exon((10, len), (30, 30 + len - 10 - 1), Strand::Plus),
            ],
        ))
        .transcript(support::transcript(
            "NR_MS.1",
            "NC_S.1",
            Strand::Minus,
            None,
            vec![
                support::exon((0, len - 10), (30, 30 + len - 10 - 1), Strand::Minus),
                support::exon((len - 10, len), (10, 19), Strand::Minus),
            ],
        ))
        .transcript(support::transcript(
            "NR_XR.1",
            "NC_X.1",
            Strand::Minus,
            None,
            vec![support::exon((0, len), (10, 10 + len - 1), Strand::Minus)],
        ))
}

fn read(hdp: &Provider, hgvs: &str) -> Result<CisReading, HgvsError> {
    VariantMapper::new(hdp).read_cis(&cis(hgvs))
}

/// What the allele makes of the reference bases `[start, end)` (0-based) of
/// the sequence it is read on.
fn makes(hdp: &Provider, hgvs: &str, start: usize, end: usize) -> String {
    let reading = read(hdp, hgvs).unwrap_or_else(|e| panic!("{hgvs}: {e}"));
    let mapper = VariantMapper::new(hdp);
    let reference = mapper.refs.reference(reading.accession(), reading.kind());
    reading.sequence(&reference, start, end).unwrap()
}

/// c.4_12 of `NM_X.1`, `AAACTGGCC`, is transcript indices 8..17.
const C4: usize = 8;
const C12_END: usize = 17;

/// Whether the two readings make the same sequence.
fn same(hdp: &Provider, a: &CisReading, b: &CisReading) -> bool {
    let mapper = VariantMapper::new(hdp);
    let reference = mapper.refs.reference(a.accession(), a.kind());
    a.same_sequence(b, &reference).unwrap()
}

/// Whether the two alleles, read on one sequence, make the same sequence.
fn same_sequence(hdp: &Provider, a: &str, b: &str) -> bool {
    let ra = read(hdp, a).unwrap_or_else(|e| panic!("{a}: {e}"));
    let rb = read(hdp, b).unwrap_or_else(|e| panic!("{b}: {e}"));
    same(hdp, &ra, &rb)
}

/// The readings the allele is ambiguous between, as written.
fn readings(hdp: &Provider, hgvs: &str) -> BTreeSet<String> {
    match read(hdp, hgvs) {
        Err(HgvsError::Ambiguous { readings, .. }) => {
            readings.iter().map(ToString::to_string).collect()
        }
        other => panic!("{hgvs} is not ambiguous: {other:?}"),
    }
}

fn set(readings: &[&str]) -> BTreeSet<String> {
    readings.iter().map(|r| r.to_string()).collect()
}

fn assert_contradictory(hdp: &Provider, hgvs: &str) {
    let err = read(hdp, hgvs).unwrap_err();
    assert!(
        matches!(&err, HgvsError::ValidationError(m) if m.contains("has no reading")),
        "{hgvs}: {err}"
    );
}

#[test]
fn members_meeting_at_a_junction_are_accepted() {
    let hdp = reading_provider();
    assert_eq!(
        makes(&hdp, "NM_X.1:c.[7_8del;9_10del]", C4, C12_END),
        "AAACC"
    );
}

#[test]
fn an_insertion_beside_a_changed_base_is_accepted() {
    let hdp = reading_provider();
    assert_eq!(
        makes(&hdp, "NM_X.1:c.[7C>T;7_8insTA]", C4, C12_END),
        "AAATTATGGCC"
    );
}

#[test]
fn an_insertion_before_a_changed_base_is_accepted() {
    let hdp = reading_provider();
    assert_eq!(
        makes(&hdp, "NM_X.1:c.[7_8insA;8T>G]", C4, C12_END),
        "AAACAGGGCC"
    );
}

#[test]
fn an_insertion_before_what_a_duplication_copies_is_accepted() {
    let hdp = reading_provider();
    assert_eq!(
        makes(&hdp, "NM_X.1:c.[7_9dup;6_7insA]", C4, C12_END),
        "AAAACTGCTGGCC"
    );
}

#[test]
fn members_that_agree_on_a_base_are_one_change() {
    let hdp = reading_provider();
    assert_eq!(
        makes(&hdp, "NM_X.1:c.[7_8delinsTG;7C>T]", C4, C12_END),
        "AAATGGGCC"
    );
}

#[test]
fn overlapping_deletions_delete_each_base_once() {
    let hdp = reading_provider();
    assert_eq!(
        makes(&hdp, "NM_X.1:c.[7_9del;8_10del]", C4, C12_END),
        "AAACC"
    );
}

#[test]
fn an_insertion_inside_a_same_length_delins_is_ambiguous() {
    let hdp = reading_provider();
    assert_eq!(
        readings(&hdp, "NM_X.1:c.[7_8delinsTG;7_8insA]"),
        set(&["NM_X.1:c.7_8delinsTG", "NM_X.1:c.7_8delinsTAG"])
    );
}

#[test]
fn an_insertion_inside_a_deletion_is_ambiguous() {
    let hdp = reading_provider();
    assert_eq!(
        readings(&hdp, "NM_X.1:c.[7_9del;7_8insA]"),
        set(&["NM_X.1:c.7_9del", "NM_X.1:c.7_9delinsA"])
    );
}

#[test]
fn an_insertion_at_the_flank_of_a_deletion_is_accepted() {
    let hdp = reading_provider();
    assert_eq!(
        makes(&hdp, "NM_X.1:c.[7_9del;6_7insA]", C4, C12_END),
        "AAAAGCC"
    );
}

#[test]
fn an_insertion_at_the_flank_of_a_delins_is_accepted() {
    let hdp = reading_provider();
    assert_eq!(
        makes(&hdp, "NM_X.1:c.[7_8delinsTG;8_9insA]", C4, C12_END),
        "AAATGAGGCC"
    );
}

#[test]
fn an_identity_beside_an_agreeing_change_stays() {
    let hdp = reading_provider();
    // c.7 is transcript index 11
    let reading = read(&hdp, "NM_X.1:c.[7=;7_8delinsCG]").unwrap();
    assert_eq!(reading.stated(), &[(11, 12)]);
}

#[test]
fn a_change_agreeing_with_an_identity_inside_it_is_accepted() {
    let hdp = reading_provider();
    // c.8 is transcript index 12
    let reading = read(&hdp, "NM_X.1:c.[8=;7_9delinsTTA]").unwrap();
    assert_eq!(reading.stated(), &[(12, 13)]);
}

#[test]
fn a_change_agreeing_with_an_identity_makes_what_the_change_makes() {
    let hdp = reading_provider();
    assert_eq!(
        makes(&hdp, "NM_X.1:c.[8=;7_9delinsTTA]", C4, C12_END),
        "AAATTAGCC"
    );
}

#[test]
fn an_edit_written_twice_is_one_edit() {
    let hdp = reading_provider();
    assert!(same_sequence(
        &hdp,
        "NM_X.1:c.[7C>T;7C>T]",
        "NM_X.1:c.[7C>T]"
    ));
}

#[test]
fn a_substitution_and_its_delins_spelling_are_one_edit() {
    let hdp = reading_provider();
    assert!(same_sequence(
        &hdp,
        "NM_X.1:c.[7C>T;7delinsT]",
        "NM_X.1:c.[7C>T]"
    ));
}

#[test]
fn a_length_changing_delins_written_twice_is_one_edit() {
    let hdp = reading_provider();
    assert_eq!(
        makes(&hdp, "NM_X.1:c.[7_9delinsA;7_9delinsA]", C4, C12_END),
        "AAAAGCC"
    );
}

#[test]
fn insertions_whose_orders_make_one_sequence_are_accepted() {
    let hdp = reading_provider();
    assert_eq!(
        makes(&hdp, "NM_X.1:c.[7_8insA;7_8insAA]", C4, C12_END),
        "AAACAAATGGCC"
    );
}

#[test]
fn two_changes_to_one_base_are_contradictory() {
    assert_contradictory(&reading_provider(), "NM_X.1:c.[7C>T;7C>G]");
}

#[test]
fn a_change_inside_a_deletion_is_contradictory() {
    assert_contradictory(&reading_provider(), "NM_X.1:c.[7_9del;8T>G]");
}

#[test]
fn a_deletion_inside_a_same_length_delins_is_contradictory() {
    assert_contradictory(&reading_provider(), "NM_X.1:c.[7_8delinsTG;8del]");
}

#[test]
fn a_change_under_an_identity_is_contradictory() {
    assert_contradictory(&reading_provider(), "NM_X.1:c.[7_9=;8T>G]");
}

#[test]
fn an_insertion_inside_an_identity_is_contradictory() {
    assert_contradictory(&reading_provider(), "NM_X.1:c.[7_8=;7_8insA]");
}

#[test]
fn a_delins_that_changes_the_length_under_another_member_is_contradictory() {
    assert_contradictory(&reading_provider(), "NM_X.1:c.[7_9delinsA;8T>G]");
}

#[test]
fn an_insertion_inside_a_delins_that_changes_the_length_is_contradictory() {
    assert_contradictory(&reading_provider(), "NM_X.1:c.[7_9delinsTT;8_9insA]");
}

#[test]
fn a_member_without_a_sequence_is_held_to_its_outer_bounds() {
    assert_contradictory(&reading_provider(), "NC_X.1:g.[22C>T;20_30copy3]");
}

#[test]
fn an_unknown_outer_bound_reaches_the_end_of_the_sequence() {
    assert_contradictory(&reading_provider(), "NC_X.1:g.[(?_5)_(10_?)del;7T>A]");
}

#[test]
fn a_member_without_a_sequence_apart_from_the_rest_is_carried_unapplied() {
    let reading = read(&reading_provider(), "NC_X.1:g.[22C>T;30_35copy3]").unwrap();
    assert_eq!(
        reading.unapplied(),
        &[parse_hgvs_variant("NC_X.1:g.30_35copy3").unwrap()]
    );
}

#[test]
fn an_insertion_of_unstated_bases_is_held_to_its_junction() {
    // the bases either side of it are not inside it
    let reading = read(&reading_provider(), "NM_X.1:c.[7_8insN[3];8T>G]").unwrap();
    assert_eq!(
        reading.unapplied(),
        &[parse_hgvs_variant("NM_X.1:c.7_8insN[3]").unwrap()]
    );
}

#[test]
fn an_unknown_base_is_not_applied() {
    let reading = read(&reading_provider(), "NM_X.1:c.[7_8insN;13T>G]").unwrap();
    assert_eq!(
        reading.unapplied(),
        &[parse_hgvs_variant("NM_X.1:c.7_8insN").unwrap()]
    );
}

#[test]
fn two_insertions_at_one_junction_are_ambiguous() {
    let hdp = reading_provider();
    assert_eq!(
        readings(&hdp, "NM_X.1:c.[7_8insA;7_8insG]"),
        set(&["NM_X.1:c.7_8insAG", "NM_X.1:c.7_8insGA"])
    );
}

#[test]
fn identical_insertions_at_one_junction_are_ambiguous() {
    let hdp = reading_provider();
    assert_eq!(
        readings(&hdp, "NM_X.1:c.[7_8insA;7_8insA]"),
        set(&["NM_X.1:c.7_8insA", "NM_X.1:c.7_8insAA"])
    );
}

#[test]
fn an_insertion_of_what_a_duplication_copies_where_it_inserts_is_ambiguous() {
    let hdp = reading_provider();
    assert_eq!(
        readings(&hdp, "NM_X.1:c.[7_9dup;9_10insCTG]"),
        set(&["NM_X.1:c.7_9dup", "NM_X.1:c.9_10insCTGCTG"])
    );
}

#[test]
fn a_change_to_what_a_duplication_copies_is_ambiguous() {
    let hdp = reading_provider();
    assert_eq!(
        readings(&hdp, "NM_X.1:c.[7_9dup;8T>G]"),
        set(&["NM_X.1:c.7_8insGGC", "NM_X.1:c.8delinsGGCG"])
    );
}

#[test]
fn a_change_to_what_an_inversion_reads_is_ambiguous() {
    let hdp = reading_provider();
    assert_eq!(
        readings(&hdp, "NM_X.1:c.[7_9inv;8T>A]"),
        set(&["NM_X.1:c.8T>A", "NM_X.1:c.7_9="])
    );
}

#[test]
fn a_change_to_what_an_inversion_reads_lists_both_inversions() {
    let hdp = reading_provider();
    assert_eq!(
        readings(&hdp, "NM_X.1:c.[7_9inv;8T>G]"),
        set(&["NM_X.1:c.8T>G", "NM_X.1:c.8T>C"])
    );
}

#[test]
fn an_identity_inside_an_inversion_is_contradictory() {
    assert_contradictory(&reading_provider(), "NM_X.1:c.[7_9inv;8=]");
}

#[test]
fn inversions_that_change_nothing_may_overlap() {
    // c.9_12 is GGCC and c.10_11 GC, each its own reverse complement
    let reading = read(&reading_provider(), "NM_X.1:c.[9_12inv;10_11inv]").unwrap();
    assert_eq!(reading.stated(), &[(13, 17)]);
}

#[test]
fn an_insertion_where_a_duplication_inserts_is_ambiguous() {
    let hdp = reading_provider();
    assert_eq!(
        readings(&hdp, "NM_X.1:c.[7_9dup;9_10insA]"),
        set(&["NM_X.1:c.6_9dup", "NM_X.1:c.9_10insCTGA"])
    );
}

#[test]
fn the_error_names_both_members() {
    let err = read(&reading_provider(), "NM_X.1:c.[13T>G;7C>T;7C>G]").unwrap_err();
    let message = err.to_string();
    assert!(message.contains("NM_X.1:c.7C>T"), "{message}");
    assert!(message.contains("NM_X.1:c.7C>G"), "{message}");
}

#[test]
fn the_ambiguity_lists_every_reading() {
    // two junctions, each with two orders: four readings
    let hdp = reading_provider();
    assert_eq!(
        readings(&hdp, "NM_X.1:c.[7_8insA;7_8insG;15_16insG;15_16insT]"),
        set(&[
            "NM_X.1:c.[7_8insAG;15_16insGT]",
            "NM_X.1:c.[7_8insAG;15_16insTG]",
            "NM_X.1:c.[7_8insGA;15_16insGT]",
            "NM_X.1:c.[7_8insGA;15_16insTG]",
        ])
    );
}

#[test]
fn an_ambiguity_is_not_a_validation_error() {
    let hdp = reading_provider();
    let mapper = VariantMapper::new(&hdp);
    let ambiguous = parse_hgvs_variant("NM_X.1:c.[7_8insA;7_8insG]").unwrap();
    let err = mapper.validate(&ambiguous).unwrap_err();
    assert!(matches!(err, HgvsError::Ambiguous { .. }), "{err}");
}

#[test]
fn an_ambiguity_carries_its_readings() {
    // the deletion cancels either insertion; which pair cancels decides the
    // bases stated
    let hdp = reading_provider();
    assert_eq!(
        readings(&hdp, "NM_X.1:c.[3_4insA;5del;6_7insA]"),
        set(&["NM_X.1:c.[3_5=;6dup]", "NM_X.1:c.[4dup;5_7=]"])
    );
}

#[test]
fn each_reading_of_an_ambiguity_has_one_reading() {
    let hdp = reading_provider();
    let err = read(&hdp, "NM_X.1:c.[3_4insA;5del;6_7insA]").unwrap_err();
    let HgvsError::Ambiguous { readings, .. } = err else {
        panic!("{err}");
    };
    for reading in &readings {
        let SequenceVariant::CisPhased(c) = reading else {
            panic!("{reading} is not an allele");
        };
        assert!(VariantMapper::new(&hdp).read_cis(c).is_ok(), "{reading}");
    }
}

#[test]
fn cancelling_sets_that_share_a_member_stay_ambiguous_beside_many_other_statements() {
    // seventeen cancelling pairs apart, then the two sets sharing g.122del
    let mut members: Vec<String> = (0..17)
        .flat_map(|i| {
            [
                format!("{}_{}insA", 7 * i + 2, 7 * i + 3),
                format!("{}del", 7 * i + 3),
            ]
        })
        .collect();
    members.extend(["120_121insA", "122del", "123_124insA"].map(String::from));
    let allele = format!("NC_P.1:g.[{}]", members.join(";"));
    let err = read(&reading_provider(), &allele).unwrap_err();
    assert!(matches!(err, HgvsError::Ambiguous { .. }), "{err}");
}

#[test]
fn validate_refuses_a_contradictory_allele_rather_than_answering_false() {
    let hdp = reading_provider();
    let mapper = VariantMapper::new(&hdp);
    let contradictory = parse_hgvs_variant("NM_X.1:c.[7C>T;7C>G]").unwrap();
    assert!(matches!(
        mapper.validate(&contradictory),
        Err(HgvsError::ValidationError(_))
    ));
}

#[test]
fn validate_answers_false_for_a_wrong_base_stated_over_an_uncertain_interval() {
    let hdp = reading_provider();
    let mapper = VariantMapper::new(&hdp);
    // g.21_23 are CCC
    let wrong = parse_hgvs_variant("NC_X.1:g.[(21_23)delAAA;40T>A]").unwrap();
    assert!(!mapper.validate(&wrong).unwrap());
}

#[test]
fn a_cancelling_pair_states_the_bases_it_names() {
    // c.4_5 are transcript indices 8 and 9
    let reading = read(&reading_provider(), "NM_X.1:c.[4_5insA;5del]").unwrap();
    assert_eq!(reading.stated(), &[(8, 10)]);
}

#[test]
fn a_cancelling_pair_beside_a_change_still_states_its_bases() {
    let reading = read(&reading_provider(), "NM_X.1:c.[4_5insA;5del;13T>G]").unwrap();
    assert_eq!(reading.stated(), &[(8, 10)]);
}

#[test]
fn where_a_cancelling_pair_is_written_is_what_it_states() {
    // one place along the run of As: c.5 and c.6
    let reading = read(&reading_provider(), "NM_X.1:c.[5_6insA;6del]").unwrap();
    assert_eq!(reading.stated(), &[(9, 11)]);
}

#[test]
fn whether_members_cancel_is_read_from_the_reference() {
    // on ten Ts the pair changes g.5 to an A, and states nothing
    let reading = read(&reading_provider(), "NC_T.1:g.[4_5insA;5del]").unwrap();
    assert!(reading.stated().is_empty());
}

#[test]
fn a_cancelling_pair_around_a_change_states_nothing() {
    let reading = read(&reading_provider(), "NM_X.1:c.[4del;5A>T;6_7insA]").unwrap();
    assert!(reading.stated().is_empty());
}

#[test]
fn a_member_that_makes_the_reference_beside_a_change_states_its_bases() {
    // c.9_12 is GGCC, its own reverse complement: transcript indices 13..17
    let reading = read(&reading_provider(), "NM_X.1:c.[9_12inv;13T>G]").unwrap();
    assert_eq!(reading.stated(), &[(13, 17)]);
}

#[test]
fn many_changes_beside_one_deletion_are_searched_for_statements() {
    // c.1_17 of ATGAAACTGGCCTATCG each changed, and c.18_19 deleted
    let members: Vec<String> = "ATGAAACTGGCCTATCG"
        .chars()
        .enumerate()
        .map(|(i, b)| format!("{}{b}>{}", i + 1, if b == 'A' { 'C' } else { 'A' }))
        .chain(["18_19del".to_string()])
        .collect();
    let allele = format!("NM_X.1:c.[{}]", members.join(";"));
    assert!(read(&reading_provider(), &allele).is_ok());
}

#[test]
fn too_many_orders_to_list_is_refused_rather_than_overflowing() {
    // six insertions at each of seven junctions: 720 orders at each
    let members: Vec<String> = (0..7)
        .flat_map(|j| {
            ["A", "C", "G", "T", "AC", "AG"].map(|b| format!("{}_{}insA{b}", 2 * j + 2, 2 * j + 3))
        })
        .collect();
    let allele = format!("NM_X.1:c.[{}]", members.join(";"));
    let err = read(&reading_provider(), &allele).unwrap_err();
    assert!(matches!(err, HgvsError::UnsupportedOperation(_)), "{err}");
}

#[test]
fn a_duplication_of_an_exons_last_base_is_read() {
    // c.5 is the last base of exon 1 of NM_S.1; c.4_9 are indices 8..14
    let hdp = reading_provider();
    assert_eq!(makes(&hdp, "NM_S.1:c.[5dup;6A>T]", C4, 14), "AAATCTG");
}

#[test]
fn a_minus_strand_allele_is_read_on_its_transcript() {
    // the A goes before the copied AGG, as in c.[7_9dup;6_7insA]
    let hdp = reading_provider();
    assert!(same_sequence(
        &hdp,
        "NR_XR.1:n.[23_25dup;22_23insA]",
        "NR_XR.1:n.[22_23insAAGG]"
    ));
}

#[test]
fn an_allele_with_an_intronic_member_is_read_on_the_genome() {
    let reading = read(&reading_provider(), "NM_S.1:c.[5+2T>A;7C>T]").unwrap();
    assert_eq!(reading.accession(), "NC_S.1");
}

#[test]
fn an_intronic_member_changes_no_other_members_reading_on_the_minus_strand() {
    // the copy goes 3' of n.25_27 on the transcript; n.25 is kept
    let hdp = reading_provider();
    assert!(read(&hdp, "NR_MS.1:n.[25_27dup;24_25delinsTG]").is_ok());
    assert!(read(&hdp, "NR_MS.1:n.[25_27dup;24_25delinsTG;30+2=]").is_ok());
}

#[test]
fn a_transcript_allele_read_on_the_genome_makes_its_projection() {
    let hdp = reading_provider();
    let mapper = VariantMapper::new(&hdp);
    let on_genome = mapper
        .read_cis_on_genome(&cis("NM_X.1:c.[7C>T;13T>G]"))
        .unwrap();
    let projection = read(&hdp, "NC_X.1:g.[22C>T;28T>G]").unwrap();
    assert!(same(&hdp, &on_genome, &projection));
}

#[test]
fn readings_on_different_sequences_are_not_compared() {
    let hdp = reading_provider();
    let on_transcript = read(&hdp, "NM_X.1:c.[7C>T;13T>G]").unwrap();
    let on_genome = read(&hdp, "NC_X.1:g.[22C>T;28T>G]").unwrap();
    let mapper = VariantMapper::new(&hdp);
    let reference = mapper
        .refs
        .reference(on_transcript.accession(), on_transcript.kind());
    assert!(matches!(
        on_transcript.same_sequence(&on_genome, &reference),
        Err(HgvsError::ValidationError(_))
    ));
}

#[test]
fn an_rna_allele_is_read_as_its_coding_spelling() {
    assert_contradictory(&reading_provider(), "NM_X.1:r.[7c>u;7c>g]");
}

#[test]
fn the_readings_of_an_rna_allele_are_written_in_rna() {
    let hdp = reading_provider();
    assert_eq!(
        readings(&hdp, "NM_X.1:r.[7_8insa;7_8insg]"),
        set(&["NM_X.1:r.7_8insag", "NM_X.1:r.7_8insga"])
    );
}

#[test]
fn a_protein_allele_is_not_read_into_a_sequence() {
    let err = read(&reading_provider(), "NP_X.1:p.[Lys2Leu;Ala4del]").unwrap_err();
    assert!(matches!(err, HgvsError::UnsupportedOperation(_)), "{err}");
}

#[test]
fn to_vrs_refuses_an_ambiguous_allele() {
    let hdp = reading_provider();
    let mapper = VariantMapper::with_refget(&hdp, &hdp);
    let ambiguous = parse_hgvs_variant("NM_X.1:c.[7_8insA;7_8insG]").unwrap();
    assert!(matches!(
        mapper.to_vrs_cis_phased(&ambiguous),
        Err(HgvsError::Ambiguous { .. })
    ));
}

#[test]
fn identical_insertions_of_unstated_bases_are_refused_as_unlisted() {
    let err = read(
        &reading_provider(),
        "NC_X.1:g.[25_26insN[20];25_26insN[20]]",
    )
    .unwrap_err();
    assert!(matches!(err, HgvsError::UnsupportedOperation(_)), "{err}");
}

#[test]
fn an_uncertain_duplication_reaches_the_junction_its_copy_may_go_to() {
    // (20_22)dup may copy up to g.22, and so insert between g.22 and g.23
    let err = read(&reading_provider(), "NC_X.1:g.[(20_22)dup;22_23insA]").unwrap_err();
    assert!(matches!(err, HgvsError::UnsupportedOperation(_)), "{err}");
}

#[test]
fn a_contradiction_is_reported_before_unlisted_readings_in_any_order() {
    let hdp = reading_provider();
    assert_contradictory(&hdp, "NC_X.1:g.[25_26insN;24_27del;40_45copy3;42T>A]");
    assert_contradictory(&hdp, "NC_X.1:g.[40_45copy3;42T>A;25_26insN;24_27del]");
}

#[test]
fn three_identical_insertions_are_one_two_or_three() {
    let hdp = reading_provider();
    assert_eq!(
        readings(&hdp, "NM_X.1:c.[7_8insA;7_8insA;7_8insA]"),
        set(&[
            "NM_X.1:c.7_8insA",
            "NM_X.1:c.7_8insAA",
            "NM_X.1:c.7_8insAAA"
        ])
    );
}

#[test]
fn identical_copy_numbers_count_once() {
    let reading = read(&reading_provider(), "NC_X.1:g.[30_35copy3;30_35copy3]").unwrap();
    assert_eq!(
        reading.unapplied(),
        &[parse_hgvs_variant("NC_X.1:g.30_35copy3").unwrap()]
    );
}

#[test]
fn members_combined_in_too_many_ways_to_search_are_refused() {
    let mut members: Vec<String> = (1..=4).map(|p| format!("{p}del")).collect();
    members.extend((5..=20).map(|p| {
        let base = &"ACGTTGCAAGCTTCGA"[(p - 1) % 16..(p - 1) % 16 + 1];
        format!("{p}{base}>{}", if base == "A" { "C" } else { "A" })
    }));
    members.extend((17..=20).map(|p| format!("{p}_{}insA", p + 1)));
    let allele = format!("NC_Q.1:g.[{}]", members.join(";"));
    let err = read(&reading_provider(), &allele).unwrap_err();
    assert!(matches!(err, HgvsError::UnsupportedOperation(_)), "{err}");
}

#[test]
fn a_reading_is_written_on_each_side_of_a_statement() {
    let hdp = reading_provider();
    assert_eq!(
        readings(&hdp, "NM_X.1:c.[4_5insC;4_5insG;5=;6del]"),
        set(&["NM_X.1:c.[4_5insCG;5=;6del]", "NM_X.1:c.[4_5insGC;5=;6del]"])
    );
}

#[test]
fn every_written_reading_has_one_reading_where_changes_shift_into_each_other() {
    let hdp = reading_provider();
    for allele in [
        "NC_U.1:g.[7_8insA;7_8insTA;7_8insT;11_13delinsTA]",
        "NC_V.1:g.[8_9insAA;3_4inv;8_10dup]",
    ] {
        let err = read(&hdp, allele).unwrap_err();
        let HgvsError::Ambiguous { readings, .. } = err else {
            panic!("{allele}: {err}");
        };
        for reading in readings {
            if let SequenceVariant::CisPhased(c) = &reading {
                assert!(
                    VariantMapper::new(&hdp).read_cis(c).is_ok(),
                    "{allele}: {reading}"
                );
            }
        }
    }
}

#[test]
fn an_insertion_on_an_exon_junction_read_on_the_genome_is_refused_as_unplaced() {
    let err = read(&reading_provider(), "NM_S.1:c.[5+1G>T;5_6insA]").unwrap_err();
    assert!(matches!(err, HgvsError::UnsupportedOperation(_)), "{err}");
}

#[test]
fn a_copy_number_reaches_the_junctions_its_copies_go_to() {
    let err = read(&reading_provider(), "NC_X.1:g.[20_30copy2;30_31insA]").unwrap_err();
    assert!(matches!(err, HgvsError::UnsupportedOperation(_)), "{err}");
}

#[test]
fn unstated_bases_inserted_inside_an_identity_are_contradictory() {
    assert_contradictory(&reading_provider(), "NM_X.1:c.[7_9=;8_9insN]");
}

#[test]
fn an_identity_inside_a_copy_number_is_accepted() {
    assert!(read(&reading_provider(), "NC_X.1:g.[20_30copy3;25=]").is_ok());
}

#[test]
fn identical_insertions_at_two_junctions_list_their_readings_rather_than_too_many() {
    // five As at each of two places in a run of As
    let members: Vec<String> = [30, 40]
        .iter()
        .flat_map(|j| (0..5).map(move |_| format!("{j}_{}insA", j + 1)))
        .collect();
    let allele = format!("NC_A.1:g.[{}]", members.join(";"));
    let hdp = reading_provider().sequence("NC_A.1", &format!("C{}G", "A".repeat(60)));
    let err = read(&hdp, &allele).unwrap_err();
    assert!(matches!(err, HgvsError::Ambiguous { .. }), "{err}");
}

#[test]
fn a_member_without_a_sequence_is_in_every_reading() {
    let hdp = reading_provider();
    assert_eq!(
        readings(&hdp, "NC_W.1:g.[15_16copy3;7_8insA;7_8insT]"),
        set(&[
            "NC_W.1:g.[6_7dup;15_16copy3]",
            "NC_W.1:g.[8_9insAT;15_16copy3]"
        ])
    );
}

#[test]
fn a_reading_that_changes_an_identity_is_no_reading() {
    // inverted from the reference, g.7_8 AT stays AT and then g.7 becomes T,
    // which keeps g.8; inverting the changed TT would make g.8 an A
    let reading = read(&reading_provider(), "NC_Y.1:g.[7A>T;8=;7_8inv]").unwrap();
    assert_eq!(reading.stated(), &[(7, 8)]);
}

#[test]
fn statements_either_side_of_an_insertion_are_written_apart() {
    let hdp = reading_provider();
    for reading in readings(&hdp, "NC_Z.1:g.[2_3insT;2_3insA;2=;3=]") {
        assert!(
            reading.contains("2=") && reading.contains(";3="),
            "{reading}"
        );
    }
}

#[test]
fn an_identity_over_the_edge_of_a_copy_number_is_accepted() {
    assert!(read(&reading_provider(), "NC_Z.1:g.[2_3=;3_4copy3]").is_ok());
}

#[test]
fn overlapping_inversions_that_change_each_others_bases_are_not_yet_listed() {
    // c.7_9 is CTG and c.8_10 TGG: each inversion changes c.8
    let err = read(&reading_provider(), "NM_X.1:c.[7_9inv;8_10inv]").unwrap_err();
    assert!(matches!(err, HgvsError::UnsupportedOperation(_)), "{err}");
}

#[test]
fn a_duplication_and_an_inversion_reading_each_others_changes_are_not_yet_listed() {
    let err = read(&reading_provider(), "NM_X.1:c.[7_9inv;7_8dup]").unwrap_err();
    assert!(matches!(err, HgvsError::UnsupportedOperation(_)), "{err}");
}

// Normalising an allele: choices.md, "Members of a cis allele that meet are
// written as one change".

/// `reading_provider` with short genomes whose runs the rule is about:
/// `NC_D.1` is `TCTCT`, `NC_E.1` `CTCTT`, `NC_F.1` `TTCAGCAGTT`, `NC_G.1`
/// `AAATAAAA`, `NC_H.1` `ACACCCCGC`, `NC_K.1` `GAAAAAG` and `NC_L.1`
/// nineteen As.
fn normalising_provider() -> Provider {
    reading_provider()
        .sequence("NC_D.1", "TCTCT")
        .sequence("NC_E.1", "CTCTT")
        .sequence("NC_F.1", "TTCAGCAGTT")
        .sequence("NC_G.1", "AAATAAAA")
        .sequence("NC_H.1", "ACACCCCGC")
        .sequence("NC_K.1", "GAAAAAG")
        .sequence("NC_L.1", &"A".repeat(19))
}

fn normalise(hdp: &Provider, hgvs: &str) -> Result<SequenceVariant, HgvsError> {
    VariantMapper::new(hdp).normalize_variant(parse_hgvs_variant(hgvs).unwrap())
}

fn normalised(hdp: &Provider, hgvs: &str) -> String {
    normalise(hdp, hgvs)
        .unwrap_or_else(|e| panic!("{hgvs}: {e}"))
        .to_string()
}

#[test]
fn members_in_one_run_are_written_as_one_change() {
    let hdp = normalising_provider();
    for (allele, expected) in [
        ("NC_D.1:g.[1del;2del;5del]", "NC_D.1:g.3_5del"),
        ("NC_E.1:g.[1del;2del;4del]", "NC_E.1:g.3_5del"),
        ("NC_F.1:g.[3_5del;7A>T]", "NC_F.1:g.4_7delinsT"),
        ("NC_F.1:g.[4A>T;6_8del]", "NC_F.1:g.4_7delinsT"),
        ("NC_G.1:g.[1_3del;4del;6A>T]", "NC_G.1:g.2_6delinsT"),
    ] {
        assert_eq!(normalised(&hdp, allele), expected, "{allele}");
    }
}

#[test]
fn members_justified_over_one_run_are_merged_before_either_is_shifted() {
    // both deletions justify over g.2_6 of GAAAAAG
    let hdp = normalising_provider();
    assert_eq!(
        normalised(&hdp, "NC_K.1:g.[2del;4A>T]"),
        "NC_K.1:g.3_4delinsT"
    );
    assert_eq!(
        normalised(&hdp, "NC_K.1:g.[6del;4A>T]"),
        "NC_K.1:g.4_5delinsT"
    );
}

#[test]
fn an_extent_holds_the_bases_a_member_is_written_over() {
    // the delins changes only g.9, but is written over g.7_9, which the run
    // of Cs the deletion is in reaches
    assert_eq!(
        normalised(&normalising_provider(), "NC_H.1:g.[4del;7_9delinsCGG]"),
        "NC_H.1:g.7_9delinsGG"
    );
}

#[test]
fn touching_members_are_one_delins() {
    assert_eq!(
        normalised(&normalising_provider(), "NM_X.1:c.[7C>T;8T>G;13T>G]"),
        "NM_X.1:c.[7_8delinsTG;13T>G]"
    );
}

#[test]
fn an_allele_merged_whole_is_a_plain_variant() {
    let merged = normalise(&normalising_provider(), "NM_X.1:c.[7C>T;8T>G]").unwrap();
    assert!(
        matches!(&merged, SequenceVariant::Coding(_)),
        "{merged} is not a plain variant"
    );
}

#[test]
fn a_group_can_reduce_to_a_deletion() {
    // shifted one at a time, the A deleted would stop at c.7, which the
    // substitution makes an A
    assert_eq!(
        normalised(&normalising_provider(), "NM_X.1:c.[4del;7C>A]"),
        "NM_X.1:c.7del"
    );
}

#[test]
fn a_group_can_reduce_to_a_duplication() {
    assert_eq!(
        normalised(&normalising_provider(), "NM_X.1:c.[7C>A;7_8insC]"),
        "NM_X.1:c.6dup"
    );
}

#[test]
fn members_with_an_unchanged_base_between_stay_apart() {
    assert_eq!(
        normalised(&normalising_provider(), "NM_X.1:c.[5A>G;9G>T]"),
        "NM_X.1:c.[5A>G;9G>T]"
    );
}

#[test]
fn members_in_one_codon_are_one_delins() {
    assert_eq!(
        normalised(&normalising_provider(), "NM_X.1:c.[4A>T;6A>T]"),
        "NM_X.1:c.4_6delinsTAT"
    );
}

#[test]
fn members_one_base_apart_in_two_codons_stay_apart() {
    assert_eq!(
        normalised(&normalising_provider(), "NM_X.1:c.[6A>T;8T>G]"),
        "NM_X.1:c.[6A>T;8T>G]"
    );
}

#[test]
fn grouping_runs_across_touch_and_codon() {
    assert_eq!(
        normalised(&normalising_provider(), "NM_X.1:c.[4A>T;6A>T;7C>G]"),
        "NM_X.1:c.4_7delinsTATG"
    );
}

#[test]
fn a_genomic_allele_has_no_codon() {
    // g.19 and g.21 are c.4 and c.6
    assert_eq!(
        normalised(&normalising_provider(), "NC_X.1:g.[19A>T;21A>T]"),
        "NC_X.1:g.[19A>T;21A>T]"
    );
}

#[test]
fn members_outside_the_cds_have_no_codon() {
    assert_eq!(
        normalised(&normalising_provider(), "NM_X.1:c.[*1C>A;*3G>T]"),
        "NM_X.1:c.[*1C>A;*3G>T]"
    );
}

#[test]
fn an_insertion_inside_a_codon_joins_it() {
    assert_eq!(
        normalised(&normalising_provider(), "NM_X.1:c.[4A>T;5_6insG]"),
        "NM_X.1:c.4_5delinsTAG"
    );
}

#[test]
fn an_insertion_between_two_codons_joins_neither() {
    assert_eq!(
        normalised(&normalising_provider(), "NM_X.1:c.[5A>T;6_7insG]"),
        "NM_X.1:c.[5A>T;6_7insG]"
    );
}

#[test]
fn a_deletion_across_two_codons_groups_by_its_region() {
    // c.9 and c.10 are GG: either deletion is justified over both codons
    let hdp = normalising_provider();
    assert_eq!(
        normalised(&hdp, "NM_X.1:c.[10del;12C>T]"),
        "NM_X.1:c.10_12delinsCT"
    );
    assert_eq!(
        normalised(&hdp, "NM_X.1:c.[9del;12C>T]"),
        "NM_X.1:c.10_12delinsCT"
    );
}

#[test]
fn members_in_a_codon_split_by_an_intron_stay_apart() {
    assert_eq!(
        normalised(&normalising_provider(), "NM_S.1:c.[4A>T;6A>T]"),
        "NM_S.1:c.[4A>T;6A>T]"
    );
}

#[test]
fn members_touching_across_a_junction_stay_apart() {
    assert_eq!(
        normalised(&normalising_provider(), "NM_S.1:c.[5A>T;6A>T]"),
        "NM_S.1:c.[5A>T;6A>T]"
    );
}

#[test]
fn members_touching_in_one_exon_beside_a_junction_are_one_change() {
    assert_eq!(
        normalised(&normalising_provider(), "NM_S.1:c.[6A>T;7C>G]"),
        "NM_S.1:c.6_7delinsTG"
    );
}

#[test]
fn a_deletion_in_an_allele_stops_at_its_exons_last_base() {
    let hdp = normalising_provider();
    assert_eq!(
        normalised(&hdp, "NM_S.1:c.[4del;7C>A]"),
        "NM_S.1:c.[5del;7C>A]"
    );
    assert_eq!(
        normalised(&hdp, "NM_S.1:c.[4del;6A>T]"),
        "NM_S.1:c.[5del;6A>T]"
    );
}

#[test]
fn a_duplication_in_an_allele_stops_at_its_exons_last_base() {
    assert_eq!(
        normalised(&normalising_provider(), "NM_S.1:c.[4dup;6A>T]"),
        "NM_S.1:c.[5dup;6A>T]"
    );
}

#[test]
fn a_shift_stops_short_of_an_identity() {
    assert_eq!(
        normalised(&normalising_provider(), "NM_X.1:c.[4del;6=]"),
        "NM_X.1:c.[5del;6=]"
    );
}

#[test]
fn a_change_against_an_identity_does_not_move() {
    assert_eq!(
        normalised(&normalising_provider(), "NM_X.1:c.[4del;5=]"),
        "NM_X.1:c.[4del;5=]"
    );
}

#[test]
fn members_are_written_in_order_of_position() {
    let hdp = normalising_provider();
    assert_eq!(
        normalised(&hdp, "NM_X.1:c.[13T>G;7C>T]"),
        "NM_X.1:c.[7C>T;13T>G]"
    );
    // along the transcript, whichever strand it is on
    assert_eq!(
        normalised(&hdp, "NR_XR.1:n.[29G>A;23A>C]"),
        "NR_XR.1:n.[23A>C;29G>A]"
    );
}

#[test]
fn a_non_coding_allele_is_normalised_on_its_transcript() {
    // n.3_4 of NR_XR.1 are TT
    assert_eq!(
        normalised(&normalising_provider(), "NR_XR.1:n.[3T>A;4T>G]"),
        "NR_XR.1:n.3_4delinsAG"
    );
}

#[test]
fn an_rna_allele_is_normalised_as_its_coding_spelling() {
    assert_eq!(
        normalised(&normalising_provider(), "NM_X.1:r.[4a>u;6a>u]"),
        "NM_X.1:r.4_6delinsuau"
    );
}

#[test]
fn a_cancelling_allele_is_the_identity_of_the_bases_its_members_name() {
    let hdp = normalising_provider();
    assert_eq!(normalised(&hdp, "NM_X.1:c.[4_5insA;5del]"), "NM_X.1:c.4_5=");
    // one place along the run of As: another statement
    assert_eq!(normalised(&hdp, "NM_X.1:c.[5_6insA;6del]"), "NM_X.1:c.5_6=");
}

#[test]
fn members_that_cancel_on_one_reference_change_a_base_on_another() {
    assert_eq!(
        normalised(&normalising_provider(), "NC_T.1:g.[4_5insA;5del]"),
        "NC_T.1:g.5T>A"
    );
}

#[test]
fn a_cancelling_group_beside_a_change_leaves_an_identity() {
    assert_eq!(
        normalised(&normalising_provider(), "NM_X.1:c.[4_5insA;5del;13T>G]"),
        "NM_X.1:c.[4_5=;13T>G]"
    );
}

#[test]
fn cancelling_members_that_would_meet_a_change_are_still_a_statement() {
    assert_eq!(
        normalised(&normalising_provider(), "NM_X.1:c.[4_5insA;5del;6A>T]"),
        "NM_X.1:c.[4_5=;6A>T]"
    );
}

#[test]
fn cancelling_members_around_a_change_are_part_of_it() {
    assert_eq!(
        normalised(&normalising_provider(), "NM_X.1:c.[4del;5A>T;6_7insA]"),
        "NM_X.1:c.4A>T"
    );
}

#[test]
fn a_single_member_that_makes_the_reference_beside_a_change_is_a_statement() {
    // c.9_12 is GGCC, its own reverse complement
    assert_eq!(
        normalised(&normalising_provider(), "NM_X.1:c.[9_12inv;13T>G]"),
        "NM_X.1:c.[9_12=;13T>G]"
    );
}

#[test]
fn a_change_is_cut_at_an_agreeing_identity() {
    assert_eq!(
        normalised(&normalising_provider(), "NM_X.1:c.[7=;7_8delinsCG]"),
        "NM_X.1:c.[7=;8T>G]"
    );
}

#[test]
fn a_change_agreeing_with_an_identity_inside_it_is_cut_there() {
    assert_eq!(
        normalised(&normalising_provider(), "NM_X.1:c.[8=;7_9delinsTTA]"),
        "NM_X.1:c.[7C>T;8=;9G>A]"
    );
}

#[test]
fn identities_that_meet_are_written_as_one() {
    let hdp = normalising_provider();
    assert_eq!(normalised(&hdp, "NM_X.1:c.[7_9=;8=]"), "NM_X.1:c.7_9=");
    assert_eq!(normalised(&hdp, "NM_X.1:c.[7_9=;10=]"), "NM_X.1:c.7_10=");
}

#[test]
fn a_cancelling_group_beside_an_identity_is_written_as_one_statement() {
    assert_eq!(
        normalised(&normalising_provider(), "NM_X.1:c.[4_5insA;5del;6=]"),
        "NM_X.1:c.4_6="
    );
}

#[test]
fn a_part_cut_at_an_identity_that_changes_nothing_is_dropped() {
    // c.7_9 is CTG: the delins changes only c.9
    assert_eq!(
        normalised(&normalising_provider(), "NM_X.1:c.[8=;7_9delinsCTA]"),
        "NM_X.1:c.[8=;9G>A]"
    );
}

#[test]
fn a_change_every_part_of_which_agrees_with_the_reference_is_a_statement() {
    assert_eq!(
        normalised(&normalising_provider(), "NM_X.1:c.[8=;7_9delinsCTG]"),
        "NM_X.1:c.7_9="
    );
}

#[test]
fn members_sharing_bases_across_an_exon_junction_are_refused() {
    // c.5_6 of NM_S.1 spans the junction; one change there deletes the intron
    let err = normalise(&normalising_provider(), "NM_S.1:c.[5_6delinsTT;5A>T]").unwrap_err();
    assert!(
        matches!(&err, HgvsError::UnsupportedOperation(m) if m.contains("across an exon junction")),
        "{err}"
    );
}

#[test]
fn a_member_across_an_exon_junction_stops_a_neighbours_shift() {
    // c.4del would shift to c.5, inside the inversion of c.5_6
    let hdp = normalising_provider();
    assert_eq!(
        normalised(&hdp, "NM_S.1:c.[4del;5_6inv]"),
        "NM_S.1:c.[4del;5_6inv]"
    );
}

#[test]
fn an_ambiguity_the_rule_cannot_write_still_lists_its_readings() {
    let hdp = normalising_provider();
    for allele in [
        "NR_MS.1:n.[30_32del;32dup]",
        "NM_S.1:c.[5_7dup;4_5insC;7_9dup;6_7=]",
    ] {
        let err = read(&hdp, allele).unwrap_err();
        assert!(
            matches!(&err, HgvsError::Ambiguous { readings, .. } if readings.len() == 2),
            "{allele}: {err}"
        );
    }
}

#[test]
fn members_beside_a_member_across_an_exon_junction_do_not_shift_onto_it() {
    let hdp = normalising_provider();
    for allele in [
        "NM_S.1:c.[4dup;5_6delinsAC]",
        "NM_S.1:c.[4del;5_6dup]",
        "NM_S.1:c.[3_4insA;4_6dup]",
    ] {
        assert_eq!(normalised(&hdp, allele), allele);
    }
}

#[test]
fn a_redundant_member_agreeing_with_a_cut_change_does_not_stop_the_cut() {
    // c.7C>T says again what the delins says of c.7
    let hdp = normalising_provider();
    let once = normalised(&hdp, "NM_X.1:c.[7_9delinsTTA;8=;7C>T]");
    assert_eq!(once, "NM_X.1:c.[7C>T;8=;9G>A]");
}

#[test]
fn a_change_cut_at_an_identity_states_nothing_a_restating_member_changes() {
    // g.22 is C: the delins changes it, as the substitution says again
    let hdp = normalising_provider();
    let reading = read(&hdp, "NC_X.1:g.[22C>T;22_23delinsTT;23=]").unwrap();
    assert_eq!(reading.stated(), &[(22, 23)]);
    assert_eq!(
        normalised(&hdp, "NC_X.1:g.[22C>T;22_23delinsTT;23=]"),
        "NC_X.1:g.[22C>T;23=]"
    );
}

#[test]
fn a_transcript_allele_is_not_written_as_another_genome() {
    // read on the transcript, the delins takes in the copy; on the genome
    // the copy stays, so no transcript form says what the allele does
    let err = normalise(
        &normalising_provider(),
        "NM_S.1:c.[5_6delinsAC;3_5dup;3_5=]",
    )
    .unwrap_err();
    assert!(
        matches!(&err, HgvsError::UnsupportedOperation(m) if m.contains("on the genome")),
        "{err}"
    );
}

#[test]
fn an_intronic_member_written_twice_is_one_member() {
    assert_eq!(
        normalised(
            &normalising_provider(),
            "NM_S.1:c.[5delinsAA;5+1G>A;5+1G>A]"
        ),
        "NM_S.1:c.5+1delinsAA"
    );
}

#[test]
fn an_intronic_edit_written_two_ways_is_written_once_whatever_the_order() {
    let hdp = normalising_provider();
    assert_eq!(
        normalised(&hdp, "NM_S.1:c.[5+5delinsT;5+5G>T]"),
        "NM_S.1:c.5+5G>T"
    );
    assert_eq!(
        normalised(&hdp, "NM_S.1:c.[5+5G>T;5+5delinsT]"),
        "NM_S.1:c.5+5G>T"
    );
}

#[test]
fn an_exonic_member_sharing_bases_with_an_intronic_group_joins_it() {
    // c.7_9del shares c.9 with the inversion and meets the duplication,
    // whose copy reaches into the intron
    let hdp = normalising_provider();
    for allele in [
        "NM_S.1:c.[9_10inv;3_6-2dup;7_9del]",
        "NM_S.1:c.[7_9del;9_10inv;3_6-2dup]",
    ] {
        assert_eq!(
            normalised(&hdp, allele),
            normalised(&hdp, "NM_S.1:c.[3_6-2dup;7_9del;9_10inv]"),
            "{allele}"
        );
    }
}

#[test]
fn an_identity_joins_no_intronic_group() {
    assert_eq!(
        normalised(&normalising_provider(), "NM_S.1:c.[5A>T;5+1=]"),
        "NM_S.1:c.[5A>T;5+1=]"
    );
}

#[test]
fn an_intronic_group_can_be_written_in_the_exon() {
    // the copy of c.5+1 sits after c.5, still in exon 1 on the genome
    assert_eq!(
        normalised(&normalising_provider(), "NM_S.1:c.[5+1dup;5A>T]"),
        "NM_S.1:c.5delinsTG"
    );
}

#[test]
fn members_sharing_bases_across_a_statement_are_refused() {
    let err = normalise(
        &normalising_provider(),
        "NM_X.1:c.[6=;6_7delinsAG;5_6delinsTA]",
    )
    .unwrap_err();
    assert!(
        matches!(&err, HgvsError::UnsupportedOperation(m) if m.contains("across a statement")),
        "{err}"
    );
}

#[test]
fn an_intronic_member_a_merged_change_reaches_joins_it() {
    // c.3_5 of NM_S.1 is GAA: deleted as one, it justifies up to c.5+1,
    // beside c.5+2, which neither deletion reached alone
    assert_eq!(
        normalised(&normalising_provider(), "NM_S.1:c.[3del;4_5del;5+2T>A]"),
        "NM_S.1:c.5_5+2del"
    );
}

#[test]
fn members_sharing_bases_across_an_exon_junction_still_read() {
    assert!(read(&normalising_provider(), "NM_S.1:c.[4_6del;6del]").is_ok());
}

#[test]
fn an_ambiguity_across_an_exon_junction_lists_its_readings() {
    let hdp = normalising_provider();
    for allele in ["NM_S.1:c.[4_6inv;5A>T]", "NM_S.1:c.[5_6del;5_6insAA]"] {
        let err = read(&hdp, allele).unwrap_err();
        assert!(
            matches!(&err, HgvsError::Ambiguous { readings, .. } if readings.len() == 2),
            "{allele}: {err}"
        );
    }
}

#[test]
fn a_cut_member_that_makes_the_reference_states_its_bases_beside_a_change() {
    // c.6_8 is ACT: the delins changes nothing; c.5 touches it
    assert_eq!(
        normalised(&normalising_provider(), "NM_X.1:c.[7=;6_8delinsACT;5A>T]"),
        "NM_X.1:c.[5A>T;6_8=]"
    );
}

#[test]
fn a_merged_change_is_not_written_where_unstated_bases_are_inserted() {
    // trimmed, the duplication and the substitution would insert at g.1_2,
    // where the N[2] goes
    let hdp = normalising_provider();
    let normalised = normalised(&hdp, "NC_L.1:g.[3_6dup;2A>G;1_2insN[2]]");
    assert!(
        VariantMapper::new(&hdp).read_cis(&cis(&normalised)).is_ok(),
        "{normalised}"
    );
}

#[test]
fn identities_apart_stay_apart() {
    assert_eq!(
        normalised(&normalising_provider(), "NM_X.1:c.[10=;20=]"),
        "NM_X.1:c.[10=;20=]"
    );
}

#[test]
fn an_identity_across_a_junction_is_not_merged_with_one_inside_an_exon() {
    // c.4_6= across the junction projects over the intron; c.4_5= does not
    assert_eq!(
        normalised(&normalising_provider(), "NM_S.1:c.[4_6=;4_5=]"),
        "NM_S.1:c.[4_5=;4_6=]"
    );
}

#[test]
fn identities_in_two_exons_stay_apart() {
    assert_eq!(
        normalised(&normalising_provider(), "NM_S.1:c.[5=;6=]"),
        "NM_S.1:c.[5=;6=]"
    );
}

#[test]
fn a_member_without_a_sequence_is_kept_as_written_in_its_place() {
    assert_eq!(
        normalised(&normalising_provider(), "NC_X.1:g.[30_35copy3;19del]"),
        "NC_X.1:g.[21del;30_35copy3]"
    );
}

#[test]
fn an_intronic_member_apart_from_the_exons_is_kept_as_written() {
    // c.5+2 is g.21, inside the intron; c.4del shifts to c.5
    assert_eq!(
        normalised(&normalising_provider(), "NM_S.1:c.[7C>T;5+2T>A;4del]"),
        "NM_S.1:c.[5del;5+2T>A;7C>T]"
    );
}

#[test]
fn an_intronic_member_and_an_exonic_member_it_meets_are_one_change() {
    // c.5 is g.20 and c.5+1 g.21: the two touch on the genome
    assert_eq!(
        normalised(&normalising_provider(), "NM_S.1:c.[5A>T;5+1G>C]"),
        "NM_S.1:c.5_5+1delinsTC"
    );
}

#[test]
fn an_intronic_change_is_not_shifted() {
    // c.5+3 and c.5+4 are AA: the deletion of c.5+3 is written where it is,
    // with the substitution it touches
    assert_eq!(
        normalised(&normalising_provider(), "NM_S.1:c.[5+2T>C;5+3del]"),
        "NM_S.1:c.5+2_5+3delinsC"
    );
}

#[test]
fn intronic_members_on_the_minus_strand_are_one_change() {
    assert_eq!(
        normalised(&normalising_provider(), "NR_MS.1:n.[30+1C>A;30T>G]"),
        "NR_MS.1:n.30_30+1inv"
    );
}

#[test]
fn an_rna_allele_with_intronic_members_is_written_in_rna() {
    assert_eq!(
        normalised(&normalising_provider(), "NM_S.1:r.[5a>u;5+1g>c]"),
        "NM_S.1:r.5_5+1delinsuc"
    );
}

#[test]
fn intronic_members_that_meet_are_one_change() {
    assert_eq!(
        normalised(&normalising_provider(), "NM_S.1:c.[5+2T>A;5+3A>G]"),
        "NM_S.1:c.5+2_5+3delinsAG"
    );
}

#[test]
fn a_splice_site_member_that_meets_nothing_is_kept_as_written() {
    assert_eq!(
        normalised(&normalising_provider(), "NM_S.1:c.[5+1G>C;13T>G]"),
        "NM_S.1:c.[5+1G>C;13T>G]"
    );
}

#[test]
fn an_allele_with_a_splice_site_member_renders_as_a_block() {
    let hdp = normalising_provider();
    assert!(matches!(
        vrs(&hdp, "NM_S.1:c.[5+1G>C;13T>G]"),
        VrsVariation::CisPhasedBlock(_)
    ));
}

#[test]
fn normalising_an_ambiguous_allele_raises_with_its_readings() {
    let err = normalise(&normalising_provider(), "NM_X.1:c.[3_4insA;5del;6_7insA]").unwrap_err();
    let HgvsError::Ambiguous { readings, .. } = err else {
        panic!("{err}");
    };
    let readings: BTreeSet<String> = readings.iter().map(ToString::to_string).collect();
    assert_eq!(
        readings,
        set(&["NM_X.1:c.[3_5=;6dup]", "NM_X.1:c.[4dup;5_7=]"])
    );
}

#[test]
fn normalising_a_contradictory_allele_raises() {
    let err = normalise(&normalising_provider(), "NM_X.1:c.[7C>T;7C>G]").unwrap_err();
    assert!(
        matches!(&err, HgvsError::ValidationError(m) if m.contains("has no reading")),
        "{err}"
    );
}

#[test]
fn normalising_an_allele_whose_readings_are_not_listed_raises() {
    let err = normalise(&normalising_provider(), "NM_X.1:c.[7_9inv;8_10inv]").unwrap_err();
    assert!(
        matches!(&err, HgvsError::UnsupportedOperation(m) if m.contains("readings are not listed")),
        "{err}"
    );
}

// to_vrs renders the normalised allele: choices.md, "A cis allele is the
// sequence its members make" and "CisPhasedBlock members are sorted before
// digesting".

fn vrs(hdp: &Provider, hgvs: &str) -> VrsVariation {
    VariantMapper::with_refget(hdp, hdp)
        .to_vrs_variation(&parse_hgvs_variant(hgvs).unwrap())
        .unwrap_or_else(|e| panic!("{hgvs}: {e}"))
}

#[test]
fn an_allele_split_into_other_alleles_has_another_vrs_identifier() {
    let hdp = normalising_provider();
    let touching = vrs(&hdp, "NM_X.1:c.[7C>T;8T>G]");
    assert!(matches!(touching, VrsVariation::Allele(_)), "{touching:?}");
    assert_eq!(touching.id(), vrs(&hdp, "NM_X.1:c.7_8delinsTG").id());
    let apart = vrs(&hdp, "NM_X.1:c.[5A>G;9G>T]");
    assert!(
        matches!(apart, VrsVariation::CisPhasedBlock(_)),
        "{apart:?}"
    );
    // the same change written as one delins: another identifier
    assert_ne!(apart.id(), vrs(&hdp, "NM_X.1:c.5_9delinsGACTT").id());
}

#[test]
fn members_justified_over_one_run_render_as_the_allele_they_make() {
    // a block of the members justified one at a time would be one block for
    // two different alleles
    let hdp = normalising_provider();
    let first = vrs(&hdp, "NC_K.1:g.[2del;4A>T]");
    let second = vrs(&hdp, "NC_K.1:g.[6del;4A>T]");
    assert_eq!(first.id(), vrs(&hdp, "NC_K.1:g.3_4delinsT").id());
    assert_eq!(second.id(), vrs(&hdp, "NC_K.1:g.4_5delinsT").id());
    assert_ne!(first.id(), second.id());
}

#[test]
fn an_allele_rendered_as_one_allele_carries_the_allele_as_written() {
    let VrsVariation::Allele(allele) = vrs(&normalising_provider(), "NM_X.1:c.[7C>T;8T>G]") else {
        panic!("not an Allele");
    };
    assert_eq!(
        allele.expressions,
        vec![VrsExpression {
            syntax: "hgvs.c".into(),
            value: "NM_X.1:c.[7C>T;8T>G]".into(),
        }]
    );
}

#[test]
fn a_one_member_allele_renders_as_its_members_allele() {
    let hdp = normalising_provider();
    assert_eq!(
        vrs(&hdp, "NM_X.1:c.[7C>T]").id(),
        vrs(&hdp, "NM_X.1:c.7C>T").id()
    );
}

#[test]
fn a_cancelling_allele_renders_as_the_reference_over_the_bases_it_names() {
    let hdp = normalising_provider();
    assert_eq!(
        vrs(&hdp, "NM_X.1:c.[4_5insA;5del]").id(),
        vrs(&hdp, "NM_X.1:c.4_5=").id()
    );
}

#[test]
fn a_block_carries_its_members_in_the_normalised_order() {
    let hdp = normalising_provider();
    let (VrsVariation::CisPhasedBlock(forward), VrsVariation::CisPhasedBlock(reversed)) = (
        vrs(&hdp, "NM_X.1:c.[7C>T;13T>G]"),
        vrs(&hdp, "NM_X.1:c.[13T>G;7C>T]"),
    ) else {
        panic!("not blocks");
    };
    assert_eq!(forward.members, reversed.members);
    let mapper = VariantMapper::with_refget(&hdp, &hdp);
    assert_eq!(
        forward.members[0].id,
        mapper
            .to_vrs(&parse_hgvs_variant("NM_X.1:c.7C>T").unwrap())
            .unwrap()
            .id
    );
}

#[test]
fn a_coding_allele_and_its_genomic_projection_share_a_vrs_identifier() {
    // c.4 and c.6 share codon 2, which groups them in the c. text but not
    // in VRS, where both spellings are the same two Alleles
    let hdp = normalising_provider();
    let coding = vrs(&hdp, "NM_X.1:c.[4A>T;6A>T]");
    assert!(
        matches!(coding, VrsVariation::CisPhasedBlock(_)),
        "{coding:?}"
    );
    assert_eq!(coding.id(), vrs(&hdp, "NC_X.1:g.[19A>T;21A>T]").id());
}

#[test]
fn a_delins_is_one_allele_whatever_it_leaves_unchanged_inside() {
    // c.4_6 is AAA, g.19_21: the delins keeps c.5
    let VrsVariation::CisPhasedBlock(block) =
        vrs(&normalising_provider(), "NM_X.1:c.[4_6delinsTAT;13T>G]")
    else {
        panic!("not a block");
    };
    assert_eq!(block.members.len(), 2);
    let location = &block.members[0].location;
    assert_eq!((location.start, location.end), (Exact(18), Exact(21)));
}

#[test]
fn members_that_touch_on_the_transcript_but_not_the_genome_are_two_alleles() {
    // c.5 and c.6 of NM_S.1 are either side of the intron
    let rendered = vrs(&normalising_provider(), "NM_S.1:c.[5del;6del]");
    assert!(
        matches!(&rendered, VrsVariation::CisPhasedBlock(b) if b.members.len() == 2),
        "{rendered:?}"
    );
}

#[test]
fn a_minus_strand_duplication_beside_a_delins_renders_as_its_genome_reading() {
    let hdp = normalising_provider();
    assert_eq!(
        vrs(&hdp, "NR_XR.1:n.[27_28delinsCA;28dup]").id(),
        vrs(&hdp, "NR_XR.1:n.28dup").id()
    );
}

#[test]
fn a_block_is_only_for_an_allele_that_normalises_to_several_members() {
    let hdp = normalising_provider();
    let mapper = VariantMapper::with_refget(&hdp, &hdp);
    let err = mapper
        .to_vrs_cis_phased(&parse_hgvs_variant("NM_X.1:c.[7C>T;8T>G]").unwrap())
        .unwrap_err();
    assert!(
        matches!(&err, HgvsError::UnsupportedOperation(m) if m.contains("to_vrs_variation")),
        "{err}"
    );
}
