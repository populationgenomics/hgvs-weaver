//! Alleles in cis, `c.[145C>T;147C>G]`: parsed, formatted, normalised and
//! validated member by member, rendered as VRS `CisPhasedBlock` objects and
//! read back.

mod support;

use hgvs_weaver::data::Strand;
use hgvs_weaver::equivalence::{EquivalenceLevel, VariantEquivalence};
use hgvs_weaver::error::HgvsError;
use hgvs_weaver::mapper::VariantMapper;
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

#[test]
fn a_cis_allele_is_a_cis_phased_block_of_its_members_alleles() {
    let hdp = provider();
    let mapper = VariantMapper::with_refget(&hdp, &hdp);
    let var = parse_hgvs_variant("NM_X.1:c.[7C>T;13T>G]").unwrap();
    let block = mapper.to_vrs_cis_phased(&var).unwrap();
    assert_eq!(block.type_, "CisPhasedBlock");
    assert!(block.id.starts_with("ga4gh:CPB."), "{}", block.id);
    assert_eq!(block.id, format!("ga4gh:CPB.{}", block.digest));

    // Each member is the Allele the member variant renders as on its own.
    let c = cis("NM_X.1:c.[7C>T;13T>G]");
    assert_eq!(block.members.len(), 2);
    for (allele, member) in block.members.iter().zip(&c.members) {
        assert_eq!(*allele, mapper.to_vrs(member).unwrap());
        assert_eq!(allele.expressions[0].value, member.to_string());
    }
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
    // The members array keeps the order written, only the digest sorts.
    assert_eq!(forward.members[0], reversed.members[1]);
    assert_ne!(forward.to_json(), reversed.to_json());
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
    let with_copy = parse_hgvs_variant("NC_X.1:g.[22C>T;20_30copy3]").unwrap();
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
        ("NM_X.1:c.[7C>T]", "NC_X.1:g.[22C>T]"),
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
    // The same members in the other order, spelled on the genome, or
    // written unnormalised, name the same molecule.
    assert_eq!(
        level("NM_X.1:c.[7C>T;13T>G]", "NM_X.1:c.[13T>G;7C>T]"),
        Analogous
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
        Analogous
    );
    assert_eq!(
        level("NP_X.1:p.[Lys2Leu;Ala4del]", "NP_X.1:p.[K2L;A4del]"),
        Identity
    );
    assert_eq!(
        level("NP_X.1:p.[Lys2Leu;Ala4del]", "NP_X.1:p.[Lys2Leu;Tyr5del]"),
        Different
    );
    assert_eq!(
        level("NM_X.1:c.[7C>T;13T>G]", "NP_X.1:p.Leu3Phe"),
        Different
    );
    // A member without a canonical allele compares Different, not an error.
    assert_eq!(
        level("NC_X.1:g.[22C>T;20_30copy3]", "NC_X.1:g.[22C>T;20_30copy3]"),
        Identity
    );
    assert_eq!(
        level("NC_X.1:g.[22C>T;20_30copy3]", "NC_X.1:g.[20_30copy3;22C>T]"),
        Different
    );
}
