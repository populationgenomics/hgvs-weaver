//! Canonical alleles, SPDI and VRS for protein variants: the same
//! normalisation as nucleotide alleles, on the protein sequence.

mod support;

use hgvs_weaver::error::HgvsError;
use hgvs_weaver::mapper::VariantMapper;
use hgvs_weaver::parse_hgvs_variant;
use hgvs_weaver::vrs::{refget_accession, VrsBound};
use support::Provider;

/// M K L A A A Y R Q
const PROTEIN: &str = "MKLAAAYRQ";
const NP: &str = "NP_TEST.1";

fn provider() -> Provider {
    Provider::new().sequence(NP, PROTEIN)
}

fn spdi(mapper: &VariantMapper, hgvs: &str) -> String {
    let var = parse_hgvs_variant(&format!("{NP}:{hgvs}")).unwrap();
    mapper
        .to_spdi_unambiguous(&var)
        .unwrap_or_else(|e| format!("ERR:{e}"))
}

#[test]
fn a_deletion_anywhere_in_a_run_is_the_same_allele() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    for hgvs in ["p.Ala4del", "p.Ala5del", "p.Ala6del"] {
        assert_eq!(spdi(&mapper, hgvs), "NP_TEST.1:3:AAA:AA", "{hgvs}");
    }
    let var = parse_hgvs_variant("NP_TEST.1:p.Ala4del").unwrap();
    assert_eq!(
        mapper.canonical_allele(&var).unwrap().repeat_subunit,
        Some(1)
    );
}

#[test]
fn substitution_insertion_delins_dup_and_repeat_resolve_on_the_protein() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_eq!(spdi(&mapper, "p.Lys2Leu"), "NP_TEST.1:1:K:L");
    assert_eq!(spdi(&mapper, "p.K2L"), "NP_TEST.1:1:K:L");
    assert_eq!(spdi(&mapper, "p.Lys2_Leu3insGlyGly"), "NP_TEST.1:2::GG");
    assert_eq!(spdi(&mapper, "p.Leu3_Ala4delinsTrp"), "NP_TEST.1:2:LA:W");
    assert_eq!(spdi(&mapper, "p.Ala4dup"), "NP_TEST.1:3:AAA:AAAA");
    assert_eq!(spdi(&mapper, "p.Ala4[5]"), "NP_TEST.1:3:AAA:AAAAA");
    assert_eq!(spdi(&mapper, "p.Lys2="), "NP_TEST.1:1:K:K");
    // The plain to_spdi has no other form for a protein.
    let var = parse_hgvs_variant("NP_TEST.1:p.Ala4dup").unwrap();
    assert_eq!(mapper.to_spdi(&var, false).unwrap(), "NP_TEST.1:3:AAA:AAAA");
}

#[test]
fn a_stop_among_the_new_residues_ends_the_protein() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    // Nonsense, a stop inserted before the residue, and a delins ending in a
    // stop all leave the same protein, so they are one allele.
    assert_eq!(spdi(&mapper, "p.Lys2Ter"), "NP_TEST.1:1:KLAAAYRQ:");
    assert_eq!(spdi(&mapper, "p.Met1_Lys2insTer"), "NP_TEST.1:1:KLAAAYRQ:");
    assert_eq!(
        spdi(&mapper, "p.Lys2_Leu3delinsTer"),
        "NP_TEST.1:1:KLAAAYRQ:"
    );
    assert_eq!(
        spdi(&mapper, "p.Leu3_Ala4delinsTrpTer"),
        "NP_TEST.1:2:LAAAYRQ:W"
    );
}

#[test]
fn consequences_have_no_allele() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    for hgvs in [
        "p.Leu3fs",
        "p.Leu3ArgfsTer4",
        "p.Gln9Terext*5",
        "p.Met1?",
        "p.?",
        "p.0?",
    ] {
        let out = spdi(&mapper, hgvs);
        assert!(out.starts_with("ERR:"), "{hgvs} gave {out}");
    }
}

/// `NP_TEST.1` as a provider may serve it, with its stop.
const NP_WITH_STOP: &str = "NP_STOP.1";

fn provider_with_stop() -> Provider {
    provider().sequence(NP_WITH_STOP, &format!("{PROTEIN}*"))
}

/// The canonical allele of `hgvs` on `ac`, as SPDI, or the error.
fn spdi_on(mapper: &VariantMapper, ac: &str, hgvs: &str) -> Result<String, HgvsError> {
    mapper.to_spdi_unambiguous(&parse_hgvs_variant(&format!("{ac}:{hgvs}")).unwrap())
}

#[test]
fn a_lost_stop_has_no_allele_whether_the_protein_is_served_with_its_stop_or_not() {
    let hdp = provider_with_stop();
    let mapper = VariantMapper::new(&hdp);
    for ac in [NP, NP_WITH_STOP] {
        for hgvs in ["p.Ter10del", "p.Ter10Gly", "p.Gln9_Ter10del"] {
            let out = spdi_on(&mapper, ac, hgvs);
            assert!(
                matches!(&out, Err(HgvsError::UnsupportedOperation(m)) if m.contains("stop")),
                "{ac}:{hgvs} gave {out:?}"
            );
        }
    }
}

#[test]
fn a_change_at_the_stop_that_keeps_it_has_an_allele() {
    let hdp = provider_with_stop();
    let mapper = VariantMapper::new(&hdp);
    for ac in [NP, NP_WITH_STOP] {
        for hgvs in ["p.Gln9_Ter10insGly", "p.Gln9_Ter10delinsGlnGlyTer"] {
            assert_eq!(
                spdi_on(&mapper, ac, hgvs).unwrap(),
                format!("{ac}:9::G"),
                "{hgvs}"
            );
        }
    }
}

#[test]
fn an_allele_past_the_stop_is_refused() {
    let hdp = provider_with_stop();
    let mapper = VariantMapper::new(&hdp);
    let out = mapper.from_spdi(&format!("{NP_WITH_STOP}:10::G"));
    assert!(
        matches!(&out, Err(HgvsError::ValidationError(m)) if m.contains("stop")),
        "{out:?}"
    );
}

#[test]
fn an_allele_that_loses_the_stop_is_refused() {
    let hdp = provider_with_stop();
    let mapper = VariantMapper::new(&hdp);
    for spdi in ["9:*:", "9:*:G", "8:Q*:"] {
        let out = mapper.from_spdi(&format!("{NP_WITH_STOP}:{spdi}"));
        assert!(
            matches!(&out, Err(HgvsError::UnsupportedOperation(m)) if m.contains("stop")),
            "{spdi} gave {out:?}"
        );
    }
}

#[test]
fn an_allele_that_replaces_the_stop_with_residues_and_a_stop_is_an_insertion() {
    let hdp = provider_with_stop();
    let mapper = VariantMapper::new(&hdp);
    assert_eq!(
        mapper
            .from_spdi(&format!("{NP_WITH_STOP}:9:*:G*"))
            .unwrap()
            .to_string(),
        "NP_STOP.1:p.Gln9_Ter10insGly"
    );
}

#[test]
fn the_vrs_allele_is_on_the_protein() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let var = parse_hgvs_variant("NP_TEST.1:p.Lys2Leu").unwrap();
    let vrs = mapper.to_vrs(&var).unwrap();
    let reference = &vrs.location.sequence_reference;
    assert_eq!(reference.residue_alphabet, "aa");
    assert_eq!(reference.molecule_type, "protein");
    assert_eq!(reference.refget_accession, refget_accession(PROTEIN));
    assert_eq!(
        (vrs.location.start, vrs.location.end),
        (VrsBound::Exact(1), VrsBound::Exact(2))
    );
    assert!(vrs.id.starts_with("ga4gh:VA."), "{}", vrs.id);
    assert_eq!(vrs.expressions[0].syntax, "hgvs.p");
    assert_eq!(vrs.expressions[0].value, "NP_TEST.1:p.Lys2Leu");
    assert!(vrs.to_json().contains(r#""sequence":"L""#));
}

#[test]
fn validation_checks_the_named_and_stated_residues() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let check = |hgvs: &str| {
        let var = parse_hgvs_variant(&format!("{NP}:{hgvs}")).unwrap();
        mapper.validate(&var).unwrap()
    };
    assert!(check("p.Lys2Leu"));
    assert!(!check("p.Arg2Leu"));
    assert!(check("p.Lys2_Ala4del"));
    assert!(!check("p.Lys2_Tyr4del"));
    assert!(check("p.Gln9del"));
    assert!(!check("p.Gln10del"), "past the end of the protein");
    assert!(
        check("p.Leu3fs"),
        "a frameshift's named residue is still checked"
    );
    assert!(!check("p.Trp3fs"));
    assert!(check("p.Ala4_Ala5delinsGly"));
    assert!(
        !check("p.Ala4_Tyr5delinsGly"),
        "the end residue is checked too"
    );
}
