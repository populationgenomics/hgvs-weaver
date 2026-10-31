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
fn spdi_on(ac: &str, hgvs: &str) -> Result<String, HgvsError> {
    let hdp = provider_with_stop();
    let mapper = VariantMapper::new(&hdp);
    mapper.to_spdi_unambiguous(&parse_hgvs_variant(&format!("{ac}:{hgvs}")).unwrap())
}

/// The variant an SPDI reads back as, printed, or the error.
fn read_back(spdi: &str) -> Result<String, HgvsError> {
    let hdp = provider_with_stop();
    let mapper = VariantMapper::new(&hdp);
    mapper.from_spdi(spdi).map(|v| v.to_string())
}

/// One test per case: on `$ac`, `$hgvs` loses the stop, so it has no allele.
macro_rules! a_lost_stop_has_no_allele {
    ($($name:ident: $ac:expr, $hgvs:literal;)*) => {$(
        #[test]
        fn $name() {
            let out = spdi_on($ac, $hgvs);
            assert!(
                matches!(&out, Err(HgvsError::UnsupportedOperation(m)) if m.contains("loses the stop")),
                "{out:?}"
            );
        }
    )*};
}

a_lost_stop_has_no_allele! {
    a_deleted_stop_has_no_allele: NP, "p.Ter10del";
    a_deleted_stop_served_has_no_allele: NP_WITH_STOP, "p.Ter10del";
    a_stop_made_a_residue_has_no_allele: NP, "p.Ter10Gly";
    a_stop_made_a_residue_served_has_no_allele: NP_WITH_STOP, "p.Ter10Gly";
    a_stop_deleted_with_the_residue_before_it_has_no_allele: NP, "p.Gln9_Ter10del";
    a_stop_deleted_with_the_residue_before_it_served_has_no_allele: NP_WITH_STOP, "p.Gln9_Ter10del";
}

/// One test per case: on `$ac`, `$hgvs` has the allele `$ac:$spdi`.
macro_rules! the_allele_is {
    ($($name:ident: $ac:expr, $hgvs:literal => $spdi:literal;)*) => {$(
        #[test]
        fn $name() {
            assert_eq!(spdi_on($ac, $hgvs).unwrap(), format!("{}:{}", $ac, $spdi));
        }
    )*};
}

the_allele_is! {
    an_insertion_before_the_stop_is_an_allele: NP, "p.Gln9_Ter10insGly" => "9::G";
    an_insertion_before_the_stop_served_is_an_allele: NP_WITH_STOP, "p.Gln9_Ter10insGly" => "9::G";
    a_delins_that_keeps_the_stop_is_the_insertion: NP, "p.Gln9_Ter10delinsGlnGlyTer" => "9::G";
    a_delins_that_keeps_the_stop_served_is_the_insertion: NP_WITH_STOP, "p.Gln9_Ter10delinsGlnGlyTer" => "9::G";
    a_duplication_across_the_stop_is_no_change: NP, "p.Gln9_Ter10dup" => "9::";
    a_duplication_across_the_stop_served_is_no_change: NP_WITH_STOP, "p.Gln9_Ter10dup" => "9::";
    a_repeat_across_the_stop_is_no_change: NP, "p.Gln9_Ter10[2]" => "9::";
    a_duplicated_stop_is_no_change: NP, "p.Ter10dup" => "9::";
    the_stop_unchanged_is_no_change: NP, "p.Ter10=" => "9::";
}

/// One test per case: the SPDI `$spdi` reads back as `$hgvs`.
macro_rules! reads_back_as {
    ($($name:ident: $spdi:literal => $hgvs:literal;)*) => {$(
        #[test]
        fn $name() {
            assert_eq!(read_back($spdi).unwrap(), $hgvs);
        }
    )*};
}

reads_back_as! {
    an_insertion_after_the_last_residue_is_one_before_the_stop:
        "NP_TEST.1:9::G" => "NP_TEST.1:p.Gln9_Ter10insGly";
    an_insertion_after_the_last_residue_served_is_one_before_the_stop:
        "NP_STOP.1:9::G" => "NP_STOP.1:p.Gln9_Ter10insGly";
    a_stop_replaced_by_residues_and_a_stop_is_an_insertion_before_it:
        "NP_STOP.1:9:*:G*" => "NP_STOP.1:p.Gln9_Ter10insGly";
    a_stop_inserted_after_the_last_residue_is_the_stop_unchanged:
        "NP_TEST.1:9::*" => "NP_TEST.1:p.Ter10=";
    a_residue_after_the_stop_is_the_stop_unchanged:
        "NP_STOP.1:9:*:*G" => "NP_STOP.1:p.Ter10=";
    no_change_at_the_stop_is_the_stop_unchanged:
        "NP_TEST.1:9::" => "NP_TEST.1:p.Ter10=";
}

/// One test per case: the SPDI `$spdi` removes the stop and makes none.
macro_rules! an_allele_that_loses_the_stop_is_refused {
    ($($name:ident: $spdi:literal;)*) => {$(
        #[test]
        fn $name() {
            let out = read_back($spdi);
            assert!(
                matches!(&out, Err(HgvsError::UnsupportedOperation(m)) if m.contains("loses the stop")),
                "{out:?}"
            );
        }
    )*};
}

an_allele_that_loses_the_stop_is_refused! {
    an_allele_deleting_the_stop_is_refused: "NP_STOP.1:9:*:";
    an_allele_replacing_the_stop_with_a_residue_is_refused: "NP_STOP.1:9:*:G";
    an_allele_deleting_the_stop_and_the_residue_before_it_is_refused: "NP_STOP.1:8:Q*:";
}

#[test]
fn an_allele_past_the_stop_is_refused() {
    let out = read_back("NP_STOP.1:10::G");
    assert!(
        matches!(&out, Err(HgvsError::ValidationError(m)) if m.contains("past it")),
        "{out:?}"
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
