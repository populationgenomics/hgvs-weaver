//! Canonical alleles, SPDI and VRS for protein variants: the same
//! normalisation as nucleotide alleles, on the protein sequence.

use hgvs_weaver::data::{DataProvider, IdentifierKind, IdentifierType, TranscriptData};
use hgvs_weaver::error::HgvsError;
use hgvs_weaver::mapper::VariantMapper;
use hgvs_weaver::parse_hgvs_variant;
use hgvs_weaver::vrs::refget_accession;

/// M K L A A A Y R Q
const PROTEIN: &str = "MKLAAAYRQ";
const NP: &str = "NP_TEST.1";

struct Provider;

impl DataProvider for Provider {
    fn get_transcript(&self, ac: &str, _: Option<&str>) -> Result<TranscriptData, HgvsError> {
        Err(HgvsError::DataProviderError(format!("no transcript {ac}")))
    }

    fn get_seq(
        &self,
        ac: &str,
        start: i32,
        end: Option<i32>,
        _kind: IdentifierType,
    ) -> Result<String, HgvsError> {
        if ac != NP {
            return Err(HgvsError::DataProviderError(format!("no sequence {ac}")));
        }
        let start = (start.max(0) as usize).min(PROTEIN.len());
        let end = end.map_or(PROTEIN.len(), |e| (e.max(0) as usize).min(PROTEIN.len()));
        Ok(PROTEIN[start..end.max(start)].to_string())
    }

    fn get_symbol_accessions(
        &self,
        _: &str,
        _: IdentifierKind,
        _: IdentifierKind,
    ) -> Result<Vec<(IdentifierType, String)>, HgvsError> {
        Ok(vec![])
    }

    fn get_identifier_type(&self, id: &str) -> Result<IdentifierType, HgvsError> {
        Ok(if id.starts_with("NP_") {
            IdentifierType::ProteinAccession
        } else {
            IdentifierType::GenomicAccession
        })
    }
}

fn spdi(mapper: &VariantMapper, hgvs: &str) -> String {
    let var = parse_hgvs_variant(&format!("{NP}:{hgvs}")).unwrap();
    mapper
        .to_spdi_unambiguous(&var)
        .unwrap_or_else(|e| format!("ERR:{e}"))
}

#[test]
fn a_deletion_anywhere_in_a_run_is_the_same_allele() {
    let hdp = Provider;
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
    let hdp = Provider;
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
fn consequences_have_no_allele() {
    let hdp = Provider;
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

#[test]
fn the_vrs_allele_is_on_the_protein() {
    let hdp = Provider;
    let mapper = VariantMapper::new(&hdp);
    let var = parse_hgvs_variant("NP_TEST.1:p.Lys2Leu").unwrap();
    let vrs = mapper.to_vrs(&var).unwrap();
    let reference = &vrs.location.sequence_reference;
    assert_eq!(reference.residue_alphabet, "aa");
    assert_eq!(reference.molecule_type, "protein");
    assert_eq!(reference.refget_accession, refget_accession(PROTEIN));
    assert_eq!((vrs.location.start, vrs.location.end), (1, 2));
    assert!(vrs.id.starts_with("ga4gh:VA."), "{}", vrs.id);
    assert_eq!(vrs.expressions[0].syntax, "hgvs.p");
    assert_eq!(vrs.expressions[0].value, "NP_TEST.1:p.Lys2Leu");
    assert!(vrs.to_json().contains(r#""sequence":"L""#));
}

#[test]
fn validation_checks_the_named_and_stated_residues() {
    let hdp = Provider;
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
