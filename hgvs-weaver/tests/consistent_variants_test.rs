//! A description on a gene symbol or an unversioned accession, as the variants
//! it can be on the versioned accessions the store holds; choices.md, "A name
//! that is not one sequence is read on each the store holds".

mod support;

use hgvs_weaver::data::{
    AccessionCatalogue, AccessionStatus, AccessionTag, IdentifierType, Strand,
};
use hgvs_weaver::equivalence::{EquivalenceLevel, VariantEquivalence};
use hgvs_weaver::error::HgvsError;
use hgvs_weaver::mapper::VariantMapper;
use hgvs_weaver::{parse_hgvs_variant, SequenceVariant};
use support::{single_exon_transcript, Provider};

const UTR5: &str = "GGGGG";
/// M K L A Y R *: c.4_6 is AAA.
const CDS: &str = "ATGAAACTGGCCTATCGCTAA";
/// M E L A Y R *: c.4_6 is GAA.
const CDS_GLU2: &str = "ATGGAACTGGCCTATCGCTAA";
const UTR3: &str = "CCGTATAAGTAAGG";

/// Gene `CG` has three transcripts, all one exon at genomic index 10 of
/// `NC_C.1`: `NM_C.1` and `NM_D.1` read `CDS`, and `NM_C.2`, the newer version
/// of `NM_C.1`, reads `CDS_GLU2`. Their proteins are `NP_C.1`, `NP_D.1` and
/// `NP_C.2`. Gene `NG` is a symbol with no transcript.
fn provider() -> Provider {
    let record = |cds: &str| format!("{UTR5}{cds}{UTR3}");
    let len = record(CDS).len() as i32;
    let cds = (UTR5.len() as i32, (UTR5.len() + CDS.len()) as i32 - 1);
    let on_gene = |ac: &str| {
        let mut t = single_exon_transcript(ac, "NC_C.1", 10, Strand::Plus, cds.0, cds.1, len);
        t.gene = "CG".to_string();
        t
    };
    Provider::new()
        .sequence("NC_C.1", &format!("TTTTTTTTTT{}CCCCCCCCCC", record(CDS)))
        .sequence("NM_C.1", &record(CDS))
        .sequence("NM_C.2", &record(CDS_GLU2))
        .sequence("NM_D.1", &record(CDS))
        .sequence("NP_C.1", "MKLAYR")
        .sequence("NP_C.2", "MELAYR")
        .sequence("NP_D.1", "MKLAYR")
        .transcript(on_gene("NM_C.1"))
        .transcript(on_gene("NM_C.2"))
        .transcript(on_gene("NM_D.1"))
        .protein_for("NM_C.1", "NP_C.1")
        .protein_for("NM_C.2", "NP_C.2")
        .protein_for("NM_D.1", "NP_D.1")
        .identifier_type("CG", IdentifierType::GeneSymbol)
        .identifier_type("NG", IdentifierType::GeneSymbol)
}

fn parse(s: &str) -> SequenceVariant {
    parse_hgvs_variant(s).unwrap_or_else(|e| panic!("{s}: {e}"))
}

/// The variants `s` can be, as text, read with `catalogue`.
fn consistent(
    hdp: &Provider,
    s: &str,
    catalogue: Option<&dyn AccessionCatalogue>,
) -> Result<Vec<String>, HgvsError> {
    let mapper = VariantMapper::new(hdp);
    Ok(mapper
        .consistent_variants(&parse(s), catalogue)?
        .into_iter()
        .map(|c| c.variant.to_string())
        .collect())
}

#[test]
fn a_gene_symbol_gives_the_description_on_each_transcript_that_holds_it() {
    let hdp = provider();
    assert_eq!(
        consistent(&hdp, "CG:c.4A>G", None).unwrap(),
        ["NM_C.1:c.4A>G", "NM_D.1:c.4A>G"]
    );
}

#[test]
fn a_protein_description_on_a_gene_symbol_gives_each_protein_that_holds_it() {
    let hdp = provider();
    assert_eq!(
        consistent(&hdp, "CG:p.Lys2Glu", None).unwrap(),
        ["NP_C.1:p.Lys2Glu", "NP_D.1:p.Lys2Glu"]
    );
}

#[test]
fn an_r_description_on_a_gene_symbol_gives_each_transcript_that_holds_it() {
    let hdp = provider();
    assert_eq!(
        consistent(&hdp, "CG:r.4a>g", None).unwrap(),
        ["NM_C.1:r.4a>g", "NM_D.1:r.4a>g"]
    );
}

#[test]
fn an_unversioned_accession_gives_each_version_that_holds_it_newest_first() {
    let hdp = provider();
    assert_eq!(
        consistent(&hdp, "NM_C:c.6A>G", Some(&hdp)).unwrap(),
        ["NM_C.2:c.6A>G", "NM_C.1:c.6A>G"]
    );
}

#[test]
fn an_unversioned_accession_gives_only_the_versions_that_hold_it() {
    let hdp = provider();
    assert_eq!(
        consistent(&hdp, "NM_C:c.4A>G", Some(&hdp)).unwrap(),
        ["NM_C.1:c.4A>G"]
    );
}

#[test]
fn a_versioned_description_that_holds_gives_itself() {
    let hdp = provider();
    assert_eq!(
        consistent(&hdp, "NM_C.1:c.4A>G", Some(&hdp)).unwrap(),
        ["NM_C.1:c.4A>G"]
    );
}

#[test]
fn a_versioned_description_that_does_not_hold_gives_nothing() {
    let hdp = provider();
    assert!(consistent(&hdp, "NM_C.2:c.4A>G", Some(&hdp))
        .unwrap()
        .is_empty());
}

#[test]
fn a_description_no_version_holds_gives_nothing() {
    let hdp = provider();
    assert!(consistent(&hdp, "NM_C:c.4T>G", Some(&hdp))
        .unwrap()
        .is_empty());
}

#[test]
fn each_variant_returned_gives_itself_alone() {
    let hdp = provider();
    assert_eq!(
        consistent(&hdp, "NM_C.1:c.6A>G", Some(&hdp)).unwrap(),
        ["NM_C.1:c.6A>G"]
    );
    assert_eq!(
        consistent(&hdp, "NM_D.1:c.6A>G", Some(&hdp)).unwrap(),
        ["NM_D.1:c.6A>G"]
    );
    assert_eq!(
        consistent(&hdp, "CG:c.6A>G", Some(&hdp)).unwrap(),
        ["NM_C.1:c.6A>G", "NM_C.2:c.6A>G", "NM_D.1:c.6A>G"]
    );
}

#[test]
fn a_description_stating_bases_is_marked_as_stating_them() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let found = mapper
        .consistent_variants(&parse("NM_C.1:c.6A>G"), None)
        .unwrap();
    assert!(found[0].states_bases);
}

#[test]
fn a_description_stating_no_bases_is_marked_as_stating_none() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let found = mapper
        .consistent_variants(&parse("NM_C:c.6del"), Some(&hdp))
        .unwrap();
    assert_eq!(
        found.len(),
        2,
        "every version holds a deletion of no stated base"
    );
    assert!(!found[0].states_bases);
    assert!(!found[1].states_bases);
}

#[test]
fn each_variant_carries_what_the_catalogue_records_about_its_accession() {
    let hdp = provider()
        .accession_info(
            "NM_C.2",
            AccessionStatus::Current,
            &[AccessionTag::ManeSelect],
        )
        .accession_info("NM_C.1", AccessionStatus::Superseded, &[]);
    let mapper = VariantMapper::new(&hdp);
    let found = mapper
        .consistent_variants(&parse("NM_C:c.6A>G"), Some(&hdp))
        .unwrap();
    let info = |i: usize| found[i].info.clone().expect("the catalogue records it");
    assert_eq!(found[0].variant.to_string(), "NM_C.2:c.6A>G");
    assert_eq!(info(0).status, AccessionStatus::Current);
    assert_eq!(info(0).tags, [AccessionTag::ManeSelect]);
    assert_eq!(found[1].variant.to_string(), "NM_C.1:c.6A>G");
    assert_eq!(info(1).status, AccessionStatus::Superseded);
    assert!(info(1).tags.is_empty());
}

#[test]
fn without_a_catalogue_a_gene_symbol_gives_its_variants_with_no_record() {
    let hdp = provider().accession_info(
        "NM_C.1",
        AccessionStatus::Current,
        &[AccessionTag::ManeSelect],
    );
    let mapper = VariantMapper::new(&hdp);
    let found = mapper
        .consistent_variants(&parse("CG:c.4A>G"), None)
        .unwrap();
    assert_eq!(found[0].variant.to_string(), "NM_C.1:c.4A>G");
    assert_eq!(found[0].info, None);
}

#[test]
fn without_a_catalogue_an_unversioned_accession_is_an_error() {
    let hdp = provider();
    let err = consistent(&hdp, "NM_C:c.4A>G", None).unwrap_err();
    assert!(
        matches!(&err, HgvsError::DataProviderError(m) if m.contains("NM_C") && m.contains("catalogue")),
        "{err}"
    );
}

#[test]
fn an_unversioned_accession_with_no_version_held_is_an_error() {
    let hdp = provider();
    let err = consistent(&hdp, "NM_Z:c.4A>G", Some(&hdp)).unwrap_err();
    assert!(
        matches!(&err, HgvsError::DataProviderError(m) if m.contains("NM_Z")),
        "{err}"
    );
}

#[test]
fn a_gene_symbol_with_no_accession_of_the_descriptions_kind_is_an_error() {
    let hdp = provider();
    let err = consistent(&hdp, "NG:c.4A>G", None).unwrap_err();
    assert!(
        matches!(&err, HgvsError::DataProviderError(m) if m.contains("NG")),
        "{err}"
    );
}

#[test]
fn an_r_description_on_a_gene_symbol_is_its_own_spelling() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let r = parse("CG:r.4a>g");
    assert_eq!(
        VariantEquivalence::new(&mapper, &hdp)
            .equivalent_level(&r, &r)
            .unwrap(),
        EquivalenceLevel::Identity
    );
}

#[test]
fn an_r_description_on_a_gene_symbol_is_its_c_spelling_on_a_transcript_of_the_gene() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert!(VariantEquivalence::new(&mapper, &hdp)
        .equivalent_level(&parse("CG:r.4a>g"), &parse("NM_C.1:c.4A>G"))
        .unwrap()
        .is_equivalent());
}
