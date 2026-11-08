//! A description on a gene symbol or an unversioned accession, as the variants
//! it can be on the versioned accessions the store holds; choices.md, "A name
//! that is not one sequence is read on each the store holds".

mod support;

use hgvs_weaver::data::{
    AccessionCatalogue, AccessionInfo, AccessionStatus, AccessionTag, IdentifierType, Strand,
};
use hgvs_weaver::equivalence::{EquivalenceLevel, VariantEquivalence};
use hgvs_weaver::error::HgvsError;
use hgvs_weaver::mapper::VariantMapper;
use hgvs_weaver::{parse_hgvs_variant, SequenceVariant};
use support::{exon, single_exon_transcript, transcript, Provider};

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
    let mapper = VariantMapper::new(&hdp);
    let again = |c: &hgvs_weaver::consistent::ConsistentVariant| {
        mapper
            .consistent_variants(&c.variant, Some(&hdp))
            .unwrap()
            .into_iter()
            .map(|c| c.variant)
            .collect::<Vec<_>>()
    };
    let found = mapper
        .consistent_variants(&parse("CG:c.6A>G"), Some(&hdp))
        .unwrap();
    assert_eq!(
        found.len(),
        3,
        "NM_C.1, NM_C.2 and NM_D.1 all hold A at c.6"
    );
    assert_eq!(again(&found[0]), [found[0].variant.clone()]);
    assert_eq!(again(&found[1]), [found[1].variant.clone()]);
    assert_eq!(again(&found[2]), [found[2].variant.clone()]);
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
        "every version holds a deletion of no written base"
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

#[test]
fn an_allele_on_a_gene_symbol_gives_the_allele_on_each_transcript_that_holds_it() {
    let hdp = provider();
    assert_eq!(
        consistent(&hdp, "CG:c.[4A>G;6A>G]", None).unwrap(),
        ["NM_C.1:c.[4A>G;6A>G]", "NM_D.1:c.[4A>G;6A>G]"]
    );
}

#[test]
fn a_record_the_catalogue_lists_as_its_own_version_gives_itself() {
    // QX1 is held under its bare name, as a UniProt accession is.
    let hdp = provider().sequence("QX1", "MKLAYR");
    assert_eq!(
        consistent(&hdp, "QX1:p.Lys2Glu", Some(&hdp)).unwrap(),
        ["QX1:p.Lys2Glu"]
    );
}

/// A catalogue that lists the accessions it was given as NM_C's versions.
struct Listing(&'static [&'static str]);

impl AccessionCatalogue for Listing {
    fn get_accession_versions(&self, _ac: &str) -> Result<Vec<String>, HgvsError> {
        Ok(self.0.iter().map(|ac| ac.to_string()).collect())
    }
    fn get_accession_info(&self, _ac: &str) -> Result<Option<AccessionInfo>, HgvsError> {
        Ok(None)
    }
}

/// The error `NM_C:c.6A>G` gives with a catalogue listing `versions`.
fn listing_error(versions: &'static [&'static str]) -> HgvsError {
    let hdp = provider();
    consistent(&hdp, "NM_C:c.6A>G", Some(&Listing(versions))).unwrap_err()
}

#[test]
fn a_catalogue_listing_another_accession_as_a_version_is_an_error() {
    let err = listing_error(&["NM_D.1"]);
    assert!(
        matches!(&err, HgvsError::DataProviderError(m) if m.contains("NM_D.1")),
        "{err}"
    );
}

#[test]
fn a_catalogue_listing_a_version_with_a_leading_zero_is_an_error() {
    let err = listing_error(&["NM_C.01"]);
    assert!(
        matches!(&err, HgvsError::DataProviderError(m) if m.contains("NM_C.01")),
        "{err}"
    );
}

#[test]
fn a_catalogue_listing_a_version_twice_is_an_error() {
    let err = listing_error(&["NM_C.2", "NM_C.2"]);
    assert!(
        matches!(&err, HgvsError::DataProviderError(m) if m.contains("twice")),
        "{err}"
    );
}

#[test]
fn a_catalogue_listing_the_bare_name_beside_its_versions_is_an_error() {
    let err = listing_error(&["NM_C.2", "NM_C"]);
    assert!(
        matches!(&err, HgvsError::DataProviderError(m) if m.contains("lists NM_C as")),
        "{err}"
    );
}

#[test]
fn a_genomic_description_on_a_gene_symbol_with_no_genomic_accession_is_an_error() {
    let hdp = provider();
    let err = consistent(&hdp, "CG:g.16A>G", None).unwrap_err();
    assert!(
        matches!(&err, HgvsError::DataProviderError(m) if m.contains("CG")),
        "{err}"
    );
}

/// Gene `MX` has the coding `NM_M.1` and the non-coding `NR_M.1`, one record
/// on `NC_M.1`.
fn mixed_provider() -> Provider {
    let record = format!("{UTR5}{CDS}{UTR3}");
    let len = record.len() as i32;
    let cds = (UTR5.len() as i32, (UTR5.len() + CDS.len()) as i32 - 1);
    let mut coding =
        single_exon_transcript("NM_M.1", "NC_M.1", 10, Strand::Plus, cds.0, cds.1, len);
    coding.gene = "MX".to_string();
    let mut non_coding = transcript(
        "NR_M.1",
        "NC_M.1",
        Strand::Plus,
        None,
        vec![exon((0, len), (10, 10 + len - 1), Strand::Plus)],
    );
    non_coding.gene = "MX".to_string();
    Provider::new()
        .sequence("NC_M.1", &format!("TTTTTTTTTT{record}CCCCCCCCCC"))
        .sequence("NM_M.1", &record)
        .sequence("NR_M.1", &record)
        .transcript(coding)
        .transcript(non_coding)
        .identifier_type("MX", IdentifierType::GeneSymbol)
}

#[test]
fn a_coding_description_on_a_gene_symbol_is_read_on_its_coding_transcripts() {
    let hdp = mixed_provider();
    assert_eq!(
        consistent(&hdp, "MX:c.4A>G", None).unwrap(),
        ["NM_M.1:c.4A>G"]
    );
}

#[test]
fn a_non_coding_description_on_a_gene_symbol_is_read_on_every_transcript() {
    // n.9 is the A of c.4 on NM_M.1, and the same base of NR_M.1.
    let hdp = mixed_provider();
    assert_eq!(
        consistent(&hdp, "MX:n.9A>G", None).unwrap(),
        ["NM_M.1:n.9A>G", "NR_M.1:n.9A>G"]
    );
}

#[test]
fn an_r_description_on_a_gene_symbol_is_read_on_coding_and_non_coding_transcripts() {
    // r.9 is c.9 on NM_M.1, a G, and n.9 on NR_M.1, an A.
    let hdp = mixed_provider();
    assert_eq!(
        consistent(&hdp, "MX:r.9a>g", None).unwrap(),
        ["NR_M.1:r.9a>g"]
    );
    assert_eq!(
        consistent(&hdp, "MX:r.9g>a", None).unwrap(),
        ["NM_M.1:r.9g>a"]
    );
}

#[test]
fn an_r_description_counting_from_the_cds_on_a_gene_symbol_is_read_on_its_coding_transcripts() {
    // r.-2 is the G two bases before the CDS on NM_M.1; NR_M.1 has no CDS to count from.
    let hdp = mixed_provider();
    assert_eq!(
        consistent(&hdp, "MX:r.-2g>a", None).unwrap(),
        ["NM_M.1:r.-2g>a"]
    );
}

/// Gene `LG` has `NM_S.1`, the `CG` transcript record, and `NM_L.1`, the same
/// with ten more bases of 3'UTR, so `c.*16` lies past the end of `NM_S.1`
/// only; and `NM_J.1`, the `NM_S.1` record split after c.5 by an intron on
/// `NC_J.1`, so c.4_6 crosses its junction and lies in one exon of the
/// others.
fn disagreeing_provider() -> Provider {
    let short = format!("{UTR5}{CDS}{UTR3}");
    let long = format!("{short}AAAAAAAAAA");
    let cds = (UTR5.len() as i32, (UTR5.len() + CDS.len()) as i32 - 1);
    let of_lg = |mut t: hgvs_weaver::data::TranscriptData| {
        t.gene = "LG".to_string();
        t
    };
    let split = 10;
    let intron = "GTAAGTCCAG";
    Provider::new()
        .sequence("NC_S.1", &format!("TTTTTTTTTT{long}CCCCCCCCCC"))
        .sequence(
            "NC_J.1",
            &format!(
                "TTTTTTTTTT{}{intron}{}CCCCCCCCCC",
                &short[..split],
                &short[split..]
            ),
        )
        .sequence("NM_S.1", &short)
        .sequence("NM_L.1", &long)
        .sequence("NM_J.1", &short)
        .transcript(of_lg(single_exon_transcript(
            "NM_S.1",
            "NC_S.1",
            10,
            Strand::Plus,
            cds.0,
            cds.1,
            short.len() as i32,
        )))
        .transcript(of_lg(single_exon_transcript(
            "NM_L.1",
            "NC_S.1",
            10,
            Strand::Plus,
            cds.0,
            cds.1,
            long.len() as i32,
        )))
        .transcript(of_lg(transcript(
            "NM_J.1",
            "NC_J.1",
            Strand::Plus,
            Some(cds),
            vec![
                exon((0, split as i32), (10, 10 + split as i32 - 1), Strand::Plus),
                exon(
                    (split as i32, short.len() as i32),
                    (
                        10 + (split + intron.len()) as i32,
                        10 + (short.len() + intron.len()) as i32 - 1,
                    ),
                    Strand::Plus,
                ),
            ],
        )))
        .identifier_type("LG", IdentifierType::GeneSymbol)
}

#[test]
fn an_r_description_on_a_gene_symbol_is_compared_on_the_transcripts_that_hold_it() {
    let hdp = disagreeing_provider();
    let mapper = VariantMapper::new(&hdp);
    assert!(VariantEquivalence::new(&mapper, &hdp)
        .equivalent_level(&parse("LG:r.*16del"), &parse("NM_L.1:c.*16del"))
        .unwrap()
        .is_equivalent());
}

#[test]
fn an_r_description_on_a_gene_symbol_crossing_one_transcripts_junction_is_compared_on_the_others() {
    let hdp = disagreeing_provider();
    let mapper = VariantMapper::new(&hdp);
    assert!(VariantEquivalence::new(&mapper, &hdp)
        .equivalent_level(&parse("LG:r.4_6del"), &parse("NM_S.1:c.4_6del"))
        .unwrap()
        .is_equivalent());
}

#[test]
fn an_r_description_on_a_gene_symbol_crossing_a_junction_is_refused_where_no_other_reading_matches()
{
    let hdp = disagreeing_provider();
    let mapper = VariantMapper::new(&hdp);
    let err = VariantEquivalence::new(&mapper, &hdp)
        .equivalent_level(&parse("LG:r.4_6del"), &parse("NM_S.1:c.10del"))
        .unwrap_err();
    assert!(
        matches!(&err, HgvsError::UnsupportedOperation(m) if m.contains("NM_J.1") && m.contains("exon junction")),
        "{err}"
    );
}

#[test]
fn an_r_description_on_a_gene_symbol_is_its_own_spelling_on_a_transcript_whose_junction_it_crosses()
{
    let hdp = disagreeing_provider();
    let mapper = VariantMapper::new(&hdp);
    let eq = VariantEquivalence::new(&mapper, &hdp);
    let (a, b) = (parse("LG:r.4_6del"), parse("NM_J.1:r.4_6del"));
    assert_eq!(
        eq.equivalent_level(&a, &b).unwrap(),
        EquivalenceLevel::Identity
    );
    assert_eq!(
        eq.equivalent_level(&b, &a).unwrap(),
        EquivalenceLevel::Identity
    );
}

#[test]
fn a_non_coding_description_on_a_gene_symbol_is_its_coding_spelling() {
    // n.9 of NM_C.1 is its c.4.
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert!(VariantEquivalence::new(&mapper, &hdp)
        .equivalent_level(&parse("CG:n.9A>G"), &parse("CG:c.4A>G"))
        .unwrap()
        .is_equivalent());
}

#[test]
fn a_gene_symbol_is_as_equivalent_as_its_closest_reading() {
    // The NM_C.1 reading is Analogous to the NM_D.1 allele; the NM_D.1
    // reading is the same text.
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    assert_eq!(
        VariantEquivalence::new(&mapper, &hdp)
            .equivalent_level(&parse("CG:c.[4A>G;6A>G]"), &parse("NM_D.1:c.[4A>G;6A>G]"))
            .unwrap(),
        EquivalenceLevel::Identity
    );
}

#[test]
fn a_gene_symbol_refused_at_a_junction_names_the_same_reading_in_either_order() {
    let hdp = disagreeing_provider();
    let mapper = VariantMapper::new(&hdp);
    let eq = VariantEquivalence::new(&mapper, &hdp);
    let (a, b) = (parse("LG:r.4_6del"), parse("NM_S.1:c.10del"));
    assert_eq!(
        eq.equivalent_level(&a, &b).unwrap_err().to_string(),
        eq.equivalent_level(&b, &a).unwrap_err().to_string()
    );
}
