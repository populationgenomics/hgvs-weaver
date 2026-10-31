//! A description is read on a reference of its own kind: a `p.` description
//! on a protein, a `g.`, `m.`, `c.`, `n.` or `r.` one on a nucleotide
//! sequence. Where the provider says the accession is of the other kind, the
//! description is refused as malformed, not read: a transcript's bases read as
//! residues would validate whenever the base is A, C, G or T and the residue
//! named Ala, Cys, Gly or Thr. An accession the provider cannot classify is
//! read as written. Issue #123.

mod support;

use hgvs_weaver::data::{IdentifierType, Strand};
use hgvs_weaver::error::{HgvsError, ReferenceUse};
use hgvs_weaver::mapper::VariantMapper;
use hgvs_weaver::{parse_hgvs_variant, CVariant, GVariant, SequenceVariant};
use support::{exon, transcript, Provider};

/// A 5'UTR of five, then `ATG AAA CTG GCC TAT CGC TAA`: c.1 is the A at
/// index 5, which read as a residue is Ala6.
const RECORD: &str = "GGGGGATGAAACTGGCCTATCGCTAA";
const PROTEIN: &str = "MKLAYR";

/// `NM_K.1` is one exon on `NC_K.1`, the record itself; `NP_K.1` its protein.
/// `XQ_K.1`, which the provider cannot classify, serves the protein too.
fn provider() -> Provider {
    Provider::new()
        .sequence("NC_K.1", RECORD)
        .sequence("NM_K.1", RECORD)
        .sequence("NP_K.1", PROTEIN)
        .sequence("XQ_K.1", PROTEIN)
        .transcript(transcript(
            "NM_K.1",
            "NC_K.1",
            Strand::Plus,
            Some((5, 25)),
            vec![exon((0, 26), (0, 25), Strand::Plus)],
        ))
        .protein_for("NM_K.1", "NP_K.1")
        .identifier_type("XQ_K.1", IdentifierType::Unknown)
}

fn parse(s: &str) -> SequenceVariant {
    parse_hgvs_variant(s).unwrap()
}

/// `result` refuses `accession` as the other kind of sequence from the one a
/// description of coordinate type `coordinate_type` is written on.
fn assert_refused_description<T: std::fmt::Debug>(
    result: Result<T, HgvsError>,
    accession: &str,
    coordinate_type: &str,
) {
    match result {
        Err(HgvsError::WrongReferenceKind {
            accession: ac,
            read_as,
            ..
        }) => assert_eq!(
            (ac.as_str(), read_as),
            (
                accession,
                ReferenceUse::Description(coordinate_type.to_string())
            )
        ),
        other => panic!(
            "expected {accession} refused for a {coordinate_type}. description, got {other:?}"
        ),
    }
}

/// `result` refuses `accession`, given to an operation as a protein if
/// `as_protein`, or else as a nucleotide sequence, as the other kind.
fn assert_refused_argument<T: std::fmt::Debug>(
    result: Result<T, HgvsError>,
    accession: &str,
    as_protein: bool,
) {
    match result {
        Err(HgvsError::WrongReferenceKind {
            accession: ac,
            read_as: read_as @ (ReferenceUse::Protein(_) | ReferenceUse::Nucleotide(_)),
            ..
        }) => assert_eq!(
            (ac.as_str(), read_as.needs_protein()),
            (accession, as_protein)
        ),
        other => panic!("expected the argument {accession} refused, got {other:?}"),
    }
}

#[test]
fn validate_refuses_a_protein_description_on_a_transcript() {
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    assert_refused_description(mapper.validate(&parse("NM_K.1:p.Ala6Pro")), "NM_K.1", "p");
}

#[test]
fn validate_refuses_a_protein_description_on_a_genome() {
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    assert_refused_description(mapper.validate(&parse("NC_K.1:p.Gly1Pro")), "NC_K.1", "p");
}

#[test]
fn normalize_refuses_a_protein_description_on_a_transcript() {
    // Read as residues, base 6 is A, so this would be a mismatch, not L.
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    assert_refused_description(
        mapper.normalize_variant(parse("NM_K.1:p.Leu6Pro")),
        "NM_K.1",
        "p",
    );
}

#[test]
fn to_spdi_refuses_a_protein_description_on_a_transcript() {
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    assert_refused_description(
        mapper.to_spdi(&parse("NM_K.1:p.Ala6Pro"), false),
        "NM_K.1",
        "p",
    );
}

#[test]
fn validate_refuses_a_nucleotide_description_on_a_protein() {
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    for (description, coordinate_type) in [
        ("NP_K.1:c.2T>A", "c"),
        ("NP_K.1:n.2T>A", "n"),
        ("NP_K.1:r.2u>a", "r"),
        ("NP_K.1:g.4A>T", "g"),
        ("NP_K.1:m.4A>T", "m"),
    ] {
        assert_refused_description(
            mapper.validate(&parse(description)),
            "NP_K.1",
            coordinate_type,
        );
    }
}

#[test]
fn normalize_refuses_a_genomic_description_on_a_protein() {
    // Residue 4 of NP_K.1 is A, so read as a base the description would hold.
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    assert_refused_description(
        mapper.normalize_variant(parse("NP_K.1:g.4A>T")),
        "NP_K.1",
        "g",
    );
}

#[test]
fn a_protein_description_on_its_protein_validates() {
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    assert!(mapper.validate(&parse("NP_K.1:p.Lys2Arg")).unwrap());
}

#[test]
fn a_protein_description_on_an_accession_of_unknown_kind_is_read_as_written() {
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    assert!(mapper.validate(&parse("XQ_K.1:p.Lys2Arg")).unwrap());
}

fn coding(s: &str) -> CVariant {
    match parse(s) {
        SequenceVariant::Coding(c) => c,
        other => panic!("{other} is not c."),
    }
}

fn genomic(s: &str) -> GVariant {
    match parse(s) {
        SequenceVariant::Genomic(g) => g,
        other => panic!("{other} is not g."),
    }
}

#[test]
fn c_to_g_refuses_a_coding_description_on_a_protein() {
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    assert_refused_description(mapper.c_to_g(&coding("NP_K.1:c.2T>A"), None), "NP_K.1", "c");
}

#[test]
fn c_to_g_refuses_a_protein_as_the_genome_to_project_onto() {
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    assert_refused_argument(
        mapper.c_to_g(&coding("NM_K.1:c.2T>A"), Some("NP_K.1")),
        "NP_K.1",
        false,
    );
}

#[test]
fn n_to_g_refuses_a_non_coding_description_on_a_protein() {
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    let n = match parse("NP_K.1:n.2G>A") {
        SequenceVariant::NonCoding(n) => n,
        other => panic!("{other} is not n."),
    };
    assert_refused_description(mapper.n_to_g(&n, None), "NP_K.1", "n");
}

#[test]
fn r_to_g_refuses_an_rna_description_on_a_protein() {
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    let r = match parse("NP_K.1:r.2u>a") {
        SequenceVariant::Rna(r) => r,
        other => panic!("{other} is not r."),
    };
    assert_refused_description(mapper.r_to_g(&r, None), "NP_K.1", "r");
}

#[test]
fn c_to_p_refuses_a_coding_description_on_a_protein() {
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    assert_refused_description(mapper.c_to_p(&coding("NP_K.1:c.2T>A"), None), "NP_K.1", "c");
}

#[test]
fn to_vrs_refuses_a_coding_description_on_a_protein() {
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    assert_refused_description(mapper.to_vrs(&parse("NP_K.1:c.2T>A")), "NP_K.1", "c");
}

#[test]
fn g_to_c_refuses_a_genomic_description_on_a_protein() {
    // Residue 4 of NP_K.1 is A, so read as a base the description would hold.
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    assert_refused_description(
        mapper.g_to_c(&genomic("NP_K.1:g.4A>T"), "NM_K.1"),
        "NP_K.1",
        "g",
    );
}

#[test]
fn g_to_c_refuses_a_protein_as_the_transcript_to_project_onto() {
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    assert_refused_argument(
        mapper.g_to_c(&genomic("NC_K.1:g.4A>T"), "NP_K.1"),
        "NP_K.1",
        false,
    );
}

#[test]
fn g_to_c_all_refuses_a_genomic_description_on_a_protein() {
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    assert_refused_description(
        mapper.g_to_c_all(&genomic("NP_K.1:g.4A>T"), &provider),
        "NP_K.1",
        "g",
    );
}

#[test]
fn c_to_p_refuses_a_genome_as_the_protein_to_predict_on() {
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    assert_refused_argument(
        mapper.c_to_p(&coding("NM_K.1:c.4A>G"), Some("NC_K.1")),
        "NC_K.1",
        true,
    );
}

#[test]
fn c_to_p_refuses_a_transcript_as_the_protein_to_predict_on() {
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    assert_refused_argument(
        mapper.c_to_p(&coding("NM_K.1:c.4A>G"), Some("NM_K.1")),
        "NM_K.1",
        true,
    );
}

#[test]
fn protein_allele_refuses_a_genome_as_the_protein_to_predict_on() {
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    assert_refused_argument(
        mapper.protein_allele(&coding("NM_K.1:c.4A>G"), Some("NC_K.1")),
        "NC_K.1",
        true,
    );
}

#[test]
fn r_to_p_refuses_a_genome_as_the_protein_to_predict_on() {
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    let r = match parse("NM_K.1:r.4a>g") {
        SequenceVariant::Rna(r) => r,
        other => panic!("{other} is not r."),
    };
    assert_refused_argument(mapper.r_to_p(&r, Some("NC_K.1")), "NC_K.1", true);
}

#[test]
fn a_refusal_names_a_description_by_the_name_of_its_letter() {
    // HGVS letters are read by name: an m., an n., an r., but a c., a g.
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    let err = mapper.validate(&parse("NP_K.1:m.4A>T")).unwrap_err();
    assert!(err.to_string().contains("an m. description"), "{err}");
}

#[test]
fn equivalence_refuses_a_coding_allele_on_a_protein_before_comparing_its_spelling() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let allele = parse("NP_K.1:c.[4A>T;6A>G]");
    assert_refused_description(
        hgvs_weaver::equivalence::VariantEquivalence::new(&mapper, &hdp)
            .equivalent_level(&allele, &allele),
        "NP_K.1",
        "c",
    );
}
