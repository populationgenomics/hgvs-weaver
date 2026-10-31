//! A description is read on a reference of its own kind: a `p.` description
//! on a protein, a `g.`, `m.`, `c.`, `n.` or `r.` one on a nucleotide
//! sequence. Where the provider says the accession is of the other kind, the
//! description is refused as malformed, not read: a transcript's bases read as
//! residues would validate whenever the base is A, C, G or T and the residue
//! named Ala, Cys, Gly or Thr. An accession the provider cannot classify is
//! read as written. Issue #123.

mod support;

use hgvs_weaver::data::{IdentifierType, Strand};
use hgvs_weaver::error::HgvsError;
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

/// `result` refuses the description for the kind of reference it names.
fn assert_wrong_kind<T: std::fmt::Debug>(result: Result<T, HgvsError>, needs: &str) {
    match result {
        Err(HgvsError::ValidationError(msg)) => assert!(msg.contains(needs), "{msg}"),
        other => panic!("expected a ValidationError saying it {needs}, got {other:?}"),
    }
}

#[test]
fn validate_refuses_a_protein_description_on_a_transcript() {
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    assert_wrong_kind(
        mapper.validate(&parse("NM_K.1:p.Ala6Pro")),
        "needs a protein reference",
    );
}

#[test]
fn validate_refuses_a_protein_description_on_a_genome() {
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    assert_wrong_kind(
        mapper.validate(&parse("NC_K.1:p.Gly1Pro")),
        "needs a protein reference",
    );
}

#[test]
fn normalize_refuses_a_protein_description_on_a_transcript() {
    // Read as residues, base 6 is A, so this would be a mismatch, not L.
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    assert_wrong_kind(
        mapper.normalize_variant(parse("NM_K.1:p.Leu6Pro")),
        "needs a protein reference",
    );
}

#[test]
fn to_spdi_refuses_a_protein_description_on_a_transcript() {
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    assert_wrong_kind(
        mapper.to_spdi(&parse("NM_K.1:p.Ala6Pro"), false),
        "needs a protein reference",
    );
}

#[test]
fn validate_refuses_a_nucleotide_description_on_a_protein() {
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    for (description, says) in [
        (
            "NP_K.1:c.2T>A",
            "a c. description needs a nucleotide reference",
        ),
        (
            "NP_K.1:n.2T>A",
            "an n. description needs a nucleotide reference",
        ),
        (
            "NP_K.1:r.2u>a",
            "an r. description needs a nucleotide reference",
        ),
        (
            "NP_K.1:g.4A>T",
            "a g. description needs a nucleotide reference",
        ),
        (
            "NP_K.1:m.4A>T",
            "an m. description needs a nucleotide reference",
        ),
    ] {
        match mapper.validate(&parse(description)) {
            Err(HgvsError::ValidationError(msg)) => {
                assert!(msg.contains(says), "{description}: {msg}")
            }
            other => panic!("{description}: expected a ValidationError, got {other:?}"),
        }
    }
}

#[test]
fn normalize_refuses_a_genomic_description_on_a_protein() {
    // Residue 4 of NP_K.1 is A, so read as a base the description would hold.
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    assert_wrong_kind(
        mapper.normalize_variant(parse("NP_K.1:g.4A>T")),
        "needs a nucleotide reference",
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
    assert_wrong_kind(
        mapper.c_to_g(&coding("NP_K.1:c.2T>A"), None),
        "needs a nucleotide reference",
    );
}

#[test]
fn c_to_g_refuses_a_protein_as_the_genome_to_project_onto() {
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    assert_wrong_kind(
        mapper.c_to_g(&coding("NM_K.1:c.2T>A"), Some("NP_K.1")),
        "the genome a transcript variant is projected onto needs a nucleotide reference",
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
    assert_wrong_kind(mapper.n_to_g(&n, None), "needs a nucleotide reference");
}

#[test]
fn r_to_g_refuses_an_rna_description_on_a_protein() {
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    let r = match parse("NP_K.1:r.2u>a") {
        SequenceVariant::Rna(r) => r,
        other => panic!("{other} is not r."),
    };
    assert_wrong_kind(mapper.r_to_g(&r, None), "needs a nucleotide reference");
}

#[test]
fn c_to_p_refuses_a_coding_description_on_a_protein() {
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    assert_wrong_kind(
        mapper.c_to_p(&coding("NP_K.1:c.2T>A"), None),
        "needs a nucleotide reference",
    );
}

#[test]
fn to_vrs_refuses_a_coding_description_on_a_protein() {
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    assert_wrong_kind(
        mapper.to_vrs(&parse("NP_K.1:c.2T>A")),
        "a c. description needs a nucleotide reference",
    );
}

#[test]
fn g_to_c_refuses_a_genomic_description_on_a_protein() {
    // Residue 4 of NP_K.1 is A, so read as a base the description would hold.
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    assert_wrong_kind(
        mapper.g_to_c(&genomic("NP_K.1:g.4A>T"), "NM_K.1"),
        "needs a nucleotide reference",
    );
}

#[test]
fn g_to_c_refuses_a_protein_as_the_transcript_to_project_onto() {
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    assert_wrong_kind(
        mapper.g_to_c(&genomic("NC_K.1:g.4A>T"), "NP_K.1"),
        "the transcript g_to_c projects onto needs a nucleotide reference",
    );
}

#[test]
fn g_to_c_all_refuses_a_genomic_description_on_a_protein() {
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    assert_wrong_kind(
        mapper.g_to_c_all(&genomic("NP_K.1:g.4A>T"), &provider),
        "needs a nucleotide reference",
    );
}

#[test]
fn to_spdi_names_the_kind_of_description_it_refuses() {
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    assert_wrong_kind(
        mapper.to_spdi(&parse("NP_K.1:c.2T>A"), false),
        "a c. description needs a nucleotide reference",
    );
}

#[test]
fn g_to_c_names_the_genomic_description_it_refuses() {
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    assert_wrong_kind(
        mapper.g_to_c(&genomic("NP_K.1:g.4A>T"), "NM_K.1"),
        "a g. description needs a nucleotide reference",
    );
}

#[test]
fn c_to_p_refuses_a_genome_as_the_protein_to_predict_on() {
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    assert_wrong_kind(
        mapper.c_to_p(&coding("NM_K.1:c.4A>G"), Some("NC_K.1")),
        "the protein a prediction is written on needs a protein reference",
    );
}

#[test]
fn c_to_p_refuses_a_transcript_as_the_protein_to_predict_on() {
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    assert_wrong_kind(
        mapper.c_to_p(&coding("NM_K.1:c.4A>G"), Some("NM_K.1")),
        "the protein a prediction is written on needs a protein reference",
    );
}

#[test]
fn protein_allele_refuses_a_genome_as_the_protein_to_predict_on() {
    let provider = provider();
    let mapper = VariantMapper::new(&provider);
    assert_wrong_kind(
        mapper.protein_allele(&coding("NM_K.1:c.4A>G"), Some("NC_K.1")),
        "the protein a prediction is written on needs a protein reference",
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
    assert_wrong_kind(
        mapper.r_to_p(&r, Some("NC_K.1")),
        "the protein a prediction is written on needs a protein reference",
    );
}
