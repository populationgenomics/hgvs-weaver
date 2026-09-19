//! Protein alleles computed from a coding variant: the residues from the
//! first change to the end of the protein become what the edited transcript
//! encodes. Covers what a p. description cannot name as a sequence.

use hgvs_weaver::coords::{GenomicPos, TranscriptPos};
use hgvs_weaver::data::{
    DataProvider, ExonData, IdentifierKind, IdentifierType, Strand, TranscriptData,
};
use hgvs_weaver::error::HgvsError;
use hgvs_weaver::mapper::VariantMapper;
use hgvs_weaver::utils::translate;
use hgvs_weaver::vrs::VrsBound;
use hgvs_weaver::{parse_hgvs_variant, CVariant, SequenceVariant};

const UTR5: &str = "GGGGG";
/// M K L A Y R *
const CDS: &str = "ATGAAACTGGCCTATCGCTAA";
/// Read through the stop: P Y K *
const UTR3: &str = "CCGTATAAGTAAGG";
const PROTEIN: &str = "MKLAYR";

fn transcript() -> String {
    format!("{UTR5}{CDS}{UTR3}")
}

/// Serves the transcript and a protein; `protein` may disagree with the CDS.
struct Provider {
    protein: &'static str,
}

impl DataProvider for Provider {
    fn get_transcript(&self, ac: &str, _: Option<&str>) -> Result<TranscriptData, HgvsError> {
        if ac != "NM_X.1" {
            return Err(HgvsError::DataProviderError(format!("no transcript {ac}")));
        }
        let len = transcript().len() as i32;
        Ok(TranscriptData {
            ac: ac.to_string(),
            gene: "X".to_string(),
            cds_start_index: Some(TranscriptPos(UTR5.len() as i32)),
            cds_end_index: Some(TranscriptPos((UTR5.len() + CDS.len()) as i32 - 1)),
            strand: Strand::Plus,
            reference_accession: "NC_X.1".to_string(),
            exons: vec![ExonData {
                transcript_start: TranscriptPos(0),
                transcript_end: TranscriptPos(len),
                reference_start: GenomicPos(100),
                reference_end: GenomicPos(100 + len - 1),
                alt_strand: Strand::Plus,
                cigar: format!("{len}M"),
            }],
        })
    }

    fn get_seq(
        &self,
        ac: &str,
        start: i32,
        end: Option<i32>,
        _kind: IdentifierType,
    ) -> Result<String, HgvsError> {
        let seq = match ac {
            "NM_X.1" => transcript(),
            "NP_X.1" => self.protein.to_string(),
            _ => return Err(HgvsError::DataProviderError(format!("no sequence {ac}"))),
        };
        let start = (start.max(0) as usize).min(seq.len());
        let end = end.map_or(seq.len(), |e| (e.max(0) as usize).min(seq.len()));
        Ok(seq[start..end.max(start)].to_string())
    }

    fn get_symbol_accessions(
        &self,
        symbol: &str,
        _: IdentifierKind,
        target: IdentifierKind,
    ) -> Result<Vec<(IdentifierType, String)>, HgvsError> {
        Ok(match (symbol, target) {
            ("NM_X.1", IdentifierKind::Protein) => {
                vec![(IdentifierType::ProteinAccession, "NP_X.1".to_string())]
            }
            _ => vec![],
        })
    }

    fn get_identifier_type(&self, id: &str) -> Result<IdentifierType, HgvsError> {
        Ok(match &id[..3] {
            "NP_" => IdentifierType::ProteinAccession,
            "NC_" => IdentifierType::GenomicAccession,
            _ => IdentifierType::TranscriptAccession,
        })
    }
}

fn coding(s: &str) -> CVariant {
    match parse_hgvs_variant(s).unwrap() {
        SequenceVariant::Coding(c) => c,
        other => panic!("{other} is not c."),
    }
}

/// The protein the edited transcript encodes, to its first stop: the oracle.
fn translated(edited_cds_onwards: &str) -> String {
    let aa = translate(edited_cds_onwards);
    aa.split('*').next().unwrap_or("").to_string()
}

/// `allele` applied to the protein.
fn applied(allele: &hgvs_weaver::allele::CanonicalAllele) -> String {
    format!(
        "{}{}{}",
        &PROTEIN[..allele.start],
        allele.alternate,
        &PROTEIN[allele.end..]
    )
}

#[test]
fn definite_changes_give_the_same_allele_by_either_route() {
    let hdp = Provider { protein: PROTEIN };
    let mapper = VariantMapper::new(&hdp);
    for c in [
        "NM_X.1:c.4A>C",          // Lys2Gln
        "NM_X.1:c.4_6del",        // Lys2del
        "NM_X.1:c.6_7insGGG",     // Lys2_Leu3insGly
        "NM_X.1:c.7_9dup",        // Leu3dup
        "NM_X.1:c.7_12delinsTGG", // Leu3_Ala4delinsTrp
        "NM_X.1:c.6A>G",          // Lys2= (AAA -> AAG)
    ] {
        let vc = coding(c);
        let from_coding = mapper.protein_allele(&vc, None).unwrap();
        let vp = mapper.c_to_p(&vc, None).unwrap();
        let from_p = mapper
            .canonical_allele(&SequenceVariant::Protein(vp.clone()))
            .unwrap();
        assert_eq!(from_coding, from_p, "{c} ({vp})");
    }
    assert_eq!(
        mapper
            .protein_allele(&coding("NM_X.1:c.4A>C"), None)
            .unwrap()
            .spdi(),
        "NP_X.1:1:K:Q"
    );
}

#[test]
fn consequences_without_a_p_sequence_still_have_an_allele() {
    let hdp = Provider { protein: PROTEIN };
    let mapper = VariantMapper::new(&hdp);
    let allele = |c: &str| mapper.protein_allele(&coding(c), None).unwrap();

    // Frameshift: deleting the A of c.5 (AAA -> AA CTG ...) shifts the frame.
    let fs = allele("NM_X.1:c.5del");
    let edited = format!("{}{}", &CDS[..4], &CDS[5..]) + UTR3;
    assert_eq!(applied(&fs), translated(&edited));
    assert_eq!(fs.start, 1, "first changed residue is Lys2");
    assert_eq!(fs.end, PROTEIN.len(), "to the end of the protein");
    assert!(mapper
        .c_to_p(&coding("NM_X.1:c.5del"), None)
        .unwrap()
        .to_string()
        .contains("fs"));

    // Nonsense: AAA -> TAA truncates after Met1.
    let stop = allele("NM_X.1:c.4A>T");
    assert_eq!(stop.spdi(), "NP_X.1:1:KLAYR:");

    // Stop loss: TAA -> CAA reads Gln then P Y K into the UTR.
    let ext = allele("NM_X.1:c.19T>C");
    assert_eq!(applied(&ext), "MKLAYRQPYK");
    assert_eq!(ext.spdi(), "NP_X.1:6::QPYK");

    // The whole CDS deleted: no protein.
    let none = allele("NM_X.1:c.-5_*14del");
    assert_eq!(none.spdi(), "NP_X.1:0:MKLAYR:");

    // A stop inside the inserted bases.
    let ins_stop = allele("NM_X.1:c.6_7insTAAGGG");
    assert_eq!(ins_stop.spdi(), "NP_X.1:2:LAYR:");
}

#[test]
fn statements_have_no_allele_and_a_wrong_protein_is_an_error() {
    let hdp = Provider { protein: PROTEIN };
    let mapper = VariantMapper::new(&hdp);
    let err = mapper
        .protein_allele(&coding("NM_X.1:c.-3G>A"), None)
        .unwrap_err();
    assert!(matches!(err, HgvsError::UnsupportedOperation(_)), "{err}");
    let err = mapper
        .protein_allele(&coding("NM_X.1:c.-3_2del"), None)
        .unwrap_err();
    assert!(matches!(err, HgvsError::UnsupportedOperation(_)), "{err}");

    // The annotation pairs the transcript with a protein its CDS does not encode.
    let wrong = Provider { protein: "MKLAYQ" };
    let mapper = VariantMapper::new(&wrong);
    let err = mapper
        .protein_allele(&coding("NM_X.1:c.4A>C"), None)
        .unwrap_err();
    assert!(
        matches!(&err, HgvsError::ValidationError(m) if m.contains("residue 6")),
        "{err}"
    );
}

#[test]
fn the_vrs_allele_sits_on_the_protein_with_the_p_description_as_expression() {
    let hdp = Provider { protein: PROTEIN };
    let mapper = VariantMapper::new(&hdp);
    let vrs = mapper.protein_vrs(&coding("NM_X.1:c.5del"), None).unwrap();
    assert_eq!(vrs.location.sequence_reference.residue_alphabet, "aa");
    assert_eq!(vrs.location.sequence_reference.molecule_type, "protein");
    assert_eq!(
        vrs.location.sequence_reference.refget_accession,
        hgvs_weaver::vrs::refget_accession(PROTEIN)
    );
    assert_eq!(vrs.location.start, VrsBound::Exact(1));
    assert_eq!(vrs.expressions[0].syntax, "hgvs.p");
    assert!(
        vrs.expressions[0].value.contains("fsTer"),
        "{}",
        vrs.expressions[0].value
    );
    assert!(vrs.id.starts_with("ga4gh:VA."));
    // The same change through a p. variant that can name it gives the same id.
    let vc = coding("NM_X.1:c.4A>C");
    let vp = SequenceVariant::Protein(mapper.c_to_p(&vc, None).unwrap());
    assert_eq!(
        mapper.protein_vrs(&vc, None).unwrap().id,
        mapper.to_vrs(&vp).unwrap().id
    );
}
