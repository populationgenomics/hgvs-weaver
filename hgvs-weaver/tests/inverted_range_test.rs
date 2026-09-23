//! An HGVS range runs from its start to its end. A range written backwards,
//! `c.100_50del`, names nothing; papers produce them (a typo, an OCR slip) and
//! one used to panic inside `validate`. It is refused when parsed, and a
//! range built in code that runs backwards is an error where it is resolved,
//! never a panic.

mod support;

use hgvs_weaver::coords::{Anchor, HgvsGenomicPos, HgvsTranscriptPos};
use hgvs_weaver::data::Strand;
use hgvs_weaver::edits::NaEdit;
use hgvs_weaver::error::HgvsError;
use hgvs_weaver::mapper::VariantMapper;
use hgvs_weaver::structs::{
    BaseOffsetInterval, BaseOffsetPosition, CVariant, GVariant, LinearVariant, PosEdit,
    SimpleInterval, SimplePosition, TranscriptVariant,
};
use hgvs_weaver::{parse_hgvs_variant, SequenceVariant};
use support::{single_exon_transcript, Provider};

#[test]
fn a_range_written_backwards_is_refused_when_parsed() {
    for s in [
        // The two the issue found in papers, and its synthetic one.
        "NM_206933.4:c.8559_2A>G",
        "NM_025137.4:c.6331_6232insG",
        "NM_025137.4:c.100_50del",
        "NC_000001.11:g.1100_1050del",
        "NP_000001.1:p.Lys10_Leu5del",
        "NM_025137.4:c.88-1_87+1del",
        "NM_025137.4:c.*3_10del",
        "NM_025137.4:r.100_50del",
    ] {
        let err = parse_hgvs_variant(s).expect_err(s);
        assert!(
            matches!(&err, HgvsError::PestError(m) if m.contains("runs backwards")),
            "{s}: {err}"
        );
    }
}

#[test]
fn ordered_ranges_still_parse() {
    for s in [
        "NM_025137.4:c.50_100del",
        "NM_025137.4:c.50_50del",
        "NM_025137.4:c.-5_10del",
        "NM_025137.4:c.10_*3del",
        "NM_025137.4:c.*3_*10del",
        "NM_025137.4:c.87+1_88-1del",
        "NM_025137.4:c.88-1_88del",
        "NC_000001.11:g.1050_1100del",
        "NC_000001.11:g.(?_100)_(200_?)del",
        "NP_000001.1:p.Leu5_Lys10del",
    ] {
        assert_eq!(parse_hgvs_variant(s).expect(s).to_string(), s);
    }
}

fn provider() -> Provider {
    let genome = "ACGT".repeat(30);
    Provider::new()
        .sequence("NC_D.1", &genome)
        .sequence("NM_D.1", &genome[..100])
        .transcript(single_exon_transcript(
            "NM_D.1",
            "NC_D.1",
            0,
            Strand::Plus,
            10,
            39,
            100,
        ))
        .protein_for("NM_D.1", "NP_D.1")
}

fn coding_position(base: i32) -> BaseOffsetPosition {
    BaseOffsetPosition {
        base: HgvsTranscriptPos(base),
        offset: None,
        anchor: Anchor::CdsStart,
        uncertain: false,
    }
}

#[test]
fn a_range_built_backwards_is_an_error_wherever_it_is_resolved() {
    let hdp = provider();
    let mapper = VariantMapper::new(&hdp);
    let del = NaEdit::Del {
        ref_: None,
        uncertain: false,
    };
    let c = SequenceVariant::Coding(CVariant::from_parts(
        "NM_D.1".into(),
        None,
        PosEdit {
            pos: Some(BaseOffsetInterval {
                start: coding_position(20),
                end: Some(coding_position(5)),
                uncertain: false,
            }),
            edit: del.clone(),
            uncertain: false,
            predicted: false,
        },
    ));
    let g = SequenceVariant::Genomic(GVariant::from_parts(
        "NC_D.1".into(),
        None,
        PosEdit {
            pos: Some(SimpleInterval {
                start: SimplePosition {
                    base: HgvsGenomicPos(40),
                    end: None,
                    uncertain: false,
                },
                end: Some(SimplePosition {
                    base: HgvsGenomicPos(20),
                    end: None,
                    uncertain: false,
                }),
                uncertain: false,
            }),
            edit: del,
            uncertain: false,
            predicted: false,
        },
    ));
    let backwards = |r: Result<String, HgvsError>| {
        assert!(
            matches!(&r, Err(HgvsError::ValidationError(m)) if m.contains("runs backwards")),
            "{r:?}"
        );
    };
    backwards(mapper.validate(&c).map(|b| b.to_string()));
    backwards(mapper.to_spdi_unambiguous(&c));
    backwards(mapper.as_genomic(&c).unwrap().map(|v| v.to_string()));
    if let SequenceVariant::Coding(cv) = &c {
        backwards(mapper.c_to_p(cv, None).map(|p| p.to_string()));
    }
    backwards(mapper.validate(&g).map(|b| b.to_string()));
    backwards(mapper.to_spdi_unambiguous(&g));
}
