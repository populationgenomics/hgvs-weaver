//! Real RefSeq transcripts whose records differ from GRCh38: MUC2
//! (NM_002457.5) and SHANK3 (NM_001372044.2), with NCBI's own
//! transcript-to-genome alignments cut into per-exon CIGARs and the genome
//! served as a window per chromosome under its own accession. Each case
//! records what VariantValidator writes and what weaver 0.5.0 wrote; the
//! difference is the projected reference, which must be the target's.

mod support;

use hgvs_weaver::mapper::VariantMapper;
use hgvs_weaver::{parse_hgvs_variant, SequenceVariant};
use serde::Deserialize;
use support::Provider;

const FIXTURE: &str = "../tests/data/projection_reference_real.json";

#[derive(Deserialize)]
struct Fixture {
    proteins: std::collections::BTreeMap<String, String>,
    cases: Vec<Case>,
}

#[derive(Deserialize)]
struct Case {
    op: String,
    input: String,
    #[serde(default)]
    reference: Option<String>,
    #[serde(default)]
    transcript: Option<String>,
    expected: serde_json::Value,
    #[serde(default)]
    weaver_0_5_0: Option<String>,
    #[serde(default)]
    note: String,
}

#[test]
fn real_records_that_differ_from_the_genome_project_to_the_targets_bases() {
    let text = std::fs::read_to_string(FIXTURE).expect("fixture");
    let fixture: Fixture = serde_json::from_str(&text).expect("fixture shape");
    let mut hdp = Provider::from_json_file(FIXTURE);
    for (tx, np) in &fixture.proteins {
        hdp = hdp.protein_for(tx, np);
    }
    let mapper = VariantMapper::new(&hdp);

    let mut failures = Vec::new();
    for case in &fixture.cases {
        let var = parse_hgvs_variant(&case.input).expect("case input parses");
        let actual: serde_json::Value = match (case.op.as_str(), &var) {
            ("c_to_g", SequenceVariant::Coding(c)) => mapper
                .c_to_g(c, case.reference.as_deref())
                .map(|g| serde_json::Value::String(g.to_string())),
            ("g_to_c", SequenceVariant::Genomic(g)) => mapper
                .g_to_c(g, case.transcript.as_deref().expect("transcript"))
                .map(|c| serde_json::Value::String(c.to_string())),
            ("validate", _) => mapper.validate(&var).map(serde_json::Value::Bool),
            (op, _) => panic!("unknown op {op} for {}", case.input),
        }
        .unwrap_or_else(|e| serde_json::Value::String(format!("error: {e}")));
        if actual != case.expected {
            failures.push(format!(
                "{} {}\n  expected {}\n  actual   {}\n  0.5.0    {}\n  {}",
                case.op,
                case.input,
                case.expected,
                actual,
                case.weaver_0_5_0.as_deref().unwrap_or("-"),
                case.note
            ));
        }
    }
    assert!(failures.is_empty(), "\n{}", failures.join("\n\n"));
}
