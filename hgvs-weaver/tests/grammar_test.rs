//! The biocommons `hgvs` grammar test table (`tests/data/grammar_test.tsv`),
//! run against our pest grammar rule by rule. A row names a rule, a `|`-separated
//! list of inputs, and whether they must all parse (`Valid`). Rules the two
//! grammars do not share are skipped. Rows the grammars disagree on are listed
//! in `KNOWN_DIFFERENCES` so that new disagreements fail the test.

use hgvs_weaver::{HgvsParser, Rule};
use pest::Parser;
use std::fs;

fn rule(name: &str) -> Option<Rule> {
    match name {
        "aa1" => Some(Rule::aa1),
        "aa13" => Some(Rule::aa13),
        "aa13_ext" => Some(Rule::aa13_ext),
        "aa13_fs" => Some(Rule::aa13_fs),
        "aa3" => Some(Rule::aa3),
        "aat1" => Some(Rule::aat1),
        "aat13" => Some(Rule::aat13),
        "aat13_seq" => Some(Rule::aat13_seq),
        "aat3" => Some(Rule::aat3),
        "accn" => Some(Rule::accn),
        "base" => Some(Rule::base),
        "c_hgvs_position" => Some(Rule::c_hgvs_position),
        "c_interval" => Some(Rule::c_interval),
        "c_pos" => Some(Rule::c_pos),
        "c_posedit" => Some(Rule::c_posedit),
        "c_variant" => Some(Rule::c_variant),
        "def_c_interval" => Some(Rule::def_c_interval),
        "def_g_interval" => Some(Rule::def_g_interval),
        "def_m_interval" => Some(Rule::def_m_interval),
        "def_n_interval" => Some(Rule::def_n_interval),
        "def_p_interval" => Some(Rule::def_p_interval),
        "def_r_interval" => Some(Rule::def_r_interval),
        "dna" => Some(Rule::dna),
        "dna_con" => Some(Rule::dna_con),
        "dna_copy" => Some(Rule::dna_copy),
        "dna_del" => Some(Rule::dna_del),
        "dna_delins" => Some(Rule::dna_delins),
        "dna_dup" => Some(Rule::dna_dup),
        "dna_edit" => Some(Rule::dna_edit),
        "dna_ident" => Some(Rule::dna_ident),
        "dna_ins" => Some(Rule::dna_ins),
        "dna_inv" => Some(Rule::dna_inv),
        "dna_subst" => Some(Rule::dna_subst),
        "ext" => Some(Rule::ext),
        "fs" => Some(Rule::fs),
        "fsext_offset" => Some(Rule::fsext_offset),
        "g_hgvs_position" => Some(Rule::g_hgvs_position),
        "g_interval" => Some(Rule::g_interval),
        "g_pos" => Some(Rule::g_pos),
        "g_posedit" => Some(Rule::g_posedit),
        "g_variant" => Some(Rule::g_variant),
        "gene_symbol" => Some(Rule::gene_symbol),
        "hgvs_position" => Some(Rule::hgvs_position),
        "hgvs_variant" => Some(Rule::hgvs_variant),
        "m_hgvs_position" => Some(Rule::m_hgvs_position),
        "m_interval" => Some(Rule::m_interval),
        "m_pos" => Some(Rule::m_pos),
        "m_posedit" => Some(Rule::m_posedit),
        "m_variant" => Some(Rule::m_variant),
        "n_hgvs_position" => Some(Rule::n_hgvs_position),
        "n_interval" => Some(Rule::n_interval),
        "n_pos" => Some(Rule::n_pos),
        "n_posedit" => Some(Rule::n_posedit),
        "n_variant" => Some(Rule::n_variant),
        "num" => Some(Rule::num),
        "offset" => Some(Rule::offset),
        "opt_gene_expr" => Some(Rule::opt_gene_expr),
        "p_hgvs_position" => Some(Rule::p_hgvs_position),
        "p_interval" => Some(Rule::p_interval),
        "p_pos" => Some(Rule::p_pos),
        "p_posedit" => Some(Rule::p_posedit),
        "p_posedit_special" => Some(Rule::p_posedit_special),
        "p_variant" => Some(Rule::p_variant),
        "pro_del" => Some(Rule::pro_del),
        "pro_delins" => Some(Rule::pro_delins),
        "pro_dup" => Some(Rule::pro_dup),
        "pro_edit" => Some(Rule::pro_edit),
        "pro_ext" => Some(Rule::pro_ext),
        "pro_fs" => Some(Rule::pro_fs),
        "pro_ident" => Some(Rule::pro_ident),
        "pro_ins" => Some(Rule::pro_ins),
        "pro_subst" => Some(Rule::pro_subst),
        "r_hgvs_position" => Some(Rule::r_hgvs_position),
        "r_interval" => Some(Rule::r_interval),
        "r_pos" => Some(Rule::r_pos),
        "r_posedit" => Some(Rule::r_posedit),
        "r_variant" => Some(Rule::r_variant),
        "rna" => Some(Rule::rna),
        "rna_con" => Some(Rule::rna_con),
        "rna_del" => Some(Rule::rna_del),
        "rna_delins" => Some(Rule::rna_delins),
        "rna_dup" => Some(Rule::rna_dup),
        "rna_edit" => Some(Rule::rna_edit),
        "rna_ident" => Some(Rule::rna_ident),
        "rna_ins" => Some(Rule::rna_ins),
        "rna_inv" => Some(Rule::rna_inv),
        "rna_subst" => Some(Rule::rna_subst),
        "snum" => Some(Rule::snum),
        "term1" => Some(Rule::term1),
        "term13" => Some(Rule::term13),
        "term3" => Some(Rule::term3),
        "uncertain_g_interval" => Some(Rule::uncertain_g_interval),
        _ => None,
    }
}

/// Whether `input` is entirely consumed by `rule`.
fn accepts(rule: Rule, input: &str) -> bool {
    HgvsParser::parse(rule, input)
        .map(|mut pairs| {
            pairs
                .next()
                .is_some_and(|p| p.as_str().len() == input.len())
        })
        .unwrap_or(false)
}

/// `(rule, input)` pairs where our grammar deliberately or currently differs
/// from biocommons hgvs. Remove an entry when the difference is resolved.
const KNOWN_DIFFERENCES: &[(&str, &str)] = &[
    // ClinVar writes stop codons mid-insertion (p.Glu26_Glu27insTerGlu); we
    // accept them so that such descriptions can be compared.
    ("aat13_seq", "TerGly"),
];

#[test]
fn biocommons_grammar_table() {
    let text = fs::read_to_string("../tests/data/grammar_test.tsv").expect("grammar table");
    let mut checked = 0;
    let mut unexpected: Vec<String> = Vec::new();
    for line in text.lines().skip(1) {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 3 || cols[0].is_empty() {
            continue;
        }
        let (func, tests, valid) = (cols[0], cols[1], cols[2] == "True");
        let Some(r) = rule(func) else { continue };
        // `InType` says how the Test column is split: a `string` is one input
        // per character, a `list` one per `|`.
        let inputs: Vec<String> = if cols.get(3).copied() == Some("string") {
            tests.chars().map(|c| c.to_string()).collect()
        } else {
            tests.split('|').map(str::to_string).collect()
        };
        for input in inputs.iter().map(String::as_str) {
            checked += 1;
            let ok = accepts(r, input);
            let known = KNOWN_DIFFERENCES.contains(&(func, input));
            if ok != valid && !known {
                unexpected.push(format!(
                    "{func}\t{input}\texpected {} got {}",
                    if valid { "accept" } else { "reject" },
                    if ok { "accept" } else { "reject" }
                ));
            }
            if ok == valid && known {
                unexpected.push(format!(
                    "{func}\t{input}\tnow agrees; remove from KNOWN_DIFFERENCES"
                ));
            }
        }
    }
    assert!(checked > 100, "only {checked} inputs checked");
    assert!(
        unexpected.is_empty(),
        "{} unexpected grammar disagreements:\n{}",
        unexpected.len(),
        unexpected.join("\n")
    );
}
