mod support;

use hgvs_weaver::*;
use support::Provider;

fn run_c_to_p(hgvsc: &str, expected_p: &str) {
    let hdp = Provider::from_json_file("../tests/data/toy_data.json");
    let mapper = VariantMapper::new(&hdp);
    let var_c = parse_hgvs_variant(hgvsc).unwrap();
    if let SequenceVariant::Coding(v) = var_c {
        let var_p = mapper.c_to_p(&v, Some("MOCK")).unwrap();
        assert_eq!(
            var_p.to_string(),
            expected_p,
            "Conversion failed for {}",
            hgvsc
        );
    } else {
        panic!("Expected coding variant for {}", hgvsc);
    }
}

#[test]
fn test_parity_substitutions() {
    run_c_to_p("NM_999999.1:c.6A>G", "MOCK:p.(Lys2=)");
    run_c_to_p("NM_999999.1:c.6A>T", "MOCK:p.(Lys2Asn)");
    run_c_to_p("NM_999996.1:c.8C>A", "MOCK:p.(Ser3Ter)");
}

#[test]
fn test_parity_insertions() {
    run_c_to_p("NM_999999.1:c.6_7insGGG", "MOCK:p.(Lys2_Ala3insGly)");
    // Rust currently produces p.(Ala8ValfsTer?) for frameshifts with unknown stop.
    // The Python test says p.(Ala8ValfsTer?).
    run_c_to_p("NM_999999.1:c.22_23insT", "MOCK:p.(Ala8ValfsTer?)");
    run_c_to_p("NM_999999.1:c.8_9insTT", "MOCK:p.(Lys4Ter)");
}

#[test]
fn test_parity_deletions() {
    run_c_to_p("NM_999999.1:c.10_12del", "MOCK:p.(Lys4del)");
    run_c_to_p("NM_999999.1:c.4_15del", "MOCK:p.(Lys2_Ala5del)");
    run_c_to_p("NM_999995.1:c.4_6del", "MOCK:p.(Lys3del)");
    run_c_to_p("NM_999994.1:c.4_9del", "MOCK:p.(Lys3_Lys4del)");
    run_c_to_p("NM_999999.1:c.5_7del", "MOCK:p.(Lys2_Ala3delinsThr)");
    run_c_to_p("NM_999993.1:c.13_24del", "MOCK:p.(Arg5_Ala8del)");
}

#[test]
fn test_parity_frameshifts() {
    run_c_to_p("NM_999999.1:c.11_12del", "MOCK:p.(Lys4SerfsTer?)");
    run_c_to_p("NM_999997.1:c.7del", "MOCK:p.(Ala3ArgfsTer6)");
}

#[test]
fn test_parity_indels() {
    run_c_to_p(
        "NM_999999.1:c.11_12delinsTCCCA",
        "MOCK:p.(Lys4delinsIlePro)",
    );
    run_c_to_p(
        "NM_999999.1:c.11_18delinsTCCCA",
        "MOCK:p.(Lys4_Phe6delinsIlePro)",
    );
}
