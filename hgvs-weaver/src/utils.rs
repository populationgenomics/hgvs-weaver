use crate::error::HgvsError;
use std::fmt::{self, Display, Formatter};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Residue {
    Ala,
    Arg,
    Asn,
    Asp,
    Cys,
    Gln,
    Glu,
    Gly,
    His,
    Ile,
    Leu,
    Lys,
    Met,
    Phe,
    Pro,
    Ser,
    Thr,
    Trp,
    Tyr,
    Val,
    // Ambiguous / Special
    Asx, // B
    Glx, // Z
    Xaa, // X
    Sec, // U
    Ter, // *
}

impl Display for Residue {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let s = match self {
            Residue::Ala => "A",
            Residue::Arg => "R",
            Residue::Asn => "N",
            Residue::Asp => "D",
            Residue::Cys => "C",
            Residue::Gln => "Q",
            Residue::Glu => "E",
            Residue::Gly => "G",
            Residue::His => "H",
            Residue::Ile => "I",
            Residue::Leu => "L",
            Residue::Lys => "K",
            Residue::Met => "M",
            Residue::Phe => "F",
            Residue::Pro => "P",
            Residue::Ser => "S",
            Residue::Thr => "T",
            Residue::Trp => "W",
            Residue::Tyr => "Y",
            Residue::Val => "V",
            Residue::Asx => "B",
            Residue::Glx => "Z",
            Residue::Xaa => "X",
            Residue::Sec => "U",
            Residue::Ter => "*",
        };
        write!(f, "{}", s)
    }
}

pub fn reverse_complement(seq: &str) -> String {
    seq.chars().rev().map(complement_dna_char).collect()
}

pub fn complement_dna_char(c: char) -> char {
    match c {
        'A' => 'T',
        'T' => 'A',
        'C' => 'G',
        'G' => 'C',
        'N' => 'N',
        'a' => 't',
        't' => 'a',
        'c' => 'g',
        'g' => 'c',
        'n' => 'n',
        'U' => 'A',
        'u' => 'a',
        _ => c,
    }
}

pub fn aa1_to_aa3(aa1: char) -> &'static str {
    match aa1.to_ascii_uppercase() {
        'A' => "Ala",
        'R' => "Arg",
        'N' => "Asn",
        'D' => "Asp",
        'C' => "Cys",
        'E' => "Glu",
        'Q' => "Gln",
        'G' => "Gly",
        'H' => "His",
        'I' => "Ile",
        'L' => "Leu",
        'K' => "Lys",
        'M' => "Met",
        'F' => "Phe",
        'P' => "Pro",
        'S' => "Ser",
        'T' => "Thr",
        'W' => "Trp",
        'Y' => "Tyr",
        'V' => "Val",
        'B' => "Asx",
        'Z' => "Glx",
        'U' => "Sec",
        '*' => "Ter",
        'X' => "Xaa",
        _ => "Xaa",
    }
}

/// A residue as weaver stores and prints it: a one-letter code becomes its
/// three-letter code (`*` is `Ter`, `X` is `Xaa`); anything longer is kept.
pub fn residue_3(aa: &str) -> &str {
    let mut chars = aa.chars();
    match (chars.next(), chars.next()) {
        (Some(one), None) => aa1_to_aa3(one),
        _ => aa,
    }
}

/// [`residue_3`] where the grammar reads `X` as a stop: after `fs` or `ext`.
pub fn terminator_3(aa: &str) -> &str {
    match aa {
        "X" | "*" => "Ter",
        other => residue_3(other),
    }
}

pub fn seq1_to_aa3(seq1: &str) -> String {
    seq1.chars().map(aa1_to_aa3).collect()
}

/// Translates every complete codon of `seq`, without stopping at a stop codon.
/// Unknown codons become `X`. Trailing bases that do not fill a codon are dropped.
pub fn translate(seq: &str) -> String {
    let bytes = seq.as_bytes();
    let mut aa = String::with_capacity(bytes.len() / 3);
    for &codon in bytes.as_chunks::<3>().0 {
        aa.push(codon_residue(codon).unwrap_or('X'));
    }
    aa
}

/// Translates a single 3-base codon (DNA or RNA, any case) to a 1-letter amino acid;
/// `None` for anything but three bases.
pub fn translate_codon(codon: &str) -> Option<char> {
    codon_residue(codon.as_bytes().try_into().ok()?)
}

/// The residue `codon` encodes, read in either case and with `U` as `T`;
/// `None` when a byte is not a base.
fn codon_residue(codon: [u8; 3]) -> Option<char> {
    let codon = codon.map(|b| match b.to_ascii_uppercase() {
        b'U' => b'T',
        b => b,
    });
    match &codon {
        b"TTT" | b"TTC" => Some('F'),
        b"TTA" | b"TTG" => Some('L'),
        b"CTT" | b"CTC" | b"CTA" | b"CTG" => Some('L'),
        b"ATT" | b"ATC" | b"ATA" => Some('I'),
        b"ATG" => Some('M'),
        b"GTT" | b"GTC" | b"GTA" | b"GTG" => Some('V'),
        b"TCT" | b"TCC" | b"TCA" | b"TCG" => Some('S'),
        b"CCT" | b"CCC" | b"CCA" | b"CCG" => Some('P'),
        b"ACT" | b"ACC" | b"ACA" | b"ACG" => Some('T'),
        b"GCT" | b"GCC" | b"GCA" | b"GCG" => Some('A'),
        b"TAT" | b"TAC" => Some('Y'),
        b"TAA" | b"TAG" | b"TGA" => Some('*'),
        b"CAT" | b"CAC" => Some('H'),
        b"CAA" | b"CAG" => Some('Q'),
        b"AAT" | b"AAC" => Some('N'),
        b"AAA" | b"AAG" => Some('K'),
        b"GAT" | b"GAC" => Some('D'),
        b"GAA" | b"GAG" => Some('E'),
        b"TGT" | b"TGC" => Some('C'),
        b"TGG" => Some('W'),
        b"CGT" | b"CGC" | b"CGA" | b"CGG" => Some('R'),
        b"AGT" | b"AGC" => Some('S'),
        b"AGA" | b"AGG" => Some('R'),
        b"GGT" | b"GGC" | b"GGA" | b"GGG" => Some('G'),
        _ => None,
    }
}

/// Returns all DNA codons that encode the given 1-letter amino acid.
pub fn codons_for_aa(aa: char) -> Vec<&'static str> {
    match aa {
        'F' => vec!["TTT", "TTC"],
        'L' => vec!["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
        'I' => vec!["ATT", "ATC", "ATA"],
        'M' => vec!["ATG"],
        'V' => vec!["GTT", "GTC", "GTA", "GTG"],
        'S' => vec!["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
        'P' => vec!["CCT", "CCC", "CCA", "CCG"],
        'T' => vec!["ACT", "ACC", "ACA", "ACG"],
        'A' => vec!["GCT", "GCC", "GCA", "GCG"],
        'Y' => vec!["TAT", "TAC"],
        '*' => vec!["TAA", "TAG", "TGA"],
        'H' => vec!["CAT", "CAC"],
        'Q' => vec!["CAA", "CAG"],
        'N' => vec!["AAT", "AAC"],
        'K' => vec!["AAA", "AAG"],
        'D' => vec!["GAT", "GAC"],
        'E' => vec!["GAA", "GAG"],
        'C' => vec!["TGT", "TGC"],
        'W' => vec!["TGG"],
        'R' => vec!["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
        'G' => vec!["GGT", "GGC", "GGA", "GGG"],
        _ => vec![],
    }
}

pub fn aa3_to_aa1(aa3: &str) -> String {
    match aa3.to_lowercase().as_str() {
        "ala" => "A".to_string(),
        "arg" => "R".to_string(),
        "asn" => "N".to_string(),
        "asp" => "D".to_string(),
        "cys" => "C".to_string(),
        "gln" => "Q".to_string(),
        "glu" => "E".to_string(),
        "gly" => "G".to_string(),
        "his" => "H".to_string(),
        "ile" => "I".to_string(),
        "leu" => "L".to_string(),
        "lys" => "K".to_string(),
        "met" => "M".to_string(),
        "phe" => "F".to_string(),
        "pro" => "P".to_string(),
        "ser" => "S".to_string(),
        "thr" => "T".to_string(),
        "trp" => "W".to_string(),
        "tyr" => "Y".to_string(),
        "val" => "V".to_string(),
        "asx" => "B".to_string(),
        "glx" => "Z".to_string(),
        "xaa" => "X".to_string(),
        "ter" | "stop" | "*" => "*".to_string(),
        // 1-letter codes are returned as-is (uppercase)
        s if s.len() == 1 => s.to_uppercase(),
        _ => "X".to_string(),
    }
}

pub fn normalize_aa(s: &str) -> String {
    aa3_to_aa1(s)
}

fn aa3_chunk_to_residue(s: &str) -> Option<Residue> {
    match s.to_lowercase().as_str() {
        "ala" => Some(Residue::Ala),
        "arg" => Some(Residue::Arg),
        "asn" => Some(Residue::Asn),
        "asp" => Some(Residue::Asp),
        "cys" => Some(Residue::Cys),
        "gln" => Some(Residue::Gln),
        "glu" => Some(Residue::Glu),
        "gly" => Some(Residue::Gly),
        "his" => Some(Residue::His),
        "ile" => Some(Residue::Ile),
        "leu" => Some(Residue::Leu),
        "lys" => Some(Residue::Lys),
        "met" => Some(Residue::Met),
        "phe" => Some(Residue::Phe),
        "pro" => Some(Residue::Pro),
        "ser" => Some(Residue::Ser),
        "thr" => Some(Residue::Thr),
        "trp" => Some(Residue::Trp),
        "tyr" => Some(Residue::Tyr),
        "val" => Some(Residue::Val),
        "asx" => Some(Residue::Asx),
        "glx" => Some(Residue::Glx),
        "xaa" => Some(Residue::Xaa),
        "ter" | "stop" => Some(Residue::Ter),
        "sec" => Some(Residue::Sec),
        _ => None,
    }
}

/// The 1-letter form of an HGVS amino acid string (`GlyGly`, `GG`, `Ter`).
pub fn residues_1(s: &str) -> Result<String, HgvsError> {
    Ok(decompose_aa(s)?.iter().map(ToString::to_string).collect())
}

pub fn decompose_aa(s: &str) -> Result<Vec<Residue>, HgvsError> {
    if s.is_empty() {
        return Ok(Vec::new());
    }

    // Attempt 1: Try to parse as only 3-letter codes
    let mut res_3 = Vec::new();
    let chars: Vec<char> = s.chars().collect();
    let mut i = 0;
    let mut all_3 = true;

    while i < chars.len() {
        if i + 3 <= chars.len() {
            let chunk: String = chars[i..i + 3].iter().collect();
            let res = aa3_chunk_to_residue(&chunk);

            if let Some(r) = res {
                res_3.push(r);
                i += 3;
                continue;
            }
        }
        all_3 = false;
        break;
    }

    if all_3 && i == chars.len() {
        return Ok(res_3);
    }

    // Attempt 2: Try to parse as only 1-letter codes
    let mut res_1 = Vec::new();
    let mut all_1 = true;
    for c in s.chars() {
        let res = match c {
            'A' => Some(Residue::Ala),
            'R' => Some(Residue::Arg),
            'N' => Some(Residue::Asn),
            'D' => Some(Residue::Asp),
            'C' => Some(Residue::Cys),
            'Q' => Some(Residue::Gln),
            'E' => Some(Residue::Glu),
            'G' => Some(Residue::Gly),
            'H' => Some(Residue::His),
            'I' => Some(Residue::Ile),
            'L' => Some(Residue::Leu),
            'K' => Some(Residue::Lys),
            'M' => Some(Residue::Met),
            'F' => Some(Residue::Phe),
            'P' => Some(Residue::Pro),
            'S' => Some(Residue::Ser),
            'T' => Some(Residue::Thr),
            'W' => Some(Residue::Trp),
            'Y' => Some(Residue::Tyr),
            'V' => Some(Residue::Val),
            'B' => Some(Residue::Asx),
            'Z' => Some(Residue::Glx),
            'X' => Some(Residue::Xaa),
            'U' => Some(Residue::Sec),
            '*' => Some(Residue::Ter),
            _ => None,
        };

        if let Some(r) = res {
            res_1.push(r);
        } else {
            all_1 = false;
            break;
        }
    }

    if all_1 {
        return Ok(res_1);
    }

    Err(HgvsError::Other(format!(
        "Could not decompose protein sequence '{}' into exclusively 3-letter or 1-letter codes",
        s
    )))
}

#[cfg(test)]
mod tests {
    use super::{translate, translate_codon};

    /// One test per `name: input => expected` row, so each row fails alone.
    macro_rules! cases {
        ($function:ident { $($name:ident: $input:expr => $expected:expr,)* }) => {
            $(
                #[test]
                fn $name() {
                    assert_eq!($function($input), $expected);
                }
            )*
        };
    }

    // A codon reads alike in either case and as RNA.
    cases!(translate_codon {
        an_upper_case_dna_codon_is_read: "TGT" => Some('C'),
        a_lower_case_dna_codon_is_read: "tgt" => Some('C'),
        an_upper_case_rna_codon_is_read: "UGU" => Some('C'),
        a_lower_case_rna_codon_is_read: "ugu" => Some('C'),
        a_mixed_case_rna_codon_is_read: "tGu" => Some('C'),
    });

    // Anything but three bases is no codon.
    cases!(translate_codon {
        an_empty_string_is_no_codon: "" => None,
        two_bases_are_no_codon: "TG" => None,
        four_bases_are_no_codon: "TGTA" => None,
        an_n_is_no_base: "TGN" => None,
        a_gap_is_no_base: "TG-" => None,
        a_non_ascii_character_is_no_base: "Té" => None,
    });

    cases!(translate {
        an_unknown_codon_is_x_and_a_partial_codon_is_dropped: "atgNNNtaaGG" => "MX*",
        a_codon_holding_a_non_ascii_byte_is_x: "ATGéTAA" => "MX",
    });
}
