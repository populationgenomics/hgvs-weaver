use super::Rule;
use crate::error::HgvsError;
use crate::structs::*;
use pest::iterators::Pair;

pub fn parse_g_posedit(pair: Pair<Rule>) -> Result<PosEdit<SimpleInterval, NaEdit>, HgvsError> {
    let mut inner = pair.into_inner();
    let pos = parse_simple_interval(
        inner
            .next()
            .ok_or_else(|| HgvsError::PestError("Missing interval".into()))?,
    )?;
    let edit = parse_na_edit(
        inner
            .next()
            .ok_or_else(|| HgvsError::PestError("Missing edit".into()))?,
    )?;
    Ok(PosEdit {
        pos: Some(pos),
        edit,
        uncertain: false,
        predicted: false,
    })
}

/// Parses a transcript-space posedit (`c.` or `n.`); positions without an
/// explicit anchor get `default_anchor`.
pub fn parse_tx_posedit(
    pair: Pair<Rule>,
    default_anchor: Anchor,
) -> Result<PosEdit<BaseOffsetInterval, NaEdit>, HgvsError> {
    let text = pair.as_str();
    let mut inner = pair.into_inner();
    let first = inner
        .next()
        .ok_or_else(|| HgvsError::PestError("Missing interval".into()))?;
    if first.as_rule() == Rule::r_posedit_special {
        return Ok(PosEdit {
            pos: None,
            edit: NaEdit::Special {
                value: first.as_str().replace(['(', ')'], ""),
                uncertain: false,
            },
            uncertain: false,
            predicted: text.starts_with('('),
        });
    }
    let pos = parse_base_offset_interval(first, default_anchor)?;
    let edit = parse_na_edit(
        inner
            .next()
            .ok_or_else(|| HgvsError::PestError("Missing edit".into()))?,
    )?;
    Ok(PosEdit {
        pos: Some(pos),
        edit,
        uncertain: false,
        predicted: false,
    })
}

pub fn parse_p_posedit(pair: Pair<Rule>) -> Result<PosEdit<AaInterval, AaEdit>, HgvsError> {
    let mut predicted = false;
    let s = pair.as_str();
    if s.starts_with('(') && s.ends_with(')') {
        predicted = true;
    }

    let mut inner = pair.into_inner();
    let inner_pair = inner
        .next()
        .ok_or_else(|| HgvsError::PestError("Empty p_posedit".into()))?;
    if inner_pair.as_rule() == Rule::p_posedit_special {
        let special = inner_pair.as_str();
        let edit = AaEdit::Special {
            value: special.replace(['(', ')'], ""),
            uncertain: false,
        };
        return Ok(PosEdit {
            pos: None,
            edit,
            uncertain: false,
            predicted,
        });
    }

    let mut pos = None;
    let mut edit = AaEdit::None;

    if inner_pair.as_rule() == Rule::p_interval {
        pos = Some(parse_aa_interval(inner_pair)?);
        if let Some(e) = inner.next() {
            edit = parse_pro_edit(e)?;
        }
    } else if inner_pair.as_rule() == Rule::pro_edit {
        edit = parse_pro_edit(inner_pair)?;
    }

    Ok(PosEdit {
        pos,
        edit,
        uncertain: false,
        predicted,
    })
}

pub fn parse_simple_interval(pair: Pair<Rule>) -> Result<SimpleInterval, HgvsError> {
    let s = pair.as_str();
    let mut uncertain = false;
    if s.starts_with('(') && s.ends_with(')') && !s.contains('_') {
        uncertain = true;
    }

    let mut inner = pair.into_inner();
    let p = inner
        .next()
        .ok_or_else(|| HgvsError::PestError("Empty interval".into()))?;
    match p.as_rule() {
        Rule::def_g_interval | Rule::def_m_interval => {
            let mut parts = p.into_inner();
            let start = parse_simple_pos(
                parts
                    .next()
                    .ok_or_else(|| HgvsError::PestError("Missing start position".into()))?,
            )?;
            let end = parts.next().map(parse_simple_pos).transpose()?;
            Ok(SimpleInterval {
                start,
                end,
                uncertain,
            })
        }
        Rule::uncertain_g_interval => {
            let mut start = None;
            let mut end = None;
            // A bound is uncertain when it is parenthesised, `(a_b)` or `(a)`;
            // the grammar does not keep the parentheses, so look at the text.
            let text = p.as_str();
            let base = p.as_span().start();
            for sub in p.into_inner() {
                if sub.as_rule() == Rule::def_g_interval {
                    let at = sub.as_span().start() - base;
                    let parenthesised = at > 0 && text.as_bytes()[at - 1] == b'(';
                    let mut parts = sub.into_inner();
                    let s =
                        parse_simple_pos(parts.next().ok_or_else(|| {
                            HgvsError::PestError("Missing start position".into())
                        })?)?;
                    let e = parts.next().map(parse_simple_pos).transpose()?;

                    let pos = SimplePosition {
                        base: s.base,
                        end: e.map(|x| x.base),
                        uncertain: parenthesised,
                    };

                    if start.is_none() {
                        start = Some(pos);
                    } else {
                        end = Some(pos);
                    }
                }
            }
            Ok(SimpleInterval {
                start: start.ok_or_else(|| {
                    HgvsError::PestError("Missing start position in uncertain interval".into())
                })?,
                end,
                uncertain: false,
            })
        }
        _ => Err(HgvsError::PestError(format!(
            "Unexpected interval rule: {:?}",
            p.as_rule()
        ))),
    }
}

pub fn parse_simple_pos(pair: Pair<Rule>) -> Result<SimplePosition, HgvsError> {
    let s = pair.as_str();
    if s == "?" {
        // Unknown is carried by the base itself; `?` is not parenthesised.
        return Ok(SimplePosition {
            base: HgvsGenomicPos::UNKNOWN,
            end: None,
            uncertain: false,
        });
    }
    let hgvs_base = s
        .parse::<i32>()
        .map_err(|_| HgvsError::PestError("Invalid position".into()))?;
    Ok(SimplePosition {
        base: HgvsGenomicPos(hgvs_base),
        end: None,
        uncertain: false,
    })
}

pub fn parse_c_posedit(pair: Pair<Rule>) -> Result<PosEdit<BaseOffsetInterval, NaEdit>, HgvsError> {
    parse_tx_posedit(pair, Anchor::CdsStart)
}

pub fn parse_n_posedit(pair: Pair<Rule>) -> Result<PosEdit<BaseOffsetInterval, NaEdit>, HgvsError> {
    parse_tx_posedit(pair, Anchor::TranscriptStart)
}

pub fn parse_base_offset_interval(
    pair: Pair<Rule>,
    default_anchor: Anchor,
) -> Result<BaseOffsetInterval, HgvsError> {
    let mut uncertain = false;
    let s = pair.as_str();
    if s.starts_with('(') && s.ends_with(')') {
        uncertain = true;
    }

    let mut inner = pair.into_inner();
    let p = inner
        .next()
        .ok_or_else(|| HgvsError::PestError("Empty base offset interval".into()))?;
    let mut p_inner = p.into_inner();

    let start = parse_base_offset_pos_with_default(
        p_inner
            .next()
            .ok_or_else(|| HgvsError::PestError("Missing start position".into()))?,
        default_anchor,
    )?;
    let end = p_inner
        .next()
        .map(|p| parse_base_offset_pos_with_default(p, default_anchor))
        .transpose()?;
    Ok(BaseOffsetInterval {
        start,
        end,
        uncertain,
    })
}

pub fn parse_base_offset_pos(pair: Pair<Rule>) -> Result<BaseOffsetPosition, HgvsError> {
    parse_base_offset_pos_with_default(pair, Anchor::CdsStart)
}

pub fn parse_base_offset_pos_with_default(
    pair: Pair<Rule>,
    default_anchor: Anchor,
) -> Result<BaseOffsetPosition, HgvsError> {
    let mut anchor = default_anchor;
    if pair.as_str().starts_with('*') {
        anchor = Anchor::CdsEnd;
    }

    let mut hgvs_base = 0;
    let mut hgvs_offset: Option<IntronicOffset> = None;

    for p in pair.into_inner() {
        match p.as_rule() {
            Rule::num | Rule::base => {
                hgvs_base = p.as_str().parse().unwrap_or(0);
            }
            Rule::offset => {
                let off_str = p.as_str();
                if !off_str.is_empty() {
                    hgvs_offset = Some(IntronicOffset(
                        off_str.replace('+', "").parse().unwrap_or(0),
                    ));
                }
            }
            _ => {}
        }
    }

    Ok(BaseOffsetPosition {
        base: HgvsTranscriptPos(hgvs_base),
        offset: hgvs_offset,
        anchor,
        uncertain: false,
    })
}

pub fn parse_aa_interval(pair: Pair<Rule>) -> Result<AaInterval, HgvsError> {
    let s = pair.as_str();
    let mut uncertain = false;
    if s.starts_with('(') && s.ends_with(')') {
        uncertain = true;
    }

    let mut inner = pair.into_inner();
    let p = inner
        .next()
        .ok_or_else(|| HgvsError::PestError("Empty AA interval".into()))?;
    let mut p_inner = p.into_inner();

    let start = parse_aa_pos(
        p_inner
            .next()
            .ok_or_else(|| HgvsError::PestError("Missing start AA position".into()))?,
    )?;
    let end = p_inner.next().map(parse_aa_pos).transpose()?;
    Ok(AaInterval {
        start,
        end,
        uncertain,
    })
}

pub fn parse_aa_pos(pair: Pair<Rule>) -> Result<AAPosition, HgvsError> {
    let mut aa = String::new();
    let mut pos = 0;

    for p in pair.into_inner() {
        match p.as_rule() {
            Rule::aa13 | Rule::term13 | Rule::aa3 | Rule::aa1 | Rule::term3 | Rule::term1 => {
                aa = p.as_str().to_string();
            }
            Rule::num => {
                pos = p.as_str().parse().unwrap_or(0);
            }
            _ => {}
        }
    }

    Ok(AAPosition {
        base: HgvsProteinPos(pos),
        aa,
        uncertain: false,
    })
}

/// The text of a pair.
fn text(p: Pair<Rule>) -> String {
    p.as_str().to_string()
}

/// The text of a pair's first child, if it has one.
fn first_child_text(pair: Pair<Rule>) -> Option<String> {
    pair.into_inner().next().map(text)
}

/// `unit[n]` or `unit(min_max)`: the stated unit (a child for which `is_unit`
/// holds) and the copy count or range. Shared by nucleotide and protein repeats.
fn repeat_parts(pair: Pair<Rule>, is_unit: fn(Rule) -> bool) -> (Option<String>, i32, i32) {
    let mut unit = None;
    let mut counts: Vec<i32> = Vec::new();
    for p in pair.into_inner() {
        if is_unit(p.as_rule()) {
            unit = Some(text(p));
        } else if p.as_rule() == Rule::num {
            counts.push(p.as_str().parse().unwrap_or(0));
        }
    }
    let min = counts.first().copied().unwrap_or(0);
    let max = counts.get(1).copied().unwrap_or(min);
    (unit, min, max)
}

/// `delins` carries the inserted bases and, optionally, the deleted bases or
/// their count before them.
fn na_delins(pair: Pair<Rule>) -> Result<NaEdit, HgvsError> {
    let parts: Vec<String> = pair.into_inner().map(text).collect();
    let (ref_, alt) = match parts.as_slice() {
        [alt] => (String::new(), alt.clone()),
        [ref_, alt] => (ref_.clone(), alt.clone()),
        _ => return Err(HgvsError::PestError("Malformed delins".into())),
    };
    Ok(NaEdit::RefAlt {
        ref_: Some(ref_),
        alt: Some(alt),
        uncertain: false,
    })
}

/// `=`, optionally preceded by the bases that are unchanged.
fn na_ident(pair: Pair<Rule>) -> NaEdit {
    let stated = pair
        .into_inner()
        .find(|p| matches!(p.as_rule(), Rule::dna | Rule::rna))
        .map(text);
    NaEdit::RefAlt {
        ref_: stated.clone(),
        alt: stated,
        uncertain: false,
    }
}

pub fn parse_na_edit(pair: Pair<Rule>) -> Result<NaEdit, HgvsError> {
    let edit = pair
        .into_inner()
        .next()
        .ok_or_else(|| HgvsError::PestError("Empty na_edit".into()))?;
    let uncertain = false;
    Ok(match edit.as_rule() {
        Rule::dna_subst | Rule::rna_subst => {
            let mut parts = edit.into_inner().map(text);
            NaEdit::RefAlt {
                ref_: parts.next(),
                alt: parts.next(),
                uncertain,
            }
        }
        Rule::dna_del | Rule::rna_del => NaEdit::Del {
            ref_: first_child_text(edit),
            uncertain,
        },
        Rule::dna_ins | Rule::rna_ins => NaEdit::Ins {
            alt: first_child_text(edit),
            uncertain,
        },
        Rule::dna_delins | Rule::rna_delins => na_delins(edit)?,
        Rule::dna_dup | Rule::rna_dup => NaEdit::Dup {
            ref_: first_child_text(edit),
            uncertain,
        },
        Rule::dna_inv | Rule::rna_inv => NaEdit::Inv {
            ref_: first_child_text(edit),
            uncertain,
        },
        Rule::dna_ident | Rule::rna_ident => na_ident(edit),
        Rule::dna_repeat | Rule::rna_repeat => {
            let (ref_, min, max) = repeat_parts(edit, |r| matches!(r, Rule::dna | Rule::rna));
            NaEdit::Repeat {
                ref_,
                min,
                max,
                uncertain,
            }
        }
        Rule::dna_copy => NaEdit::NACopy {
            copy: first_child_text(edit)
                .and_then(|n| n.parse().ok())
                .unwrap_or(0),
            uncertain,
        },
        // Conversions are parsed but not modelled.
        _ => NaEdit::None,
    })
}

/// The parts of `fs...` or `ext...`: the terminator or residue named, and the
/// count or offset, wherever the grammar nests them.
fn fs_ext_parts(pair: Pair<Rule>) -> (Option<String>, Option<String>) {
    let mut named = None;
    let mut count = None;
    for p in pair.into_inner().flatten() {
        match p.as_rule() {
            Rule::term13 | Rule::aa13 => named = Some(text(p)),
            Rule::fsext_offset | Rule::snum => count = Some(text(p)),
            _ => {}
        }
    }
    (named, count)
}

/// `Xxx#Yyyfs*N`: the new residue, then the frameshift's terminator and distance.
fn pro_fs(pair: Pair<Rule>) -> AaEdit {
    let mut alt = String::new();
    let mut term = None;
    let mut length = None;
    for p in pair.into_inner() {
        match p.as_rule() {
            Rule::aat13 => alt = text(p),
            Rule::fs => (term, length) = fs_ext_parts(p),
            _ => {}
        }
    }
    AaEdit::Fs {
        ref_: String::new(),
        alt,
        term,
        length,
        uncertain: false,
    }
}

/// `Ter#Xxxext*N`: the residue read through the stop, then the extension's
/// terminator and distance (or a residue and signed offset).
fn pro_ext(pair: Pair<Rule>) -> AaEdit {
    let mut alt = String::new();
    let mut aaterm = None;
    let mut length = None;
    for p in pair.into_inner() {
        match p.as_rule() {
            Rule::aat13 => alt = text(p),
            Rule::ext => (aaterm, length) = fs_ext_parts(p),
            _ => {}
        }
    }
    AaEdit::Ext {
        ref_: String::new(),
        alt,
        aaterm,
        length,
        uncertain: false,
    }
}

pub fn parse_pro_edit(pair: Pair<Rule>) -> Result<AaEdit, HgvsError> {
    let edit = pair
        .into_inner()
        .next()
        .ok_or_else(|| HgvsError::PestError("Empty pro_edit".into()))?;
    let uncertain = false;
    Ok(match edit.as_rule() {
        Rule::pro_ident => AaEdit::Identity { uncertain },
        Rule::pro_subst => AaEdit::Subst {
            ref_: String::new(),
            alt: text(edit),
            uncertain,
        },
        Rule::pro_del => AaEdit::Del {
            ref_: String::new(),
            uncertain,
        },
        Rule::pro_ins => AaEdit::Ins {
            alt: first_child_text(edit).unwrap_or_default(),
            uncertain,
        },
        Rule::pro_dup => AaEdit::Dup {
            ref_: None,
            uncertain,
        },
        Rule::pro_delins => AaEdit::DelIns {
            ref_: String::new(),
            alt: first_child_text(edit).unwrap_or_default(),
            uncertain,
        },
        Rule::pro_fs => pro_fs(edit),
        Rule::pro_ext => pro_ext(edit),
        Rule::pro_repeat => {
            let (ref_, min, max) = repeat_parts(edit, |r| r == Rule::aat13_seq);
            AaEdit::Repeat {
                ref_,
                min,
                max,
                uncertain,
            }
        }
        _ => AaEdit::None,
    })
}

#[cfg(test)]
mod tests {
    use crate::parse_hgvs_variant;

    #[test]
    fn test_unsupported_circular_no_panic() {
        // 'o.' is mentioned in docs but not supported in SequenceVariant.
        // Ensure it doesn't panic, just returns an error.
        let result = parse_hgvs_variant("NC_000001.11:o.123A>G");
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_repeats() {
        let v_c = parse_hgvs_variant("NM_001291285.3:c.7035TGGAAC[3]").unwrap();
        match v_c {
            crate::coords::SequenceVariant::Coding(v) => match v.posedit.edit {
                crate::edits::NaEdit::Repeat { ref_, min, max, .. } => {
                    assert_eq!(ref_, Some("TGGAAC".to_string()));
                    assert_eq!(min, 3);
                    assert_eq!(max, 3);
                }
                _ => panic!("Expected Repeat edit"),
            },
            _ => panic!("Expected Coding variant"),
        }

        let v_p = parse_hgvs_variant("NP_001278214.1:p.2346GT[3]").unwrap();
        match v_p {
            crate::coords::SequenceVariant::Protein(v) => match v.posedit.edit {
                crate::edits::AaEdit::Repeat { ref_, min, max, .. } => {
                    assert_eq!(ref_, Some("GT".to_string()));
                    assert_eq!(min, 3);
                    assert_eq!(max, 3);
                }
                _ => panic!("Expected Repeat edit for protein"),
            },
            _ => panic!("Expected Protein variant"),
        }
    }

    #[test]
    fn test_parse_extension() {
        let v_p = parse_hgvs_variant("NP_001116078.1:p.Ter312Argext*5").unwrap();
        match v_p {
            crate::coords::SequenceVariant::Protein(v) => match v.posedit.edit {
                crate::edits::AaEdit::Ext {
                    ref_,
                    alt,
                    aaterm,
                    length,
                    ..
                } => {
                    assert_eq!(ref_, "");
                    assert_eq!(alt, "Arg");
                    assert_eq!(aaterm, Some("*".to_string()));
                    assert_eq!(length, Some("5".to_string()));
                }
                _ => panic!("Expected Extension edit"),
            },
            _ => panic!("Expected Protein variant"),
        }
    }
}
