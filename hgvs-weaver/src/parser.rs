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

/// Parses a transcript-space posedit (`c.`, `n.` or `r.`); positions without
/// an explicit anchor get `default_anchor`. Only the r. grammar has a
/// predicted form, `r.(123a>g)`, its own rule. A statement about the whole
/// transcript, `c.=` or `r.spl`, has no position.
pub fn parse_tx_posedit(
    pair: Pair<Rule>,
    default_anchor: Anchor,
) -> Result<PosEdit<BaseOffsetInterval, NaEdit>, HgvsError> {
    let text = pair.as_str();
    let mut inner = pair.into_inner();
    let mut first = inner
        .next()
        .ok_or_else(|| HgvsError::PestError("Missing interval".into()))?;
    // The parentheses of r.(=) mark a prediction; those of r.(?) belong to
    // its statement.
    let statement = match first.as_rule() {
        Rule::whole_identity => Some((NaEdit::WholeIdentity, text.starts_with('('))),
        Rule::r_posedit_special => Some((
            NaEdit::Statement(parse_rna_statement(first.clone())?),
            false,
        )),
        _ => None,
    };
    if let Some((edit, predicted)) = statement {
        return Ok(PosEdit {
            pos: None,
            edit,
            uncertain: false,
            predicted,
        });
    }
    let predicted = first.as_rule() == Rule::r_posedit_predicted;
    if predicted {
        inner = first.into_inner();
        first = inner
            .next()
            .ok_or_else(|| HgvsError::PestError("Missing interval".into()))?;
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
        predicted,
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
    if inner_pair.as_rule() == Rule::p_unknown_from {
        let mut at = parse_aa_pos(
            inner_pair
                .into_inner()
                .next()
                .ok_or_else(|| HgvsError::PestError("Missing residue".into()))?,
        )?;
        at.aa = crate::utils::residue_3(&at.aa).to_string();
        return Ok(PosEdit {
            pos: Some(AaInterval {
                start: at,
                end: None,
                uncertain: false,
            }),
            edit: AaEdit::Statement(ProteinStatement::Unknown),
            uncertain: false,
            predicted,
        });
    }
    if inner_pair.as_rule() == Rule::p_posedit_special {
        return Ok(PosEdit {
            pos: None,
            edit: AaEdit::Statement(parse_protein_statement(inner_pair)?),
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
        if let Some(pos) = pos.as_mut() {
            pos.start.aa = edit.canonical_start(&pos.start.aa);
            if let Some(end) = pos.end.as_mut() {
                end.aa = crate::utils::residue_3(&end.aa).to_string();
            }
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

/// The statement an `r_posedit_special` names.
fn parse_rna_statement(pair: Pair<Rule>) -> Result<RnaStatement, HgvsError> {
    match statement_parts(pair)? {
        (Rule::r_no_rna, uncertain) => Ok(RnaStatement::NoRna { uncertain }),
        (Rule::r_splicing_affected, uncertain) => Ok(RnaStatement::SplicingAffected { uncertain }),
        (Rule::r_unknown, _) => Ok(RnaStatement::Unknown),
        (Rule::r_follows_dna, _) => Ok(RnaStatement::FollowsDna),
        (r, _) => Err(HgvsError::PestError(format!(
            "Unexpected RNA statement {r:?}"
        ))),
    }
}

/// The statement a `p_posedit_special` names.
fn parse_protein_statement(pair: Pair<Rule>) -> Result<ProteinStatement, HgvsError> {
    match statement_parts(pair)? {
        (Rule::p_unchanged, _) => Ok(ProteinStatement::Unchanged),
        (Rule::p_no_protein, uncertain) => Ok(ProteinStatement::NoProtein { uncertain }),
        (Rule::p_unknown, _) => Ok(ProteinStatement::Unknown),
        (r, _) => Err(HgvsError::PestError(format!(
            "Unexpected protein statement {r:?}"
        ))),
    }
}

/// The rule of a statement, and whether it carries `statement_uncertain`.
fn statement_parts(pair: Pair<Rule>) -> Result<(Rule, bool), HgvsError> {
    let statement = pair
        .into_inner()
        .next()
        .ok_or_else(|| HgvsError::PestError("Missing statement".into()))?;
    let rule = statement.as_rule();
    let uncertain = statement
        .into_inner()
        .any(|p| p.as_rule() == Rule::statement_uncertain);
    Ok((rule, uncertain))
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
            if let Some(end) = &end {
                if let (Some(a), Some(b)) = (start.base.known(), end.base.known()) {
                    ordered(a.0, b.0, s)?;
                }
            }
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
    let base = if s == "?" {
        GenomicBound::Unknown
    } else {
        GenomicBound::At(HgvsGenomicPos::try_from(number::<i32>(s)?).map_err(HgvsError::PestError)?)
    };
    Ok(SimplePosition {
        base,
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
    if let Some(end) = &end {
        ordered(start.order_key(), end.order_key(), s)?;
    }
    Ok(BaseOffsetInterval {
        start,
        end,
        uncertain,
    })
}

/// An HGVS range runs from its start to its end; `c.100_50del` names nothing.
/// Papers write such ranges (a typo, an OCR slip), and letting one through
/// used to panic when its bases were sliced.
fn ordered<K: Ord>(start: K, end: K, text: &str) -> Result<(), HgvsError> {
    if end < start {
        return Err(HgvsError::PestError(format!(
            "range {text} runs backwards: its start is after its end"
        )));
    }
    Ok(())
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

    let text = pair.as_str();
    let mut hgvs_base = None;
    let mut hgvs_offset: Option<IntronicOffset> = None;

    for p in pair.into_inner() {
        match p.as_rule() {
            Rule::num | Rule::base => {
                hgvs_base = Some(
                    HgvsTranscriptPos::try_from(number::<i32>(p.as_str())?)
                        .map_err(HgvsError::PestError)?,
                );
            }
            Rule::offset => {
                let off_str = p.as_str();
                if !off_str.is_empty() {
                    hgvs_offset = Some(IntronicOffset(number(&off_str.replace('+', ""))?));
                }
            }
            _ => {}
        }
    }

    Ok(BaseOffsetPosition {
        base: hgvs_base
            .ok_or_else(|| HgvsError::PestError(format!("position {text} has no base")))?,
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
    if let Some(end) = &end {
        ordered(start.base.0, end.base.0, s)?;
    }
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
            Rule::residue_num => {
                pos = p.as_str().parse().map_err(|_| {
                    HgvsError::PestError(format!("residue number {} is out of range", p.as_str()))
                })?;
            }
            _ => {}
        }
    }

    Ok(AAPosition {
        base: HgvsProteinPos(pos),
        aa,
    })
}

/// The text of a pair.
fn text(p: Pair<Rule>) -> String {
    p.as_str().to_string()
}

/// A number the grammar matched, in the type that holds it.
fn number<T: std::str::FromStr>(s: &str) -> Result<T, HgvsError> {
    s.parse()
        .map_err(|_| HgvsError::PestError(format!("number {s} is out of range")))
}

/// The text of a pair's first child, if it has one.
fn first_child_text(pair: Pair<Rule>) -> Option<String> {
    pair.into_inner().next().map(text)
}

/// `unit[n]` or `unit(min_max)`: the stated unit (a child for which `is_unit`
/// holds) and the copy count or range. Shared by nucleotide and protein repeats.
fn repeat_parts(
    pair: Pair<Rule>,
    is_unit: fn(Rule) -> bool,
) -> Result<(Option<String>, i32, i32), HgvsError> {
    let mut unit = None;
    let mut counts: Vec<i32> = Vec::new();
    for p in pair.into_inner() {
        if is_unit(p.as_rule()) {
            unit = Some(text(p));
        } else if p.as_rule() == Rule::num {
            counts.push(number(p.as_str())?);
        }
    }
    let (min, max) = min_max(&counts)?;
    Ok((unit, min, max))
}

/// The first count and the second, or the first again when there is one.
fn min_max<T: Copy>(counts: &[T]) -> Result<(T, T), HgvsError> {
    let min = *counts
        .first()
        .ok_or_else(|| HgvsError::PestError("a count range with no count".into()))?;
    Ok((min, counts.get(1).copied().unwrap_or(min)))
}

/// `N[20]`, `N[(20_30)]`, `(20)` or `(20_30)`: the number of inserted bases,
/// as `(min, max)`, equal when exact.
fn ins_length(pair: Pair<Rule>) -> Result<(usize, usize), HgvsError> {
    let counts = pair
        .into_inner()
        .flatten()
        .filter(|p| p.as_rule() == Rule::num)
        .map(|p| number(p.as_str()))
        .collect::<Result<Vec<usize>, _>>()?;
    min_max(&counts)
}

/// `ins` carries the inserted bases, or their number alone.
fn na_ins(pair: Pair<Rule>) -> Result<NaEdit, HgvsError> {
    let uncertain = false;
    Ok(match pair.into_inner().next() {
        Some(p) if p.as_rule() == Rule::ins_length => {
            let (min, max) = ins_length(p)?;
            NaEdit::InsLength {
                min,
                max,
                uncertain,
            }
        }
        p => NaEdit::Ins {
            alt: p.map(text),
            uncertain,
        },
    })
}

/// `delins` carries the inserted bases (or their number alone) and,
/// optionally, the deleted bases or their count before them.
fn na_delins(pair: Pair<Rule>) -> Result<NaEdit, HgvsError> {
    let mut parts: Vec<Pair<Rule>> = pair.into_inner().collect();
    let inserted = parts
        .pop()
        .ok_or_else(|| HgvsError::PestError("Malformed delins".into()))?;
    if parts.len() > 1 {
        return Err(HgvsError::PestError("Malformed delins".into()));
    }
    let deleted = parts.pop().map(text);
    if inserted.as_rule() == Rule::ins_length {
        let (min, max) = ins_length(inserted)?;
        return Ok(NaEdit::DelInsLength {
            ref_: deleted,
            min,
            max,
            uncertain: false,
        });
    }
    Ok(NaEdit::RefAlt {
        ref_: Some(deleted.unwrap_or_default()),
        alt: Some(text(inserted)),
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
        Rule::dna_ins | Rule::rna_ins => na_ins(edit)?,
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
            let (ref_, min, max) = repeat_parts(edit, |r| matches!(r, Rule::dna | Rule::rna))?;
            NaEdit::Repeat {
                ref_,
                min,
                max,
                uncertain,
            }
        }
        Rule::dna_copy => NaEdit::NACopy {
            copy: number(
                &first_child_text(edit)
                    .ok_or_else(|| HgvsError::PestError("copy with no count".into()))?,
            )?,
            uncertain,
        },
        // Conversions are parsed but not modelled.
        _ => NaEdit::None,
    })
}

/// `Xxx#Yyyfs*N`: the new residue, then the frameshift's stop and distance.
fn pro_fs(pair: Pair<Rule>) -> Result<AaEdit, HgvsError> {
    let mut alt = String::new();
    let mut stop = None;
    for p in pair.into_inner() {
        match p.as_rule() {
            Rule::aat13 => alt = text(p),
            Rule::fs => {
                let offset = p
                    .into_inner()
                    .flatten()
                    .find(|q| q.as_rule() == Rule::fsext_offset);
                if let Some(offset) = offset {
                    stop = Some(stop_distance(offset.as_str())?);
                }
            }
            _ => {}
        }
    }
    Ok(AaEdit::Fs {
        ref_: String::new(),
        alt,
        stop,
        uncertain: false,
    })
}

/// `fsext_offset`: the distance to a new stop, `12`, `?` or nothing.
fn stop_distance(s: &str) -> Result<StopDistance, HgvsError> {
    Ok(match s {
        "" => StopDistance::Unstated,
        "?" => StopDistance::Unknown,
        n => StopDistance::Residues { count: number(n)? },
    })
}

/// `Ter#Xxxext*N`, `Met1ext-5`: the residue that replaces the stop or start,
/// then what follows `ext`.
fn pro_ext(pair: Pair<Rule>) -> Result<AaEdit, HgvsError> {
    let mut alt = String::new();
    let mut end = ExtensionEnd::Unstated;
    for p in pair.into_inner() {
        match p.as_rule() {
            Rule::aat13 => alt = text(p),
            Rule::ext => {
                if let Some(tail) = p.into_inner().next() {
                    end = extension_end(tail)?;
                }
            }
            _ => {}
        }
    }
    Ok(AaEdit::Ext { alt, end })
}

/// `aa13_ext`: a stop and its distance, or an optional residue and a count,
/// signed or not.
fn extension_end(pair: Pair<Rule>) -> Result<ExtensionEnd, HgvsError> {
    let mut residue = None;
    for p in pair.into_inner() {
        match p.as_rule() {
            Rule::aa13 => residue = Some(text(p)),
            Rule::fsext_offset => {
                return Ok(ExtensionEnd::Stop {
                    distance: stop_distance(p.as_str())?,
                });
            }
            Rule::snum => {
                let n = p.as_str();
                return Ok(match n.strip_prefix('-') {
                    Some(n) => ExtensionEnd::Upstream {
                        residue,
                        distance: number(n)?,
                    },
                    None => ExtensionEnd::Count {
                        residue,
                        count: number(n.trim_start_matches('+'))?,
                    },
                });
            }
            _ => {}
        }
    }
    Err(HgvsError::PestError(
        "an extension's end has neither a stop distance nor a count".into(),
    ))
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
        Rule::pro_fs => pro_fs(edit)?,
        Rule::pro_ext => pro_ext(edit)?,
        Rule::pro_repeat => {
            let (ref_, min, max) = repeat_parts(edit, |r| r == Rule::aat13_seq)?;
            AaEdit::Repeat {
                ref_,
                min,
                max,
                uncertain,
            }
        }
        _ => AaEdit::None,
    }
    .into_canonical())
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
    fn insertions_of_a_stated_length_keep_the_length_not_letters() {
        use crate::edits::NaEdit;
        // Both the recommended and the older spellings parse to the same
        // edit and print in the recommended spelling.
        for (input, expected) in [
            (
                "NC_000001.11:g.123_124insN[20]",
                "NC_000001.11:g.123_124insN[20]",
            ),
            (
                "NC_000001.11:g.123_124ins(20)",
                "NC_000001.11:g.123_124insN[20]",
            ),
            (
                "NC_000001.11:g.123_124insN[(20_30)]",
                "NC_000001.11:g.123_124insN[(20_30)]",
            ),
            (
                "NC_000001.11:g.123_124insN[20_30]",
                "NC_000001.11:g.123_124insN[(20_30)]",
            ),
            (
                "NC_000001.11:g.123_124ins(20_30)",
                "NC_000001.11:g.123_124insN[(20_30)]",
            ),
            (
                "NM_004006.2:c.761_762insN[5]",
                "NM_004006.2:c.761_762insN[5]",
            ),
            (
                "NM_004006.2:r.761_762insn[5]",
                "NM_004006.2:r.761_762insN[5]",
            ),
        ] {
            let v = parse_hgvs_variant(input).unwrap();
            assert_eq!(v.to_string(), expected, "{input}");
        }
        let v = parse_hgvs_variant("NC_000001.11:g.123_124insN[(20_30)]").unwrap();
        let crate::coords::SequenceVariant::Genomic(g) = v else {
            panic!("expected a g. variant");
        };
        assert_eq!(
            g.posedit.edit,
            NaEdit::InsLength {
                min: 20,
                max: 30,
                uncertain: false
            }
        );
        // A single N is still one inserted base.
        let v = parse_hgvs_variant("NC_000001.11:g.123_124insN").unwrap();
        let crate::coords::SequenceVariant::Genomic(g) = v else {
            panic!("expected a g. variant");
        };
        assert_eq!(
            g.posedit.edit,
            NaEdit::Ins {
                alt: Some("N".into()),
                uncertain: false
            }
        );
    }

    #[test]
    fn an_identity_prints_the_bases_it_states() {
        for hgvs in [
            "NM_004006.2:c.123A=",
            "NC_000001.11:g.5_7CAG=",
            "NM_004006.2:r.123a=",
            "NM_004006.2:c.123=",
        ] {
            assert_eq!(parse_hgvs_variant(hgvs).unwrap().to_string(), hgvs);
        }
    }

    #[test]
    fn delins_of_a_stated_length_keeps_the_deleted_part() {
        use crate::edits::NaEdit;
        for (input, expected) in [
            (
                "NM_004006.2:c.812_829delinsN[12]",
                "NM_004006.2:c.812_829delinsN[12]",
            ),
            (
                "NM_004006.2:c.812_829del18insN[12]",
                "NM_004006.2:c.812_829del18insN[12]",
            ),
            (
                "NC_000001.11:g.812_814delACGins(2)",
                "NC_000001.11:g.812_814delACGinsN[2]",
            ),
            (
                "NC_000001.11:g.812_814delinsN[(2_4)]",
                "NC_000001.11:g.812_814delinsN[(2_4)]",
            ),
        ] {
            let v = parse_hgvs_variant(input).unwrap();
            assert_eq!(v.to_string(), expected, "{input}");
        }
        let v = parse_hgvs_variant("NC_000001.11:g.812_814delACGinsN[2]").unwrap();
        let crate::coords::SequenceVariant::Genomic(g) = v else {
            panic!("expected a g. variant");
        };
        assert_eq!(
            g.posedit.edit,
            NaEdit::DelInsLength {
                ref_: Some("ACG".into()),
                min: 2,
                max: 2,
                uncertain: false
            }
        );
        assert_eq!(g.posedit.edit.stated_ref(), Some("ACG"));
    }

    #[test]
    fn test_parse_extension() {
        let v_p = parse_hgvs_variant("NP_001116078.1:p.Ter312Argext*5").unwrap();
        match v_p {
            crate::coords::SequenceVariant::Protein(v) => match v.posedit.edit {
                crate::edits::AaEdit::Ext { alt, end } => {
                    assert_eq!(alt, "Arg");
                    assert_eq!(
                        end,
                        crate::edits::ExtensionEnd::Stop {
                            distance: crate::edits::StopDistance::Residues { count: 5 }
                        }
                    );
                }
                _ => panic!("Expected Extension edit"),
            },
            _ => panic!("Expected Protein variant"),
        }
    }

    /// One test per case: a number past its type's range in `$hgvs` is a
    /// parse error naming it, not a value read as 0.
    macro_rules! out_of_range_is_a_parse_error {
        ($($name:ident: $hgvs:literal, $number:literal;)*) => {$(
            #[test]
            fn $name() {
                let err = parse_hgvs_variant($hgvs).unwrap_err();
                assert!(
                    matches!(&err, crate::error::HgvsError::PestError(m) if m.contains($number)),
                    "{err}"
                );
            }
        )*};
    }

    out_of_range_is_a_parse_error! {
        a_genomic_position_past_i32: "NC_W.1:g.99999999999A>G", "99999999999";
        a_residue_number_past_i32: "NP_W.1:p.Arg99999999999Gly", "99999999999";
        a_transcript_position_past_i32: "NM_W.1:c.99999999999A>G", "99999999999";
        an_intronic_offset_past_i32: "NM_W.1:c.10+99999999999A>G", "99999999999";
        a_repeat_count_past_i32: "NM_W.1:c.10_12CAG[99999999999]", "99999999999";
        a_repeat_range_end_past_i32: "NM_W.1:c.10_12CAG[(5_99999999999)]", "99999999999";
        an_inserted_length_past_usize: "NC_W.1:g.10_11insN[99999999999999999999999]", "99999999999999999999999";
        a_copy_count_past_i32: "NC_W.1:g.10_20copy99999999999", "99999999999";
    }

    #[test]
    fn a_protein_repeat_with_a_number_before_its_count_is_a_parse_error() {
        assert!(
            matches!(
                parse_hgvs_variant("NP_W.1:p.(Gln18)5[(6_7)]"),
                Err(crate::error::HgvsError::PestError(_))
            ),
            "the 5 is neither the unit nor the count"
        );
    }

    #[test]
    fn a_frameshift_count_past_u32_is_a_parse_error() {
        let err = parse_hgvs_variant("NP_W.1:p.Arg97ProfsTer4294967296").unwrap_err();
        assert!(
            matches!(&err, crate::error::HgvsError::PestError(m) if m.contains("4294967296")),
            "{err}"
        );
    }

    #[test]
    fn an_extension_count_past_u32_is_a_parse_error() {
        let err = parse_hgvs_variant("NP_W.1:p.Ter30GlnextTer4294967296").unwrap_err();
        assert!(
            matches!(&err, crate::error::HgvsError::PestError(m) if m.contains("4294967296")),
            "{err}"
        );
    }
}
