"""Statements about a whole transcript through the bindings: dicts, old and new."""

from __future__ import annotations

from typing import Any

import pytest

import weaver


def _special(variant_type: str, ac: str, value: str, *, predicted: bool = False) -> dict[str, Any]:
    """A dict in the shape statements were written in before they had types: a `Special` value."""
    return {
        "variant_type": variant_type,
        "ac": ac,
        "gene": None,
        "posedit": {
            "pos": None,
            "edit": {"type": "Special", "value": value, "uncertain": False},
            "uncertain": False,
            "predicted": predicted,
        },
    }


@pytest.mark.parametrize(
    "hgvs",
    ["NM_W.1:r.0?", "NM_W.1:r.spl?", "NM_W.1:r.?"],
)
def test_a_statement_round_trips_through_a_dict(hgvs: str) -> None:
    assert str(weaver.Variant.from_dict(weaver.parse(hgvs).to_dict())) == hgvs


def test_a_dict_names_an_rna_statement_and_its_uncertainty() -> None:
    edit = weaver.parse("NM_W.1:r.spl?").to_dict()["posedit"]["edit"]
    assert edit == {"type": "Statement", "statement": "SplicingAffected", "uncertain": True}


@pytest.mark.parametrize(
    ("variant_type", "ac", "value", "predicted", "expected"),
    [
        ("Rna", "NM_W.1", "spl?", False, "NM_W.1:r.spl?"),
        ("Rna", "NM_W.1", "0", False, "NM_W.1:r.0"),
    ],
)
def test_a_dict_in_the_special_shape_still_loads(
    variant_type: str,
    ac: str,
    value: str,
    predicted: bool,
    expected: str,
) -> None:
    loaded = weaver.Variant.from_dict(_special(variant_type, ac, value, predicted=predicted))
    assert str(loaded) == expected


def test_a_special_dict_for_a_statement_the_grammar_does_not_write_is_refused() -> None:
    with pytest.raises(ValueError, match="not a statement about RNA"):
        weaver.Variant.from_dict(_special("Rna", "NM_W.1", "??"))
