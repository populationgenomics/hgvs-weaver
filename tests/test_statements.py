"""Statements about a whole transcript or protein through the bindings: dicts, old and new."""

from __future__ import annotations

from typing import Any

import pytest

import weaver


def _special(variant_type: str, ac: str, value: str, *, predicted: bool = False) -> dict[str, Any]:
    """A dict in the `Special` shape: the statement as its HGVS text."""
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
    ["NM_W.1:r.spl?", "NP_W.1:p.0", "NP_W.1:p.(=)"],
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
        ("Protein", "NP_W.1", "=", True, "NP_W.1:p.(=)"),
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
