"""c.=, the whole transcript unchanged, through the bindings: dicts and pickling."""

from __future__ import annotations

import pickle

import pytest

import weaver


def test_c_identity_round_trips_through_a_dict() -> None:
    var = weaver.parse("NM_W.1:c.=")
    assert str(weaver.Variant.from_dict(var.to_dict())) == "NM_W.1:c.="


def test_c_identity_round_trips_through_pickle() -> None:
    var = weaver.parse("NM_W.1:c.=")
    assert str(pickle.loads(pickle.dumps(var))) == "NM_W.1:c.="  # noqa: S301 -- the bytes are our own


def test_a_dict_naming_a_c_splicing_statement_is_refused() -> None:
    d = {
        "variant_type": "Coding",
        "ac": "NM_W.1",
        "gene": None,
        "posedit": {
            "pos": None,
            "edit": {"type": "Special", "value": "spl", "uncertain": False},
            "uncertain": False,
            "predicted": False,
        },
    }
    with pytest.raises(ValueError, match=r"c\. cannot write c\.spl"):
        weaver.Variant.from_dict(d)
