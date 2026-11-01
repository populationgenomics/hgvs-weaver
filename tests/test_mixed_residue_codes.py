"""A protein description that mixes the amino acid codes is a ParseError a caller can catch."""

from __future__ import annotations

import pytest

import weaver


@pytest.mark.parametrize(
    "hgvs",
    ["NP_000021.1:p.Arg111X", "NP_000081.1:p.Ala356_Gly357insDPA", "NP_W.1:p.[Arg2Lys;L3P]"],
)
def test_a_description_mixing_the_codes_is_a_parse_error(hgvs: str) -> None:
    with pytest.raises(weaver.ParseError, match=r"in the three-letter amino acid code and .* in the one-letter code"):
        weaver.parse(hgvs)


def test_the_parse_error_names_the_residue_written_in_each_code() -> None:
    with pytest.raises(weaver.ParseError, match=r"writes Arg in the three-letter .* and X in the one-letter"):
        weaver.parse("NP_000021.1:p.Arg111X")
