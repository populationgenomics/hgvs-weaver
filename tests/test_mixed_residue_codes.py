"""A protein description that mixes the amino acid codes raises an error a caller can recover from."""

from __future__ import annotations

import pytest

import weaver


@pytest.mark.parametrize(
    "hgvs",
    ["NP_W.1:p.Arg2K", "NP_000081.1:p.Ala356_Gly357insDPA", "NP_W.1:p.[Arg2Lys;L3P]"],
)
def test_a_description_mixing_the_codes_is_a_mixed_residue_codes_error(hgvs: str) -> None:
    with pytest.raises(weaver.MixedResidueCodesError):
        weaver.parse(hgvs)


def test_a_mixed_residue_codes_error_is_a_parse_error() -> None:
    with pytest.raises(weaver.ParseError):
        weaver.parse("NP_W.1:p.Arg2K")


def test_a_mixed_residue_codes_error_is_an_hgvs_error() -> None:
    with pytest.raises(weaver.HGVSError):
        weaver.parse("NP_W.1:p.Arg2K")


def test_the_error_carries_the_description() -> None:
    with pytest.raises(weaver.MixedResidueCodesError) as raised:
        weaver.parse("NP_W.1:p.Arg2K")
    assert raised.value.description == "NP_W.1:p.Arg2K"


def test_the_error_carries_each_one_letter_token_where_written() -> None:
    with pytest.raises(weaver.MixedResidueCodesError) as raised:
        weaver.parse("NP_000081.1:p.Ala356_Gly357insDPA")
    assert raised.value.one_letter == [("DPA", 30, 33)]


def test_the_error_carries_each_three_letter_token_where_written() -> None:
    with pytest.raises(weaver.MixedResidueCodesError) as raised:
        weaver.parse("NP_000081.1:p.Ala356_Gly357insDPA")
    assert raised.value.three_letter == [("Ala", 14, 17), ("Gly", 21, 24)]


def test_the_message_names_a_residue_written_in_each_code() -> None:
    with pytest.raises(
        weaver.MixedResidueCodesError,
        match=r"writes Arg in the three-letter .* K in the one-letter",
    ):
        weaver.parse("NP_W.1:p.Arg2K")
