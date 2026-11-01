"""A protein description that mixes the amino acid codes raises an error a caller can recover from."""

from __future__ import annotations

import pytest

import weaver


@pytest.mark.parametrize(
    "hgvs",
    ["NP_000021.1:p.Arg111X", "NP_000081.1:p.Ala356_Gly357insDPA", "NP_W.1:p.[Arg2Lys;L3P]"],
)
def test_a_description_mixing_the_codes_is_a_mixed_residue_codes_error(hgvs: str) -> None:
    with pytest.raises(weaver.MixedResidueCodesError):
        weaver.parse(hgvs)


def test_a_mixed_residue_codes_error_is_a_parse_error() -> None:
    with pytest.raises(weaver.ParseError):
        weaver.parse("NP_000021.1:p.Arg111X")


def test_a_mixed_residue_codes_error_is_an_hgvs_error() -> None:
    with pytest.raises(weaver.HGVSError):
        weaver.parse("NP_000021.1:p.Arg111X")


def test_the_error_carries_the_description() -> None:
    with pytest.raises(weaver.MixedResidueCodesError) as raised:
        weaver.parse("NP_000021.1:p.Arg111X")
    assert raised.value.description == "NP_000021.1:p.Arg111X"


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
        match=r"writes Arg in the three-letter .* X in the one-letter",
    ):
        weaver.parse("NP_000021.1:p.Arg111X")


def test_a_caller_can_rewrite_a_lone_x_as_ter_and_parse_again() -> None:
    with pytest.raises(weaver.MixedResidueCodesError) as raised:
        weaver.parse("NP_000021.1:p.Arg111X")
    ((text, start, end),) = raised.value.one_letter
    assert text == "X"
    description = raised.value.description
    assert str(weaver.parse(description[:start] + "Ter" + description[end:])) == "NP_000021.1:p.Arg111Ter"
