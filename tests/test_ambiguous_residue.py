"""A one-letter X that could be a stop or any residue raises an error a caller can recover from."""

from __future__ import annotations

import pytest

import weaver


def test_a_one_letter_substitution_to_x_is_an_ambiguous_residue_error() -> None:
    with pytest.raises(weaver.AmbiguousResidueError):
        weaver.parse("NP_000026.2:p.R60X")


def test_an_ambiguous_residue_error_is_a_parse_error() -> None:
    with pytest.raises(weaver.ParseError):
        weaver.parse("NP_000026.2:p.R60X")


def test_the_error_carries_the_description() -> None:
    with pytest.raises(weaver.AmbiguousResidueError) as raised:
        weaver.parse("NP_000026.2:p.R60X")
    assert raised.value.description == "NP_000026.2:p.R60X"


def test_the_error_carries_each_x_where_written() -> None:
    with pytest.raises(weaver.AmbiguousResidueError) as raised:
        weaver.parse("NP_W.1:p.[R2X;A4_G5insKXX]")
    assert raised.value.ambiguous == [("X", 12, 13), ("X", 23, 24), ("X", 24, 25)]


def test_a_caller_can_rewrite_an_x_as_a_stop_and_parse_again() -> None:
    with pytest.raises(weaver.AmbiguousResidueError) as raised:
        weaver.parse("NP_000026.2:p.R60X")
    ((_, start, end),) = raised.value.ambiguous
    description = raised.value.description
    assert str(weaver.parse(description[:start] + "*" + description[end:])) == "NP_000026.2:p.Arg60Ter"
