"""An X that could be a stop or any residue raises an error, unless the caller says how to read it."""

from __future__ import annotations

import pickle

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


def test_a_three_letter_substitution_to_x_is_an_ambiguous_residue_error() -> None:
    with pytest.raises(weaver.AmbiguousResidueError):
        weaver.parse("NP_000021.1:p.Arg111X")


def test_parse_refuses_an_ambiguous_x_when_asked() -> None:
    with pytest.raises(weaver.AmbiguousResidueError):
        weaver.parse("NP_000026.2:p.R60X", x_treatment=weaver.XTreatment.Error)


def test_parse_reads_an_ambiguous_x_as_a_stop_when_asked() -> None:
    variant = weaver.parse("NP_000026.2:p.R60X", x_treatment=weaver.XTreatment.Stop)
    assert str(variant) == "NP_000026.2:p.Arg60Ter"


def test_parse_reads_an_ambiguous_x_as_any_residue_when_asked() -> None:
    variant = weaver.parse("NP_000026.2:p.R60X", x_treatment=weaver.XTreatment.Unknown)
    assert str(variant) == "NP_000026.2:p.Arg60Xaa"


def test_an_x_treatment_round_trips_through_pickle() -> None:
    # A worker process receives the treatment it parses with pickled.
    restored = pickle.loads(pickle.dumps(weaver.XTreatment.Stop))  # noqa: S301 -- the bytes are our own
    assert restored == weaver.XTreatment.Stop
