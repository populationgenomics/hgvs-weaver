"""weaver's exceptions behave as Python exceptions do: constructed, subclassed, referenced and read."""

from __future__ import annotations

import weakref
from collections.abc import Callable

import pytest

import weaver


def _mixed_residue_codes() -> weaver.MixedResidueCodesError:
    with pytest.raises(weaver.MixedResidueCodesError) as raised:
        weaver.parse("NP_W.1:p.Arg2K")
    return raised.value


def test_a_base_error_takes_any_arguments_as_exception_does() -> None:
    assert weaver.HGVSError().args == ()
    assert weaver.ValidationError("a", 2).args == ("a", 2)


def test_a_python_subclass_of_a_base_error_shapes_its_own_constructor() -> None:
    class ProviderFailure(weaver.DataProviderError):
        def __init__(self, accession: str, attempts: int) -> None:
            super().__init__(f"{accession} failed after {attempts} attempts")
            self.attempts = attempts

    error = ProviderFailure("NM_X.1", 3)

    assert isinstance(error, weaver.HGVSError)
    assert error.attempts == 3


def test_an_error_with_attributes_says_its_message_alone() -> None:
    error = _mixed_residue_codes()

    assert str(error) == error.args[0]
    assert "NP_W.1:p.Arg2K" in str(error)


def test_an_errors_attributes_are_read_only() -> None:
    error = _mixed_residue_codes()

    with pytest.raises(AttributeError):
        error.description = "NP_X.1:p.Arg1Ter"  # type: ignore[misc]


@pytest.mark.parametrize(
    "make",
    [lambda: weaver.HGVSError("m"), lambda: weaver.ParseError("m"), _mixed_residue_codes],
    ids=["HGVSError", "ParseError", "MixedResidueCodesError"],
)
def test_an_error_can_be_weakly_referenced(make: Callable[[], BaseException]) -> None:
    error = make()

    assert weakref.ref(error)() is error


def test_from_vrs_refuses_an_object_that_is_neither_json_nor_a_dict() -> None:
    class Provider:
        def get_identifier_type(self, identifier: str) -> weaver.IdentifierType:  # noqa: ARG002
            return weaver.IdentifierType.GenomicAccession

    with pytest.raises(TypeError, match=r"str \| dict"):
        weaver.VariantMapper(Provider()).from_vrs(["not", "an", "allele"])  # type: ignore[arg-type]
