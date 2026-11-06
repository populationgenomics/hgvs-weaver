"""weaver's exceptions survive pickling, as a worker process's error must to reach its parent."""

from __future__ import annotations

import concurrent.futures
import copy
import pickle
from collections.abc import Callable

import pytest

import weaver
from weaver import _weaver

EXCEPTION_CLASSES = sorted(
    (value for value in vars(_weaver).values() if isinstance(value, type) and issubclass(value, BaseException)),
    key=lambda cls: cls.__name__,
)
_MISMATCH = _weaver.ReferenceMismatchError("expected C, found T", "NM_X.1", "C", "T", 9, 10)

# Each class that carries attributes, made as the core makes it, with the attributes it carries.
WITH_ATTRIBUTES: dict[str, tuple[tuple[object, ...], set[str]]] = {
    "AmbiguousDescriptionError": (
        ("m", "NM_X.1:c.[7_8insA;7_8insG]", [weaver.parse("NM_X.1:c.7_8insAG")]),
        {"description", "readings"},
    ),
    "AmbiguousResidueError": (("m", "NP_W.1:p.R2X", [("X", 10, 11)]), {"description", "ambiguous"}),
    "CisMemberError": (
        ("m", "NM_X.1:c.[1A>G;10C>T]", 1, "NM_X.1:c.10C>T", _MISMATCH),
        {"allele", "index", "member", "source"},
    ),
    "MixedResidueCodesError": (
        ("m", "NP_W.1:p.Arg2K", [("K", 13, 14)], [("Arg", 9, 12)]),
        {"description", "one_letter", "three_letter"},
    ),
    "NotAnExonBorderError": (("m", "NM_X.1", "c.4+1"), {"accession", "position"}),
    "OffsetPastIntronError": (("m", "NM_X.1", "c.5+6", 5), {"accession", "position", "intron_length"}),
    "PositionOutsideSequenceError": (("m", "NM_X.1"), {"accession"}),
    "PositionPastEndError": (("m", "NM_X.1", 16, 17, 12), {"accession", "start", "end", "length"}),
    "ReferenceMismatchError": (
        ("m", "NM_X.1", "C", "T", 9, 10),
        {"accession", "expected", "found", "start", "end"},
    ),
    "UnplacedError": (("m", "NM_U.1"), {"accession"}),
    "WrongReferenceKindError": (
        (
            "m",
            "NM_X.1",
            weaver.ReferenceUse.Description(coordinate_type="p"),
            weaver.IdentifierType.TranscriptAccession,
        ),
        {"accession", "read_as", "kind"},
    ),
}
MESSAGE_ONLY = [cls for cls in EXCEPTION_CLASSES if cls.__name__ not in WITH_ATTRIBUTES]


def _comparable(value: object) -> object:
    """An attribute's value with each Variant, which has no equality, as its dict, and each error as its args."""
    if isinstance(value, list):
        return [_comparable(item) for item in value]
    if isinstance(value, weaver.Variant):
        return value.to_dict()
    if isinstance(value, BaseException):
        return type(value), _comparable(list(value.args))
    return value


def _through_pickle(error: BaseException) -> BaseException:
    return pickle.loads(pickle.dumps(error))  # noqa: S301 -- the bytes are our own


def _parse_in_worker(description: str) -> weaver.Variant:
    return weaver.parse(description)


def _raise_in_worker(name: str) -> None:
    raise getattr(_weaver, name)(*WITH_ATTRIBUTES[name][0])


def test_the_extension_declares_exception_classes() -> None:
    assert {cls.__name__ for cls in EXCEPTION_CLASSES} >= {"HGVSError", "ParseError", "AmbiguousDescriptionError"}


@pytest.mark.parametrize("cls", MESSAGE_ONLY, ids=lambda cls: cls.__name__)
def test_every_exception_class_round_trips_through_pickle(cls: type[BaseException]) -> None:
    restored = _through_pickle(cls("the message"))

    assert type(restored) is cls
    assert str(restored) == "the message"


@pytest.mark.parametrize("name", sorted(WITH_ATTRIBUTES))
def test_an_error_with_attributes_is_not_made_from_a_message_alone(name: str) -> None:
    with pytest.raises(TypeError, match="missing"):
        getattr(_weaver, name)("the message")


@pytest.mark.parametrize("name", sorted(WITH_ATTRIBUTES))
@pytest.mark.parametrize(
    "duplicate",
    [_through_pickle, copy.copy, copy.deepcopy],
    ids=["pickle", "copy", "deepcopy"],
)
def test_an_error_keeps_its_attributes_through_a_copy(
    name: str,
    duplicate: Callable[[BaseException], BaseException],
) -> None:
    args, attributes = WITH_ATTRIBUTES[name]
    original = getattr(_weaver, name)(*args)

    restored = duplicate(original)

    assert type(restored) is type(original)
    assert str(restored) == "m"
    assert {attribute: _comparable(getattr(restored, attribute)) for attribute in attributes} == {
        attribute: _comparable(getattr(original, attribute)) for attribute in attributes
    }


def test_an_error_raised_in_a_worker_process_reaches_the_parent() -> None:
    with concurrent.futures.ProcessPoolExecutor(max_workers=1) as pool:
        future = pool.submit(_parse_in_worker, "NP_W.1:p.Arg2K")
        with pytest.raises(weaver.MixedResidueCodesError) as raised:
            future.result()

    assert isinstance(raised.value, weaver.ParseError)
    assert raised.value.description == "NP_W.1:p.Arg2K"
    assert raised.value.three_letter == [("Arg", 9, 12)]
    assert raised.value.one_letter == [("K", 13, 14)]


def test_an_error_holding_a_members_error_reaches_the_parent_with_it() -> None:
    with concurrent.futures.ProcessPoolExecutor(max_workers=1) as pool:
        future = pool.submit(_raise_in_worker, "CisMemberError")
        with pytest.raises(weaver.CisMemberError) as raised:
            future.result()

    source = raised.value.source
    assert isinstance(source, weaver.ReferenceMismatchError)
    assert (source.accession, source.expected, source.found, source.start, source.end) == ("NM_X.1", "C", "T", 9, 10)
