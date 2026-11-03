"""weaver's exceptions survive pickling, as a worker process's error must to reach its parent."""

from __future__ import annotations

import concurrent.futures
import copy
import pickle
import typing
from collections.abc import Callable

import pytest

import weaver
from weaver import _weaver

EXCEPTION_CLASSES = sorted(
    (value for value in vars(_weaver).values() if isinstance(value, type) and issubclass(value, BaseException)),
    key=lambda cls: cls.__name__,
)
# These carry attributes, so a message alone does not make one.
WITH_ATTRIBUTES = {"AmbiguousDescriptionError", "AmbiguousResidueError", "MixedResidueCodesError"}
MESSAGE_ONLY = [cls for cls in EXCEPTION_CLASSES if cls.__name__ not in WITH_ATTRIBUTES]


def _raised(call: Callable[[], object]) -> weaver.HGVSError:
    with pytest.raises(weaver.HGVSError) as raised:
        call()
    return raised.value


def _mixed_residue_codes() -> weaver.HGVSError:
    return _raised(lambda: weaver.parse("NP_000021.1:p.Arg111X"))


def _ambiguous_residue() -> weaver.HGVSError:
    return _raised(lambda: weaver.parse("NP_W.1:p.[R2X;A4_G5insKXX]"))


TRANSCRIPT = "GGGGGATGAAACTGGCCTATCGCTAACCGTATAAGTAAGG"


class Provider:
    """One single-exon transcript, NM_X.1, laid on its genome NC_X.1 at index 10."""

    def get_transcript(self, ac: str, _ref: str | None) -> dict[str, typing.Any]:
        return {
            "ac": ac,
            "gene": "TEST",
            "cds_start_index": 5,
            "cds_end_index": 25,
            "strand": 1,
            "reference_accession": "NC_X.1",
            "exons": [
                {
                    "transcript_start": 0,
                    "transcript_end": len(TRANSCRIPT),
                    "reference_start": 10,
                    "reference_end": 10 + len(TRANSCRIPT) - 1,
                    "alt_strand": 1,
                    "cigar": f"{len(TRANSCRIPT)}=",
                },
            ],
        }

    def get_seq(self, ac: str, start: int, end: int | None, _kind: str) -> str:
        sequences = {"NM_X.1": TRANSCRIPT, "NC_X.1": "T" * 10 + TRANSCRIPT + "C" * 10}
        return sequences[ac][start:end]

    def get_identifier_type(self, identifier: str) -> weaver.IdentifierType:
        if identifier.startswith("NC_"):
            return weaver.IdentifierType.GenomicAccession
        return weaver.IdentifierType.TranscriptAccession


def _ambiguous_description() -> weaver.HGVSError:
    return _raised(lambda: weaver.parse("NM_X.1:c.[7_8insA;7_8insG]").validate(Provider()))


def _comparable(value: object) -> object:
    """An attribute's value with each Variant, which has no equality, as its dict."""
    if isinstance(value, list):
        return [_comparable(item) for item in value]
    if isinstance(value, weaver.Variant):
        return value.to_dict()
    return value


def _through_pickle(error: BaseException) -> BaseException:
    return pickle.loads(pickle.dumps(error))  # noqa: S301 -- the bytes are our own


def _parse_in_worker(description: str) -> weaver.Variant:
    return weaver.parse(description)


def test_the_extension_declares_exception_classes() -> None:
    assert {cls.__name__ for cls in EXCEPTION_CLASSES} >= {"HGVSError", "ParseError", "AmbiguousDescriptionError"}


@pytest.mark.parametrize("cls", MESSAGE_ONLY, ids=lambda cls: cls.__name__)
def test_every_exception_class_round_trips_through_pickle(cls: type[BaseException]) -> None:
    restored = _through_pickle(cls("the message"))

    assert type(restored) is cls
    assert str(restored) == "the message"


@pytest.mark.parametrize("name", sorted(WITH_ATTRIBUTES))
def test_an_error_with_attributes_is_not_made_from_a_message_alone(name: str) -> None:
    with pytest.raises(TypeError):
        getattr(_weaver, name)("the message")


@pytest.mark.parametrize(
    ("make", "attributes"),
    [
        (_mixed_residue_codes, {"description", "one_letter", "three_letter"}),
        (_ambiguous_residue, {"description", "ambiguous"}),
        (_ambiguous_description, {"readings"}),
    ],
    ids=["MixedResidueCodesError", "AmbiguousResidueError", "AmbiguousDescriptionError"],
)
@pytest.mark.parametrize(
    "duplicate",
    [_through_pickle, copy.copy, copy.deepcopy],
    ids=["pickle", "copy", "deepcopy"],
)
def test_an_error_keeps_its_attributes_through_a_copy(
    make: Callable[[], weaver.HGVSError],
    attributes: set[str],
    duplicate: Callable[[BaseException], BaseException],
) -> None:
    original = make()

    restored = duplicate(original)

    assert type(restored) is type(original)
    assert str(restored) == str(original)
    assert {name: _comparable(getattr(restored, name)) for name in attributes} == {
        name: _comparable(getattr(original, name)) for name in attributes
    }


def test_an_error_raised_in_a_worker_process_reaches_the_parent() -> None:
    with concurrent.futures.ProcessPoolExecutor(max_workers=1) as pool:
        future = pool.submit(_parse_in_worker, "NP_000021.1:p.Arg111X")
        with pytest.raises(weaver.MixedResidueCodesError) as raised:
            future.result()

    assert isinstance(raised.value, weaver.ParseError)
    assert raised.value.description == "NP_000021.1:p.Arg111X"
    assert raised.value.three_letter == [("Arg", 14, 17)]
    assert raised.value.one_letter == [("X", 20, 21)]
