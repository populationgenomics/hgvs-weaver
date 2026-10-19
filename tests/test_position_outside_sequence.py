"""The binding of a position past the end of its reference: the error it raises."""

from __future__ import annotations

import typing

import pytest

import weaver


class AmpliconProvider:
    """One 9-base amplicon, `AMP1`, served with the bases that exist over any window."""

    def get_transcript(self, ac: str, _ref: str | None) -> dict[str, typing.Any]:
        raise weaver.DataProviderError(f"no transcript {ac}")

    def get_seq(self, ac: str, start: int, end: int | None, _kind: str) -> str:
        if ac != "AMP1":
            raise weaver.DataProviderError(f"no sequence {ac}")
        return "TTGGGGCTT"[start:end]

    def get_symbol_accessions(self, _symbol: str, _s: str, _t: str) -> list[typing.Any]:
        return []

    def get_identifier_type(self, _identifier: str) -> weaver.IdentifierType:
        return weaver.IdentifierType.GenomicAccession


def test_a_deletion_past_the_end_raises_an_error_naming_the_accession_and_its_length() -> None:
    mapper = weaver.VariantMapper(AmpliconProvider())
    with pytest.raises(weaver.PositionOutsideSequenceError, match=r"AMP1 is 9 long"):
        mapper.to_spdi(weaver.parse("AMP1:g.8_12del"))
