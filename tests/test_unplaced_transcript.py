"""A provider serves a transcript record with no genome placement by leaving the placement keys out."""

from __future__ import annotations

import pytest

import weaver

RECORD = "ATGAAAGGCCTGTAACC"  # Met Lys Gly Leu Ter


class Provider:
    """One transcript, TX_U.1, that no assembly places."""

    def get_transcript(self, transcript_ac: str, _reference_ac: str | None) -> weaver.TranscriptData:
        """Returns the model with no strand, reference_accession or exons."""
        return {"ac": transcript_ac, "gene": "UNPLACED", "cds_start_index": 0, "cds_end_index": 14}

    def get_seq(self, _ac: str, start: int, end: int | None, _kind: str) -> str:
        """Returns a slice of the record."""
        return RECORD[max(start, 0) : end]

    def get_identifier_type(self, _identifier: str) -> weaver.IdentifierType:
        """Every accession here is the transcript."""
        return weaver.IdentifierType.TranscriptAccession


def test_a_record_without_placement_keys_is_read_on_its_bases() -> None:
    assert weaver.parse("TX_U.1:c.7G>T").validate(Provider())


def test_projecting_an_unplaced_record_is_unsupported() -> None:
    mapper = weaver.VariantMapper(Provider())

    with pytest.raises(weaver.UnsupportedOperationError, match=r"TX_U\.1 has no genome placement"):
        mapper.c_to_g(weaver.parse("TX_U.1:c.7G>T"))
