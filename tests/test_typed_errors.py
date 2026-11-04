"""Each error the core returns is raised as its own class, carrying the core's fields, under the class it had before."""

from __future__ import annotations

import typing

import pytest

import weaver

TRANSCRIPT = "ATGAAGGCCTAA"  # Met Lys Ala Ter; c.1 to c.5 in exon 1, c.6 to c.*0 in exon 2
INTRON = "GTCAG"
GENOME = "C" * 10 + TRANSCRIPT[:5] + INTRON + TRANSCRIPT[5:] + "C" * 10


class Provider:
    """NM_T.1, two exons on NC_T.1 at g.11 to g.15 and g.21 to g.27; NM_U.1, its sequence placed on no genome."""

    def get_transcript(self, ac: str, _ref: str | None) -> dict[str, typing.Any]:
        record = {"ac": ac, "gene": "TEST", "cds_start_index": 0, "cds_end_index": 11}
        if ac == "NM_U.1":
            return record
        return record | {
            "strand": 1,
            "reference_accession": "NC_T.1",
            "exons": [
                {
                    "transcript_start": 0,
                    "transcript_end": 5,
                    "reference_start": 10,
                    "reference_end": 14,
                    "alt_strand": 1,
                    "cigar": "5=",
                },
                {
                    "transcript_start": 5,
                    "transcript_end": 12,
                    "reference_start": 20,
                    "reference_end": 26,
                    "alt_strand": 1,
                    "cigar": "7=",
                },
            ],
        }

    def get_seq(self, ac: str, start: int, end: int | None, _kind: str) -> str:
        sequence = GENOME if ac == "NC_T.1" else TRANSCRIPT
        return sequence[max(start, 0) : end]

    def get_identifier_type(self, identifier: str) -> weaver.IdentifierType:
        if identifier.startswith("NC_"):
            return weaver.IdentifierType.GenomicAccession
        return weaver.IdentifierType.TranscriptAccession


@pytest.fixture
def mapper() -> weaver.VariantMapper:
    return weaver.VariantMapper(Provider())


def test_a_stated_base_the_record_does_not_hold_carries_both(mapper: weaver.VariantMapper) -> None:
    with pytest.raises(weaver.ReferenceMismatchError) as raised:
        mapper.c_to_g(weaver.parse("NM_T.1:c.4C>T"))

    error = raised.value
    assert (error.accession, error.expected, error.found, error.start, error.end) == ("NM_T.1", "C", "A", 3, 4)


def test_a_range_past_the_records_end_carries_the_range_and_the_length(mapper: weaver.VariantMapper) -> None:
    with pytest.raises(weaver.PositionOutsideSequenceError) as raised:
        mapper.c_to_g(weaver.parse("NM_T.1:c.*5del"))

    error = raised.value
    assert isinstance(error, weaver.PositionPastEndError)
    assert (error.accession, error.start, error.end, error.length) == ("NM_T.1", 16, 17, 12)


def test_an_offset_from_inside_an_exon_carries_the_position_as_written(mapper: weaver.VariantMapper) -> None:
    with pytest.raises(weaver.PositionOutsideSequenceError) as raised:
        mapper.c_to_g(weaver.parse("NM_T.1:c.4+1del"))

    error = raised.value
    assert isinstance(error, weaver.NotAnExonBorderError)
    assert (error.accession, error.position) == ("NM_T.1", "c.4+1")


def test_an_offset_past_its_intron_carries_the_introns_length(mapper: weaver.VariantMapper) -> None:
    with pytest.raises(weaver.PositionOutsideSequenceError) as raised:
        mapper.c_to_g(weaver.parse("NM_T.1:c.5+6del"))

    error = raised.value
    assert isinstance(error, weaver.OffsetPastIntronError)
    assert (error.accession, error.position, error.intron_length) == ("NM_T.1", "c.5+6", 5)


def test_a_projection_of_an_unplaced_record_carries_its_accession(mapper: weaver.VariantMapper) -> None:
    with pytest.raises(weaver.UnsupportedOperationError) as raised:
        mapper.c_to_g(weaver.parse("NM_U.1:c.4A>G"))

    error = raised.value
    assert isinstance(error, weaver.UnplacedError)
    assert error.accession == "NM_U.1"


def test_a_protein_description_on_a_transcript_carries_both_kinds() -> None:
    with pytest.raises(weaver.ValidationError) as raised:
        weaver.parse("NM_T.1:p.Met1Val").validate(Provider())

    error = raised.value
    assert isinstance(error, weaver.WrongReferenceKindError)
    assert error.accession == "NM_T.1"
    assert error.read_as == weaver.ReferenceUse.Description(coordinate_type="p")
    assert error.kind == weaver.IdentifierType.TranscriptAccession


def test_a_failing_member_of_an_allele_is_named_by_its_place(mapper: weaver.VariantMapper) -> None:
    with pytest.raises(weaver.CisMemberError) as raised:
        mapper.c_to_g(weaver.parse("NM_T.1:c.[1A>G;10C>T]"))

    error = raised.value
    assert (error.allele, error.index, error.member) == ("NM_T.1:c.[1A>G;10C>T]", 1, "NM_T.1:c.10C>T")


def test_a_failing_member_of_an_allele_holds_the_members_own_error(mapper: weaver.VariantMapper) -> None:
    with pytest.raises(weaver.CisMemberError) as raised:
        mapper.c_to_g(weaver.parse("NM_T.1:c.[1A>G;10C>T]"))

    source = raised.value.source
    assert isinstance(source, weaver.ReferenceMismatchError)
    assert (source.expected, source.found, source.start, source.end) == ("C", "T", 9, 10)
    assert str(raised.value).endswith(str(source))


def test_an_ambiguous_allele_carries_the_description_read() -> None:
    with pytest.raises(weaver.AmbiguousDescriptionError) as raised:
        weaver.parse("NM_T.1:c.[7_8insA;7_8insG]").validate(Provider())

    assert raised.value.description == "NM_T.1:c.[7_8insA;7_8insG]"
