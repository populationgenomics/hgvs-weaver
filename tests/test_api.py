"""Tests for the weaver Python API."""

from __future__ import annotations

import typing

import pytest

import weaver


class MockProvider:
    """Mock DataProvider for testing."""

    def get_transcript(self, ac: str, _ref: str | None) -> dict[str, typing.Any]:
        """Returns a mock transcript model."""
        return {
            "ac": ac,
            "gene": "TEST",
            "cds_start_index": 10,
            "cds_end_index": 20,
            "strand": 1,
            "reference_accession": "NC_TEST.1",
            "reference_alignment_method": "splign",
            "exons": [
                {
                    "transcript_start": 0,
                    "transcript_end": 20,
                    "reference_start": 1000,
                    "reference_end": 1019,
                    "alt_strand": 1,
                    "cigar": "20=",
                },
                {
                    "transcript_start": 20,
                    "transcript_end": 100,
                    "reference_start": 1030,
                    "reference_end": 1109,
                    "alt_strand": 1,
                    "cigar": "80=",
                },
            ],
        }

    def get_seq(self, ac: str, start: int, end: int | None, _kind: str) -> str:
        """Returns a mock sequence: the transcript, and the genome carrying its first exon at index 1000."""
        # Index 10 is c.1: ATG GGG CCC AAA ...
        seq = "A" * 10 + "ATGGGGCCCAAA" + "A" * 100
        if ac.startswith("NC_"):
            seq = "N" * 1000 + seq
        return seq[start:end]

    def get_symbol_accessions(self, symbol: str, _s: str, t: str) -> list[tuple[weaver.IdentifierType, str]]:
        """Looks up mock symbols."""
        if t == "p":
            return [(weaver.IdentifierType.ProteinAccession, "NP_TEST.1")]
        return [(weaver.IdentifierType.GeneSymbol, symbol)]

    def get_identifier_type(self, identifier: str) -> weaver.IdentifierType:
        """Genomic for NC_, protein for NP_, transcript otherwise."""
        if identifier.startswith("NC_"):
            return weaver.IdentifierType.GenomicAccession
        if identifier.startswith("NP_"):
            return weaver.IdentifierType.ProteinAccession
        return weaver.IdentifierType.TranscriptAccession


def test_parsing() -> None:
    """Tests description parsing."""
    v = weaver.parse("NM_0001.1:c.123A>G")
    assert v.ac == "NM_0001.1"
    assert v.coordinate_type == "c"
    assert v.format() == "NM_0001.1:c.123A>G"


def test_normalization() -> None:
    """Tests description normalisation (3' shift)."""
    provider = MockProvider()
    mapper = weaver.VariantMapper(provider)

    # c.4_5del (GG) in GGG
    v = weaver.parse("NM_TEST:c.4_5del")
    v_norm = mapper.normalize_variant(v)
    # Should shift to 5_6del
    assert v_norm.format() == "NM_TEST:c.5_6del"


def test_c_to_p() -> None:
    """Tests c. to p. projection for a start-codon substitution."""
    provider = MockProvider()
    mapper = weaver.VariantMapper(provider)

    # c.1A>G changes the start codon ATG: what translation does is not predicted.
    v_c = weaver.parse("NM_TEST:c.1A>G")
    v_p = mapper.c_to_p(v_c)
    assert v_p.format().endswith(":p.(Met1?)")


def test_c_to_p_missense() -> None:
    """Tests that a missense substitution at a non-start codon produces a amino-acid change."""
    provider = MockProvider()
    mapper = weaver.VariantMapper(provider)

    # Mock CDS: ATG(Met1) GGG(Gly2) CCC(Pro3)...
    # c.4G>C changes codon 2 GGG(Gly) → CGG(Arg): p.(Gly2Arg)
    v_c = weaver.parse("NM_TEST:c.4G>C")
    v_p = mapper.c_to_p(v_c)
    assert "p.(Gly2Arg)" in v_p.format()


def test_c_to_g() -> None:
    """Tests c. to g. projection."""
    provider = MockProvider()
    mapper = weaver.VariantMapper(provider)

    # c.1A>G. cds_start=10. Index 10.
    # reference_start=1000 (index 0).
    # Index 10 -> g.1011
    v_c = weaver.parse("NM_TEST:c.1A>G")
    v_g = mapper.c_to_g(v_c, "NC_TEST.1")
    assert v_g.format() == "NC_TEST.1:g.1011A>G"


class MockMinusProvider:
    """Mock DataProvider for minus strand testing."""

    def get_transcript(self, ac: str, _ref: str | None) -> dict[str, typing.Any]:
        """Returns a mock transcript model on the minus strand."""
        return {
            "ac": ac,
            "gene": "TEST",
            "cds_start_index": 10,
            "cds_end_index": 20,
            "strand": -1,
            "reference_accession": "NC_TEST.1",
            "reference_alignment_method": "splign",
            "exons": [
                {
                    "transcript_start": 0,
                    "transcript_end": 100,
                    "reference_start": 1000,
                    "reference_end": 1100,
                    "alt_strand": -1,
                    "cigar": "100=",
                },
            ],
        }

    def get_seq(self, ac: str, start: int, end: int | None, _kind: str) -> str:
        """Returns a mock sequence: an all-A transcript, read from an all-T genome on the minus strand."""
        return (("T" if ac.startswith("NC_") else "A") * 1200)[start:end]

    def get_symbol_accessions(self, symbol: str, _s: str, _t: str) -> list[tuple[weaver.IdentifierType, str]]:
        """Looks up mock symbols."""
        return [(weaver.IdentifierType.GeneSymbol, symbol)]

    def get_identifier_type(self, identifier: str) -> weaver.IdentifierType:
        """Genomic for NC_, transcript otherwise."""
        if identifier.startswith("NC_"):
            return weaver.IdentifierType.GenomicAccession
        return weaver.IdentifierType.TranscriptAccession


def test_minus_strand() -> None:
    """Tests projection on the minus strand."""
    provider = MockMinusProvider()
    mapper = weaver.VariantMapper(provider)

    v_c = weaver.parse("NM_TEST:c.1A>G")
    v_g = mapper.c_to_g(v_c, "NC_TEST.1")
    # Coordinate 1090 is HgvsGenomicPos(1091).
    assert v_g.format() == "NC_TEST.1:g.1091T>C"


def test_intronic() -> None:
    """Tests projection of intronic descriptions."""
    provider = MockProvider()
    mapper = weaver.VariantMapper(provider)

    # c.10+1A>G.
    v_c = weaver.parse("NM_TEST:c.10+1A>G")
    v_g = mapper.c_to_g(v_c, "NC_TEST.1")
    assert v_g.format() == "NC_TEST.1:g.1021A>G"


def test_an_offset_from_a_base_inside_an_exon_raises_position_outside_sequence_error() -> None:
    """c.10 is the last base of exon 1, so c.9 has no base after it in the intron."""
    mapper = weaver.VariantMapper(MockProvider())
    with pytest.raises(weaver.PositionOutsideSequenceError, match=r"NM_TEST:c\.9\+1 "):
        mapper.c_to_g(weaver.parse("NM_TEST:c.9+1del"), "NC_TEST.1")


def test_validate_answers_false_for_an_offset_from_a_base_inside_an_exon() -> None:
    assert not weaver.parse("NM_TEST:c.9+1del").validate(MockProvider())


def test_validate_refuses_a_protein_description_on_a_transcript() -> None:
    """Base 11 of NM_TEST is A, which read as a residue would be Ala11."""
    with pytest.raises(weaver.ValidationError, match="NM_TEST"):
        weaver.parse("NM_TEST:p.Ala11Pro").validate(MockProvider())
