"""The binding of stated-base checks: validate's answer, and the error a mismatch raises."""

from __future__ import annotations

import typing

import pytest

import weaver

TRANSCRIPT = "GCTAGCTAGCATGGCTGGATCCAAGTTCCTGCACGATGAAGTCATCAACCCTCGATACTGGAGCACATAAACGTTGCAAGGTCCATGACC"


class Provider:
    """A single-exon transcript with the CDS at indices 10..69; c.4_5 is GC."""

    def get_transcript(self, ac: str, _ref: str | None) -> dict[str, typing.Any]:
        return {
            "ac": ac,
            "gene": "TEST",
            "cds_start_index": 10,
            "cds_end_index": 69,
            "strand": 1,
            "reference_accession": "NC_TEST.1",
            "reference_alignment_method": "splign",
            "exons": [
                {
                    "transcript_start": 0,
                    "transcript_end": 90,
                    "reference_start": 1000,
                    "reference_end": 1089,
                    "alt_strand": 1,
                    "cigar": "90=",
                },
            ],
        }

    def get_seq(self, ac: str, start: int, end: int | None, _kind: str) -> str:
        seq = "N" * 1000 + TRANSCRIPT if ac.startswith("NC_") else TRANSCRIPT
        return seq[start:end]

    def get_symbol_accessions(self, symbol: str, _s: str, _t: str) -> list[typing.Any]:
        return [(weaver.IdentifierType.GeneSymbol, symbol)]


@pytest.mark.parametrize(
    ("hgvs", "expected"),
    [
        ("NM_TEST:c.4G>A", True),
        ("NM_TEST:c.4T>A", False),
        ("NM_TEST:c.4delG", True),
        ("NM_TEST:c.4delC", False),
        ("NM_TEST:c.4_5delGC", True),
        ("NM_TEST:c.4_5delAA", False),
        ("NM_TEST:c.4dupG", True),
        ("NM_TEST:c.4dupT", False),
        ("NM_TEST:c.4_5delGCinsTT", True),
        ("NM_TEST:c.4_5delAAinsTT", False),
        ("NM_TEST:c.4_5del2insTT", True),
        ("NM_TEST:c.4del", True),
    ],
)
def test_bases_stated_after_del_dup_and_delins_are_checked(hgvs: str, expected: bool) -> None:
    assert weaver.parse(hgvs).validate(Provider()) is expected


INTRON = "GTAAGTCCAG"


class SplicedProvider(Provider):
    """TRANSCRIPT with INTRON after index 50: c.40 ends exon 1, and c.40+1 is the intron's G."""

    def get_transcript(self, ac: str, _ref: str | None) -> dict[str, typing.Any]:
        model = super().get_transcript(ac, _ref)
        model["exons"] = [
            {
                "transcript_start": 0,
                "transcript_end": 50,
                "reference_start": 1000,
                "reference_end": 1049,
                "alt_strand": 1,
                "cigar": "50=",
            },
            {
                "transcript_start": 50,
                "transcript_end": 90,
                "reference_start": 1060,
                "reference_end": 1099,
                "alt_strand": 1,
                "cigar": "40=",
            },
        ]
        return model

    def get_seq(self, ac: str, start: int, end: int | None, _kind: str) -> str:
        if ac.startswith("NC_"):
            return ("N" * 1000 + TRANSCRIPT[:50] + INTRON + TRANSCRIPT[50:])[start:end]
        return TRANSCRIPT[start:end]


def test_a_base_the_genome_lacks_is_an_error_naming_the_genome() -> None:
    mapper = weaver.VariantMapper(SplicedProvider())
    with pytest.raises(weaver.ReferenceMismatchError, match=r"expected C, found G on NC_TEST\.1"):
        mapper.c_to_g(weaver.parse("NM_TEST:c.40+1C>A"))
