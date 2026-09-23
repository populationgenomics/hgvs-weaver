"""Issue #38: validate checks the bases stated after del, dup and delins."""

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
