"""Insertions of a stated length, insN[20], as VRS Alleles with a LengthExpression state."""

import json
import typing

import pytest

import weaver

#         1234567890123456789012345678901234567890
GENOME = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"


class Provider:
    """One genomic sequence, with the refget lookups implemented."""

    def get_transcript(self, ac: str, _ref: str | None) -> dict[str, typing.Any]:
        raise weaver.DataProviderError(f"no transcript {ac}")

    def get_seq(self, _ac: str, start: int, end: int | None, _kind: str) -> str:
        return GENOME[start:end]

    def get_symbol_accessions(self, _symbol: str, _s: str, _t: str) -> list[tuple[weaver.IdentifierType, str]]:
        return []

    def get_identifier_type(self, _identifier: str) -> weaver.IdentifierType:
        return weaver.IdentifierType.GenomicAccession

    def get_refget_accession(self, _ac: str) -> str | None:
        return None  # let weaver compute it

    def get_accession_for_refget(self, _refget: str) -> str | None:
        return "NC_TEST.1"


@pytest.fixture
def mapper() -> weaver.VariantMapper:
    provider = Provider()
    return weaver.VariantMapper(provider, refget=provider)


def test_to_vrs_gives_a_length_expression_at_the_insertion_point(mapper: weaver.VariantMapper) -> None:
    allele = mapper.to_vrs(weaver.parse("NC_TEST.1:g.10_11insN[20]"))
    assert allele["type"] == "Allele"
    assert allele["id"].startswith("ga4gh:VA.")
    assert allele["state"] == {"type": "LengthExpression", "length": 20}
    assert allele["location"]["start"] == 10
    assert allele["location"]["end"] == 10
    assert allele["expressions"] == [{"syntax": "hgvs.g", "value": "NC_TEST.1:g.10_11insN[20]"}]
    assert mapper.vrs_id(weaver.parse("NC_TEST.1:g.10_11insN[20]")) == allele["id"]


def test_both_spellings_give_the_same_identifier(mapper: weaver.VariantMapper) -> None:
    recommended = mapper.vrs_id(weaver.parse("NC_TEST.1:g.10_11insN[20]"))
    assert mapper.vrs_id(weaver.parse("NC_TEST.1:g.10_11ins(20)")) == recommended
    assert str(weaver.parse("NC_TEST.1:g.10_11ins(20)")) == "NC_TEST.1:g.10_11insN[20]"


def test_an_uncertain_length_is_a_range(mapper: weaver.VariantMapper) -> None:
    allele = mapper.to_vrs(weaver.parse("NC_TEST.1:g.10_11insN[(20_30)]"))
    assert allele["state"] == {"type": "LengthExpression", "length": [20, 30]}
    assert mapper.vrs_id(weaver.parse("NC_TEST.1:g.10_11ins(20_30)")) == allele["id"]


@pytest.mark.parametrize(
    ("hgvs", "expected"),
    [
        ("NC_TEST.1:g.10_11insN[20]", "NC_TEST.1:g.10_11insN[20]"),
        ("NC_TEST.1:g.10_11ins(20)", "NC_TEST.1:g.10_11insN[20]"),
        ("NC_TEST.1:g.10_11insN[(20_30)]", "NC_TEST.1:g.10_11insN[(20_30)]"),
        ("NC_TEST.1:g.10_12delinsN[5]", "NC_TEST.1:g.10_12delinsN[5]"),
        ("NC_TEST.1:m.10_11insN[20]", "NC_TEST.1:g.10_11insN[20]"),
    ],
)
def test_length_expressions_round_trip(mapper: weaver.VariantMapper, hgvs: str, expected: str) -> None:
    allele = mapper.to_vrs(weaver.parse(hgvs))
    back = mapper.from_vrs(allele)
    assert str(back) == expected
    assert mapper.vrs_id(back) == allele["id"]
    assert str(mapper.from_vrs(json.dumps(allele), "NC_TEST.1")) == expected


def test_a_delins_of_a_stated_length_covers_the_deleted_range(mapper: weaver.VariantMapper) -> None:
    allele = mapper.to_vrs(weaver.parse("NC_TEST.1:g.10_12delinsN[5]"))
    assert allele["location"]["start"] == 9
    assert allele["location"]["end"] == 12
    assert allele["state"] == {"type": "LengthExpression", "length": 5}


def test_an_open_range_of_lengths_has_no_hgvs(mapper: weaver.VariantMapper) -> None:
    allele = mapper.to_vrs(weaver.parse("NC_TEST.1:g.10_11insN[20]"))
    allele["state"]["length"] = [20, None]
    with pytest.raises(weaver.UnsupportedOperationError):
        mapper.from_vrs(allele)
