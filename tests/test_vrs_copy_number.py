"""Copy-number edits rendered as VRS CopyNumberCount dicts and read back."""

import json
import typing

import pytest

import weaver

GENOME = "ACGTTTGCAAGGCTAGCTAGCTTTTAACGGGATCGATCGA"


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


def test_to_vrs_returns_a_copy_number_count(mapper: weaver.VariantMapper) -> None:
    count = mapper.to_vrs(weaver.parse("NC_TEST.1:g.5_12copy3"))
    assert count["type"] == "CopyNumberCount"
    assert count["id"].startswith("ga4gh:CN.")
    assert count["copies"] == 3
    assert count["location"]["start"] == 4
    assert count["location"]["end"] == 12
    assert count["expressions"] == [{"syntax": "hgvs.g", "value": "NC_TEST.1:g.5_12copy3"}]
    assert "state" not in count
    assert mapper.vrs_id(weaver.parse("NC_TEST.1:g.5_12copy3")) == count["id"]


@pytest.mark.parametrize(
    ("hgvs", "expected"),
    [
        ("NC_TEST.1:g.5_12copy3", "NC_TEST.1:g.5_12copy3"),
        ("NC_TEST.1:g.5copy2", "NC_TEST.1:g.5copy2"),
        ("NC_TEST.1:g.(3_5)_(10_12)copy4", "NC_TEST.1:g.(3_5)_(10_12)copy4"),
        ("NC_TEST.1:m.5_12copy3", "NC_TEST.1:g.5_12copy3"),
    ],
)
def test_copy_number_count_round_trip(mapper: weaver.VariantMapper, hgvs: str, expected: str) -> None:
    count = mapper.to_vrs(weaver.parse(hgvs))
    back = mapper.from_vrs(count)
    assert str(back) == expected
    assert mapper.vrs_id(back) == count["id"]
    assert str(mapper.from_vrs(json.dumps(count), "NC_TEST.1")) == expected


def test_alleles_are_still_alleles(mapper: weaver.VariantMapper) -> None:
    allele = mapper.to_vrs(weaver.parse("NC_TEST.1:g.7G>C"))
    assert allele["type"] == "Allele"
    assert mapper.vrs_id(weaver.parse("NC_TEST.1:g.7G>C")).startswith("ga4gh:VA.")


def test_a_range_of_copies_has_no_hgvs(mapper: weaver.VariantMapper) -> None:
    count = mapper.to_vrs(weaver.parse("NC_TEST.1:g.5_12copy3"))
    count["copies"] = [3, None]
    with pytest.raises(weaver.UnsupportedOperationError):
        mapper.from_vrs(count)
