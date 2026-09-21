"""Duplications with uncertain breakpoints as VRS CopyNumberChange dicts, and read back."""

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


def test_an_imprecise_duplication_is_a_copy_number_change(mapper: weaver.VariantMapper) -> None:
    change = mapper.to_vrs(weaver.parse("NC_TEST.1:g.(3_5)_(10_12)dup"))
    assert change["type"] == "CopyNumberChange"
    assert change["id"].startswith("ga4gh:CX.")
    assert change["copyChange"] == "gain"
    assert change["location"]["start"] == [2, 4]
    assert change["location"]["end"] == [10, 12]
    assert change["expressions"] == [{"syntax": "hgvs.g", "value": "NC_TEST.1:g.(3_5)_(10_12)dup"}]
    assert "state" not in change
    assert "copies" not in change
    assert mapper.vrs_id(weaver.parse("NC_TEST.1:g.(3_5)_(10_12)dup")) == change["id"]


def test_an_imprecise_deletion_is_still_an_allele(mapper: weaver.VariantMapper) -> None:
    allele = mapper.to_vrs(weaver.parse("NC_TEST.1:g.(3_5)_(10_12)del"))
    assert allele["type"] == "Allele"
    assert allele["state"] == {"type": "LiteralSequenceExpression", "sequence": ""}
    # And an exact duplication is a normalised Allele.
    assert mapper.to_vrs(weaver.parse("NC_TEST.1:g.5_12dup"))["type"] == "Allele"


@pytest.mark.parametrize(
    ("hgvs", "expected"),
    [
        ("NC_TEST.1:g.(3_5)_(10_12)dup", "NC_TEST.1:g.(3_5)_(10_12)dup"),
        ("NC_TEST.1:g.(?_5)_(10_?)dup", "NC_TEST.1:g.(?_5)_(10_?)dup"),
        ("NC_TEST.1:g.5_(10_12)dup", "NC_TEST.1:g.5_(10_12)dup"),
        ("NC_TEST.1:m.(3_5)_(10_12)dup", "NC_TEST.1:g.(3_5)_(10_12)dup"),
    ],
)
def test_copy_number_changes_round_trip(mapper: weaver.VariantMapper, hgvs: str, expected: str) -> None:
    change = mapper.to_vrs(weaver.parse(hgvs))
    back = mapper.from_vrs(change)
    assert str(back) == expected
    assert mapper.vrs_id(back) == change["id"]
    assert str(mapper.from_vrs(json.dumps(change), "NC_TEST.1")) == expected


@pytest.mark.parametrize(
    ("term", "expected"),
    [
        ("low-level gain", "NC_TEST.1:g.(3_5)_(10_12)dup"),
        ("EFO:0030072", "NC_TEST.1:g.(3_5)_(10_12)dup"),
        ("loss", "NC_TEST.1:g.(3_5)_(10_12)del"),
        ("complete genomic loss", "NC_TEST.1:g.(3_5)_(10_12)del"),
        ("EFO:0030068", "NC_TEST.1:g.(3_5)_(10_12)del"),
        (
            {"primaryCoding": {"code": "EFO:0030067", "system": "https://www.ebi.ac.uk/efo/"}},
            "NC_TEST.1:g.(3_5)_(10_12)del",
        ),
    ],
)
def test_gains_read_back_as_dup_and_losses_as_del(
    mapper: weaver.VariantMapper,
    term: str | dict[str, typing.Any],
    expected: str,
) -> None:
    change = mapper.to_vrs(weaver.parse("NC_TEST.1:g.(3_5)_(10_12)dup"))
    change["copyChange"] = term
    assert str(mapper.from_vrs(change)) == expected


def test_other_copy_changes_have_no_hgvs(mapper: weaver.VariantMapper) -> None:
    change = mapper.to_vrs(weaver.parse("NC_TEST.1:g.(3_5)_(10_12)dup"))
    change["copyChange"] = "regional base ploidy"
    with pytest.raises(weaver.UnsupportedOperationError):
        mapper.from_vrs(change)
