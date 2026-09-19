"""VRS Alleles and SPDI strings read back into HGVS through the bindings."""

import json
import typing

import pytest

import weaver

GENOME = "ACGTTTGCAAGGCTAGCTAGCTTTTAACGGGATCGATCGA"
PROTEIN = "MKLAAAYRQ"
SEQUENCES = {"NC_TEST.1": GENOME, "NP_TEST.1": PROTEIN}


class Provider:
    """Two sequences, with the optional refget lookups implemented."""

    def __init__(self, *, lookup: bool) -> None:
        self.lookup = lookup

    def get_transcript(self, ac: str, _ref: str | None) -> dict[str, typing.Any]:
        raise weaver.DataProviderError(f"no transcript {ac}")

    def get_seq(self, ac: str, start: int, end: int | None, _kind: str) -> str:
        return SEQUENCES[ac][start:end]

    def get_symbol_accessions(self, _symbol: str, _s: str, _t: str) -> list[tuple[weaver.IdentifierType, str]]:
        return []

    def get_identifier_type(self, identifier: str) -> weaver.IdentifierType:
        if identifier.startswith("NP_"):
            return weaver.IdentifierType.ProteinAccession
        return weaver.IdentifierType.GenomicAccession

    def get_refget_accession(self, _ac: str) -> str | None:
        return None  # let weaver compute it

    def get_accession_for_refget(self, refget: str) -> str | None:
        if not self.lookup:
            return None
        mapper = weaver.VariantMapper(self)
        for ac in SEQUENCES:
            probe = "p.Lys2Leu" if ac.startswith("NP_") else "g.7G>C"
            allele = mapper.to_vrs(weaver.parse(f"{ac}:{probe}"))
            if allele["location"]["sequenceReference"]["refgetAccession"] == refget:
                return ac
        return None


@pytest.mark.parametrize(
    ("hgvs", "expected"),
    [
        ("NC_TEST.1:g.7G>C", "NC_TEST.1:g.7G>C"),
        ("NC_TEST.1:g.4del", "NC_TEST.1:g.6del"),
        ("NC_TEST.1:g.3_4insT", "NC_TEST.1:g.6dup"),
        ("NC_TEST.1:g.(?_5)_(10_?)del", "NC_TEST.1:g.(?_5)_(10_?)del"),
        ("NP_TEST.1:p.Ala4del", "NP_TEST.1:p.Ala6del"),
    ],
)
def test_vrs_round_trip(hgvs: str, expected: str) -> None:
    provider = Provider(lookup=True)
    mapper = weaver.VariantMapper(provider, refget=provider)
    allele = mapper.to_vrs(weaver.parse(hgvs))
    back = mapper.from_vrs(allele)
    assert str(back) == expected
    assert mapper.vrs_id(back) == allele["id"]


def test_accession_is_passed_when_the_provider_cannot_look_it_up() -> None:
    mapper = weaver.VariantMapper(Provider(lookup=False))
    allele = mapper.to_vrs(weaver.parse("NC_TEST.1:g.7G>C"))
    with pytest.raises(weaver.DataProviderError):
        mapper.from_vrs(allele)
    assert str(mapper.from_vrs(allele, accession="NC_TEST.1")) == "NC_TEST.1:g.7G>C"
    assert str(mapper.from_vrs(json.dumps(allele), "NC_TEST.1")) == "NC_TEST.1:g.7G>C"


def test_spdi_reads_back() -> None:
    mapper = weaver.VariantMapper(Provider(lookup=False))
    assert str(mapper.from_spdi("NC_TEST.1:6:G:C")) == "NC_TEST.1:g.7G>C"
    assert str(mapper.from_spdi("NC_TEST.1:6:1:C")) == "NC_TEST.1:g.7G>C"
    assert str(mapper.from_spdi("NP_TEST.1:3:AAA:AA")) == "NP_TEST.1:p.Ala6del"
    with pytest.raises(weaver.ValidationError):
        mapper.from_spdi("NC_TEST.1:6:A:C")
