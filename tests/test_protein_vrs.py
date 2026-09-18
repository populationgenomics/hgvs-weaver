"""Protein variants have canonical alleles, SPDI and VRS on their protein sequence."""

import typing

import weaver

PROTEIN = "MKLAAAYRQ"


class ProteinProvider:
    """Serves one protein sequence."""

    def get_transcript(self, ac: str, _ref: str | None) -> dict[str, typing.Any]:
        raise weaver.DataProviderError(f"no transcript {ac}")

    def get_seq(self, ac: str, start: int, end: int | None, _kind: str) -> str:
        assert ac == "NP_TEST.1"
        return PROTEIN[start:end]

    def get_symbol_accessions(self, _symbol: str, _s: str, _t: str) -> list[tuple[weaver.IdentifierType, str]]:
        return []

    def get_identifier_type(self, identifier: str) -> weaver.IdentifierType:
        if identifier.startswith("NP_"):
            return weaver.IdentifierType.ProteinAccession
        return weaver.IdentifierType.GenomicAccession


def test_protein_allele_spdi_and_vrs() -> None:
    provider = ProteinProvider()
    mapper = weaver.VariantMapper(provider)
    var = weaver.parse("NP_TEST.1:p.Ala5del")
    assert mapper.to_spdi_unambiguous(var) == "NP_TEST.1:3:AAA:AA"
    vrs = mapper.to_vrs(var)
    assert vrs["location"]["sequenceReference"]["residueAlphabet"] == "aa"
    assert vrs["location"]["sequenceReference"]["moleculeType"] == "protein"
    assert vrs["state"]["type"] == "ReferenceLengthExpression"
    assert mapper.vrs_id(var) == vrs["id"]
    assert var.validate(provider)
    assert not weaver.parse("NP_TEST.1:p.Tyr5del").validate(provider)
