"""The binding of consistent_variants and the AccessionCatalogue protocol."""

from __future__ import annotations

import typing
from collections.abc import Sequence

import pytest

import weaver

# c.1 is index 5; c.4_6 is AAA on NM_C.1 and GAA on NM_C.2.
_RECORDS = {
    "NM_C.1": "GGGGGATGAAACTGGCCTATCGCTAACCGTATAAGTAAGG",
    "NM_C.2": "GGGGGATGGAACTGGCCTATCGCTAACCGTATAAGTAAGG",
}


class Provider:
    """Two versions of one transcript, NM_C.1 and NM_C.2, of gene CG, one exon each on NC_C.1."""

    def get_transcript(self, transcript_ac: str, reference_ac: str | None) -> dict[str, typing.Any]:  # noqa: ARG002
        """One exon at genomic index 10, the CDS at transcript indices 5 to 25."""
        if transcript_ac not in _RECORDS:
            raise weaver.DataProviderError(f"no transcript {transcript_ac}")
        length = len(_RECORDS[transcript_ac])
        return {
            "ac": transcript_ac,
            "gene": "CG",
            "cds_start_index": 5,
            "cds_end_index": 25,
            "strand": 1,
            "reference_accession": "NC_C.1",
            "reference_alignment_method": "splign",
            "exons": [
                {
                    "transcript_start": 0,
                    "transcript_end": length,
                    "reference_start": 10,
                    "reference_end": 10 + length - 1,
                    "alt_strand": 1,
                    "cigar": f"{length}=",
                },
            ],
        }

    def get_seq(self, ac: str, start: int, end: int | None, kind: str | weaver.IdentifierType) -> str:  # noqa: ARG002
        """A transcript's record, or the genome holding NM_C.1 at index 10."""
        sequence = "T" * 10 + _RECORDS["NM_C.1"] + "C" * 10 if ac == "NC_C.1" else _RECORDS[ac]
        return sequence[start:end]

    def get_symbol_accessions(
        self,
        symbol: str,
        source_kind: str,  # noqa: ARG002
        target_kind: str,
    ) -> list[tuple[weaver.IdentifierType, str]]:
        """CG's transcripts."""
        if symbol == "CG" and target_kind == "c":
            return [(weaver.IdentifierType.TranscriptAccession, ac) for ac in _RECORDS]
        return []

    def get_identifier_type(self, identifier: str) -> weaver.IdentifierType:
        """CG is a gene symbol; NC_ is genomic; anything else a transcript."""
        if identifier == "CG":
            return weaver.IdentifierType.GeneSymbol
        if identifier.startswith("NC_"):
            return weaver.IdentifierType.GenomicAccession
        return weaver.IdentifierType.TranscriptAccession


class Catalogue:
    """Lists NM_C's two versions and answers each one's info as it was given."""

    def __init__(
        self,
        info: dict[str, tuple[str | weaver.AccessionStatus, Sequence[str | weaver.AccessionTag]]],
    ) -> None:
        self.info = info

    def get_accession_versions(self, ac: str) -> list[str]:
        """NM_C's versions, newest first."""
        return ["NM_C.2", "NM_C.1"] if ac == "NM_C" else []

    def get_accession_info(
        self,
        ac: str,
    ) -> tuple[str | weaver.AccessionStatus, Sequence[str | weaver.AccessionTag]] | None:
        """The info this catalogue was given for ac."""
        return self.info.get(ac)


def test_a_catalogue_answering_names_gives_each_variant_its_status_and_tags() -> None:
    catalogue = Catalogue({"NM_C.2": ("current", ["mane_select"]), "NM_C.1": ("superseded", [])})
    found = weaver.VariantMapper(Provider()).consistent_variants(weaver.parse("NM_C:c.6A>G"), catalogue)
    assert [(c.variant.format(), c.status, c.tags) for c in found] == [
        ("NM_C.2:c.6A>G", weaver.AccessionStatus.Current, [weaver.AccessionTag.ManeSelect]),
        ("NM_C.1:c.6A>G", weaver.AccessionStatus.Superseded, []),
    ]


def test_a_catalogue_answering_members_is_read_as_one_answering_names() -> None:
    catalogue = Catalogue(
        {
            "NM_C.1": (
                weaver.AccessionStatus.Current,
                [weaver.AccessionTag.ManeSelect, weaver.AccessionTag.RefSeqSelect],
            ),
        },
    )
    found = weaver.VariantMapper(Provider()).consistent_variants(weaver.parse("NM_C:c.4A>G"), catalogue)
    assert [(c.variant.format(), c.status, c.tags) for c in found] == [
        (
            "NM_C.1:c.4A>G",
            weaver.AccessionStatus.Current,
            [weaver.AccessionTag.ManeSelect, weaver.AccessionTag.RefSeqSelect],
        ),
    ]


def test_a_status_the_protocol_does_not_name_is_a_provider_error() -> None:
    catalogue = Catalogue({"NM_C.1": ("retired", [])})
    with pytest.raises(weaver.DataProviderError, match="retired"):
        weaver.VariantMapper(Provider()).consistent_variants(weaver.parse("NM_C:c.4A>G"), catalogue)


def test_a_tag_the_protocol_does_not_name_is_a_provider_error() -> None:
    catalogue = Catalogue({"NM_C.1": ("current", ["mane"])})
    with pytest.raises(weaver.DataProviderError, match="mane"):
        weaver.VariantMapper(Provider()).consistent_variants(weaver.parse("NM_C:c.4A>G"), catalogue)


def test_an_unversioned_accession_without_a_catalogue_is_a_provider_error() -> None:
    with pytest.raises(weaver.DataProviderError, match="catalogue"):
        weaver.VariantMapper(Provider()).consistent_variants(weaver.parse("NM_C:c.4A>G"))


def test_a_gene_symbol_without_a_catalogue_gives_its_variants_with_no_status() -> None:
    found = weaver.VariantMapper(Provider()).consistent_variants(weaver.parse("CG:c.4A>G"))
    assert [(c.variant.format(), c.status, c.tags) for c in found] == [("NM_C.1:c.4A>G", None, [])]


def test_a_description_stating_a_base_is_marked_as_stating_one() -> None:
    found = weaver.VariantMapper(Provider()).consistent_variants(weaver.parse("NM_C.1:c.4A>G"))
    assert [(c.variant.format(), c.states_bases) for c in found] == [("NM_C.1:c.4A>G", True)]
