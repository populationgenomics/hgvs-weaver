"""weaver-validate's rs_nuc_equiv: weaver's verdict on a row's nucleotide description against ClinVar's protein one."""

from __future__ import annotations

import re

import hgvs.parser
import hgvs.variantmapper
import pytest

import weaver
from weaver.cli import provider as cli_provider
from weaver.cli import validate

# c.1 is index 10; codon 2, c.4_6, is GCT (Ala).
TRANSCRIPT = "GCTAGCTAGCATGGCTGGATCCAAGTTCCTGCACGATGAAGTCATCAACCCTCGATACTGGAGCACATAAACGTTGCAAGGTCCATGACC"
PROTEIN = "MAGSKFLHDEVINPRYWST"


class Provider:
    """NM_TEST.1, one exon at NC_TEST.1 indices 1000..1089 on the forward strand, encoding NP_TEST.1."""

    def get_transcript(self, transcript_ac: str, reference_ac: str | None) -> weaver.TranscriptData:  # noqa: ARG002
        return {
            "ac": transcript_ac,
            "gene": "TEST",
            "cds_start_index": 10,
            "cds_end_index": 69,
            "strand": 1,
            "reference_accession": "NC_TEST.1",
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

    def get_seq(self, ac: str, start: int, end: int | None, kind: str | weaver.IdentifierType) -> str:  # noqa: ARG002
        if ac.startswith("NP_"):
            return PROTEIN[start:end]
        seq = "N" * 1000 + TRANSCRIPT if ac.startswith("NC_") else TRANSCRIPT
        return seq[start:end]

    def get_symbol_accessions(
        self,
        symbol: str,
        source_kind: str,  # noqa: ARG002
        target_kind: str,
    ) -> list[tuple[weaver.IdentifierType, str]]:
        if target_kind == "p":
            return [(weaver.IdentifierType.ProteinAccession, "NP_TEST.1")]
        return [(weaver.IdentifierType.GeneSymbol, symbol)]

    def get_identifier_type(self, identifier: str) -> weaver.IdentifierType:
        if identifier.startswith("NC_"):
            return weaver.IdentifierType.GenomicAccession
        if identifier.startswith("NM_"):
            return weaver.IdentifierType.TranscriptAccession
        if identifier.startswith("NP_"):
            return weaver.IdentifierType.ProteinAccession
        return weaver.IdentifierType.GeneSymbol

    def get_transcripts_for_region(self, chrom: str, start: int, end: int) -> list[str]:  # noqa: ARG002
        return ["NM_TEST.1"]


@pytest.fixture(scope="module")
def parser() -> hgvs.parser.Parser:
    return hgvs.parser.Parser()


@pytest.fixture
def worker(monkeypatch: pytest.MonkeyPatch, parser: hgvs.parser.Parser) -> None:
    rp = Provider()
    monkeypatch.setattr(
        validate,
        "_worker",
        validate._Worker(
            rp=rp,
            mapper=weaver.VariantMapper(rp),
            ref_vm=hgvs.variantmapper.VariantMapper(cli_provider.ReferenceHgvsDataProvider(rp)),
            ref_hp=parser,
            fh_results={},
        ),
    )


def row(variant_nuc: str, variant_prot: str) -> dict[str, str]:
    return {"variant_nuc": variant_nuc, "variant_prot": variant_prot, "spdi": "NC_TEST.1:1013:G:A"}


@pytest.mark.usefixtures("worker")
@pytest.mark.parametrize(
    ("variant_prot", "verdict"),
    [("NP_TEST.1:p.(Ala2Thr)", "Identity"), ("NP_TEST.1:p.(Ala2Ser)", "Different")],
)
def test_the_column_is_weavers_verdict_on_the_two_descriptions(variant_prot: str, verdict: str) -> None:
    assert validate.process_variant(row("NM_TEST.1:c.4G>A", variant_prot))["rs_nuc_equiv"] == verdict


@pytest.mark.usefixtures("worker")
def test_a_row_without_a_predicted_protein_is_still_judged() -> None:
    result = validate.process_variant(row("NC_TEST.1:g.1014G>A", "NP_TEST.1:p.(Ala2Thr)"))
    assert result["rs_equiv"] == "Unknown", "the case needs a row the harness predicts no protein for"
    assert result["rs_nuc_equiv"] == "Identity"


@pytest.mark.usefixtures("worker")
@pytest.mark.parametrize(
    ("variant_nuc", "variant_prot", "message"),
    [
        ("NM_TEST.1:c.4T>A", "NP_TEST.1:p.(Ala2Thr)", "expected T, found G"),
        ("NM_TEST.1:c.4G>A", "NP_TEST.1:p.Garbage", "NP_TEST.1:p.Garbage"),
    ],
)
def test_a_description_that_cannot_be_judged_is_an_error_in_the_column(
    variant_nuc: str,
    variant_prot: str,
    message: str,
) -> None:
    result = validate.process_variant(row(variant_nuc, variant_prot))["rs_nuc_equiv"]
    assert result.startswith("ERR:")
    assert re.search(message, result), result


@pytest.mark.usefixtures("worker")
def test_a_row_without_a_clinvar_protein_has_no_verdict() -> None:
    assert validate.process_variant(row("NM_TEST.1:c.4G>A", ""))["rs_nuc_equiv"] == ""
