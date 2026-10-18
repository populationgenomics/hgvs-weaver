"""An identity that covers no whole codon makes no protein statement, and the binding says so by type."""

from __future__ import annotations

import typing

import pytest

import weaver

CHROM = "NC_000099.1"
RECORD = "GG" + "ATGGCTGAACAA" + "TAA" + "GGCCGGCC"  # the CDS is indices 2..16: Met Ala Glu Gln Ter
PROTEIN = "MAEQ"


class Provider:
    """One single-exon transcript, TX_W.1, encoding NP_W.1."""

    def get_transcript(self, transcript_ac: str, _reference_ac: str | None) -> dict[str, typing.Any]:
        """Returns the transcript model."""
        return {
            "ac": transcript_ac,
            "gene": "GENE",
            "cds_start_index": 2,
            "cds_end_index": 16,
            "strand": 1,
            "reference_accession": CHROM,
            "exons": [
                {
                    "transcript_start": 0,
                    "transcript_end": len(RECORD),
                    "reference_start": 0,
                    "reference_end": len(RECORD) - 1,
                    "alt_strand": 1,
                    "cigar": f"{len(RECORD)}=",
                },
            ],
        }

    def get_seq(self, ac: str, start: int, end: int | None, _kind: str) -> str:
        """Returns a slice of the record (which is also the genome) or the protein."""
        seq = PROTEIN if ac.startswith("NP_") else RECORD
        return seq[max(start, 0) : end]

    def get_symbol_accessions(self, _symbol: str, _source_kind: str, target_kind: str) -> list[typing.Any]:
        """Names the transcript's protein."""
        if target_kind == "p":
            return [(weaver.IdentifierType.ProteinAccession, "NP_W.1")]
        return []

    def get_identifier_type(self, identifier: str) -> weaver.IdentifierType:
        """Genomic for the chromosome, protein for NP_, transcript otherwise."""
        if identifier == CHROM:
            return weaver.IdentifierType.GenomicAccession
        if identifier.startswith("NP_"):
            return weaver.IdentifierType.ProteinAccession
        return weaver.IdentifierType.TranscriptAccession


def test_one_base_unchanged_raises_no_protein_statement_error() -> None:
    mapper = weaver.VariantMapper(Provider())
    with pytest.raises(weaver.NoProteinStatementError, match=r"TX_W\.1:c\.5="):
        mapper.c_to_p(weaver.parse("TX_W.1:c.5="), None)
