"""A provider says translation does not initiate at codon 1 with cds_start_initiates."""

from __future__ import annotations

import typing

import weaver

CHROM = "NC_000099.1"
RECORD = "GATGCTGAACAACCACTTTCTGAAGCTTAA" + "GGCCGGCCGGC"  # Asp Ala Glu Gln Pro Leu Ser Glu Ala Ter
PROTEIN = "DAEQPLSEA"


class Provider:
    """One transcript, TX_FRAG.1, numbered from GAT without translation initiating there."""

    def get_transcript(self, transcript_ac: str, _reference_ac: str | None) -> dict[str, typing.Any]:
        """Returns the model, with cds_start_initiates False."""
        return {
            "ac": transcript_ac,
            "gene": "GENE",
            "cds_start_index": 0,
            "cds_end_index": 29,
            "cds_start_initiates": False,
            "strand": 1,
            "reference_accession": CHROM,
            "exons": [
                {
                    "transcript_start": 0,
                    "transcript_end": 41,
                    "reference_start": 0,
                    "reference_end": 40,
                    "alt_strand": 1,
                    "cigar": "41=",
                },
            ],
        }

    def get_seq(self, ac: str, start: int, end: int | None, _kind: str) -> str:
        """Returns a slice of the record, which is also the genome, or of the protein."""
        seq = PROTEIN if ac.startswith("NP_") else RECORD
        return seq[max(start, 0) : end]

    def get_symbol_accessions(self, symbol: str, _source_kind: str, target_kind: str) -> list[typing.Any]:
        """Names the transcript's protein."""
        if target_kind == "p":
            return [(weaver.IdentifierType.ProteinAccession, symbol.replace("TX_", "NP_"))]
        return []

    def get_identifier_type(self, identifier: str) -> weaver.IdentifierType:
        """Genomic for the chromosome, protein for NP_, transcript otherwise."""
        if identifier == CHROM:
            return weaver.IdentifierType.GenomicAccession
        if identifier.startswith("NP_"):
            return weaver.IdentifierType.ProteinAccession
        return weaver.IdentifierType.TranscriptAccession


def test_a_change_in_codon_1_is_a_substitution_where_translation_does_not_initiate() -> None:
    """With cds_start_initiates False, c.1G>A on GAT is Asp1Asn, not p.(Met1?)."""
    mapper = weaver.VariantMapper(Provider())
    assert mapper.c_to_p(weaver.parse("TX_FRAG.1:c.1G>A"), None).format() == "NP_FRAG.1:p.(Asp1Asn)"
