"""A CDS that runs off an end of the record refuses the positions numbered from that end (issue #47)."""

import typing

import pytest

import weaver

GENOME = "GG" + "ATG" + "GCTGAACAACCACTTTCTGAAGCT" + "TAA" + "GGCCGGCCGGC" + "ACGTACGTACGTACGTACGT"
CHROM = "NC_000099.1"
RECORD = GENOME[:43]  # one exon; the CDS is indices 2..31
PROTEIN = "MAEQPLSEA"
OPEN = {"TX_WHOLE.1": {}, "TX_OPEN5.1": {"cds_start_open": True}, "TX_OPEN3.1": {"cds_end_open": True}}


class Provider:
    """Three copies of one transcript model, differing only in what they say about the CDS ends."""

    def get_transcript(self, transcript_ac: str, _reference_ac: str | None) -> dict[str, typing.Any]:
        """Returns the model, with the open-end flags only where the publisher states them."""
        return {
            "ac": transcript_ac,
            "gene": "GENE",
            "cds_start_index": 2,
            "cds_end_index": 31,
            "strand": 1,
            "reference_accession": CHROM,
            "exons": [
                {
                    "transcript_start": 0,
                    "transcript_end": 43,
                    "reference_start": 0,
                    "reference_end": 42,
                    "alt_strand": 1,
                    "cigar": "43=",
                },
            ],
            **OPEN[transcript_ac],
        }

    def get_seq(self, ac: str, start: int, end: int | None, _kind: str) -> str:
        """Returns a slice of the genome, the record or the protein."""
        seq = GENOME if ac == CHROM else PROTEIN if ac.startswith("NP_") else RECORD
        return seq[max(start, 0) : end]

    def get_symbol_accessions(self, symbol: str, _source_kind: str, target_kind: str) -> list[typing.Any]:
        """Names each transcript's protein."""
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


def test_a_position_numbered_from_an_open_end_is_refused_in_both_directions() -> None:
    """The flags take effect on c. to g. and g. to c. alike; a whole CDS numbers from both codons."""
    mapper = weaver.VariantMapper(Provider())
    assert mapper.c_to_g(weaver.parse("TX_WHOLE.1:c.*1G>T"), CHROM).format() == "NC_000099.1:g.33G>T"
    with pytest.raises(weaver.ValidationError, match="open at the 3' end"):
        mapper.c_to_g(weaver.parse("TX_OPEN3.1:c.*1G>T"), CHROM)
    with pytest.raises(weaver.ValidationError, match="open at the 5' end"):
        mapper.g_to_c(weaver.parse("NC_000099.1:g.3A>T"), "TX_OPEN5.1")
    assert mapper.c_to_g(weaver.parse("TX_OPEN5.1:c.*1G>T"), CHROM).format() == "NC_000099.1:g.33G>T"


def test_a_protein_consequence_that_depends_on_an_open_end_is_refused() -> None:
    """A 5'-open CDS predicts nothing; a 3'-open one predicts only an in-frame change inside what it carries."""
    mapper = weaver.VariantMapper(Provider())
    with pytest.raises(weaver.ValidationError, match="no protein consequence"):
        mapper.c_to_p(weaver.parse("TX_OPEN5.1:c.5C>T"), None)
    assert mapper.c_to_p(weaver.parse("TX_OPEN3.1:c.5C>T"), None).format() == "NP_OPEN3.1:p.(Ala2Val)"
    with pytest.raises(weaver.ValidationError, match="open at the 3' end"):
        mapper.c_to_p(weaver.parse("TX_OPEN3.1:c.5del"), None)
