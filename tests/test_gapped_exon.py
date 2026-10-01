"""A transcript base the genome lacks projects as an insertion, or fails at an exon edge (issue #43)."""

import typing

import pytest

import weaver

GENOME = "ACGTTGCAAGCTAGCTTACGGATCCATGCAAGTCGATCGGCTAGCTAAGGCCTTAAGCTT"
CHROM = "NC_000099.1"
# One exon on genome 10..39 with a T inserted after genome index 24: n.16 is the inserted base.
TX_INS = GENOME[10:25] + "T" + GENOME[25:40]
# One exon whose first five bases are absent from the genome; the rest aligns to genome 30..49.
TX_CLIP = "GGGCC" + GENOME[30:50]
# The same exon with three bases, TCG, inserted after genome index 24: n.16_18 is the run.
TX_GAP3 = GENOME[10:25] + "TCG" + GENOME[25:40]
# The same exon lacking genome index 25, the A at g.26: c.15 and c.16 flank the gap (the CDS runs from n.1).
TX_DEL = GENOME[10:25] + GENOME[26:40]
MODELS = {
    "TX_INS.1": (TX_INS, [(0, 31, 10, 39, "15=1I15=")]),
    "TX_GAP3.1": (TX_GAP3, [(0, 33, 10, 39, "15=3I15=")]),
    "TX_DEL.1": (TX_DEL, [(0, 29, 10, 39, "15=1D14=")]),
    "TX_CLIP.1": (TX_CLIP, [(0, 25, 30, 49, "5I20=")]),
}


class Provider:
    """An in-memory DataProvider with per-exon cigars carrying insertions."""

    def get_transcript(self, transcript_ac: str, _reference_ac: str | None) -> dict[str, typing.Any]:
        """Returns the transcript model."""
        _, exons = MODELS[transcript_ac]
        coding = transcript_ac == "TX_DEL.1"
        return {
            "ac": transcript_ac,
            "gene": "GENE",
            "cds_start_index": 0 if coding else None,
            "cds_end_index": 28 if coding else None,
            "strand": 1,
            "reference_accession": CHROM,
            "exons": [
                {
                    "transcript_start": t0,
                    "transcript_end": t1,
                    "reference_start": g0,
                    "reference_end": g1,
                    "alt_strand": 1,
                    "cigar": cigar,
                }
                for t0, t1, g0, g1, cigar in exons
            ],
        }

    def get_seq(self, ac: str, start: int, end: int | None, _kind: str) -> str:
        """Returns a slice of the genome or a transcript."""
        return (GENOME if ac == CHROM else MODELS[ac][0])[max(start, 0) : end]

    def get_symbol_accessions(self, _symbol: str, _source_kind: str, _target_kind: str) -> list[typing.Any]:
        """No symbols."""
        return []

    def get_identifier_type(self, identifier: str) -> weaver.IdentifierType:
        """Genomic for the chromosome, transcript for everything else."""
        if identifier == CHROM:
            return weaver.IdentifierType.GenomicAccession
        return weaver.IdentifierType.TranscriptAccession


def n_to_g(mapper: weaver.VariantMapper, name: str) -> str:
    """Projects an n. variant to the genome."""
    return mapper.n_to_g(weaver.parse(name), CHROM).format()


def test_a_base_the_genome_lacks_projects_as_an_insertion() -> None:
    """The inserted base goes between its genomic neighbours; its neighbours project as before."""
    mapper = weaver.VariantMapper(Provider())
    assert n_to_g(mapper, "TX_INS.1:n.15C>G") == "NC_000099.1:g.25C>G"
    assert n_to_g(mapper, "TX_INS.1:n.16T>G") == "NC_000099.1:g.25_26insG"
    assert n_to_g(mapper, "TX_INS.1:n.17A>G") == "NC_000099.1:g.26A>G"


def test_an_edit_inside_a_run_the_genome_lacks_carries_the_whole_run() -> None:
    """No position inside the run exists on the genome, so the run as changed is what is inserted."""
    mapper = weaver.VariantMapper(Provider())
    assert n_to_g(mapper, "TX_GAP3.1:n.17C>A") == "NC_000099.1:g.25_26insTAG"
    assert n_to_g(mapper, "TX_GAP3.1:n.17del") == "NC_000099.1:g.25_26insTG"
    assert n_to_g(mapper, "TX_GAP3.1:n.16_18del") == "NC_000099.1:g.25="


def test_a_soft_clipped_base_has_no_genomic_position() -> None:
    """A clipped base lies outside any flanking pair and raises."""
    mapper = weaver.VariantMapper(Provider())
    for name in ("TX_CLIP.1:n.1G>C", "TX_CLIP.1:n.3G>C", "TX_CLIP.1:n.5C>G"):
        with pytest.raises(weaver.ValidationError, match="edge of an exon"):
            n_to_g(mapper, name)
    assert n_to_g(mapper, "TX_CLIP.1:n.6A>G") == "NC_000099.1:g.31A>G"


def test_a_genome_base_the_transcript_lacks_projects_as_an_insertion() -> None:
    """The genome reads C A T over g.25..27 and the record C T over c.15..16; g.26 is the base it lacks."""
    mapper = weaver.VariantMapper(Provider())
    assert mapper.g_to_c(weaver.parse("NC_000099.1:g.26A>G"), "TX_DEL.1").format() == "TX_DEL.1:c.15_16insG"
    assert mapper.g_to_c(weaver.parse("NC_000099.1:g.26del"), "TX_DEL.1").format() == "TX_DEL.1:c.15="
    assert mapper.g_to_c(weaver.parse("NC_000099.1:g.25C>G"), "TX_DEL.1").format() == "TX_DEL.1:c.15C>G"
