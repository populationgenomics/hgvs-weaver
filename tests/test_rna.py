"""r. variants through the bindings: respelling, projection and protein prediction."""

import typing

import pytest

import weaver

# c.1 is index 10: ATG GGG CCC AAA ...
SEQUENCE = "A" * 10 + "ATGGGGCCCAAA" + "A" * 100


class Provider:
    """One single-exon coding transcript."""

    def get_transcript(self, ac: str, _ref: str | None) -> dict[str, typing.Any]:
        return {
            "ac": ac,
            "gene": "TEST",
            "cds_start_index": 10,
            "cds_end_index": 20,
            "strand": 1,
            "reference_accession": "NC_TEST.1",
            "exons": [
                {
                    "transcript_start": 0,
                    "transcript_end": 122,
                    "reference_start": 1000,
                    "reference_end": 1121,
                    "alt_strand": 1,
                    "cigar": "122=",
                },
            ],
        }

    def get_seq(self, ac: str, start: int, end: int | None, _kind: str) -> str:
        # The exon puts the transcript at genome index 1000.
        seq = "N" * 1000 + SEQUENCE if ac.startswith("NC_") else SEQUENCE
        return seq[start:end]

    def get_symbol_accessions(self, symbol: str, _s: str, t: str) -> list[tuple[weaver.IdentifierType, str]]:
        if t == "p":
            return [(weaver.IdentifierType.ProteinAccession, "NP_TEST.1")]
        return [(weaver.IdentifierType.GeneSymbol, symbol)]

    def get_identifier_type(self, identifier: str) -> weaver.IdentifierType:
        if identifier.startswith("NC_"):
            return weaver.IdentifierType.GenomicAccession
        return weaver.IdentifierType.TranscriptAccession


def test_r_is_c_in_rna_letters() -> None:
    mapper = weaver.VariantMapper(Provider())
    r = weaver.parse("NM_TEST.1:r.4g>a")
    c = mapper.r_to_c(r)
    assert str(c) == "NM_TEST.1:c.4G>A"
    assert str(mapper.c_to_r(c)) == "NM_TEST.1:r.4g>a"
    assert str(mapper.r_to_g(r)) == str(mapper.c_to_g(c))
    assert str(mapper.c_to_p(r)) == str(mapper.c_to_p(c)) == "NP_TEST.1:p.(Gly2Arg)"
    assert str(mapper.c_to_p(weaver.parse("NM_TEST.1:r.spl"))) == "NP_TEST.1:p.?"
    assert str(mapper.c_to_p(weaver.parse("NM_TEST.1:r.0"))) == "NP_TEST.1:p.0"
    with pytest.raises(weaver.UnsupportedOperationError):
        mapper.r_to_n(r)  # the transcript has a CDS
    with pytest.raises(ValueError, match="RNA variant"):
        mapper.r_to_c(c)
