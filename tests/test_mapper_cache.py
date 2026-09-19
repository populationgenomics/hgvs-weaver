"""A VariantMapper keeps what it fetched for as long as it lives."""

import typing

import weaver

SEQUENCE = "A" * 10 + "ATGGGGCCCAAA" + "A" * 100


class CountingProvider:
    """A single-exon transcript whose sequence fetches are counted."""

    def __init__(self) -> None:
        self.seq_calls = 0

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

    def get_seq(self, _ac: str, start: int, end: int | None, _kind: str) -> str:
        self.seq_calls += 1
        return SEQUENCE[start:end]

    def get_symbol_accessions(self, symbol: str, _s: str, t: str) -> list[tuple[weaver.IdentifierType, str]]:
        if t == "p":
            return [(weaver.IdentifierType.ProteinAccession, "NP_TEST.1")]
        return [(weaver.IdentifierType.GeneSymbol, symbol)]

    def get_identifier_type(self, identifier: str) -> weaver.IdentifierType:
        if identifier.startswith("NC_"):
            return weaver.IdentifierType.GenomicAccession
        return weaver.IdentifierType.TranscriptAccession


def test_sequences_are_fetched_once_per_mapper() -> None:
    provider = CountingProvider()
    mapper = weaver.VariantMapper(provider)
    var = weaver.parse("NM_TEST.1:c.4G>A")
    assert str(mapper.c_to_p(var)) == "NP_TEST.1:p.(Gly2Arg)"
    fetched = provider.seq_calls
    assert fetched > 0
    assert str(mapper.c_to_p(var)) == "NP_TEST.1:p.(Gly2Arg)"
    assert str(mapper.normalize_variant(weaver.parse("NM_TEST.1:c.4del"))) == "NM_TEST.1:c.6del"
    assert provider.seq_calls == fetched, "a later call refetched what the mapper had"
    # A new mapper starts cold.
    weaver.VariantMapper(provider).c_to_p(var)
    assert provider.seq_calls > fetched


class Digests:
    """A Refget lookup that answers from a table."""

    def __init__(self, table: dict[str, str]) -> None:
        self.table = table
        self.asked: list[str] = []

    def get_refget_accession(self, ac: str) -> str | None:
        self.asked.append(ac)
        return self.table.get(ac)

    def get_accession_for_refget(self, refget: str) -> str | None:
        return next((ac for ac, d in self.table.items() if d == refget), None)


def test_a_refget_lookup_is_asked_before_the_sequence_is_hashed() -> None:
    provider = CountingProvider()
    digests = Digests({"NC_TEST.1": "SQ.not-really-a-digest"})
    mapper = weaver.VariantMapper(provider, refget=digests)
    allele = mapper.to_vrs(weaver.parse("NC_TEST.1:g.14G>A"))
    assert allele["location"]["sequenceReference"]["refgetAccession"] == "SQ.not-really-a-digest"
    assert digests.asked == ["NC_TEST.1"]
    # Asked once; the answer is cached with the sequence.
    mapper.to_vrs(weaver.parse("NC_TEST.1:g.14G>A"))
    assert digests.asked == ["NC_TEST.1"]
    # And it resolves the other way.
    assert str(mapper.from_vrs(allele)) == "NC_TEST.1:g.14G>A"
