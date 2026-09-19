"""RefgetProvider against a recorded refget server."""

import base64
import hashlib
import json
import re
import typing
import urllib.parse

import pytest

import weaver
from weaver.refget import RefgetProvider, sq_digest

GENOME = "ACGTTTGCAAGGCTAGCTAGCTTTTAACGGGATCGATCGA"
DIGEST = "SQ." + base64.urlsafe_b64encode(hashlib.sha512(GENOME.encode()).digest()[:24]).decode().rstrip("=")
BASE = "http://refget.test/seqrepo/1"


class Server:
    """A SeqRepo-like server: knows the sequence as refseq:NC_000099.1 and by digest."""

    def __init__(self) -> None:
        self.requests: list[str] = []
        self.metadata: dict[str, typing.Any] = {
            "length": len(GENOME),
            "md5": hashlib.md5(GENOME.encode()).hexdigest(),  # noqa: S324 - refget's own field
            "aliases": [
                {"naming_authority": "refseq", "alias": "refseq:NC_000099.1"},
                {"naming_authority": "INSDC", "alias": "insdc:CM000000.1"},
                {"naming_authority": "ga4gh", "alias": f"ga4gh:{DIGEST}"},
            ],
        }

    def __call__(self, url: str, accept: str) -> tuple[int, bytes]:
        self.requests.append(url)
        assert url.startswith(BASE + "/sequence/")
        rest = url[len(BASE) + len("/sequence/") :]
        path, _, query = rest.partition("?")
        sequence_id = urllib.parse.unquote(path.removesuffix("/metadata"))
        if sequence_id not in ("refseq:NC_000099.1", f"ga4gh:{DIGEST}"):
            return 404, b"not found"
        if path.endswith("/metadata"):
            assert "json" in accept
            return 200, json.dumps({"metadata": self.metadata}).encode()
        assert "plain" in accept
        params = dict(urllib.parse.parse_qsl(query))
        start, end = int(params["start"]), int(params["end"])
        if end > len(GENOME):
            return 416, b"range not satisfiable"
        return 200, GENOME[start:end].encode()


def test_sequences_come_from_the_server_with_ranges_clamped() -> None:
    server = Server()
    provider = RefgetProvider(BASE, transport=server)
    assert provider.get_seq("NC_000099.1", 6, 7, "genomic_accession") == "G"
    assert provider.get_seq("NC_000099.1", 30, None, "genomic_accession") == GENOME[30:]
    assert provider.get_seq("NC_000099.1", 35, 100, "genomic_accession") == GENOME[35:]
    assert provider.get_seq("NC_000099.1", 100, 110, "genomic_accession") == ""
    # The bare accession was tried, then the refseq: alias, then remembered.
    assert [u for u in server.requests if u.endswith("/metadata")] == [
        f"{BASE}/sequence/NC_000099.1/metadata",
        f"{BASE}/sequence/refseq:NC_000099.1/metadata",
    ]


def test_refget_lookups_both_ways() -> None:
    provider = RefgetProvider(BASE, transport=Server())
    assert provider.get_refget_accession("NC_000099.1") == DIGEST
    assert provider.get_accession_for_refget(DIGEST) == "NC_000099.1"
    assert provider.get_accession_for_refget(f"ga4gh:{DIGEST}") == "NC_000099.1"
    assert provider.get_accession_for_refget("SQ.nope") is None
    with pytest.raises(weaver.DataProviderError, match="knows no sequence"):
        provider.get_seq("NC_OTHER.1", 0, 1, "genomic_accession")


def test_v1_metadata_gives_the_digest_from_trunc512_and_insdc_names() -> None:
    trunc512 = hashlib.sha512(GENOME.encode()).digest()[:24].hex()
    assert sq_digest({"trunc512": trunc512}) == DIGEST
    assert sq_digest({"ga4gh": f"ga4gh:{DIGEST}"}) == DIGEST
    assert sq_digest({"md5": "x"}) is None
    server = Server()
    server.metadata["aliases"] = [{"naming_authority": "INSDC", "alias": "CM000000.1"}]
    assert RefgetProvider(BASE, transport=server).get_accession_for_refget(DIGEST) == "CM000000.1"


def test_vrs_round_trips_through_the_server() -> None:
    provider = RefgetProvider(BASE, transport=Server())
    mapper = weaver.VariantMapper(provider, refget=provider)
    allele = mapper.to_vrs(weaver.parse("NC_000099.1:g.4del"))
    assert allele["location"]["sequenceReference"]["refgetAccession"] == DIGEST
    assert str(mapper.from_vrs(allele)) == "NC_000099.1:g.6del"
    assert mapper.to_spdi_unambiguous(weaver.parse("NC_000099.1:g.4del")) == "NC_000099.1:3:TTT:TT"


def test_transcripts_come_from_the_wrapped_provider_or_not_at_all() -> None:
    bare = RefgetProvider(BASE, transport=Server())
    with pytest.raises(weaver.DataProviderError, match="sequences only"):
        bare.get_transcript("NM_TEST.1", None)
    assert bare.get_identifier_type("NP_1.1") == weaver.IdentifierType.ProteinAccession
    assert bare.get_identifier_type("NM_1.1") == weaver.IdentifierType.TranscriptAccession
    assert bare.get_identifier_type("NC_1.1") == weaver.IdentifierType.GenomicAccession

    class Transcripts:
        def get_transcript(self, ac: str, _ref: str | None) -> weaver.TranscriptData:
            return typing.cast("weaver.TranscriptData", {"ac": ac})

        def get_seq(self, _ac: str, _start: int, _end: int | None, _kind: str) -> str:
            raise AssertionError("sequences must come from refget")

        def get_symbol_accessions(self, _symbol: str, _source: str, _target: str) -> list[tuple[str, str]]:
            return [("transcript_accession", "NM_TEST.1")]

        def get_identifier_type(self, _identifier: str) -> str:
            return "gene_symbol"

    wrapped = RefgetProvider(BASE, Transcripts(), transport=Server())
    assert wrapped.get_transcript("NM_TEST.1", None) == {"ac": "NM_TEST.1"}
    assert wrapped.get_symbol_accessions("X", "gene_symbol", "transcript_accession") == [
        ("transcript_accession", "NM_TEST.1"),
    ]
    assert wrapped.get_identifier_type("anything") == "gene_symbol"
    assert re.match(r"^SQ\.[A-Za-z0-9_-]{32}$", wrapped.get_refget_accession("NC_000099.1") or "")
