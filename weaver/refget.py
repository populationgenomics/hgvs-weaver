"""A DataProvider over the GA4GH refget protocol.

Refget servers hand out sequences and their metadata by digest or alias:
``GET {base}/sequence/{id}?start=&end=`` returns an interbase sub-range and
``GET {base}/sequence/{id}/metadata`` the digests, length and aliases. That is
exactly what weaver needs for ``get_seq`` and for the ``Refget`` lookup, so
:class:`RefgetProvider` is both a DataProvider (sequences from the server,
transcript models from another provider) and a Refget; pass it as either or
both: ``VariantMapper(provider, refget=provider)``.

Any refget v1 or v2 server works, for example the biocommons SeqRepo REST
service (``http://localhost:5000/seqrepo/1``) or EBI's
(``https://www.ebi.ac.uk/ena/cram``). Which names a server knows a sequence by
is up to the server: SeqRepo lists RefSeq accessions, EBI lists INSDC ones.
"""

from __future__ import annotations

import base64
import json
import re
import typing
import urllib.error
import urllib.parse
import urllib.request

from weaver import DataProviderError, IdentifierType, TranscriptData

#: ``(url, accept) -> (status, body)``: how the provider talks HTTP. Tests
#: substitute a recording.
Transport = typing.Callable[[str, str], tuple[int, bytes]]


class TranscriptSource(typing.Protocol):
    """The part of a DataProvider that is not about sequences."""

    def get_transcript(self, transcript_ac: str, reference_ac: str | None) -> TranscriptData: ...

    def get_symbol_accessions(
        self,
        symbol: str,
        source_kind: str,
        target_kind: str,
    ) -> list[tuple[str, str]] | list[tuple[IdentifierType, str]]: ...

    def get_identifier_type(self, identifier: str) -> str | IdentifierType: ...


_SEQUENCE_ACCEPT = "text/vnd.ga4gh.refget.v2.0.0+plain, text/vnd.ga4gh.refget.v1.0.0+plain, text/plain"
_METADATA_ACCEPT = (
    "application/vnd.ga4gh.refget.v2.0.0+json, application/vnd.ga4gh.refget.v1.0.0+json, application/json"
)
_REFSEQ = re.compile(r"^(?:[NX][CGMPRTW]|A[CP])_\d+(?:\.\d+)?$")


def urllib_transport(timeout: float = 30.0) -> Transport:
    """The default transport: ``urllib`` over http or https, with a timeout."""

    def fetch(url: str, accept: str) -> tuple[int, bytes]:
        if urllib.parse.urlparse(url).scheme not in ("http", "https"):
            raise DataProviderError(f"refget URL must be http or https: {url}")
        request = urllib.request.Request(url, headers={"Accept": accept})  # noqa: S310 - scheme checked above
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:  # noqa: S310 - scheme checked above
                return response.status, response.read()
        except urllib.error.HTTPError as e:
            return e.code, e.read()
        except urllib.error.URLError as e:
            raise DataProviderError(f"refget request to {url} failed: {e.reason}") from e

    return fetch


def sq_digest(metadata: dict[str, typing.Any]) -> str | None:
    """The ``SQ.`` refget accession in a metadata record, however the server spells it.

    Refget v2 gives it as ``ga4gh``; SeqRepo lists it as a ``ga4gh:SQ.`` alias; v1
    servers give ``trunc512``, the same 24 bytes in hex.
    """
    ga4gh = metadata.get("ga4gh")
    if isinstance(ga4gh, str) and ga4gh:
        return ga4gh.removeprefix("ga4gh:")
    for entry in metadata.get("aliases", []):
        alias = str(entry.get("alias", ""))
        if alias.startswith("ga4gh:SQ."):
            return alias.removeprefix("ga4gh:")
        if alias.startswith("SQ."):
            return alias
    trunc512 = metadata.get("trunc512")
    if isinstance(trunc512, str) and trunc512:
        return "SQ." + base64.urlsafe_b64encode(bytes.fromhex(trunc512)).decode("ascii").rstrip("=")
    return None


class RefgetProvider:
    """Sequences and refget lookups from a refget server; transcripts from ``transcripts``.

    Args:
        base_url: The server's refget root, the part before ``/sequence/``.
        transcripts: A provider answering ``get_transcript`` and
            ``get_symbol_accessions`` (and ``get_identifier_type``, if it has one).
            Without it those calls raise ``DataProviderError`` and identifier types
            come from the accession prefix.
        timeout: Seconds per request.
        transport: Replaces the HTTP layer; see :data:`Transport`.
    """

    def __init__(
        self,
        base_url: str,
        transcripts: TranscriptSource | None = None,
        *,
        timeout: float = 30.0,
        transport: Transport | None = None,
    ) -> None:
        self._base = base_url.rstrip("/")
        self._transcripts = transcripts
        self._fetch = transport or urllib_transport(timeout)
        # Per accession: the id the server answers to, and its metadata.
        self._known: dict[str, tuple[str, dict[str, typing.Any]]] = {}

    # --- the refget protocol ---

    def _get(self, path: str, accept: str) -> tuple[int, bytes]:
        return self._fetch(f"{self._base}/{path}", accept)

    def _metadata(self, sequence_id: str) -> dict[str, typing.Any] | None:
        """The metadata record for one id, or None if the server does not know it."""
        status, body = self._get(f"sequence/{urllib.parse.quote(sequence_id, safe=':.')}/metadata", _METADATA_ACCEPT)
        if status == 404:  # noqa: PLR2004 - HTTP status
            return None
        if status != 200:  # noqa: PLR2004 - HTTP status
            raise DataProviderError(f"refget metadata for {sequence_id}: HTTP {status}")
        record = json.loads(body)
        metadata = record.get("metadata", record)
        if not isinstance(metadata, dict):
            raise DataProviderError(f"refget metadata for {sequence_id} is not an object")
        return metadata

    def _resolve(self, ac: str) -> tuple[str, dict[str, typing.Any]]:
        """The id the server knows ``ac`` by, and its metadata; cached."""
        if ac in self._known:
            return self._known[ac]
        candidates = [f"ga4gh:{ac}", ac] if ac.startswith("SQ.") else [ac, f"refseq:{ac}"]
        for sequence_id in candidates:
            metadata = self._metadata(sequence_id)
            if metadata is not None:
                self._known[ac] = (sequence_id, metadata)
                return self._known[ac]
        raise DataProviderError(f"refget server at {self._base} knows no sequence {ac}")

    # --- DataProvider ---

    def get_seq(self, ac: str, start: int, end: int | None, kind: str | IdentifierType) -> str:  # noqa: ARG002
        """The bases of ``ac`` over ``[start, end)``; a range past the end returns what exists."""
        sequence_id, metadata = self._resolve(ac)
        length = int(metadata["length"])
        start = max(0, min(start, length))
        end = length if end is None else max(start, min(end, length))
        if start == end:
            return ""
        path = f"sequence/{urllib.parse.quote(sequence_id, safe=':.')}?start={start}&end={end}"
        status, body = self._get(path, _SEQUENCE_ACCEPT)
        if status != 200:  # noqa: PLR2004 - HTTP status
            raise DataProviderError(f"refget sequence {ac}[{start}:{end}]: HTTP {status}")
        return body.decode("ascii").strip()

    def get_refget_accession(self, ac: str) -> str | None:
        """The ``SQ.`` accession of ``ac`` from the server's metadata."""
        _, metadata = self._resolve(ac)
        return sq_digest(metadata)

    def get_accession_for_refget(self, refget: str) -> str | None:
        """An accession the server lists for a ``SQ.`` digest: RefSeq if it has one, else INSDC."""
        digest = refget.removeprefix("ga4gh:")
        metadata = self._metadata(f"ga4gh:{digest}")
        if metadata is None:
            metadata = self._metadata(digest)
        if metadata is None:
            return None
        aliases = metadata.get("aliases", [])
        names = [(str(a.get("naming_authority", "")), str(a.get("alias", ""))) for a in aliases]
        for _, alias in names:
            name = alias.split(":", 1)[1] if ":" in alias else alias
            if _REFSEQ.match(name):
                return name
        for authority, alias in names:
            if authority.lower() == "insdc":
                return alias.split(":", 1)[1] if ":" in alias else alias
        return None

    def get_transcript(self, transcript_ac: str, reference_ac: str | None) -> TranscriptData:
        if self._transcripts is None:
            raise DataProviderError("RefgetProvider serves sequences only; give it a transcript provider")
        return self._transcripts.get_transcript(transcript_ac, reference_ac)

    def get_symbol_accessions(
        self,
        symbol: str,
        source_kind: str,
        target_kind: str,
    ) -> list[tuple[str, str]] | list[tuple[IdentifierType, str]]:
        if self._transcripts is None:
            return []
        return self._transcripts.get_symbol_accessions(symbol, source_kind, target_kind)

    def get_identifier_type(self, identifier: str) -> str | IdentifierType:
        if self._transcripts is not None:
            return self._transcripts.get_identifier_type(identifier)
        if identifier.startswith(("NP_", "XP_", "AP_")):
            return IdentifierType.ProteinAccession
        if identifier.startswith(("NM_", "NR_", "XM_", "XR_")):
            return IdentifierType.TranscriptAccession
        if identifier.startswith(("NC_", "NT_", "NW_", "NG_", "AC_", "SQ.")):
            return IdentifierType.GenomicAccession
        return IdentifierType.GeneSymbol
