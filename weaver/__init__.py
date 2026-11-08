from __future__ import annotations

from collections.abc import Sequence
from typing import Protocol, TypedDict

from weaver._weaver import (  # type: ignore[attr-defined]
    AccessionStatus,
    AccessionTag,
    AmbiguousDescriptionError,
    AmbiguousResidueError,
    CigarError,
    CisMemberError,
    ConsistentVariant,
    DataProviderError,
    EquivalenceLevel,
    HGVSError,
    IdentifierType,
    MixedResidueCodesError,
    NoProteinStatementError,
    NotAnExonBorderError,
    OffsetPastIntronError,
    ParseError,
    PositionOutsideSequenceError,
    PositionPastEndError,
    ReferenceMismatchError,
    ReferenceUse,
    UnplacedError,
    UnsupportedOperationError,
    ValidationError,
    Variant,
    VariantMapper,
    WrongReferenceKindError,
    XTreatment,
    parse,
)

__all__ = [
    "AccessionCatalogue",
    "AccessionStatus",
    "AccessionTag",
    "AmbiguousDescriptionError",
    "AmbiguousResidueError",
    "CigarError",
    "CisMemberError",
    "ConsistentVariant",
    "DataProvider",
    "DataProviderError",
    "EquivalenceLevel",
    "ExonData",
    "HGVSError",
    "IdentifierType",
    "MixedResidueCodesError",
    "NoProteinStatementError",
    "NotAnExonBorderError",
    "OffsetPastIntronError",
    "ParseError",
    "PositionOutsideSequenceError",
    "PositionPastEndError",
    "ReferenceMismatchError",
    "ReferenceUse",
    "Refget",
    "TranscriptData",
    "TranscriptSearch",
    "UnplacedError",
    "UnsupportedOperationError",
    "ValidationError",
    "Variant",
    "VariantMapper",
    "WrongReferenceKindError",
    "XTreatment",
    "parse",
]


class ExonData(TypedDict):
    """Represents an exon's coordinates and alignment.

    Coordinates are 0-based:
    - transcript_start: inclusive start index in transcript.
    - transcript_end: exclusive end index in transcript.
    - reference_start: inclusive start index on the genome.
    - reference_end: inclusive end index on the genome.
    """

    transcript_start: int
    transcript_end: int
    reference_start: int
    reference_end: int
    alt_strand: int  # 1 for plus, -1 for minus
    cigar: str  # Extended CIGAR string (e.g., "100=")


class _TranscriptRecord(TypedDict):
    """The keys every transcript model has."""

    ac: str
    gene: str
    cds_start_index: int | None
    cds_end_index: int | None


class TranscriptData(_TranscriptRecord, total=False):
    """Represents a full transcript model.

    Coordinates are 0-based:
    - cds_start_index: inclusive index of the first base of the start codon, or
      of the codon c.1 is numbered from where cds_start_initiates is False.
    - cds_end_index: inclusive index of the last base of the stop codon.

    A CDS that runs off an end of the record has no start or stop codon there;
    say so with cds_start_open or cds_end_open, and weaver refuses the positions
    numbered from that end rather than number them from a codon the record lacks.
    Both default to False.

    cds_start_initiates, default True, says translation initiates at the codon
    cds_start_index names, so residue 1 is Met whatever the codon. A coding
    sequence numbered from a first codon translation does not initiate at, as
    MAVE-HGVS numbers a target that is part of a protein, sets it False: codon 1
    is then read by the standard table, as any other codon is. It cannot be
    True on a CDS whose start is open.

    strand, reference_accession and exons place the transcript on a genome, and
    come together or not at all. A record no assembly places, such as a GenBank
    mRNA, leaves all three out: weaver reads it on its own sequence, and refuses
    with UnsupportedOperationError whatever needs the genome, a projection, an
    intronic position or a genomic form. With no genome to answer to, a change
    on it is normalised across where its exon junctions would be.
    """

    cds_start_open: bool
    cds_end_open: bool
    cds_start_initiates: bool
    strand: int  # 1 or -1
    reference_accession: str  # The genome's accession (e.g., NC_000001.11)
    exons: list[ExonData]


class TranscriptSearch(Protocol):
    """Optional interface for regional discovery."""

    def get_transcripts_for_region(self, chrom: str, start: int, end: int) -> list[str]:
        """Return list of transcript accessions overlapping the given genomic region."""
        ...


class AccessionCatalogue(Protocol):
    """Optional lookup passed to ``VariantMapper.consistent_variants(var, catalogue)``.

    What a store holds under an accession: the versions of an unversioned one, and the
    status and tags of a versioned one. Without it, ``consistent_variants`` cannot read an
    unversioned accession, and gives no status or tags.
    """

    def get_accession_versions(self, ac: str) -> list[str]:
        """Return every versioned accession held under the unversioned ``ac``, newest first.

        Each is ``ac`` followed by ``.`` and its version; a record with no version of its
        own, such as a UniProt accession, is listed as ``ac`` itself, alone. Only the store knows
        what newest means. Empty when none is held.
        """
        ...

    def get_accession_info(
        self,
        ac: str,
    ) -> tuple[str | AccessionStatus, Sequence[str | AccessionTag]] | None:
        """Return the status and tags of the versioned ``ac``, or None if nothing is held under it.

        The status is ``'current'``, ``'superseded'`` or ``'suppressed'``, or an
        AccessionStatus; each tag is ``'mane_select'``, ``'mane_plus_clinical'``,
        ``'refseq_select'`` or ``'ensembl_canonical'``, or an AccessionTag. A protein
        carries the tags of the transcripts that encode it.
        """
        ...


class Refget(Protocol):
    """Optional lookup passed to ``VariantMapper(provider, refget=...)``.

    Without one, weaver computes refget accessions from the whole sequence (slow
    for a chromosome) and ``from_vrs`` needs the accession passed by the caller.
    ``weaver.refget.RefgetProvider`` implements this as well as DataProvider.
    """

    def get_refget_accession(self, ac: str) -> str | None:
        """Return the refget accession ("SQ." + sha512t24u of the sequence) for ``ac``, or None."""
        ...

    def get_accession_for_refget(self, refget: str) -> str | None:
        """Return the accession of the sequence whose refget accession is ``refget``, or None."""
        ...


class DataProvider(Protocol):
    """Required interface for the object passed to VariantMapper."""

    def get_transcript(self, transcript_ac: str, reference_ac: str | None) -> TranscriptData:
        """Retrieve transcript model for the given accession.

        If reference_ac is provided, returns the placement on that genome.
        A record placed on no genome is returned unplaced whatever reference_ac is, so
        weaver can say the placement is missing rather than report the lookup failing.
        """
        ...

    def get_seq(self, ac: str, start: int, end: int | None, kind: str | IdentifierType) -> str:
        """Fetch the bases of ``ac`` in the 0-based half-open range ``[start, end)``.

        ``end`` is ``None`` when the whole sequence from ``start`` is wanted; a plain
        ``seq[start:end]`` already does the right thing in that case. A range that
        extends past the end of the sequence must return the bases that exist rather
        than raising or returning an empty string. ``kind`` should be an IdentifierType.
        """
        ...

    def get_symbol_accessions(
        self,
        symbol: str,
        source_kind: str,
        target_kind: str,
    ) -> list[tuple[str, str]] | list[tuple[IdentifierType, str]]:
        """Look up an identifier's counterparts in another namespace.

        ``source_kind`` and ``target_kind`` are ``'g'``, ``'c'``, ``'p'`` or ``'gene'``:
        ``(transcript, 'c', 'p')`` asks for a transcript's protein and ``(symbol, 'gene', 'p')``
        for a gene's proteins. Returns a list of tuples (identifier_type, accession).
        identifier_type should be one of 'genomic_accession', 'transcript_accession',
        'protein_accession', 'gene_symbol', or a member of the IdentifierType enum.
        """
        ...

    def get_identifier_type(self, identifier: str) -> str | IdentifierType:
        """Identify what kind of identifier a string is.

        Should return one of:
        - 'genomic_accession'
        - 'transcript_accession'
        - 'protein_accession'
        - 'gene_symbol'
        - 'unknown'
        Or the equivalent IdentifierType enum value.
        """
        ...
