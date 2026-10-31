"""RefSeqDataProvider: a library module over a GFF3 and a FASTA, with no side effects on import."""

from __future__ import annotations

import pathlib
import re
import subprocess
import sys

import pytest

import weaver
from weaver import refseq

# NC_000099.1, 60 bases. FWD's NM_000099.1 has exons g.11-25 and g.36-50 and a CDS g.14-25 + g.36-47:
# c.1 is g.14, c.12 the exon's last base g.25, c.13 the next exon's first base g.36 (a T).
GENOME = "TTCAGTCAGT" + "CCA" + "ATGGCTGAACGT" + "GTAAGTCCAG" + "TTCAAGCTGTAA" + "GCA" + "TCGATCGATC"
GFF3 = """\
##gff-version 3
NC_000099.1\tRefSeq\tgene\t11\t50\t.\t+\t.\tID=gene-FWD;gene=FWD
NC_000099.1\tRefSeq\tmRNA\t11\t50\t.\t+\t.\tID=rna-NM_000099.1;Parent=gene-FWD;gene=FWD;transcript_id=NM_000099.1
NC_000099.1\tRefSeq\texon\t11\t25\t.\t+\t.\tID=exon-1;Parent=rna-NM_000099.1;transcript_id=NM_000099.1
NC_000099.1\tRefSeq\texon\t36\t50\t.\t+\t.\tID=exon-2;Parent=rna-NM_000099.1;transcript_id=NM_000099.1
NC_000099.1\tRefSeq\tCDS\t14\t25\t.\t+\t0\tID=cds-1;Parent=rna-NM_000099.1;protein_id=NP_000099.1
NC_000099.1\tRefSeq\tCDS\t36\t47\t.\t+\t0\tID=cds-1;Parent=rna-NM_000099.1;protein_id=NP_000099.1
"""


@pytest.fixture(autouse=True)
def live_sequences(monkeypatch: pytest.MonkeyPatch) -> None:
    """Reads the FASTA, whatever sequence-cache mode the environment running the tests asks for."""
    monkeypatch.delenv("WEAVER_SEQ_MODE", raising=False)
    monkeypatch.delenv("WEAVER_SEQ_CACHE", raising=False)


def provider(tmp_path: pathlib.Path, gff3: str = GFF3) -> refseq.RefSeqDataProvider:
    """A provider over the synthetic chromosome and FWD, and an unannotated sequence under a custom name."""
    fasta = tmp_path / "genome.fa"
    fasta.write_text(f">NC_000099.1\n{GENOME}\n>AMPLICON1\nACGTACGTAC\n")
    gff = tmp_path / "annotation.gff3"
    gff.write_text(gff3)
    return refseq.RefSeqDataProvider(str(gff), str(fasta))


def run_python(code: str) -> subprocess.CompletedProcess[str]:
    """Runs `code` in a fresh interpreter, so what an import does is not hidden by what this one already did."""
    return subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=False)  # noqa: S603


def test_importing_the_provider_leaves_the_root_logger_unconfigured() -> None:
    result = run_python("import logging, weaver.refseq; print(logging.getLogger().handlers)")
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "[]"


def test_importing_the_provider_does_not_import_biocommons_hgvs() -> None:
    result = run_python("import sys, weaver.refseq; print('hgvs' in sys.modules)")
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "False"


def test_a_missing_pysam_is_an_import_error_a_caller_can_catch() -> None:
    result = run_python(
        "import sys\n"
        "sys.modules['pysam'] = None\n"
        "try:\n"
        "    import weaver.refseq\n"
        "except ImportError:\n"
        "    print('ImportError')\n",
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "ImportError"


def test_a_coding_variant_projects_across_the_intron_the_annotation_places(tmp_path: pathlib.Path) -> None:
    mapper = weaver.VariantMapper(provider(tmp_path))
    assert mapper.c_to_g(weaver.parse("NM_000099.1:c.13T>A")).format() == "NC_000099.1:g.36T>A"


@pytest.mark.parametrize(
    ("ac", "kind"),
    [
        ("NC_000098.1", weaver.IdentifierType.GenomicAccession),
        ("NM_000098.1", weaver.IdentifierType.TranscriptAccession),
        ("NP_000098.1", weaver.IdentifierType.ProteinAccession),
    ],
)
def test_a_sequence_the_provider_does_not_hold_is_refused_by_name(
    tmp_path: pathlib.Path,
    ac: str,
    kind: weaver.IdentifierType,
) -> None:
    with pytest.raises(weaver.DataProviderError, match=re.escape(f"{ac}: not in the reference data")):
        provider(tmp_path).get_seq(ac, 0, 5, kind)


@pytest.mark.parametrize("variant", ["NC_000098.1:g.10del", "NC_000098.1:g.10A>T"])
def test_a_variant_on_a_sequence_the_reference_lacks_fails_naming_the_sequence(
    tmp_path: pathlib.Path,
    variant: str,
) -> None:
    mapper = weaver.VariantMapper(provider(tmp_path))
    with pytest.raises(weaver.DataProviderError, match=re.escape("NC_000098.1: not in the reference data")):
        mapper.normalize_variant(weaver.parse(variant))


@pytest.mark.parametrize(
    ("start", "end", "bases"),
    [
        (55, 70, "CGATC"),  # the 5 bases before the end of the 60
        (60, 70, ""),
        (100, 110, ""),
    ],
)
def test_a_range_past_the_end_of_a_held_sequence_answers_the_bases_there(
    tmp_path: pathlib.Path,
    start: int,
    end: int,
    bases: str,
) -> None:
    assert provider(tmp_path).get_seq("NC_000099.1", start, end, weaver.IdentifierType.GenomicAccession) == bases


@pytest.mark.parametrize(
    ("identifier", "kind"),
    [
        ("NC_000099.1", "genomic_accession"),
        ("AMPLICON1", "genomic_accession"),
        ("NM_000099.1", "transcript_accession"),
        ("NP_000099.1", "protein_accession"),
        ("FWD", "gene_symbol"),
        ("NC_000098.1", "unknown"),
        ("NM_000098.1", "unknown"),
        ("AY286018", "unknown"),
    ],
)
def test_an_identifier_is_classified_by_what_the_provider_holds_under_it(
    tmp_path: pathlib.Path,
    identifier: str,
    kind: str,
) -> None:
    assert provider(tmp_path).get_identifier_type(identifier) == kind


# Two versions of REV's transcript, each on a different sequence: .1 on AMPLICON1, .2 on NC_000099.1.
TWO_VERSIONS = (
    GFF3
    + "AMPLICON1\tRefSeq\tmRNA\t1\t10\t.\t-\t.\tID=rna-NM_000097.1;gene=REV;transcript_id=NM_000097.1\n"
    + "AMPLICON1\tRefSeq\texon\t1\t10\t.\t-\t.\tParent=rna-NM_000097.1;transcript_id=NM_000097.1\n"
    + "NC_000099.1\tRefSeq\tmRNA\t51\t60\t.\t-\t.\tID=rna-NM_000097.2;gene=REV;transcript_id=NM_000097.2\n"
    + "NC_000099.1\tRefSeq\texon\t51\t60\t.\t-\t.\tParent=rna-NM_000097.2;transcript_id=NM_000097.2\n"
)


@pytest.mark.parametrize("ac", ["NM_000097.3", "NM_000097"])
def test_a_transcript_version_the_annotation_lacks_has_no_sequence(tmp_path: pathlib.Path, ac: str) -> None:
    with pytest.raises(weaver.DataProviderError, match=re.escape(f"{ac}: not in the reference data")):
        provider(tmp_path, TWO_VERSIONS).get_seq(ac, 0, 5, weaver.IdentifierType.TranscriptAccession)


@pytest.mark.parametrize("ac", ["NM_000097.3", "NM_000097"])
def test_a_transcript_version_the_annotation_lacks_has_no_model(tmp_path: pathlib.Path, ac: str) -> None:
    with pytest.raises(weaver.DataProviderError, match=re.escape(f"{ac}: not in the reference data")):
        provider(tmp_path, TWO_VERSIONS).get_transcript(ac, None)


@pytest.mark.parametrize(("ac", "reference"), [("NM_000097.1", "AMPLICON1"), ("NM_000097.2", "NC_000099.1")])
def test_each_transcript_version_is_placed_where_the_annotation_places_it(
    tmp_path: pathlib.Path,
    ac: str,
    reference: str,
) -> None:
    assert provider(tmp_path, TWO_VERSIONS).get_transcript(ac, None)["reference_accession"] == reference


def test_a_transcript_version_is_read_from_its_own_placement(tmp_path: pathlib.Path) -> None:
    # NM_000097.2 is g.51-60 of NC_000099.1, TCGATCGATC, read on the minus strand.
    rev = provider(tmp_path, TWO_VERSIONS).get_seq("NM_000097.2", 0, 10, weaver.IdentifierType.TranscriptAccession)
    assert rev == "GATCGATCGA"


@pytest.mark.parametrize(("identifier", "kind"), [("NM_000097.2", "transcript_accession"), ("NM_000097.3", "unknown")])
def test_an_identifier_is_a_transcript_exactly_when_its_version_is_annotated(
    tmp_path: pathlib.Path,
    identifier: str,
    kind: str,
) -> None:
    assert provider(tmp_path, TWO_VERSIONS).get_identifier_type(identifier) == kind


def test_a_transcript_asked_for_on_a_sequence_the_annotation_does_not_place_it_on_has_no_model(
    tmp_path: pathlib.Path,
) -> None:
    with pytest.raises(weaver.DataProviderError, match=re.escape("NM_000097.1: not placed on NC_000099.1")):
        provider(tmp_path, TWO_VERSIONS).get_transcript("NM_000097.1", "NC_000099.1")


@pytest.mark.parametrize(
    ("ac", "kind", "start", "end"),
    [
        ("NC_000099.1", weaver.IdentifierType.GenomicAccession, 10, 5),
        ("NC_000099.1", weaver.IdentifierType.GenomicAccession, -5, 5),
        ("NC_000099.1", weaver.IdentifierType.GenomicAccession, 0, -1),
        ("NM_000099.1", weaver.IdentifierType.TranscriptAccession, -5, 5),
        ("NP_000099.1", weaver.IdentifierType.ProteinAccession, 3, 1),
    ],
)
def test_a_range_that_is_not_one_is_refused_naming_the_sequence_and_the_range(
    tmp_path: pathlib.Path,
    ac: str,
    kind: weaver.IdentifierType,
    start: int,
    end: int,
) -> None:
    with pytest.raises(weaver.DataProviderError, match=re.escape(f"{ac}: [{start}, {end}) is not a range")):
        provider(tmp_path).get_seq(ac, start, end, kind)


def recorded(tmp_path: pathlib.Path, start: int, end: int | None) -> refseq.SequenceProxy:
    """A proxy replaying what a recording proxy read of NC_000099.1 over ``[start, end)``."""
    fasta = tmp_path / "genome.fa"
    fasta.write_text(f">NC_000099.1\n{GENOME}\n")
    cache = tmp_path / "cache" / "sequences.json"
    recording = refseq.SequenceProxy(str(fasta), str(cache), "record")
    recording.fetch("NC_000099.1", start, end)
    recording.save_cache()
    return refseq.SequenceProxy(str(fasta), str(cache), "replay")


def test_a_replayed_read_to_the_end_answers_what_the_recorded_one_did(tmp_path: pathlib.Path) -> None:
    assert recorded(tmp_path, 50, None).fetch("NC_000099.1", 50, None) == "TCGATCGATC"


def test_a_replayed_read_of_a_range_that_is_not_one_is_refused(tmp_path: pathlib.Path) -> None:
    with pytest.raises(weaver.DataProviderError, match=re.escape("NC_000099.1: [0, -1) is not a range")):
        recorded(tmp_path, 50, None).fetch("NC_000099.1", 0, -1)


@pytest.mark.parametrize(
    ("symbol", "source", "target", "accessions"),
    [
        ("FWD", "g", "c", [(weaver.IdentifierType.TranscriptAccession, "NM_000099.1")]),
        ("FWD", "g", "p", [(weaver.IdentifierType.ProteinAccession, "NP_000099.1")]),
        ("NM_000099.1", "c", "p", [(weaver.IdentifierType.ProteinAccession, "NP_000099.1")]),
        ("NP_000099.1", "p", "c", [(weaver.IdentifierType.TranscriptAccession, "NM_000099.1")]),
        ("REV", "g", "p", []),
        ("NM_000097.2", "c", "p", []),
        ("NOSUCH", "g", "c", []),
        ("NM_000099.1", "c", "g", []),
    ],
)
def test_a_symbol_maps_to_the_accessions_of_the_kind_asked_for(
    tmp_path: pathlib.Path,
    symbol: str,
    source: str,
    target: str,
    accessions: list[tuple[weaver.IdentifierType, str]],
) -> None:
    assert provider(tmp_path, TWO_VERSIONS).get_symbol_accessions(symbol, source, target) == accessions


def test_a_variant_on_a_gene_symbol_is_equivalent_to_the_same_change_on_its_transcript(tmp_path: pathlib.Path) -> None:
    rp = provider(tmp_path)
    mapper = weaver.VariantMapper(rp)
    assert mapper.equivalent(weaver.parse("FWD:c.13T>A"), weaver.parse("NM_000099.1:c.13T>A"), rp)
