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


def provider(tmp_path: pathlib.Path) -> refseq.RefSeqDataProvider:
    """A provider over the synthetic chromosome and FWD, and an unannotated sequence under a custom name."""
    fasta = tmp_path / "genome.fa"
    fasta.write_text(f">NC_000099.1\n{GENOME}\n>AMPLICON1\nACGTACGTAC\n")
    gff = tmp_path / "annotation.gff3"
    gff.write_text(GFF3)
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
