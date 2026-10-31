"""RefSeqDataProvider: a library module over a GFF3 and a FASTA, with no side effects on import."""

from __future__ import annotations

import pathlib
import subprocess
import sys

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


def provider(tmp_path: pathlib.Path) -> refseq.RefSeqDataProvider:
    """A provider over the synthetic chromosome and FWD."""
    fasta = tmp_path / "genome.fa"
    fasta.write_text(f">NC_000099.1\n{GENOME}\n")
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
