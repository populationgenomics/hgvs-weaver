"""weaver-gate's ClinVar protein report: the rows it selects from ClinVar's HGVS report, and its comparisons."""

from __future__ import annotations

import csv
import gzip
import pathlib
import re
from collections.abc import Iterable

import pytest

from weaver.cli import gate

HGVS4VARIATION_HEADER = (
    "##Report of HGVS expressions for VariationIDs and AlleleIDs\n"
    "#GeneSymbol:           the offical symbol from HGNC\n"
    "#Symbol\tGeneID\tVariationID\tAlleleID\tType\tAssembly\tNucleotideExpression\tNucleotideChange"
    "\tProteinExpression\tProteinChange\tUsedForNaming\tSubmitted\tOnRefSeqGene\n"
)


def hgvs4variation_line(nuc: str, prot: str, kind: str = "coding") -> str:
    change = nuc.partition(":")[2] if nuc != "-" else "-"
    protein_change = prot.partition(":")[2] if prot != "-" else "-"
    return f"FWD\t99\t1\t2\t{kind}\t-\t{nuc}\t{change}\t{prot}\t{protein_change}\tYes\tCurrent\tNo\n"


def write_hgvs4variation(path: pathlib.Path, lines: Iterable[str], header: str = HGVS4VARIATION_HEADER) -> pathlib.Path:
    with gzip.open(path, "wt") as fh:
        fh.write(header)
        fh.writelines(lines)
    return path


def write_report(path: pathlib.Path, rows: Iterable[tuple[str, str, str]]) -> pathlib.Path:
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(("variant_nuc", "variant_prot", "rs_p"))
        writer.writerows(rows)
    return path


def test_a_coding_indel_listed_twice_is_selected_once(tmp_path: pathlib.Path) -> None:
    source = write_hgvs4variation(
        tmp_path / "hgvs4variation.txt.gz",
        [
            hgvs4variation_line("NM_000099.1:c.10_12del", "NP_000099.1:p.Leu4del"),
            hgvs4variation_line("NM_000099.1:c.20dup", "NP_000099.1:p.Ser8fs"),
            hgvs4variation_line("NM_000099.1:c.10_12del", "NP_000099.1:p.Leu4del"),
        ],
    )
    assert gate.clinvar_proteins(source) == [
        ("NM_000099.1:c.10_12del", "NP_000099.1:p.Leu4del"),
        ("NM_000099.1:c.20dup", "NP_000099.1:p.Ser8fs"),
    ]


def test_an_allele_of_substitutions_is_selected(tmp_path: pathlib.Path) -> None:
    source = write_hgvs4variation(
        tmp_path / "hgvs4variation.txt.gz",
        [hgvs4variation_line("NM_000099.1:c.[4G>A;6A>T]", "NP_000099.1:p.[Ala2Thr;Lys2Asn]")],
    )
    assert gate.clinvar_proteins(source) == [("NM_000099.1:c.[4G>A;6A>T]", "NP_000099.1:p.[Ala2Thr;Lys2Asn]")]


@pytest.mark.parametrize(
    ("nuc", "prot", "kind"),
    [
        ("NM_000099.1:c.4G>A", "NP_000099.1:p.Ala2Thr", "coding"),
        ("NM_000099.1:c.10_12del", "-", "coding"),
        ("NC_000099.1:g.1010_1012del", "-", "genomic"),
        ("XM_000099.1:c.10_12del", "XP_000099.1:p.Leu4del", "coding"),
        ("NR_000099.1:n.10_12del", "-", "non-coding"),
    ],
)
def test_a_row_that_is_not_a_coding_indel_on_a_named_protein_is_left_out(
    tmp_path: pathlib.Path,
    nuc: str,
    prot: str,
    kind: str,
) -> None:
    source = write_hgvs4variation(tmp_path / "hgvs4variation.txt.gz", [hgvs4variation_line(nuc, prot, kind)])
    assert gate.clinvar_proteins(source) == []


def test_a_report_without_its_header_line_is_refused(tmp_path: pathlib.Path) -> None:
    source = write_hgvs4variation(
        tmp_path / "hgvs4variation.txt.gz",
        [hgvs4variation_line("NM_000099.1:c.10_12del", "NP_000099.1:p.Leu4del")],
        header="##Report of HGVS expressions\n",
    )
    with pytest.raises(ValueError, match="no #Symbol header"):
        gate.clinvar_proteins(source)


@pytest.mark.parametrize(
    ("clinvar", "predicted", "agree"),
    [
        ("NP_000099.1:p.Leu4del", "NP_000099.1:p.(Leu4del)", 1),
        ("NP_000099.1:p.Gly5*", "NP_000099.1:p.(Gly5Ter)", 1),
        ("NP_000099.1:p.Ser8fs", "NP_000099.1:p.(Ser8LysfsTer3)", 1),
        ("NP_000099.1:p.Ser8LysfsTer3", "NP_000099.1:p.(Ser8Lysfs*3)", 1),
        ("NP_000099.1:p.Ser8fs", "NP_000099.1:p.(Ala9GlyfsTer5)", 0),
        ("NP_000099.1:p.Ser8LysfsTer3", "NP_000099.1:p.(Ser8LysfsTer4)", 0),
        ("NP_000099.1:p.Leu4del", "NP_000099.1:p.(Leu5del)", 0),
    ],
)
def test_a_prediction_spells_clinvar_up_to_parentheses_the_stop_and_a_short_frameshift(
    tmp_path: pathlib.Path,
    capsys: pytest.CaptureFixture[str],
    clinvar: str,
    predicted: str,
    agree: int,
) -> None:
    gate.compare_proteins(None, write_report(tmp_path / "new.tsv", [("NM_000099.1:c.20dup", clinvar, predicted)]))
    assert f"new.tsv: 1 rows; 1 answered, {agree} of them spelled as ClinVar spells them" in capsys.readouterr().out


def test_a_prediction_that_failed_is_not_counted_as_answered(
    tmp_path: pathlib.Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    output = write_report(
        tmp_path / "new.tsv",
        [
            ("NM_000099.1:c.10_12del", "NP_000099.1:p.Leu4del", "NP_000099.1:p.(Leu4del)"),
            ("NM_000099.1:c.30del", "NP_000099.1:p.Val10fs", "ERR:no CDS on NM_000099.1"),
            ("NM_000099.1:c.33del", "NP_000099.1:p.Ala11fs", "PANIC"),
        ],
    )
    gate.compare_proteins(None, output)
    assert "new.tsv: 3 rows; 1 answered, 1 of them spelled as ClinVar spells them" in capsys.readouterr().out


def test_agreement_with_clinvar_is_counted_by_the_kind_of_clinvar_edit(
    tmp_path: pathlib.Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    output = write_report(
        tmp_path / "new.tsv",
        [
            ("NM_000099.1:c.10_12del", "NP_000099.1:p.Leu4del", "NP_000099.1:p.(Leu4del)"),
            ("NM_000099.1:c.20dup", "NP_000099.1:p.Ser8fs", "NP_000099.1:p.(Ala9GlyfsTer5)"),
            ("NM_000099.1:c.23dup", "NP_000099.1:p.Ala9fs", "NP_000099.1:p.(Ala9Ter)"),
        ],
    )
    gate.compare_proteins(None, output)
    out = capsys.readouterr().out
    assert "ClinVar fs: 2 answered, 0 of them spelled" in out
    assert "ClinVar Ter" not in out
    assert "ClinVar del: 1 answered, 1 of them spelled" in out


def test_a_moved_prediction_is_counted_by_the_kind_of_edit_on_each_side(
    tmp_path: pathlib.Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    unchanged = ("NM_000099.1:c.20dup", "NP_000099.1:p.Ser8fs", "NP_000099.1:p.(Ser8LysfsTer3)")
    baseline = write_report(
        tmp_path / "old.tsv",
        [("NM_000099.1:c.10_15del", "NP_000099.1:p.Leu4Ter", "NP_000099.1:p.(Leu4_Gly5del)"), unchanged],
    )
    output = write_report(
        tmp_path / "new.tsv",
        [("NM_000099.1:c.10_15del", "NP_000099.1:p.Leu4Ter", "NP_000099.1:p.(Leu4Ter)"), unchanged],
    )
    gate.compare_proteins(baseline, output)
    out = capsys.readouterr().out
    assert "rs_p: 1 rows moved" in out
    assert re.search(
        r"1  del -> Ter\n +e\.g\. NM_000099\.1:c\.10_15del NP_000099\.1:p\.Leu4Ter: "
        r"'NP_000099\.1:p\.\(Leu4_Gly5del\)' -> 'NP_000099\.1:p\.\(Leu4Ter\)'",
        out,
    )
    assert "on the shared rows, old.tsv: 2 answered, 1 of them spelled" in out
    assert "on the shared rows, new.tsv: 2 answered, 2 of them spelled" in out


def test_a_row_one_report_lacks_is_not_compared(
    tmp_path: pathlib.Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    shared = ("NM_000099.1:c.10_12del", "NP_000099.1:p.Leu4del", "NP_000099.1:p.(Leu4del)")
    baseline = write_report(
        tmp_path / "old.tsv",
        [shared, ("NM_000099.1:c.20dup", "NP_000099.1:p.Ser8fs", "NP_000099.1:p.(Ser8LysfsTer3)")],
    )
    gate.compare_proteins(baseline, write_report(tmp_path / "new.tsv", [shared]))
    out = capsys.readouterr().out
    assert "old.tsv: 2 rows; 1 shared" in out
    assert "rs_p: 0 rows moved" in out


def test_a_file_without_the_report_columns_is_refused(tmp_path: pathlib.Path) -> None:
    output = tmp_path / "new.tsv"
    output.write_text("variant_nuc\tvariant_prot\nNM_000099.1:c.10_12del\tNP_000099.1:p.Leu4del\n")
    with pytest.raises(ValueError, match="no column rs_p"):
        gate.compare_proteins(None, output)
