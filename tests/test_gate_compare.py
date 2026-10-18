"""weaver-gate's comparison of two validation outputs, and its report of the rows the two verdicts disagree on."""

from __future__ import annotations

import csv
import pathlib
import re
from collections.abc import Iterable, Sequence

import pytest

pytest.importorskip("weaver_data_provider", reason="weaver-gate needs hgvs-weaver-data, on Python 3.12 or later")

from weaver.cli import gate

FIELDS = ("variant_nuc", "variant_prot", "spdi", "rs_p", "rs_spdi", "rs_equiv", "rs_nuc_equiv")
FIELDS_WITHOUT_NUC_EQUIV = tuple(f for f in FIELDS if f != "rs_nuc_equiv")


def output_row(**fields: str) -> dict[str, str]:
    row = {
        "variant_nuc": "NM_TEST.1:c.4G>A",
        "variant_prot": "NP_TEST.1:p.Ala2Thr",
        "spdi": "NC_TEST.1:1013:G:A",
        "rs_p": "p.(Ala2Thr)",
        "rs_spdi": "NC_TEST.1:1013:G:A",
        "rs_equiv": "Identity",
        "rs_nuc_equiv": "Identity",
    }
    row.update(fields)
    return row


def write(path: pathlib.Path, rows: Iterable[dict[str, str]], fields: Sequence[str] = FIELDS) -> pathlib.Path:
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return path


def test_a_baseline_without_the_column_reports_it_absent(
    tmp_path: pathlib.Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    rows = [output_row(), output_row(variant_nuc="NM_TEST.1:c.5C>T", rs_nuc_equiv="Analogous")]
    baseline = write(tmp_path / "old.tsv", rows, FIELDS_WITHOUT_NUC_EQUIV)
    gate.compare(baseline, write(tmp_path / "new.tsv", rows))
    out = capsys.readouterr().out
    assert "rs_nuc_equiv: absent from old.tsv" in out
    assert not re.search(r"rs_nuc_equiv: \d+ rows moved", out)
    assert "rs_equiv: 0 rows moved" in out


def test_a_moved_verdict_is_counted_by_its_direction(
    tmp_path: pathlib.Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    baseline = write(tmp_path / "old.tsv", [output_row(rs_nuc_equiv="Identity")])
    gate.compare(baseline, write(tmp_path / "new.tsv", [output_row(rs_nuc_equiv="Different")]))
    out = capsys.readouterr().out
    assert "rs_nuc_equiv: 1 rows moved" in out
    assert re.search(r"1  Identity -> Different", out)


def test_rows_whose_two_verdicts_disagree_are_listed_by_pattern(
    tmp_path: pathlib.Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    rows = [
        output_row(variant_nuc="NM_TEST.1:c.4G>A", rs_equiv="Identity", rs_nuc_equiv="Identity"),
        output_row(variant_nuc="NM_TEST.1:c.5C>T", rs_equiv="Identity", rs_nuc_equiv="Analogous"),
        output_row(variant_nuc="NM_TEST.1:c.6T>A", rs_equiv="Identity", rs_nuc_equiv="Analogous"),
        output_row(variant_nuc="NM_TEST.1:c.7G>C", rs_equiv="Unknown", rs_nuc_equiv="Different"),
        output_row(variant_nuc="NM_TEST.1:c.8G>C", rs_equiv="Unknown", rs_nuc_equiv="ERR:no protein"),
    ]
    output = write(tmp_path / "new.tsv", rows)
    gate.compare(write(tmp_path / "old.tsv", rows), output)
    out = capsys.readouterr().out
    report = out[out.index("rs_equiv against rs_nuc_equiv") :]
    assert re.search(r"3 rows disagree; 1 judged by neither", report)
    assert re.search(r"2  rs_equiv Identity, rs_nuc_equiv Analogous", report)
    assert re.search(r"1  rs_equiv Unknown, rs_nuc_equiv Different", report)
    assert "NM_TEST.1:c.7G>C" in report
    assert "NM_TEST.1:c.4G>A" not in report
    assert "NM_TEST.1:c.8G>C" not in report
