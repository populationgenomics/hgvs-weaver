"""weaver-gate: the ClinVar validation from NCBI's published files to a row-by-row comparison, in one command.

The validation needs a GRCh38 store and genome built by hgvs-weaver-data's `weaver-data-build`: the current
annotation release's shard with NCBI's historical set of retired transcript versions stacked under it, so
that the variants ClinVar names on retired versions project too. This fetches the inputs, builds what is
missing, runs `weaver-validate` over the result and compares the output with a baseline run, row by row.
Every step is skipped when its output is already there, so a rerun after a weaver change costs only the
validation itself.

    uv run --group gate weaver-gate --data data/gate --variants data/clinvar_variants_100k.tsv \\
        --baseline data/validation_release_0.7.0.tsv --output data/validation_<label>.tsv

The downloads are a few gigabytes, the Entrez status fetch about twenty minutes, and the builds about
fifteen; they happen once. The gate dependency group installs what this needs, Python 3.12 or later.

With --clinvar-proteins, the gate also predicts the protein of every coding indel in ClinVar's HGVS report
over the same store, and compares the predictions with ClinVar's spelling and, given
--clinvar-proteins-baseline, with an earlier report's, row by row.
"""

from __future__ import annotations

import argparse
import collections
import concurrent.futures
import csv
import dataclasses
import datetime
import functools
import gzip
import hashlib
import logging
import pathlib
import re
import sys
import urllib.request
from collections.abc import Callable, Iterable, Mapping, Sequence

import tqdm

from weaver.cli import validate

NCBI = "https://ftp.ncbi.nlm.nih.gov"
HGNC_URL = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
# ClinVar's every HGVS expression, nucleotide beside protein; republished weekly.
CLINVAR_HGVS_URL = f"{NCBI}/pub/clinvar/tab_delimited/hgvs4variation.txt.gz"
PROTEIN_FIELDS = ("variant_nuc", "variant_prot", "rs_p")
# A frameshift's new residue and stop, which its short form p.Arg97fs leaves out.
_FRAMESHIFT = re.compile(r"(\d+)(?:[A-Z][a-z]{2})?fs(?:Ter(?:\d+|\?))?$")
# The SPDI column against ClinVar's; the equivalence verdicts as they stand, biocommons' answer judged by weaver's
# equivalence among them; weaver's protein text against its own in the baseline, a moved row counted by the kind of
# edit on each side.
COMPARED = (("rs_spdi", "spdi"), ("rs_equiv", None), ("rs_nuc_equiv", None), ("ref_equiv", None), ("rs_p", None))
VERDICTS = frozenset({"rs_equiv", "rs_nuc_equiv", "ref_equiv"})
PROTEINS = frozenset({"rs_p"})
JUDGED = frozenset({"Identity", "Analogous", "Different"})


@dataclasses.dataclass(frozen=True)
class Release:
    """One NCBI annotation release of GRCh38, and the historical set published beside it."""

    release: str  # "RS_2024_08"
    accession: str  # "GCF_000001405.40"
    assembly: str  # "GRCh38.p14"
    anchor: str  # the release the historical set is anchored on, "RS_2023_03"
    mane: str  # the MANE summary's file name under refseq/MANE/MANE_human/current

    @property
    def directory(self) -> str:
        return f"{NCBI}/genomes/all/annotation_releases/9606/{self.accession}-{self.release}"

    def file(self, suffix: str) -> str:
        return f"{self.accession}_{self.assembly}_{suffix}"

    def historical(self, suffix: str) -> str:
        return f"{self.accession}-{self.anchor}_{suffix}"


def download(url: str, dest: pathlib.Path, *, md5: str | None = None) -> pathlib.Path:
    """Fetch `url` to `dest` unless it is already there and whole: matching `md5`, or the server's length."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        if md5 is None:
            with urllib.request.urlopen(urllib.request.Request(url, method="HEAD"), timeout=60) as head:  # noqa: S310
                length = head.headers.get("Content-Length")
            if length is None or int(length) == dest.stat().st_size:
                return dest
        elif _md5(dest) == md5:
            return dest
        print(f"{dest.name}: present but not whole; fetching again", file=sys.stderr)
    print(f"fetching {url}", file=sys.stderr)
    part = dest.with_name(dest.name + ".part")
    with urllib.request.urlopen(url, timeout=600) as response, part.open("wb") as out:  # noqa: S310
        while chunk := response.read(1 << 22):
            out.write(chunk)
    if md5 is not None and (found := _md5(part)) != md5:
        part.unlink()
        raise RuntimeError(f"{url}: md5 {found}, NCBI's checksums say {md5}")
    part.replace(dest)
    return dest


def _md5(path: pathlib.Path) -> str:
    digest = hashlib.md5()  # noqa: S324 — NCBI publishes md5s, so md5 is what verifies a download
    with path.open("rb") as fh:
        while chunk := fh.read(1 << 22):
            digest.update(chunk)
    return digest.hexdigest()


def _checksums(release: Release, data: pathlib.Path) -> dict[str, str]:
    """NCBI's md5 for each file of the release, by file name."""
    path = download(f"{release.directory}/md5checksums.txt", data / "md5checksums.txt")
    out = {}
    for line in path.read_text().splitlines():
        md5, _, name = line.partition("  ")
        out[name.rsplit("/", 1)[-1]] = md5
    return out


@dataclasses.dataclass(frozen=True)
class Inputs:
    """The builder's inputs on local disk."""

    annotation: pathlib.Path
    transcripts: pathlib.Path
    proteins: pathlib.Path
    records: pathlib.Path
    alignments: tuple[pathlib.Path, ...]
    genome_fasta: pathlib.Path
    historical_annotation: pathlib.Path
    historical_records: pathlib.Path
    historical_alignments: pathlib.Path
    hgnc: pathlib.Path
    mane: pathlib.Path


def fetch(release: Release, data: pathlib.Path) -> Inputs:
    """Every input the builds need, downloaded unless present and whole."""
    md5s = _checksums(release, data)

    def release_file(suffix: str, subdirectory: str = "") -> pathlib.Path:
        name = release.file(suffix)
        where = f"{release.directory}/{subdirectory}/{name}" if subdirectory else f"{release.directory}/{name}"
        return download(where, data / name, md5=md5s.get(name))

    def historical_file(suffix: str) -> pathlib.Path:
        name = release.historical(suffix)
        return download(f"{release.directory}/RefSeq_historical_alignments/{name}", data / name)

    alignments = []
    for kind in ("knownrefseq", "modelrefseq"):
        alignments.append(release_file(f"{kind}_alns.bam", "RefSeq_transcripts_alignments"))
        release_file(f"{kind}_alns.bam.bai", "RefSeq_transcripts_alignments")
    historical_alignments = historical_file("knownrefseq_alns.bam")
    historical_file("knownrefseq_alns.bam.bai")
    return Inputs(
        annotation=release_file("genomic.gff.gz"),
        transcripts=release_file("rna.fna.gz"),
        proteins=release_file("protein.faa.gz"),
        records=release_file("rna.gbff.gz"),
        alignments=tuple(alignments),
        genome_fasta=release_file("genomic.fna.gz"),
        historical_annotation=historical_file("genomic.gff.gz"),
        historical_records=historical_file("knownrefseq_rna.gbff.gz"),
        historical_alignments=historical_alignments,
        hgnc=download(HGNC_URL, data / "hgnc_complete_set.txt"),
        mane=download(f"{NCBI}/refseq/MANE/MANE_human/current/{release.mane}", data / release.mane),
    )


def _shard(shards: pathlib.Path, prefix: str) -> pathlib.Path | None:
    """The one shard under `prefix` already cut, if any; a shard's record beside it is not a shard."""
    found = [p for p in shards.glob(f"{prefix}*.bagz") if not p.name.endswith(".shard.bagz")]
    if len(found) > 1:
        raise RuntimeError(f"{shards}: more than one shard under {prefix}: {sorted(p.name for p in found)}")
    return found[0] if found else None


def build(release: Release, inputs: Inputs, data: pathlib.Path) -> tuple[pathlib.Path, pathlib.Path]:
    """The store and genome under `data`, cut by `weaver-data-build` where not already there."""
    try:
        from weaver_data_provider import genome as wd_genome  # noqa: PLC0415
        from weaver_data_provider import store as wd_store  # noqa: PLC0415
        from weaver_data_provider.build import cli as wd_build  # noqa: PLC0415
    except ImportError:
        print(
            "Error: the gate needs 'hgvs-weaver-data[build]' on Python 3.12 or later: uv run --group gate",
            file=sys.stderr,
        )
        sys.exit(1)
    shards, store, genome = data / "shards", data / "store", data / "genome"
    status = data / f"{release.anchor}_status.tsv"
    wd_build.main(["status", "--records", str(inputs.historical_records), "--out", str(status)])

    historical = _shard(shards, f"01-{release.anchor}-historical-")
    if historical is None:
        wd_build.main([
            "historical", "--assembly", "GRCh38", "--release", f"{release.anchor}-historical",
            "--annotation", str(inputs.historical_annotation), "--records", str(inputs.historical_records),
            "--alignments", str(inputs.historical_alignments), "--status", str(status),
            "--hgnc", str(inputs.hgnc), "--mane", str(inputs.mane), "--shards", str(shards), "--prefix", "01-",
        ])  # fmt: skip
        historical = _shard(shards, f"01-{release.anchor}-historical-")
    current = _shard(shards, f"02-{release.release}-")
    if current is None:
        alignments = [arg for path in inputs.alignments for arg in ("--alignments", str(path))]
        wd_build.main([
            "refseq", "--assembly", "GRCh38", "--release", release.release,
            "--annotation", str(inputs.annotation), "--transcripts", str(inputs.transcripts),
            "--proteins", str(inputs.proteins), "--records", str(inputs.records), *alignments,
            "--hgnc", str(inputs.hgnc), "--mane", str(inputs.mane), "--shards", str(shards), "--prefix", "02-",
        ])  # fmt: skip
        current = _shard(shards, f"02-{release.release}-")
    if historical is None or current is None:
        raise RuntimeError(f"{shards}: a build wrote no shard")
    if not (store / wd_store.MANIFEST).exists():
        wd_build.main(["index", "--assembly", "GRCh38", "--out", str(store), str(historical), str(current)])
    if not (genome / wd_genome.CATALOGUE).exists():
        wd_build.main(["genome", "--assembly", "GRCh38", "--fasta", str(inputs.genome_fasta), "--out", str(genome)])
    return store, genome


def _load(path: pathlib.Path) -> tuple[list[str], dict[str, dict[str, str]]]:
    """A validation output's columns, and its rows by nucleotide description."""
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = {row["variant_nuc"]: row for row in reader}
        return list(reader.fieldnames or []), rows


def _protein_kind(description: str) -> str:
    """The kind of edit a protein description writes, by which a moved protein is counted.

    One of del, dup, ins, delins, fs, ext, Ter (a nonsense change), = (no change), sub (another residue),
    ? (unknown), 0 (no protein) or allele; error and blank for a cell that holds no description. Read from the
    text, with or without an accession or parentheses, so that a description weaver cannot parse has a kind too.
    """
    if not description:
        return "blank"
    if description.startswith("ERR") or description == "PANIC":
        return "error"
    body = description.rpartition(":")[2].removeprefix("p.")
    if body.startswith("["):
        return "allele"
    body = body.removeprefix("(").removesuffix(")")
    if body in {"0", "0?"}:
        return "0"
    for edit in ("fs", "ext", "delins", "del", "dup", "ins"):
        if edit in body:
            return edit
    if body.endswith("="):
        return "="
    if body.endswith(("Ter", "*")):
        return "Ter"
    if body.endswith("?"):
        return "?"
    return "sub"


def _state(row: Mapping[str, str], column: str, truth: str | None) -> str:
    """The kind of answer a cell holds, by which a moved row is counted; a protein by the kind of its edit."""
    value = row[column]
    if not value:
        return "blank"
    if value.startswith("ERR"):
        return "error"
    if column in VERDICTS:
        return value
    if column in PROTEINS:
        return _protein_kind(value)
    if truth is None:
        return "answered"
    return "match" if value == row[truth] else "mismatch"


def compare(baseline: pathlib.Path, output: pathlib.Path) -> None:
    """Print, per compared column, the rows that moved between two validation outputs and in which direction.

    A protein whose text moved is counted by the kind of edit each run writes it as (del -> Ter). A column one of
    the two outputs lacks is reported absent rather than compared. Then the output's rows on which weaver's two
    verdicts against ClinVar's protein disagree: through its own c_to_p, and directly.
    """
    (before_fields, before), (after_fields, after) = _load(baseline), _load(output)
    shared = before.keys() & after.keys()
    print(f"\n{baseline.name}: {len(before)} rows; {output.name}: {len(after)} rows; {len(shared)} shared")

    for column, truth in COMPARED:
        absent = [
            path.name for path, fields in ((baseline, before_fields), (output, after_fields)) if column not in fields
        ]
        if absent:
            print(f"{column}: absent from {' and '.join(absent)}; not compared")
            continue
        moved = _moved(before, after, shared, column, functools.partial(_state, column=column, truth=truth))
        line = f"{column}: {moved.total()} rows moved"
        if truth is not None:
            pct = [
                100 * sum(_state(r, column, truth) == "match" for r in rows.values()) / len(rows)
                for rows in (before, after)
            ]
            line += f"; identity {pct[0]:.3f}% -> {pct[1]:.3f}%"
        print(line)
        moved.print()

    if "rs_nuc_equiv" in after_fields:
        _print_disagreements(after)


@dataclasses.dataclass
class _Moves:
    """The rows whose cell moved between two runs, counted by direction (the cell's state before and after)."""

    counts: collections.Counter[str] = dataclasses.field(default_factory=collections.Counter)
    examples: collections.defaultdict[str, list[tuple[str, str, str]]] = dataclasses.field(
        default_factory=lambda: collections.defaultdict(list),
    )

    def total(self) -> int:
        return self.counts.total()

    def print(self) -> None:
        for direction, n in self.counts.most_common():
            print(f"    {n:6d}  {direction}")
            for key, was, now in self.examples.get(direction, []):
                print(f"        e.g. {key}: {was!r} -> {now!r}")


def _moved(
    before: Mapping[str, Mapping[str, str]],
    after: Mapping[str, Mapping[str, str]],
    keys: Iterable[str],
    column: str,
    state: Callable[[Mapping[str, str]], str],
) -> _Moves:
    """The rows of `keys` whose `column` differs between the two runs, with up to three examples a direction.

    A move from one error to another gets no example: an error's text says nothing about what weaver writes.
    """
    moves = _Moves()
    for key in sorted(keys):
        a, b = before[key], after[key]
        if a[column] == b[column]:
            continue
        kinds = (state(a), state(b))
        direction = f"{kinds[0]} -> {kinds[1]}"
        moves.counts[direction] += 1
        if len(moves.examples[direction]) < 3 and kinds != ("error", "error"):
            moves.examples[direction].append((key, a[column][:60], b[column][:60]))
    return moves


def _print_disagreements(rows: dict[str, dict[str, str]]) -> None:
    """The rows whose verdict through c_to_p (rs_equiv) differs from the direct one (rs_nuc_equiv), by pattern.

    Rows both columns judge come first, a contradiction (equivalent one way, Different the other) apart from a
    difference of level; then rows only one column judges. A row with no ClinVar protein, or that neither judges,
    says nothing about transitivity and is only counted.
    """
    groups: dict[str, collections.Counter[tuple[str, str]]] = {
        "contradict": collections.Counter(),
        "differ in level": collections.Counter(),
        "judged by one column only": collections.Counter(),
    }
    examples: dict[tuple[str, str], list[dict[str, str]]] = collections.defaultdict(list)
    no_protein = unjudged = 0
    for key in sorted(rows):
        row = rows[key]
        pattern = (_state(row, "rs_equiv", None), _state(row, "rs_nuc_equiv", None))
        judged = [verdict in JUDGED for verdict in pattern]
        if pattern[1] == "blank":
            no_protein += 1
            continue
        if not any(judged):
            unjudged += 1
            continue
        if pattern[0] == pattern[1]:
            continue
        if not all(judged):
            group = "judged by one column only"
        elif "Different" in pattern:
            group = "contradict"
        else:
            group = "differ in level"
        groups[group][pattern] += 1
        if len(examples[pattern]) < 3:
            examples[pattern].append(row)
    print(f"rs_equiv against rs_nuc_equiv: {no_protein} rows with no ClinVar protein, {unjudged} judged by neither")
    for group, patterns in groups.items():
        print(f"  {sum(patterns.values())} rows {group}")
        for pattern, n in patterns.most_common():
            print(f"    {n:6d}  rs_equiv {pattern[0]}, rs_nuc_equiv {pattern[1]}")
            for row in examples[pattern]:
                print(f"        e.g. {row['variant_nuc']}  ClinVar {row['variant_prot']}  weaver {row['rs_p'][:60]!r}")


def clinvar_proteins(hgvs4variation: pathlib.Path) -> list[tuple[str, str]]:
    """The coding indels in ClinVar's HGVS report, each an NM_ c. description and the NP_ protein ClinVar names.

    Every nucleotide change but a single substitution: the 100k validation already holds thousands of
    substitutions, which weaver writes one way, and few of the deletions, duplications, insertions and delins
    whose protein spelling weaver chooses. An allele, c.[...], is kept whatever its members. A pair ClinVar lists
    more than once is returned once, in the order first listed.
    """
    pairs: dict[tuple[str, str], None] = {}
    with gzip.open(hgvs4variation, "rt", newline="") as fh:
        header: list[str] | None = None
        for line in fh:
            if line.startswith("#"):
                if line.startswith("#Symbol\t"):
                    header = line[1:].rstrip("\n").split("\t")
                continue
            if header is None:
                raise ValueError(f"{hgvs4variation}: no #Symbol header line before the first row")
            row = dict(zip(header, line.rstrip("\n").split("\t"), strict=True))
            nuc, change, prot = row["NucleotideExpression"], row["NucleotideChange"], row["ProteinExpression"]
            if not (nuc.startswith("NM_") and ":c." in nuc and prot.startswith("NP_")):
                continue
            if ">" in change and "[" not in change:
                continue
            pairs[nuc, prot] = None
    return list(pairs)


def predict_proteins(
    pairs: Sequence[tuple[str, str]],
    output: pathlib.Path,
    *,
    store: pathlib.Path,
    genome: pathlib.Path,
    workers: int,
) -> None:
    """weaver's protein for each pair, on ClinVar's protein, over the validation's worker pool, written to `output`.

    The report is written beside `output` and moved into place whole, so an interrupted run leaves no report that
    could later pass for a baseline.
    """
    part = output.with_name(output.name + ".part")
    with (
        part.open("w", newline="") as fh,
        concurrent.futures.ProcessPoolExecutor(
            max_workers=workers,
            initializer=validate.init_worker,
            initargs=(str(store), str(genome)),
        ) as pool,
    ):
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(PROTEIN_FIELDS)
        predictions = pool.map(validate.predict_protein, pairs, chunksize=1000)
        for (nuc, prot), predicted in tqdm.tqdm(zip(pairs, predictions, strict=True), total=len(pairs)):
            writer.writerow((nuc, prot, predicted))
    part.replace(output)


def _spells(predicted: str, clinvar: str) -> bool:
    """Whether weaver's protein spells ClinVar's, parentheses and a stop written * or Ter aside.

    A frameshift ClinVar writes short, p.Arg97fs, is compared in that form: HGVS allows it, and ClinVar writes
    nearly every frameshift so, where weaver writes the new residue and the stop, p.(Arg97ProfsTer23).
    """
    predicted, clinvar = (d.replace("(", "").replace(")", "").replace("*", "Ter") for d in (predicted, clinvar))
    if clinvar.endswith("fs"):
        predicted, clinvar = (_FRAMESHIFT.sub(r"\1fs", d) for d in (predicted, clinvar))
    return predicted == clinvar


def _load_proteins(path: pathlib.Path) -> dict[str, dict[str, str]]:
    """A protein report's rows, by nucleotide description and ClinVar's protein."""
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = set(PROTEIN_FIELDS) - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path}: not a protein report; no column {', '.join(sorted(missing))}")
        rows: dict[str, dict[str, str]] = {}
        for row in reader:
            key = f"{row['variant_nuc']} {row['variant_prot']}"
            if key in rows:
                raise ValueError(f"{path}: {key} is listed twice")
            rows[key] = row
        return rows


@dataclasses.dataclass
class _Agreement:
    """The rows weaver answered and those of them it spells as ClinVar does, by the kind of ClinVar's edit."""

    answered: collections.Counter[str] = dataclasses.field(default_factory=collections.Counter)
    agree: collections.Counter[str] = dataclasses.field(default_factory=collections.Counter)

    @classmethod
    def of(cls, rows: Iterable[Mapping[str, str]]) -> _Agreement:
        """The agreement over `rows`; a row weaver answered with an error is not counted."""
        agreement = cls()
        for row in rows:
            if _protein_kind(row["rs_p"]) in {"error", "blank"}:
                continue
            kind = _protein_kind(row["variant_prot"])
            agreement.answered[kind] += 1
            agreement.agree[kind] += _spells(row["rs_p"], row["variant_prot"])
        return agreement

    def __str__(self) -> str:
        return _share(self.answered.total(), self.agree.total())

    def print_by_kind(self) -> None:
        for kind, answered in self.answered.most_common():
            print(f"    ClinVar {kind}: {_share(answered, self.agree[kind])}")


def _share(answered: int, agree: int) -> str:
    share = f" ({100 * agree / answered:.3f}%)" if answered else ""
    return f"{answered} answered, {agree} of them spelled as ClinVar spells them{share}"


def compare_proteins(baseline: pathlib.Path | None, output: pathlib.Path) -> None:
    """Print how often a protein report's predictions spell ClinVar's protein, and what moved since a baseline.

    First the output's agreement with ClinVar's spelling, in all and by the kind of ClinVar's edit. Given a
    baseline, over the rows both reports hold: each run's agreement, and the rows whose prediction moved, counted
    by the kind of edit each run writes (del -> Ter).
    """
    after = _load_proteins(output)
    agreement = _Agreement.of(after.values())
    print(f"\n{output.name}: {len(after)} rows; {agreement}")
    agreement.print_by_kind()
    if baseline is None:
        return
    before = _load_proteins(baseline)
    shared = before.keys() & after.keys()
    print(f"{baseline.name}: {len(before)} rows; {len(shared)} shared")
    for path, rows in ((baseline, before), (output, after)):
        print(f"    on the shared rows, {path.name}: {_Agreement.of(rows[key] for key in shared)}")
    moved = _moved(before, after, shared, "rs_p", lambda row: _protein_kind(row["rs_p"]))
    print(f"rs_p: {moved.total()} rows moved")
    moved.print()


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data", type=pathlib.Path, required=True, help="where inputs, shards, store and genome live")
    parser.add_argument("--variants", type=pathlib.Path, required=True, help="the ClinVar variants TSV to validate")
    parser.add_argument("--output", type=pathlib.Path, required=True, help="the validation output TSV to write")
    parser.add_argument(
        "--baseline",
        type=pathlib.Path,
        default=None,
        help="an earlier output to compare with, row by row",
    )
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--release", default="RS_2024_08", help="the NCBI annotation release")
    parser.add_argument("--accession", default="GCF_000001405.40", help="the assembly's RefSeq accession")
    parser.add_argument("--assembly", default="GRCh38.p14", help="the assembly's name in NCBI's file names")
    parser.add_argument(
        "--anchor",
        default="RS_2023_03",
        help="the release the historical set beside it is anchored on",
    )
    parser.add_argument("--mane", default="MANE.GRCh38.v1.5.summary.txt.gz", help="the MANE summary's file name")
    parser.add_argument(
        "--clinvar-proteins",
        type=pathlib.Path,
        default=None,
        help="also predict the protein of every coding indel in ClinVar's HGVS report, to this TSV",
    )
    parser.add_argument(
        "--clinvar-proteins-baseline",
        type=pathlib.Path,
        default=None,
        help="an earlier --clinvar-proteins report to compare with, row by row",
    )
    args = parser.parse_args(argv)
    if args.clinvar_proteins_baseline is not None and args.clinvar_proteins is None:
        parser.error("--clinvar-proteins-baseline compares a --clinvar-proteins report; name the report to write")
    for written, read in ((args.output, args.baseline), (args.clinvar_proteins, args.clinvar_proteins_baseline)):
        if written is not None and read is not None and written.resolve() == read.resolve():
            parser.error(f"{written} is both written and the baseline it is compared with")
    logging.basicConfig(level=logging.INFO)

    release = Release(args.release, args.accession, args.assembly, args.anchor, args.mane)
    inputs = fetch(release, args.data / "ncbi")
    store, genome = build(release, inputs, args.data)
    proteins: tuple[pathlib.Path, list[tuple[str, str]]] | None = None
    if args.clinvar_proteins is not None:
        hgvs4variation = download(CLINVAR_HGVS_URL, args.data / "clinvar" / "hgvs4variation.txt.gz")
        fetched = datetime.datetime.fromtimestamp(hgvs4variation.stat().st_mtime, tz=datetime.UTC).date()
        proteins = (args.clinvar_proteins, clinvar_proteins(hgvs4variation))
        print(f"{len(proteins[1])} coding indels from ClinVar's HGVS report, fetched {fetched}")
    validate.main([
        str(args.variants), "--no-ferro", "--workers", str(args.workers),
        "--store", str(store), "--genome", str(genome), "--output-file", str(args.output),
    ])  # fmt: skip
    if args.baseline is not None:
        compare(args.baseline, args.output)
    if proteins is not None:
        report, pairs = proteins
        predict_proteins(pairs, report, store=store, genome=genome, workers=args.workers)
        compare_proteins(args.clinvar_proteins_baseline, report)


if __name__ == "__main__":
    main()
