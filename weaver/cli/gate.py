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
"""

from __future__ import annotations

import argparse
import collections
import csv
import dataclasses
import hashlib
import pathlib
import sys
import urllib.request

from weaver.cli import validate

NCBI = "https://ftp.ncbi.nlm.nih.gov"
HGNC_URL = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
# The SPDI column against ClinVar's; the equivalence verdict as it stands; the protein column only as answered,
# an error or blank, since ClinVar's protein strings carry an NP accession and a spelling weaver's do not share.
COMPARED = (("rs_spdi", "spdi"), ("rs_equiv", None), ("rs_p", None))


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


def compare(baseline: pathlib.Path, output: pathlib.Path) -> None:
    """Print, per compared column, the rows that moved between two validation outputs and in which direction."""

    def load(path: pathlib.Path) -> dict[str, dict[str, str]]:
        with path.open(newline="") as fh:
            return {row["variant_nuc"]: row for row in csv.DictReader(fh, delimiter="\t")}

    before, after = load(baseline), load(output)
    shared = before.keys() & after.keys()
    print(f"\n{baseline.name}: {len(before)} rows; {output.name}: {len(after)} rows; {len(shared)} shared")

    def state(row: dict[str, str], column: str, truth: str | None) -> str:
        value = row[column]
        if not value:
            return "blank"
        if column == "rs_equiv":
            return value
        if value.startswith("ERR"):
            return "error"
        if truth is None:
            return "answered"
        return "match" if value == row[truth] else "mismatch"

    for column, truth in COMPARED:
        moved: collections.Counter[str] = collections.Counter()
        examples: list[tuple[str, str, str]] = []
        for key in sorted(shared):
            a, b = before[key], after[key]
            if a[column] != b[column]:
                kinds = (state(a, column, truth), state(b, column, truth))
                moved[f"{kinds[0]} -> {kinds[1]}"] += 1
                if len(examples) < 5 and kinds != ("error", "error"):
                    examples.append((key, a[column][:60], b[column][:60]))
        line = f"{column}: {sum(moved.values())} rows moved"
        if truth is not None:
            pct = [
                100 * sum(state(r, column, truth) == "match" for r in rows.values()) / len(rows)
                for rows in (before, after)
            ]
            line += f"; identity {pct[0]:.3f}% -> {pct[1]:.3f}%"
        print(line)
        for kind, n in moved.most_common():
            print(f"    {n:6d}  {kind}")
        for key, was, now in examples:
            print(f"    e.g. {key}: {was!r} -> {now!r}")


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
    args = parser.parse_args(argv)

    release = Release(args.release, args.accession, args.assembly, args.anchor, args.mane)
    inputs = fetch(release, args.data / "ncbi")
    store, genome = build(release, inputs, args.data)
    validate.main([
        str(args.variants), "--no-ferro", "--workers", str(args.workers),
        "--store", str(store), "--genome", str(genome), "--output-file", str(args.output),
    ])  # fmt: skip
    if args.baseline is not None:
        compare(args.baseline, args.output)


if __name__ == "__main__":
    main()
