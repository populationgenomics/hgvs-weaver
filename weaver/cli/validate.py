# /// script
# requires-python = ">=3.12"
# dependencies = [
#   "hgvs-weaver-data",
#   "tqdm",
#   "parsley",
#   "bioutils",
# ]
# ///

"""Full validation script against ClinVar variants."""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import dataclasses
import logging
import sys
import typing

try:
    import hgvs.parser
    import hgvs.variantmapper
except ImportError:
    print(
        "Error: 'hgvs' package not found. Please install it manually (e.g. without dependencies to avoid psycopg2) with:",
    )
    print("  pip install hgvs --no-deps")
    sys.exit(1)

try:
    import tqdm
except ImportError:
    print("Error: 'tqdm' package not found. Please install it manually with: pip install tqdm")
    sys.exit(1)

import weaver

try:
    from weaver_data_provider import genome as wd_genome
    from weaver_data_provider import provider as wd_provider
    from weaver_data_provider import store as wd_store
except ImportError:
    print("Error: 'hgvs-weaver-data' package not found (it needs Python 3.12 or later). Install it with:")
    print("  pip install hgvs-weaver-data")
    sys.exit(1)

from weaver.cli import provider


@dataclasses.dataclass(frozen=True)
class _Worker:
    """What one worker process holds: the two mappers over the store, the biocommons parser, the ferro results."""

    rp: wd_provider.BundleProvider
    mapper: weaver.VariantMapper
    ref_vm: hgvs.variantmapper.VariantMapper
    ref_hp: typing.Any  # hgvs.parser.Parser builds its parse_* methods at run time, so a checker cannot see them
    fh_results: dict[str, str]  # pre-computed ferro normalize results: nuc_hgvs -> normalized string or "ERR:..."


# The pool's initializer sets this once per worker process; it is the one module-level slot a
# ProcessPoolExecutor initializer can hand its work to.
_worker: _Worker | None = None


def init_worker(store: str, genome: str, fh_results_path: str | None = None) -> None:
    """Initializes this worker process's mappers over a weaver-data-provider store and genome."""
    global _worker
    rp = wd_provider.BundleProvider(wd_store.BundleStore(store), wd_genome.Genome(genome))
    fh_results: dict[str, str] = {}
    if fh_results_path:
        import json  # noqa: PLC0415

        with open(fh_results_path) as f:
            fh_results = json.load(f)
    _worker = _Worker(
        rp=rp,
        mapper=weaver.VariantMapper(rp),
        ref_vm=hgvs.variantmapper.VariantMapper(provider.ReferenceHgvsDataProvider(rp)),
        ref_hp=hgvs.parser.Parser(),
        fh_results=fh_results,
    )


def hgvs_lib_to_spdi(v: typing.Any, data_provider: typing.Any) -> str | None:
    """Converts a standard hgvs library Variant object to SPDI string format."""
    try:
        ac = v.ac
        start_1 = v.posedit.pos.start.base
        end_1 = v.posedit.pos.end.base if v.posedit.pos.end else start_1
        if hasattr(v.posedit.edit, "ref"):
            ref = v.posedit.edit.ref or ""
            alt = v.posedit.edit.alt or ""
            if not ref or ref.isdigit():
                ref = data_provider.get_seq(ac, start_1 - 1, end_1, "g")
            if v.posedit.edit.type == "ins":
                return f"{ac}:{start_1}:{ref}:{alt}"
            return f"{ac}:{start_1 - 1}:{ref}:{alt}"
        return "UnsupportedType"
    except Exception:
        return "ERR:SPDI"


def process_variant(row: dict[str, str]) -> dict[str, str]:
    """Maps a single variant using both weaver and ref-hgvs for comparison."""
    nuc_hgvs = row["variant_nuc"]
    spdi_ac = row["spdi"].split(":")[0]

    rs_p = "ERR"
    rs_spdi = "ERR"
    ref_p = "ERR"
    ref_spdi = "ERR"

    # weaver block
    v_p = None
    rs_p = "ERR"
    rs_spdi = "ERR"
    worker = _worker
    if worker is None:
        res_row = row.copy()
        res_row.update(
            {
                "rs_p": "ERR:MapperNotInit",
                "rs_spdi": "ERR:MapperNotInit",
                "rs_equiv": "ERR:MapperNotInit",
                "ref_equiv": "ERR:MapperNotInit",
            },
        )
        return res_row
    mapper, rp = worker.mapper, worker.rp
    try:
        v_rs_raw = weaver.parse(nuc_hgvs)
        v_rs = mapper.normalize_variant(v_rs_raw)
        try:
            if v_rs.coordinate_type == "c":
                v_p = mapper.c_to_p(v_rs)
                rs_p = v_p.format().split(":")[-1]
        except Exception as e:
            rs_p = f"ERR:{e!s}"

        try:
            rs_spdi = mapper.to_spdi(v_rs_raw, unambiguous=True)
        except Exception as e:
            rs_spdi = f"ERR:{e!s}"
    except Exception as e:
        rs_p = rs_spdi = f"ERR:{e!s}"
    except BaseException:
        rs_p = rs_spdi = "PANIC"

    # ref-hgvs block
    ref_p = "ERR"
    ref_spdi = "ERR"
    try:
        v_ref = worker.ref_hp.parse_hgvs_variant(nuc_hgvs)
        if True:
            try:
                if v_ref.type == "c":
                    v_p_ref = worker.ref_vm.c_to_p(v_ref)
                    ref_p = str(v_p_ref).split(":")[-1]
            except Exception as e:
                ref_p = f"ERR:{e!s}"

            try:
                vg_ref = worker.ref_vm.c_to_g(v_ref, spdi_ac) if v_ref.type != "g" else v_ref
                ref_spdi = hgvs_lib_to_spdi(vg_ref, rp)
            except Exception as e:
                ref_spdi = f"ERR:{e!s}"
    except Exception:
        ref_p = ref_spdi = "ERR:Parse"
    except BaseException:
        ref_p = ref_spdi = "PANIC"

    # ferro-hgvs block: look up pre-computed normalize result
    fh_parse = worker.fh_results.get(nuc_hgvs, "SKIP")

    # Equivalence Checks (Using weaver to judge both)
    rs_equiv = "Unknown"
    ref_equiv = "Unknown"
    gt_p_str = row.get("variant_prot", "")

    if gt_p_str and gt_p_str != "ERR" and not gt_p_str.startswith("ERR"):
        try:
            v_gt = weaver.parse(gt_p_str)
            # RS Equivalence
            if v_p:
                rs_equiv = str(mapper.equivalent_level(v_p, v_gt, rp)).split(".")[-1]

            # REF Equivalence (judged by weaver)
            if ref_p and not ref_p.startswith("ERR"):
                try:
                    # Construct full protein string for weaver parsing
                    ref_p_val = f"{v_gt.ac}:{ref_p}" if ":" not in ref_p else ref_p
                    v_ref_p = weaver.parse(ref_p_val)
                    ref_equiv = str(mapper.equivalent_level(v_ref_p, v_gt, rp)).split(".")[-1]
                except Exception:
                    try:
                        v_ref_p = weaver.parse(ref_p)
                        ref_equiv = str(mapper.equivalent_level(v_ref_p, v_gt, rp)).split(".")[-1]
                    except Exception:
                        logging.exception("Failed to judge equivalence")
        except Exception:
            logging.exception("Failed to judge equivalence")

    res_row = row.copy()
    res_row.update(
        {
            "rs_p": rs_p or "",
            "rs_spdi": rs_spdi or "",
            "ref_p": ref_p or "",
            "ref_spdi": ref_spdi or "",
            "rs_equiv": rs_equiv,
            "ref_equiv": ref_equiv,
            "fh_parse": fh_parse,
        },
    )

    return res_row


def run_ferro_normalize(variants: list[str], reference_dir: str) -> dict[str, str]:
    """Batch-normalizes variants via the ferro CLI; returns nuc_hgvs → result mapping."""
    import shutil  # noqa: PLC0415
    import subprocess  # noqa: PLC0415
    import tempfile  # noqa: PLC0415

    ferro_bin = shutil.which("ferro")
    if not ferro_bin:
        print("Warning: 'ferro' binary not found in PATH; skipping ferro normalization.")
        return {}

    print(f"Running ferro normalize on {len(variants):,} variants (reference: {reference_dir})...")
    with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as tmp_in:
        tmp_in.write("\n".join(variants))
        tmp_in_path = tmp_in.name

    results: dict[str, str] = {}
    try:
        import json  # noqa: PLC0415

        proc = subprocess.run(  # noqa: S603
            [ferro_bin, "normalize", "--reference", reference_dir, "-i", tmp_in_path, "-f", "json"],
            capture_output=True,
            text=True,
            check=False,
        )
        # JSON mode: one JSON object per line with {input, success, output?, error?}
        for raw_line in proc.stdout.splitlines():
            line = raw_line.strip()
            if not line:
                continue
            try:
                obj = json.loads(line)
                variant = obj.get("input", "")
                if not variant:
                    continue
                if obj.get("success"):
                    results[variant] = obj.get("output") or variant
                else:
                    err = obj.get("error") or "UnknownError"
                    results[variant] = f"ERR:{err}"
            except json.JSONDecodeError:
                continue
        # Any variant with no output line (e.g. ferro crashed mid-run)
        for variant in variants:
            if variant not in results:
                results[variant] = "ERR:NoOutput"
    except Exception as e:
        print(f"Warning: ferro normalize failed: {e}")
    finally:
        import os  # noqa: PLC0415

        os.unlink(tmp_in_path)

    ok_count = sum(1 for v in results.values() if not v.startswith("ERR:"))
    print(f"ferro normalize: {ok_count:,}/{len(variants):,} succeeded.")
    return results


def main(argv: list[str] | None = None) -> None:
    """Main entry point for validation."""
    parser = argparse.ArgumentParser(description="Full validation against ClinVar variants.")
    parser.add_argument("input_file", help="Input ClinVar TSV file.")
    parser.add_argument("--max-variants", type=int, default=None, help="Maximum variants to process.")
    parser.add_argument("--output-file", default="clinvar_full_validation.tsv", help="Output validation TSV.")
    parser.add_argument(
        "--store",
        required=True,
        help="A weaver-data-provider store (directory or gs:// prefix) built by 'weaver-data-build index'.",
    )
    parser.add_argument(
        "--genome",
        required=True,
        help="The matching genome (directory or gs:// prefix) built by 'weaver-data-build genome'.",
    )
    parser.add_argument("--workers", type=int, default=4, help="Number of worker processes.")
    parser.add_argument(
        "--ferro-reference",
        default=None,
        help="Path to ferro reference directory (produced by 'ferro prepare'). "
        "When provided, ferro normalize is run in batch before validation.",
    )
    parser.add_argument(
        "--no-ferro",
        action="store_true",
        help="Disable ferro-hgvs comparison entirely (fh_parse column will be SKIP).",
    )
    args = parser.parse_args(argv)

    with open(args.input_file) as f_in:
        reader = csv.DictReader(f_in, delimiter="\t")
        base_fields = [
            f
            for f in (reader.fieldnames or [])
            if f
            not in {"rs_p", "rs_spdi", "ref_p", "ref_spdi", "rs_equiv", "ref_equiv", "equivalence_level", "fh_parse"}
        ]
        fieldnames = [*base_fields, "rs_p", "rs_spdi", "ref_p", "ref_spdi", "rs_equiv", "ref_equiv", "fh_parse"]
        rows: list[dict[str, str]] = (
            [next(reader) for _ in range(args.max_variants)] if args.max_variants else list(reader)
        )

    print(f"Processing {len(rows)} variants with ProcessPool...")

    # Pre-run ferro normalize in batch if reference directory supplied
    fh_results_path: str | None = None
    if not args.no_ferro and args.ferro_reference:
        import json  # noqa: PLC0415
        import tempfile  # noqa: PLC0415

        all_nuc = [row["variant_nuc"] for row in rows]
        fh_results = run_ferro_normalize(all_nuc, args.ferro_reference)
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as tmp:
            json.dump(fh_results, tmp)
            fh_results_path = tmp.name

    with open(args.output_file, "w", newline="") as f_out:
        writer = csv.DictWriter(f_out, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()

        with concurrent.futures.ProcessPoolExecutor(
            max_workers=args.workers,
            initializer=init_worker,
            initargs=(args.store, args.genome, fh_results_path),
        ) as executor:
            # map instead of executor.map to catch task-level errors
            results_iter = executor.map(process_variant, rows)

            pbar = tqdm.tqdm(total=len(rows))
            while True:
                try:
                    row_res = next(results_iter)
                    writer.writerow(row_res)
                    pbar.update(1)
                except StopIteration:
                    break
                except Exception as e:
                    print(f"\nWorker crashed: {e}")
                    pbar.update(1)
            pbar.close()

    if fh_results_path:
        import os  # noqa: PLC0415

        os.unlink(fh_results_path)


if __name__ == "__main__":
    main()
