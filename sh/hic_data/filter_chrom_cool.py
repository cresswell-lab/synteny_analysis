#!/usr/bin/env python3
"""
Filter an .mcool file down to a specific list of chromosomes, preserving
(or regenerating) all resolutions present in the original multi-res file.

Strategy:
  1. Filter the finest (smallest bin size) resolution to the requested
     chromosomes -> write as a standalone .cool.
  2. Re-zoomify that filtered .cool into a new .mcool with the same set
     of resolutions as the original (or a user-specified set), with
     matrix balancing re-run on each level.

Usage:
  python filter_mcool_chroms.py \
      --input /path/to/file.mcool \
      --output /path/to/file.filtered.mcool \
      --chroms chr1,chr2,chr3,...,chrX,chrY \
      [--chroms-file chroms.txt]   # one chrom name per line, alternative to --chroms
      [--resolutions 5000,10000,50000]  # defaults to original file's resolutions
      [--no-balance]               # skip cooler balance on the rebuilt pyramid
      [--nproc 4]                  # parallel workers for zoomify/balance

Requires: cooler (pip install cooler)
"""

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import cooler


def get_mcool_resolutions(mcool_path: str):
    """Return sorted list of integer resolutions present in an .mcool file."""
    entries = cooler.fileops.list_coolers(mcool_path)
    resolutions = []
    for e in entries:
        # entries look like '/resolutions/5000'
        parts = e.strip("/").split("/")
        if len(parts) == 2 and parts[0] == "resolutions":
            resolutions.append(int(parts[1]))
    if not resolutions:
        raise ValueError(
            f"No '/resolutions/<N>' groups found in {mcool_path}. "
            "Is this really a multi-resolution .mcool file?"
        )
    return sorted(resolutions)


def load_chrom_list(chroms_arg: str, chroms_file_arg: str):
    if chroms_file_arg:
        with open(chroms_file_arg) as fh:
            chroms = [line.strip() for line in fh if line.strip()]
    elif chroms_arg:
        chroms = [c.strip() for c in chroms_arg.split(",") if c.strip()]
    else:
        raise ValueError("Must supply either --chroms or --chroms-file")
    return chroms


def filter_base_resolution(input_uri: str, keep_chroms, out_cool_path: str):
    """Filter a single-resolution cooler URI down to keep_chroms, write out_cool_path."""
    c = cooler.Cooler(input_uri)

    missing = set(keep_chroms) - set(c.chromnames)
    if missing:
        raise ValueError(
            f"These requested chroms aren't in {input_uri}: {sorted(missing)}\n"
            f"Available chroms: {c.chromnames}"
        )

    dropped = [ch for ch in c.chromnames if ch not in keep_chroms]
    print(f"[filter] Keeping {len(keep_chroms)} chroms, dropping {len(dropped)}: {dropped}")

    bins = c.bins()[:]
    keep_mask = bins["chrom"].isin(keep_chroms)
    new_bins = bins[keep_mask].reset_index(drop=True)

    old_to_new = {old_id: new_id for new_id, old_id in enumerate(bins.index[keep_mask])}

    def pixel_iter(chunksize=5_000_000):
        nnz = c.info["nnz"]
        selector = c.pixels()
        for lo in range(0, nnz, chunksize):
            hi = min(lo + chunksize, nnz)
            chunk = selector[lo:hi]
            chunk = chunk[chunk["bin1_id"].isin(old_to_new) & chunk["bin2_id"].isin(old_to_new)]
            chunk = chunk.copy()
            chunk["bin1_id"] = chunk["bin1_id"].map(old_to_new)
            chunk["bin2_id"] = chunk["bin2_id"].map(old_to_new)
            yield chunk

    cooler.create_cooler(
        out_cool_path,
        bins=new_bins,
        pixels=pixel_iter(),
        dtypes={"count": c.pixels().dtypes["count"]},
        ordered=True,
        symmetric_upper=True,
    )
    print(f"[filter] Wrote filtered base-resolution cooler: {out_cool_path}")


def zoomify(
    base_cool_path: str,
    out_mcool_path: str,
    resolutions,
    do_balance: bool,
    nproc: int,
):
    """Rebuild an .mcool pyramid from a filtered base-resolution .cool."""
    res_str = ",".join(str(r) for r in resolutions)

    cmd = [
        "cooler", "zoomify",
        "--nproc", str(nproc),
        "--resolutions", res_str,
        "--out", out_mcool_path,
    ]
    if do_balance:
        cmd.append("--balance")

    cmd.append(base_cool_path)

    print(f"[zoomify] Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    print(f"[zoomify] Wrote multi-resolution file: {out_mcool_path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, help="Path to input .mcool file")
    ap.add_argument("--output", required=True, help="Path to output .mcool file")
    ap.add_argument("--chroms", default=None, help="Comma-separated list of chromosomes to keep")
    ap.add_argument("--chroms-file", default=None, help="File with one chromosome name per line (alt to --chroms)")
    ap.add_argument(
        "--resolutions",
        default=None,
        help="Comma-separated resolutions to include in the output "
             "(default: same resolutions as the input .mcool)",
    )
    ap.add_argument("--no-balance", action="store_true", help="Skip cooler balance on the rebuilt pyramid")
    ap.add_argument("--nproc", type=int, default=4, help="Parallel workers for zoomify/balance (default: 4)")
    ap.add_argument("--keep-tmp", action="store_true", help="Keep intermediate filtered base .cool file")
    args = ap.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if not input_path.exists():
        sys.exit(f"Input file not found: {input_path}")

    keep_chroms = load_chrom_list(args.chroms, args.chroms_file)

    all_resolutions = get_mcool_resolutions(str(input_path))
    finest = min(all_resolutions)
    print(f"[main] Resolutions found in input: {all_resolutions}")
    print(f"[main] Using finest resolution as base for filtering: {finest}")

    target_resolutions = (
        [int(r) for r in args.resolutions.split(",")] if args.resolutions else all_resolutions
    )
    if finest not in target_resolutions:
        # zoomify needs the base resolution included in its pyramid; if user
        # asked for a custom resolution set that excludes the finest one,
        # we still zoom from `finest`, cooler will just also materialize it.
        print(
            f"[main] Note: base resolution {finest} not in requested output "
            f"resolutions {target_resolutions}; it will still be produced by zoomify."
        )

    with tempfile.TemporaryDirectory() as tmpdir:
        base_cool = str(Path(tmpdir) / "filtered_base.cool")
        input_uri = f"{input_path}::/resolutions/{finest}"

        filter_base_resolution(input_uri, keep_chroms, base_cool)

        zoomify(
            base_cool_path=base_cool,
            out_mcool_path=str(output_path),
            resolutions=target_resolutions,
            do_balance=not args.no_balance,
            nproc=args.nproc,
        )

        if args.keep_tmp:
            kept_path = output_path.with_suffix("").as_posix() + f".base_{finest}.cool"
            shutil.copy(base_cool, kept_path)
            print(f"[main] Kept intermediate base cooler at: {kept_path}")

    print("[main] Done.")


if __name__ == "__main__":
    main()