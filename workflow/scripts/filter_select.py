#!/usr/bin/env python3
"""
filter_select.py — SNP Filter Gate (Segment 2)

Reads pass1.mafs.gz and shows a grid of SNP counts across MAF and minInd
thresholds. User selects parameters; script writes filter_params.yaml so
Snakemake can proceed with pass 2 (angsd_gl_snps).

Usage:
    python workflow/scripts/filter_select.py
    python workflow/scripts/filter_select.py --mafs results/angsd/pass1.mafs.gz
    python workflow/scripts/filter_select.py --maf 0.05 --min-ind-frac 0.90  # non-interactive
"""

import argparse
import gzip
import math
import os
import sys

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
MAFS_FILE    = "results/angsd/pass1.mafs.gz"
PARAMS_FILE  = "results/angsd/filter_params.yaml"
APPROVED_FILE = "results/qc/samples_approved.txt"

MAF_THRESHOLDS      = [0.01, 0.02, 0.05, 0.10, 0.15, 0.20]
MIN_IND_FRACTIONS   = [0.60, 0.70, 0.80, 0.90, 0.95, 1.00]


def count_snps(mafs_file, n_approved):
    """
    Scan pass1.mafs.gz once. Return a dict:
        (maf_threshold, min_ind_frac) -> count
    """
    print(f"Scanning {mafs_file} ...", flush=True)

    # Pre-compute absolute minInd cutoffs for each fraction
    min_ind_abs = {frac: math.ceil(n_approved * frac) for frac in MIN_IND_FRACTIONS}

    # counts[(maf, frac)] = int
    counts = {(m, f): 0 for m in MAF_THRESHOLDS for f in MIN_IND_FRACTIONS}
    total_sites = 0

    with gzip.open(mafs_file, "rt") as fh:
        next(fh)  # skip header: chromo position major minor ref unknownEM pK-EM nInd
        for line in fh:
            parts = line.split("\t")
            try:
                freq  = float(parts[5])   # unknownEM = estimated minor allele freq
                n_ind = int(parts[7])
            except (ValueError, IndexError):
                continue

            maf = min(freq, 1.0 - freq)   # minor allele frequency
            total_sites += 1

            for m in MAF_THRESHOLDS:
                if maf < m:
                    break               # MAF thresholds are sorted ascending
                for f in MIN_IND_FRACTIONS:
                    if n_ind >= min_ind_abs[f]:
                        counts[(m, f)] += 1

    print(f"Done. {total_sites:,} sites scanned.\n")
    return counts, min_ind_abs


def print_table(counts, n_approved, min_ind_abs):
    """Print SNP count grid with minInd fractions as columns, MAF as rows."""
    col_width = 13
    header = f"{'MAF \\ minInd':>12}" + "".join(
        f"  {f:.0%} (≥{min_ind_abs[f]:d})".rjust(col_width)
        for f in MIN_IND_FRACTIONS
    )
    print("=" * len(header))
    print(f"SNP counts from pass 1  (N_approved = {n_approved})")
    print("=" * len(header))
    print(header)
    print("-" * len(header))
    for m in MAF_THRESHOLDS:
        row = f"  MAF ≥ {m:.2f}" + "".join(
            f"{counts[(m, f)]:>{col_width},}"
            for f in MIN_IND_FRACTIONS
        )
        print(row)
    print("=" * len(header))
    print()


def prompt_float(prompt, default, valid):
    while True:
        raw = input(f"{prompt} [{default}]: ").strip()
        val = float(raw) if raw else default
        if val in valid:
            return val
        print(f"  Please choose from: {valid}")


def get_n_approved():
    if os.path.exists(APPROVED_FILE):
        with open(APPROVED_FILE) as f:
            n = sum(
                1 for l in f
                if l.strip() and not l.startswith('#') and not l.startswith('sample')
            )
        return n
    # Fall back to counting samples.csv — load config to find it
    try:
        import yaml
        with open("config/config.yaml") as f:
            cfg = yaml.safe_load(f)
        import pandas as pd
        df = pd.read_csv(cfg["samples_csv"])
        return len(df)
    except Exception:
        return None


def write_params(maf, frac, min_ind):
    os.makedirs(os.path.dirname(PARAMS_FILE), exist_ok=True)
    with open(PARAMS_FILE, "w") as f:
        f.write(f"# Written by filter_select.py\n")
        f.write(f"min_maf: {maf}\n")
        f.write(f"min_ind_frac: {frac}\n")
        f.write(f"min_ind: {min_ind}\n")
    print(f"\nWritten: {PARAMS_FILE}")
    print(f"  min_maf:      {maf}")
    print(f"  min_ind_frac: {frac}")
    print(f"  min_ind:      {min_ind}  (absolute sample count)")
    print("\nNow re-run snakemake to start pass 2.\n")


def main():
    parser = argparse.ArgumentParser(description="SNP filter gate — select MAF and minInd for pass 2")
    parser.add_argument("--mafs",        default=MAFS_FILE, help="pass1.mafs.gz path")
    parser.add_argument("--maf",         type=float, help="MAF threshold (non-interactive)")
    parser.add_argument("--min-ind-frac",type=float, help="minInd fraction (non-interactive)")
    args = parser.parse_args()

    if not os.path.exists(args.mafs):
        print(f"ERROR: {args.mafs} not found. Run pass 1 first.", file=sys.stderr)
        sys.exit(1)

    n_approved = get_n_approved()
    if n_approved is None:
        print("ERROR: could not determine N_approved. Is samples_approved.txt present?", file=sys.stderr)
        sys.exit(1)

    # Non-interactive mode
    if args.maf is not None and args.min_ind_frac is not None:
        maf  = args.maf
        frac = args.min_ind_frac
        min_ind = math.ceil(n_approved * frac)
        write_params(maf, frac, min_ind)
        return

    # Interactive mode: scan pass1.mafs.gz and show grid
    counts, min_ind_abs = count_snps(args.mafs, n_approved)
    print_table(counts, n_approved, min_ind_abs)

    print("Select parameters for pass 2 (angsd_gl_snps).")
    print("These control which SNPs enter the beagle file for PCA/admixture/LD.")
    print("SAF/diversity/FST uses all non-repeat sites regardless of these settings.\n")

    maf  = prompt_float("MAF threshold",      default=0.05, valid=MAF_THRESHOLDS)
    frac = prompt_float("minInd fraction",    default=0.90, valid=MIN_IND_FRACTIONS)
    min_ind = min_ind_abs[frac]

    selected = counts[(maf, frac)]
    print(f"\nSelected: MAF ≥ {maf}, minInd ≥ {frac:.0%} ({min_ind} samples)")
    print(f"Pass 2 will process {selected:,} SNPs.")

    confirm = input("\nWrite filter_params.yaml? [Y/n]: ").strip().lower()
    if confirm in ("", "y", "yes"):
        write_params(maf, frac, min_ind)
    else:
        print("Aborted — filter_params.yaml not written.")


if __name__ == "__main__":
    main()
