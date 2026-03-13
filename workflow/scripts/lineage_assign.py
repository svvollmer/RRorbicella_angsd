#!/usr/bin/env python3
"""
lineage_assign.py — Lineage Gate (Segment 3 → Segment 4)

Shows admixture log-likelihoods across all K values so you can choose the
biologically appropriate K, then assigns each unrelated sample to a lineage:

  Q[:,0] > threshold  → lineageA
  Q[:,1] > threshold  → lineageB  (for K=2; generalizes to selected component)
  otherwise           → admixed   (excluded from lineage FST)

Writes results/admixture/lineage_assignments.txt so Snakefile.diversity can
build per-lineage bamlists and run lineage-aware FST.

Usage:
    python workflow/scripts/lineage_assign.py
    python workflow/scripts/lineage_assign.py --q-threshold 0.80
    python workflow/scripts/lineage_assign.py --k 2 --q-threshold 0.80 --yes
"""

import argparse
import os
import sys

import numpy as np

UNRELATED_FILE  = "results/relatedness/unrelated_samples.txt"
ADMIXTURE_DIR   = "results/admixture"
OUTPUT_FILE     = "results/admixture/lineage_assignments.txt"


def scan_k_results():
    """
    Return list of (k, log_likelihood) for all available K runs.
    PCAngsd .log files contain a line like:
        "Log-like: -12345678.90" or "Converged" with log-like nearby.
    """
    results = []
    for k in range(2, 20):
        q_path   = os.path.join(ADMIXTURE_DIR, f"pcangsd_K{k}.Q")
        log_path = os.path.join(ADMIXTURE_DIR, f"pcangsd_K{k}.log")
        if not os.path.exists(q_path):
            break
        ll = None
        if os.path.exists(log_path):
            with open(log_path) as f:
                for line in f:
                    if "Log-like" in line or "loglike" in line.lower():
                        parts = line.strip().split()
                        for i, p in enumerate(parts):
                            try:
                                ll = float(parts[i + 1])
                                break
                            except (ValueError, IndexError):
                                pass
                        if ll is not None:
                            break
        results.append((k, ll))
    return results


def print_k_summary(k_results):
    """Print K vs log-likelihood table with delta-K hint."""
    print()
    print("=" * 52)
    print("Available K values (admixture results)")
    print("=" * 52)
    print(f"  {'K':>3}   {'Log-likelihood':>18}   {'ΔLL':>12}")
    print("-" * 52)
    prev_ll = None
    for k, ll in k_results:
        ll_str  = f"{ll:.2f}" if ll is not None else "  (log not found)"
        delta   = ""
        if ll is not None and prev_ll is not None:
            delta = f"{ll - prev_ll:+.2f}"
        print(f"  {k:>3}   {ll_str:>18}   {delta:>12}")
        if ll is not None:
            prev_ll = ll
    print("=" * 52)
    print()


def load_q(k):
    path = os.path.join(ADMIXTURE_DIR, f"pcangsd_K{k}.Q")
    if not os.path.exists(path):
        print(f"ERROR: {path} not found.", file=sys.stderr)
        sys.exit(1)
    return np.loadtxt(path)


def load_samples():
    if not os.path.exists(UNRELATED_FILE):
        print(f"ERROR: {UNRELATED_FILE} not found. Run Segment 2 first.", file=sys.stderr)
        sys.exit(1)
    with open(UNRELATED_FILE) as f:
        return [l.strip() for l in f if l.strip()]


def assign_lineages(q, samples, threshold):
    """
    Assign lineages from a K-column Q matrix.
    For K=2: col 0 = lineageA, col 1 = lineageB.
    For K>2: same logic — col 0 = lineageA, col 1 = lineageB; other cols are admixed.
    """
    if q.ndim == 1:
        q = q.reshape(-1, 1)

    assignments = []
    for i, sample in enumerate(samples):
        q0 = float(q[i, 0])
        q1 = float(q[i, 1]) if q.shape[1] > 1 else 0.0
        q_vals = [float(q[i, j]) for j in range(q.shape[1])]

        if q0 > threshold:
            lineage = "lineageA"
        elif q1 > threshold:
            lineage = "lineageB"
        else:
            lineage = "admixed"

        assignments.append((sample, lineage, q0, q1, q_vals))
    return assignments


def print_assignment_summary(assignments, threshold, k, metadata_path=None):
    pop_map = {}
    if metadata_path and os.path.exists(metadata_path):
        try:
            import pandas as pd
            df = pd.read_csv(metadata_path).set_index("sample_id")
            pop_map = df["population"].to_dict()
        except Exception:
            pass

    counts = {"lineageA": 0, "lineageB": 0, "admixed": 0}
    for _, lineage, _, _, _ in assignments:
        counts[lineage] += 1

    print()
    print("=" * 76)
    print(f"Lineage Assignment from K={k}  (Q threshold = {threshold})")
    print("=" * 76)
    print(f"{'Sample':<22} {'Population':<12} {'Q_lineageA':>12} {'Q_lineageB':>12}  Assignment")
    print("-" * 76)

    for sample, lineage, q0, q1, _ in sorted(assignments, key=lambda x: (x[1], x[0])):
        pop    = pop_map.get(sample, "")
        marker = "  ← admixed" if lineage == "admixed" else ""
        print(f"{sample:<22} {pop:<12} {q0:>12.4f} {q1:>12.4f}  {lineage}{marker}")

    print("=" * 76)
    print(f"  lineageA : {counts['lineageA']:>3} samples  (Q_col0 > {threshold})")
    print(f"  lineageB : {counts['lineageB']:>3} samples  (Q_col1 > {threshold})")
    print(f"  admixed  : {counts['admixed']:>3} samples  (excluded from lineage FST)")
    print("=" * 76)
    print()

    if counts["lineageA"] < 5:
        print(f"WARNING: lineageA has only {counts['lineageA']} samples — FST will be noisy.", file=sys.stderr)
    if counts["lineageB"] < 5:
        print(f"WARNING: lineageB has only {counts['lineageB']} samples — FST will be noisy.", file=sys.stderr)


def write_assignments(assignments, k, threshold):
    os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)
    with open(OUTPUT_FILE, "w") as f:
        f.write(f"# lineage_assignments.txt — written by lineage_assign.py\n")
        f.write(f"# K={k}, Q_threshold={threshold}\n")
        f.write("# sample_id\tlineage\tQ_lineageA\tQ_lineageB\n")
        for sample, lineage, q0, q1, _ in assignments:
            f.write(f"{sample}\t{lineage}\t{q0:.6f}\t{q1:.6f}\n")

    print(f"Written: {OUTPUT_FILE}")
    print()
    print("Next: start Segment 4:")
    print("  snakemake --snakefile workflow/Snakefile.diversity --profile profiles/aws")
    print()


def main():
    parser = argparse.ArgumentParser(
        description="Lineage gate — assign samples to lineages from admixture Q-file"
    )
    parser.add_argument("--k",           type=int,   default=None,  help="K to use for assignment (prompted if not set)")
    parser.add_argument("--q-threshold", type=float, default=None,  help="Q threshold (default from config or 0.80)")
    parser.add_argument("--metadata",    default="config/samples.csv")
    parser.add_argument("--yes", "-y",   action="store_true",       help="Skip confirmation prompts")
    args = parser.parse_args()

    # Read threshold from config if not set on command line
    threshold = args.q_threshold
    if threshold is None:
        try:
            import yaml
            with open("config/config.yaml") as f:
                cfg = yaml.safe_load(f)
            threshold = cfg.get("lineage_q_threshold", 0.80)
        except Exception:
            threshold = 0.80

    # Show available K results
    k_results = scan_k_results()
    if not k_results:
        print(f"ERROR: No admixture Q files found in {ADMIXTURE_DIR}/", file=sys.stderr)
        sys.exit(1)

    print_k_summary(k_results)

    # Select K
    available_ks = [k for k, _ in k_results]
    if args.k is not None:
        k = args.k
        if k not in available_ks:
            print(f"ERROR: K={k} not in available K values {available_ks}", file=sys.stderr)
            sys.exit(1)
    elif args.yes:
        k = 2  # default for non-interactive
        print(f"Using K={k} (default for --yes mode)")
    else:
        default_k = 2
        raw = input(f"Select K for lineage assignment [{default_k}]: ").strip()
        k = int(raw) if raw else default_k
        if k not in available_ks:
            print(f"ERROR: K={k} not available.", file=sys.stderr)
            sys.exit(1)

    # Load samples + Q matrix
    samples = load_samples()
    q       = load_q(k)

    if len(samples) != q.shape[0]:
        print(f"ERROR: {len(samples)} samples in unrelated_samples.txt "
              f"but K={k} Q matrix has {q.shape[0]} rows.", file=sys.stderr)
        sys.exit(1)

    # Assign and display
    assignments = assign_lineages(q, samples, threshold)
    print_assignment_summary(assignments, threshold, k, args.metadata)

    # Optionally adjust threshold
    if not args.yes:
        ans = input(f"Write assignments? [Y/n/t=adjust threshold (current: {threshold})]: ").strip().lower()
        if ans == 't':
            raw = input(f"New threshold [0.0–1.0]: ").strip()
            try:
                threshold = float(raw)
            except ValueError:
                print("Invalid threshold.", file=sys.stderr)
                sys.exit(1)
            assignments = assign_lineages(q, samples, threshold)
            print_assignment_summary(assignments, threshold, k, args.metadata)
            ans = input("Write lineage_assignments.txt? [Y/n]: ").strip().lower()
        if ans not in ("", "y", "yes"):
            print("Aborted — lineage_assignments.txt not written.")
            sys.exit(0)

    write_assignments(assignments, k, threshold)


if __name__ == "__main__":
    main()
