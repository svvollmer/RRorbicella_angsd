#!/usr/bin/env python3
"""
QC Approval Gate — interactive sample review for coral-angsd-pipeline.

Reads results/qc/samples_pass.txt (written by the qc_report Snakemake rule),
walks the user through flagged samples with per-metric context, collects
confirmations, and writes:

  results/qc/samples_approved.txt   — one sample_id per line (pipeline input)
  results/qc/samples_exclusions.txt — audit log with reasons, decisions, timestamp

Usage:
    python workflow/scripts/qc_approve.py [--results results/qc]
    python workflow/scripts/qc_approve.py --yes   # auto-confirm all recommendations (CI)
"""

import argparse
import csv
import os
import socket
import sys
from datetime import datetime


# ── helpers ──────────────────────────────────────────────────────────────────

def bold(s):   return f"\033[1m{s}\033[0m"
def red(s):    return f"\033[31m{s}\033[0m"
def yellow(s): return f"\033[33m{s}\033[0m"
def green(s):  return f"\033[32m{s}\033[0m"
def cyan(s):   return f"\033[36m{s}\033[0m"
def dim(s):    return f"\033[2m{s}\033[0m"

SEP = "─" * 62

def ask_yn(prompt, default="y", auto=False):
    """Prompt for yes/no. Returns True for yes."""
    hint = "[Y/n]" if default == "y" else "[y/N]"
    if auto:
        print(f"{prompt} {hint}: (auto) {'yes' if default == 'y' else 'no'}")
        return default == "y"
    while True:
        resp = input(f"{prompt} {hint}: ").strip().lower()
        if resp == "":
            return default == "y"
        if resp in ("y", "yes"):
            return True
        if resp in ("n", "no"):
            return False
        print("  Please enter y or n.")

def ask_reason(default, auto=False):
    if auto:
        return default
    resp = input(f"  Reason [{default}]: ").strip()
    return resp if resp else default


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--results", default="results/qc",
                        help="Path to QC results directory (default: results/qc)")
    parser.add_argument("--yes", action="store_true",
                        help="Non-interactive: auto-confirm all recommendations")
    args = parser.parse_args()

    pass_file  = os.path.join(args.results, "samples_pass.txt")
    dist_file  = os.path.join(args.results, "qc_distribution.txt")
    out_approved  = os.path.join(args.results, "samples_approved.txt")
    out_exclusions = os.path.join(args.results, "samples_exclusions.txt")

    if not os.path.exists(pass_file):
        sys.exit(
            red(f"\nError: {pass_file} not found.\n") +
            "Run snakemake first to complete alignment and QC:\n"
            "  snakemake --keep-going ...\n"
            "Then re-run this script when the QC report is ready.\n"
        )

    # ── load distribution ────────────────────────────────────────────────────
    dist = {}
    if os.path.exists(dist_file):
        with open(dist_file) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                dist[row["metric"]] = row

    # ── load samples ─────────────────────────────────────────────────────────
    samples = []
    with open(pass_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            samples.append(row)

    flagged  = [s for s in samples if s["flags"] != "PASS"]
    passing  = [s for s in samples if s["flags"] == "PASS"]

    # ── header ───────────────────────────────────────────────────────────────
    print()
    print(bold("=" * 62))
    print(bold("  QC SAMPLE REVIEW — coral-angsd-pipeline"))
    print(bold("=" * 62))
    print(f"  Samples total:   {len(samples)}")
    print(f"  Flagged:         {yellow(str(len(flagged)))}")
    print(f"  Passing:         {green(str(len(passing)))}")
    print()

    # ── distribution summary ──────────────────────────────────────────────────
    if dist:
        print(bold("  Distribution across all samples:"))
        for metric, row in dist.items():
            mean = float(row["mean"])
            sd   = float(row["sd"])
            mn   = row["min_threshold"]
            mx   = row["max_threshold"]
            if metric == "mean_depth":
                print(f"    Mean depth:    {mean:.1f}x  (±{sd:.1f})  "
                      f"thresholds: {float(mn):.1f}x – {float(mx):.1f}x")
            elif metric == "mapping_rate":
                print(f"    Mapping rate:  {mean:.3f}  (±{sd:.3f})  "
                      f"threshold:  >{float(mn):.3f}")
            elif metric == "dup_rate":
                print(f"    Dup rate:      {mean:.3f}  (±{sd:.3f})  "
                      f"threshold:  <{float(mx):.3f}")
        print()

    # ── flagged samples ───────────────────────────────────────────────────────
    decisions = []   # list of (sample_id, action, reason)

    if flagged:
        print(bold(f"  FLAGGED SAMPLES ({len(flagged)}) — review each:"))
        print(f"  {SEP}")

        for i, s in enumerate(flagged, 1):
            sid   = s["sample_id"]
            depth = float(s["mean_depth"])
            rate  = float(s["mapping_rate"])
            dup   = float(s["dup_rate"])
            flags = s["flags"]

            print(f"\n  [{i}/{len(flagged)}]  {bold(sid)}")
            print(f"    Mean depth:    {depth:.1f}x")
            print(f"    Mapping rate:  {rate:.3f}")
            print(f"    Dup rate:      {dup:.3f}")
            print(f"    Flags:         {yellow(flags)}")
            print()

            exclude = ask_yn(
                f"  {red('Exclude')} {bold(sid)}?",
                default="y",
                auto=args.yes,
            )

            if exclude:
                default_reason = flags
                reason = ask_reason(default_reason, auto=args.yes)
                decisions.append((sid, "excluded", reason))
                print(f"  {red('✗')} {sid} excluded — {reason}")
            else:
                decisions.append((sid, "approved_manual", f"flag_overridden: {flags}"))
                print(f"  {green('✓')} {sid} retained despite flag")

        print(f"\n  {SEP}")

    else:
        print(green("  No flagged samples — all pass QC thresholds.\n"))

    # ── passing samples preview ───────────────────────────────────────────────
    if passing:
        print(f"\n  {bold('PASSING SAMPLES')} ({len(passing)}):")
        col = "    {:<20} {:>8}  {:>7}  {:>7}  {}"
        print(dim(col.format("sample_id", "depth", "maprate", "dup_rate", "flags")))
        for s in passing:
            print(col.format(
                s["sample_id"],
                s["mean_depth"] + "x",
                s["mapping_rate"],
                s["dup_rate"],
                s["flags"],
            ))

    # ── manual exclusions ─────────────────────────────────────────────────────
    print()
    if args.yes:
        extra_ids = []
    else:
        extra_raw = input(
            "  Additional manual exclusions? "
            "(comma-separated sample IDs, or Enter to skip): "
        ).strip()
        extra_ids = [x.strip() for x in extra_raw.split(",") if x.strip()]

    known_ids = {s["sample_id"] for s in samples}
    for eid in extra_ids:
        if eid not in known_ids:
            print(yellow(f"  Warning: '{eid}' not in samples list — skipping"))
            continue
        if any(d[0] == eid for d in decisions):
            print(dim(f"  '{eid}' already in decisions — skipping"))
            continue
        if args.yes:
            reason = "manual_exclusion"
        else:
            reason = ask_reason("manual_exclusion")
        decisions.append((eid, "excluded_manual", reason))
        print(f"  {red('✗')} {eid} manually excluded — {reason}")

    # ── build approved list ───────────────────────────────────────────────────
    excluded_ids = {d[0] for d in decisions if d[1].startswith("excluded")}
    approved_ids = [s["sample_id"] for s in samples if s["sample_id"] not in excluded_ids]

    # add passing+retained flagged to decisions for audit completeness
    decided_ids = {d[0] for d in decisions}
    for s in samples:
        if s["sample_id"] not in decided_ids:
            decisions.append((s["sample_id"], "approved", "PASS"))

    # ── summary ───────────────────────────────────────────────────────────────
    print()
    print(bold("=" * 62))
    print(bold("  SUMMARY"))
    print(bold("=" * 62))
    print(f"  {green('Approved:')}  {len(approved_ids)} samples")
    print(f"  {red('Excluded:')}  {len(excluded_ids)} samples"
          + (f"  ({', '.join(sorted(excluded_ids))})" if excluded_ids else ""))
    print()
    print(f"  Will write:")
    print(f"    {cyan(out_approved)}")
    print(f"    {cyan(out_exclusions)}")
    print()

    if not args.yes:
        if not ask_yn("  Confirm and write?", default="y"):
            print(red("\n  Aborted — no files written.\n"))
            sys.exit(1)

    # ── write outputs ─────────────────────────────────────────────────────────
    timestamp = datetime.now().isoformat(timespec="seconds")
    operator  = f"{os.environ.get('USER', 'unknown')}@{socket.gethostname()}"

    # samples_approved.txt — one sample_id per line
    with open(out_approved, "w") as f:
        for sid in approved_ids:
            f.write(sid + "\n")

    # samples_exclusions.txt — full audit
    with open(out_exclusions, "w") as f:
        f.write(f"# QC exclusion audit — {timestamp}\n")
        f.write(f"# Operator: {operator}\n")
        f.write(f"# Approved: {len(approved_ids)}/{len(samples)} samples\n")
        f.write("sample_id\tdecision\treason\toperator\ttimestamp\n")
        for sid, action, reason in sorted(decisions, key=lambda x: x[0]):
            f.write(f"{sid}\t{action}\t{reason}\t{operator}\t{timestamp}\n")

    print()
    print(green(f"  ✓ Written: {out_approved}  ({len(approved_ids)} samples)"))
    print(green(f"  ✓ Written: {out_exclusions}"))
    print()
    print("  Re-run snakemake to continue the pipeline:")
    print(dim("    snakemake --snakefile workflow/Snakefile --keep-going ..."))
    print()


if __name__ == "__main__":
    main()
