#!/usr/bin/env python3
"""
Clone Approval Gate — interactive relatedness review for coral-angsd-pipeline.

Reads results/relatedness/clones_report.txt (written by the filter_clones
Snakemake rule after ngsRelate) and walks the user through:

  1. Clone pairs (KING >= threshold) — with recommendation for which to drop
  2. Close relatives (KING 0.2-0.45) — flagged only, no auto-exclusion
  3. Optional manual exclusions

Writes:
  results/relatedness/unrelated_samples.txt   — pipeline input for SAF/PCA/LD
  results/relatedness/clone_exclusions.txt    — audit log

Usage:
    python workflow/scripts/clone_approve.py [--results results/relatedness]
    python workflow/scripts/clone_approve.py --yes   # auto-confirm recommendations
"""

import argparse
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
def magenta(s):return f"\033[35m{s}\033[0m"

SEP = "─" * 62

def ask_yn(prompt, default="y", auto=False):
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

def ask_choice(prompt, choices, auto_choice=None, auto=False):
    """Ask user to pick from a list. Returns chosen value."""
    choices_str = "/".join(choices)
    if auto:
        print(f"{prompt} [{choices_str}]: (auto) {auto_choice}")
        return auto_choice
    while True:
        resp = input(f"{prompt} [{choices_str}]: ").strip()
        if resp in choices:
            return resp
        print(f"  Please enter one of: {choices_str}")

def ask_reason(default, auto=False):
    if auto:
        return default
    resp = input(f"  Reason [{default}]: ").strip()
    return resp if resp else default


# ── parse clones_report.txt ───────────────────────────────────────────────────

def parse_report(path):
    """Returns (approved, clones, relatives, thresholds)."""
    approved   = []   # (sample_id, depth)
    clones     = []   # (a, b, king, depth_a, depth_b, recommended_exclude)
    relatives  = []   # (a, b, king, relationship)
    thresholds = {}

    section = None
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("# clone_threshold:"):
                thresholds["clone"] = float(line.split(":")[1].strip())
            elif line.startswith("# kinship_threshold:"):
                thresholds["kinship"] = float(line.split(":")[1].strip())
            elif line.startswith("[approved]"):
                section = "approved"
            elif line.startswith("[clones]"):
                section = "clones"
            elif line.startswith("[relatives]"):
                section = "relatives"
            elif line.startswith("#") or line.startswith("sample") or not line.strip():
                continue
            elif section == "approved":
                parts = line.split("\t")
                approved.append((parts[0], float(parts[1]) if len(parts) > 1 else 0.0))
            elif section == "clones":
                parts = line.split("\t")
                clones.append((
                    parts[0], parts[1],
                    float(parts[2]),
                    float(parts[3]), float(parts[4]),
                    parts[5],
                ))
            elif section == "relatives":
                parts = line.split("\t")
                relatives.append((parts[0], parts[1], float(parts[2]),
                                   parts[3] if len(parts) > 3 else ""))

    return approved, clones, relatives, thresholds


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--results", default="results/relatedness",
                        help="Path to relatedness results directory (default: results/relatedness)")
    parser.add_argument("--yes", action="store_true",
                        help="Non-interactive: auto-confirm all recommendations")
    parser.add_argument("--exclude", metavar="ID", nargs="+", default=[],
                        help="Additional sample IDs to exclude (used with --yes)")
    args = parser.parse_args()

    report_file    = os.path.join(args.results, "clones_report.txt")
    out_unrelated  = os.path.join(args.results, "unrelated_samples.txt")
    out_audit      = os.path.join(args.results, "clone_exclusions.txt")

    if not os.path.exists(report_file):
        sys.exit(
            red(f"\nError: {report_file} not found.\n") +
            "Run snakemake first to complete ngsRelate + filter_clones:\n"
            "  snakemake --keep-going ...\n"
        )

    approved, clones, relatives, thresholds = parse_report(report_file)
    approved_ids = [s[0] for s in approved]
    depth_map    = {s[0]: s[1] for s in approved}

    clone_thresh   = thresholds.get("clone",   0.45)
    kinship_thresh = thresholds.get("kinship", 0.20)

    # ── header ────────────────────────────────────────────────────────────────
    print()
    print(bold("=" * 62))
    print(bold("  CLONE / RELATEDNESS REVIEW — coral-angsd-pipeline"))
    print(bold("=" * 62))
    print(f"  Approved samples:   {len(approved_ids)}")
    print(f"  Clone pairs:        {magenta(str(len(clones)))}  "
          f"  (KING ≥ {clone_thresh})")
    print(f"  Close relatives:    {yellow(str(len(relatives)))}  "
          f"  (KING {kinship_thresh}–{clone_thresh}, flagged only)")
    print()

    # ── KING reference ────────────────────────────────────────────────────────
    print(dim("  KING coefficient reference:"))
    print(dim("    ≈ 0.50  Identical / clone (asexual ramet)"))
    print(dim("    ≈ 0.25  First-degree (parent-offspring, full sibling)"))
    print(dim("    ≈ 0.125 Second-degree (half-sibling, grandparent)"))
    print(dim("    ≈ 0.0   Unrelated"))
    print()

    decisions = []   # (sample_id, action, reason)

    # ── clone pairs ───────────────────────────────────────────────────────────
    if clones:
        print(bold(f"  CLONE PAIRS ({len(clones)}) — review each:"))
        print(f"  {SEP}")

        for i, (a, b, king, da, db, rec) in enumerate(clones, 1):
            other = b if rec == a else a
            print(f"\n  [{i}/{len(clones)}]  {bold(a)}  ×  {bold(b)}")
            print(f"    KING:          {magenta(f'{king:.4f}')}")
            print(f"    {a:<20} depth: {da:.1f}x"
                  + ("  ← recommended keep" if rec != a else "  ← recommended EXCLUDE"))
            print(f"    {b:<20} depth: {db:.1f}x"
                  + ("  ← recommended keep" if rec != b else "  ← recommended EXCLUDE"))
            print(f"    Recommendation: exclude {bold(red(rec))}")
            print()

            if args.yes:
                choice = rec
                print(f"  (auto) Excluding {red(rec)}")
            else:
                choice = ask_choice(
                    f"  Which sample to exclude?",
                    choices=[a, b, "neither"],
                    auto_choice=rec,
                    auto=False,
                )

            if choice == "neither":
                decisions.append((a, "clone_retained", f"clone_pair_with_{b}_KING={king:.4f}"))
                decisions.append((b, "clone_retained", f"clone_pair_with_{a}_KING={king:.4f}"))
                print(f"  {yellow('!')} Both retained — downstream analyses may include clones")
            else:
                reason_default = f"clone_of_{other}_KING={king:.4f}"
                reason = ask_reason(reason_default, auto=args.yes)
                decisions.append((choice, "excluded_clone", reason))
                kept = b if choice == a else a
                decisions.append((kept, "clone_kept", f"higher_depth_ramet_of_{choice}"))
                print(f"  {red('✗')} {choice} excluded — {reason}")
                print(f"  {green('✓')} {kept} retained")

        print(f"\n  {SEP}")
    else:
        print(green("  No clone pairs detected.\n"))

    # ── close relatives summary ───────────────────────────────────────────────
    if relatives:
        print(f"\n  {bold('CLOSE RELATIVES')} ({len(relatives)}) — flagged, not excluded:")
        print(dim(f"  {'sample_a':<20} {'sample_b':<20} {'KING':>7}  relationship"))
        for a, b, king, rel in relatives:
            print(f"  {a:<20} {b:<20} {king:>7.4f}  {yellow(rel)}")
        print(dim("\n  Close relatives are retained in all analyses."))
        print(dim("  If you want to exclude any, add them as manual exclusions below."))

    # ── manual exclusions ─────────────────────────────────────────────────────
    print()
    if args.yes:
        extra_ids = args.exclude
    else:
        extra_raw = input(
            "  Additional manual exclusions? "
            "(comma-separated sample IDs, or Enter to skip): "
        ).strip()
        extra_ids = [x.strip() for x in extra_raw.split(",") if x.strip()]

    decided_ids = {d[0] for d in decisions}
    for eid in extra_ids:
        if eid not in {s[0] for s in approved}:
            print(yellow(f"  Warning: '{eid}' not in approved sample list — skipping"))
            continue
        if eid in decided_ids:
            print(dim(f"  '{eid}' already in decisions — skipping"))
            continue
        reason = ask_reason("manual_exclusion", auto=args.yes)
        decisions.append((eid, "excluded_manual", reason))
        decided_ids.add(eid)
        print(f"  {red('✗')} {eid} manually excluded — {reason}")

    # ── build unrelated list ──────────────────────────────────────────────────
    excluded_ids  = {d[0] for d in decisions if d[1].startswith("excluded")}
    unrelated_ids = [s for s in approved_ids if s not in excluded_ids]

    # fill in decisions for all retained samples
    for sid in approved_ids:
        if not any(d[0] == sid for d in decisions):
            decisions.append((sid, "approved", "unrelated"))

    # ── summary ───────────────────────────────────────────────────────────────
    print()
    print(bold("=" * 62))
    print(bold("  SUMMARY"))
    print(bold("=" * 62))
    print(f"  {green('Unrelated:')}  {len(unrelated_ids)} samples")
    print(f"  {red('Excluded:')}   {len(excluded_ids)} samples"
          + (f"  ({', '.join(sorted(excluded_ids))})" if excluded_ids else ""))
    print()
    print(f"  Will write:")
    print(f"    {cyan(out_unrelated)}")
    print(f"    {cyan(out_audit)}")
    print()

    if not args.yes:
        if not ask_yn("  Confirm and write?", default="y"):
            print(red("\n  Aborted — no files written.\n"))
            sys.exit(1)

    # ── write outputs ─────────────────────────────────────────────────────────
    timestamp = datetime.now().isoformat(timespec="seconds")
    operator  = f"{os.environ.get('USER', 'unknown')}@{socket.gethostname()}"

    with open(out_unrelated, "w") as f:
        for sid in unrelated_ids:
            f.write(sid + "\n")

    with open(out_audit, "w") as f:
        f.write(f"# Clone exclusion audit — {timestamp}\n")
        f.write(f"# Operator: {operator}\n")
        f.write(f"# Unrelated: {len(unrelated_ids)}/{len(approved_ids)} samples\n")
        f.write(f"# KING clone threshold: {clone_thresh}\n")
        f.write("sample_id\tdecision\treason\toperator\ttimestamp\n")
        for sid, action, reason in sorted(decisions, key=lambda x: x[0]):
            f.write(f"{sid}\t{action}\t{reason}\t{operator}\t{timestamp}\n")

    print()
    print(green(f"  ✓ Written: {out_unrelated}  ({len(unrelated_ids)} samples)"))
    print(green(f"  ✓ Written: {out_audit}"))
    print()
    print("  Re-run snakemake to continue the pipeline:")
    print(dim("    snakemake --snakefile workflow/Snakefile --keep-going ..."))
    print()


if __name__ == "__main__":
    main()
