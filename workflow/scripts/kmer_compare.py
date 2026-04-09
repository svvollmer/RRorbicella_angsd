"""
kmer_compare.py — Compare k-mer sharing between Orbicella species pairs

Uses sort + streaming merge to avoid loading 800M+ k-mer sets into memory.
Each dump file is sorted on disk, then pairs are compared via a streaming
merge-join that counts intersection without holding full sets in RAM.

Usage:
    python kmer_compare.py --oann oann_kmers.txt \
                           --ofav ofav_kmers.txt \
                           --ofra ofra_kmers.txt \
                           --out kmer_summary.txt
"""

import argparse
import subprocess
import sys
import os
import tempfile


def sort_kmers(path, tmpdir):
    """Extract k-mer column and sort. Returns path to sorted k-mer-only file.

    Jellyfish dump -c produces 'KMER COUNT' lines. comm -12 requires exact-line
    matches, so we must strip counts before sorting — otherwise k-mers present
    in both files but with different counts would not be counted as shared.
    """
    sorted_path = os.path.join(tmpdir, os.path.basename(path) + ".sorted")
    print(f"  Sorting {os.path.basename(path)}...", file=sys.stderr)
    # Extract k-mer column only, then sort
    awk = subprocess.Popen(["awk", "{print $1}", path], stdout=subprocess.PIPE)
    with open(sorted_path, "w") as fout:
        subprocess.run(
            ["sort", "--parallel=4", f"-T{tmpdir}"],
            stdin=awk.stdout, stdout=fout, check=True
        )
    awk.stdout.close()
    awk.wait()
    return sorted_path


def count_lines(path):
    """Count lines in a file efficiently."""
    result = subprocess.run(["wc", "-l", path], capture_output=True, text=True, check=True)
    return int(result.stdout.strip().split()[0])


def count_intersection(sorted_a, sorted_b):
    """
    Count k-mers present in both files using comm -12 (streaming, no RAM overhead).
    Both files must be pre-sorted on the k-mer column.
    """
    # comm -12 prints lines common to both — pipe directly to wc -l
    comm = subprocess.Popen(
        ["comm", "-12", sorted_a, sorted_b],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL
    )
    wc = subprocess.Popen(
        ["wc", "-l"],
        stdin=comm.stdout, stdout=subprocess.PIPE, text=True
    )
    comm.stdout.close()
    out, _ = wc.communicate()
    comm.wait()
    return int(out.strip())


def jaccard(n_a, n_b, n_shared):
    union = n_a + n_b - n_shared
    return n_shared / union if union > 0 else 0.0


def containment(n_a, n_shared):
    return n_shared / n_a if n_a > 0 else 0.0


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--oann", required=True)
    ap.add_argument("--ofav", required=True)
    ap.add_argument("--ofra", required=True)
    ap.add_argument("--out",  default=None)
    args = ap.parse_args()

    out = open(args.out, "w") if args.out else sys.stdout

    with tempfile.TemporaryDirectory(dir=os.path.dirname(args.oann)) as tmpdir:

        # Sort all three files
        print("Sorting k-mer files...", file=sys.stderr)
        s_oann = sort_kmers(args.oann, tmpdir)
        s_ofav = sort_kmers(args.ofav, tmpdir)
        s_ofra = sort_kmers(args.ofra, tmpdir)

        # Count totals
        print("Counting totals...", file=sys.stderr)
        n = {
            "Oannularis": count_lines(s_oann),
            "Ofaveolata":  count_lines(s_ofav),
            "Ofranksi":    count_lines(s_ofra),
        }
        for sp, cnt in n.items():
            print(f"  {sp}: {cnt:,} k-mers", file=sys.stderr)

        # Pairwise intersections
        pairs = [
            ("Oannularis", "Ofaveolata",  s_oann, s_ofav),
            ("Oannularis", "Ofranksi",    s_oann, s_ofra),
            ("Ofaveolata", "Ofranksi",    s_ofav, s_ofra),
        ]

        print("Computing pairwise intersections...", file=sys.stderr)
        results = []
        for sp1, sp2, f1, f2 in pairs:
            print(f"  {sp1} vs {sp2}...", file=sys.stderr)
            shared = count_intersection(f1, f2)
            j   = jaccard(n[sp1], n[sp2], shared)
            c12 = containment(n[sp1], shared)
            c21 = containment(n[sp2], shared)
            results.append((sp1, sp2, shared, j, c12, c21))
            print(f"    shared={shared:,}  Jaccard={j:.4f}", file=sys.stderr)

    # Write output
    out.write("K-mer sharing between Orbicella species pairs\n")
    out.write("=" * 60 + "\n\n")

    out.write(f"{'Species':30s}  {'N k-mers':>12}\n")
    out.write("-" * 45 + "\n")
    for sp, cnt in n.items():
        out.write(f"{sp:30s}  {cnt:>12,}\n")
    out.write("\n")

    out.write(f"{'Pair':40s}  {'Shared':>10}  {'Jaccard':>8}  {'A in B':>8}  {'B in A':>8}\n")
    out.write("-" * 80 + "\n")
    for sp1, sp2, shared, j, c12, c21 in results:
        label = f"{sp1} vs {sp2}"
        out.write(f"{label:40s}  {shared:>10,}  {j:>8.4f}  {c12:>8.4f}  {c21:>8.4f}\n")

    out.write("\n")
    out.write("Columns:\n")
    out.write("  Shared  = k-mers present in both species\n")
    out.write("  Jaccard = shared / union (symmetric similarity)\n")
    out.write("  A in B  = fraction of species A k-mers found in species B\n")
    out.write("  B in A  = fraction of species B k-mers found in species A\n")

    if args.out:
        out.close()
        print(f"Results written to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
