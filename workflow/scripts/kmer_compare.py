"""
kmer_compare.py — Compare k-mer sharing between Orbicella species pairs

After jellyfish dump produces per-species k-mer files, this script computes:
  - Total k-mers per species (above min-count threshold)
  - Pairwise intersection sizes
  - Jaccard similarity: |A∩B| / |A∪B|
  - Containment: |A∩B| / |A|  (fraction of A's k-mers found in B)

Usage:
    python kmer_compare.py --oann oann_kmers.txt \
                           --ofav ofav_kmers.txt \
                           --ofra ofra_kmers.txt \
                           --out kmer_summary.txt
"""

import argparse
import sys


def load_kmers(path):
    """Load k-mers from jellyfish dump output (kmer count per line). Returns set."""
    kmers = set()
    with open(path) as f:
        for line in f:
            parts = line.strip().split()
            if parts:
                kmers.add(parts[0])
    return kmers


def jaccard(a, b):
    intersection = len(a & b)
    union = len(a | b)
    return intersection / union if union > 0 else 0.0


def containment(a, b):
    """Fraction of a's k-mers found in b."""
    if not a:
        return 0.0
    return len(a & b) / len(a)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--oann", required=True, help="Oannularis jellyfish dump file")
    ap.add_argument("--ofav", required=True, help="Ofaveolata jellyfish dump file")
    ap.add_argument("--ofra", required=True, help="Ofranksi jellyfish dump file")
    ap.add_argument("--out", default=None, help="Output file (default: stdout)")
    args = ap.parse_args()

    out = open(args.out, "w") if args.out else sys.stdout

    print("Loading k-mers...", file=sys.stderr)
    species = {
        "Oannularis": load_kmers(args.oann),
        "Ofaveolata": load_kmers(args.ofav),
        "Ofranksi":   load_kmers(args.ofra),
    }

    for sp, kmers in species.items():
        print(f"  {sp}: {len(kmers):,} k-mers", file=sys.stderr)

    print("Computing pairwise comparisons...", file=sys.stderr)

    pairs = [
        ("Oannularis", "Ofaveolata"),
        ("Oannularis", "Ofranksi"),
        ("Ofaveolata", "Ofranksi"),
    ]

    out.write("K-mer sharing between Orbicella species pairs\n")
    out.write("=" * 60 + "\n\n")

    out.write(f"{'Species':30s}  {'N k-mers':>12}\n")
    out.write("-" * 45 + "\n")
    for sp, kmers in species.items():
        out.write(f"{sp:30s}  {len(kmers):>12,}\n")
    out.write("\n")

    out.write(f"{'Pair':40s}  {'Shared':>10}  {'Jaccard':>8}  {'A in B':>8}  {'B in A':>8}\n")
    out.write("-" * 80 + "\n")
    for sp1, sp2 in pairs:
        a = species[sp1]
        b = species[sp2]
        shared = len(a & b)
        j = jaccard(a, b)
        c_ab = containment(a, b)
        c_ba = containment(b, a)
        pair_label = f"{sp1} vs {sp2}"
        out.write(f"{pair_label:40s}  {shared:>10,}  {j:>8.4f}  {c_ab:>8.4f}  {c_ba:>8.4f}\n")

    out.write("\n")
    out.write("Columns:\n")
    out.write("  Shared  = k-mers present in both species\n")
    out.write("  Jaccard = shared / union (symmetric similarity)\n")
    out.write("  A in B  = fraction of species A's k-mers found in species B\n")
    out.write("  B in A  = fraction of species B's k-mers found in species A\n")

    if args.out:
        out.close()
        print(f"Results written to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
