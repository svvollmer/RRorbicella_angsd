#!/usr/bin/env python3
"""
fst_outlier_genes.py — Intersect FST outlier windows with A. palmata gene annotation.

Reads ANGSD sliding-window FST output, identifies outlier windows (top N% or
Z-score threshold), converts to BED, intersects with a GFF via bedtools, and
outputs an annotated gene list.

Usage:
    python fst_outlier_genes.py \\
        --fst results/fst/lineageB_FL_vs_lineageA_FL.fst.windows \\
        --gff /projects/vollmer/RR_heat-tolerance/Acropora/reference/jaAcrPala1.3_genomic.gff.gz \\
        --out results/selection/fst_outlier_genes \\
        --top-pct 1 \\
        --min-sites 1000
"""
import argparse
import os
import re
import subprocess
import sys
import tempfile

import numpy as np


def parse_fst_windows(path, min_sites):
    """Parse ANGSD fst.windows file. Returns list of (chr, start, end, midPos, nsites, fst)."""
    windows = []
    with open(path) as f:
        for line in f:
            if line.startswith("region"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 5:
                continue
            region, chrom, mid, nsites, fst = parts[0], parts[1], int(parts[2]), int(parts[3]), float(parts[4])
            # Parse window coords from last (start,end) group in region field
            m = re.findall(r"\((\d+),(\d+)\)", region)
            if not m:
                continue
            start, end = int(m[-1][0]), int(m[-1][1])
            if nsites < min_sites:
                continue
            windows.append((chrom, start, end, mid, nsites, fst))
    return windows


def write_bed(windows, path):
    with open(path, "w") as f:
        for chrom, start, end, mid, nsites, fst in windows:
            f.write(f"{chrom}\t{start}\t{end}\t{fst:.6f}\t{nsites}\t{mid}\n")


def run_bedtools_intersect(bed_path, gff_path, out_path):
    """Intersect BED with GFF gene features. Requires bedtools in PATH or via module."""
    # Extract gene-only lines from GFF to a temp file
    tmp_genes = out_path + ".tmp_genes.gff"
    if gff_path.endswith(".gz"):
        cat_cmd = f"zcat {gff_path}"
    else:
        cat_cmd = f"cat {gff_path}"

    os.system(f'{cat_cmd} | awk \'$3=="gene"\' > {tmp_genes}')

    cmd = (
        f"bedtools intersect -a {bed_path} -b {tmp_genes} -wa -wb "
        f"> {out_path}"
    )
    ret = os.system(cmd)
    os.remove(tmp_genes)
    if ret != 0:
        sys.exit("ERROR: bedtools intersect failed. Is bedtools in PATH? Try: module load bedtools/2.29.0")
    return out_path


def parse_gff_attr(attr_str, key):
    m = re.search(rf"{key}=([^;]+)", attr_str)
    return m.group(1) if m else ""


def parse_intersect(intersect_path):
    """Parse bedtools intersect output (BED6 + GFF columns)."""
    genes = []
    seen = set()
    with open(intersect_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 12:
                continue
            # BED cols: chr start end fst nsites mid
            # GFF cols: chr source type gff_start gff_end score strand phase attributes
            chrom = parts[0]
            fst   = float(parts[3])
            nsites = int(parts[4])
            mid   = int(parts[5])
            gff_start = int(parts[9])
            gff_end   = int(parts[10])
            attrs = parts[14] if len(parts) > 14 else parts[-1]

            gene_id   = parse_gff_attr(attrs, "ID").replace("gene-", "")
            name      = parse_gff_attr(attrs, "Name")
            desc      = parse_gff_attr(attrs, "description").replace("%2C", ",").replace("%20", " ")
            ncbi_id   = parse_gff_attr(attrs, "GeneID")
            biotype   = parse_gff_attr(attrs, "gene_biotype")

            key = (ncbi_id or gene_id, chrom, gff_start)
            if key in seen:
                continue
            seen.add(key)

            genes.append(dict(
                chr=chrom, gene_start=gff_start, gene_end=gff_end,
                gene_id=gene_id, name=name, description=desc,
                ncbi_gene_id=ncbi_id, gene_biotype=biotype,
                window_mid=mid, fst=fst, nsites=nsites,
            ))
    return genes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fst",       required=True, help="ANGSD fst.windows file")
    ap.add_argument("--gff",       required=True, help="Gene annotation GFF/GFF.gz")
    ap.add_argument("--out",       required=True, help="Output prefix")
    ap.add_argument("--top-pct",   type=float, default=1.0,
                    help="Top N%% of windows as outliers (default: 1)")
    ap.add_argument("--min-sites", type=int, default=1000,
                    help="Min sites per window (default: 1000)")
    ap.add_argument("--fst-floor", type=float, default=0.0,
                    help="Minimum absolute FST (applied after top-pct filter, default: 0)")
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)

    print(f"Loading FST windows from {args.fst}...", flush=True)
    windows = parse_fst_windows(args.fst, args.min_sites)
    print(f"  {len(windows):,} windows passing min_sites={args.min_sites}")

    fst_vals = np.array([w[5] for w in windows])
    threshold = np.percentile(fst_vals, 100 - args.top_pct)
    threshold = max(threshold, args.fst_floor)
    outliers = [w for w in windows if w[5] >= threshold]
    print(f"  Top {args.top_pct}% threshold: FST ≥ {threshold:.4f}")
    print(f"  {len(outliers):,} outlier windows")

    # Also compute per-chromosome thresholds for reference
    chroms = sorted(set(w[0] for w in windows))
    print("\nPer-chromosome FST summary:")
    print(f"  {'Chr':<15} {'N':>6}  {'Mean FST':>9}  {'99th pct':>9}  {'Max FST':>9}")
    for chrom in chroms:
        cv = np.array([w[5] for w in windows if w[0] == chrom])
        print(f"  {chrom:<15} {len(cv):>6}  {cv.mean():>9.4f}  {np.percentile(cv,99):>9.4f}  {cv.max():>9.4f}")

    # Write BED files
    all_bed     = args.out + ".all.bed"
    outlier_bed = args.out + ".outliers.bed"
    write_bed(windows, all_bed)
    write_bed(outliers, outlier_bed)
    print(f"\nWrote: {all_bed}")
    print(f"Wrote: {outlier_bed}")

    # Bedtools intersect
    print("\nIntersecting with gene annotation...", flush=True)
    intersect_out = args.out + ".intersect.txt"
    run_bedtools_intersect(outlier_bed, args.gff, intersect_out)

    genes = parse_intersect(intersect_out)
    genes.sort(key=lambda g: g["fst"], reverse=True)
    print(f"  {len(genes):,} unique genes in FST outlier windows")

    # Write gene table
    tsv_out = args.out + ".genes.tsv"
    with open(tsv_out, "w") as f:
        f.write("chr\tgene_start\tgene_end\tncbi_gene_id\tgene_name\tdescription\t"
                "gene_biotype\twindow_mid\tfst\tnsites\n")
        for g in genes:
            f.write(f"{g['chr']}\t{g['gene_start']}\t{g['gene_end']}\t"
                    f"{g['ncbi_gene_id']}\t{g['name']}\t{g['description']}\t"
                    f"{g['gene_biotype']}\t{g['window_mid']}\t{g['fst']:.6f}\t{g['nsites']}\n")
    print(f"Wrote: {tsv_out}")

    # Write gene ID list for GO enrichment
    id_out = args.out + ".gene_ids.txt"
    with open(id_out, "w") as f:
        for g in genes:
            if g["ncbi_gene_id"]:
                f.write(g["ncbi_gene_id"] + "\n")
    print(f"Wrote: {id_out}  (NCBI GeneIDs for GO enrichment)")

    # Print top 20
    print(f"\nTop 20 genes by FST:")
    print(f"  {'Chr':<15} {'GeneID':>12}  {'Name':<20}  {'FST':>7}  {'Description'}")
    print("  " + "-" * 95)
    for g in genes[:20]:
        desc = g["description"][:45] + "…" if len(g["description"]) > 45 else g["description"]
        print(f"  {g['chr']:<15} {g['ncbi_gene_id']:>12}  {g['name']:<20}  {g['fst']:>7.4f}  {desc}")


if __name__ == "__main__":
    main()
