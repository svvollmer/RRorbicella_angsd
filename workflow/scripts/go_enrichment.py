#!/usr/bin/env python3
"""
go_enrichment.py — GO term enrichment for FST outlier genes.

Downloads NCBI gene2go (if not cached) and runs Fisher's exact test
comparing outlier gene GO terms vs. the full A. palmata gene background.

Requires: goatools  (pip install goatools)
          wget or curl (for downloading gene2go)

Usage:
    python go_enrichment.py \\
        --outlier-ids  results/selection/fst_outlier_genes.gene_ids.txt \\
        --background-gff /projects/vollmer/RR_heat-tolerance/Acropora/reference/jaAcrPala1.3_genomic.gff.gz \\
        --out           results/selection/go_enrichment \\
        --taxid         6131 \\
        --gene2go       results/selection/gene2go.gz
"""
import argparse
import gzip
import os
import re
import subprocess
import sys


def load_gene_ids_from_gff(gff_path):
    """Extract all NCBI GeneIDs from gene features in GFF."""
    ids = set()
    opener = gzip.open if gff_path.endswith(".gz") else open
    with opener(gff_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            m = re.search(r"GeneID:(\d+)", parts[8])
            if m:
                ids.add(m.group(1))
    return ids


def download_gene2go(out_path):
    url = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz"
    print(f"Downloading gene2go from NCBI to {out_path}...", flush=True)
    ret = os.system(f"wget -q -O {out_path} {url} || curl -s -o {out_path} {url}")
    if ret != 0 or not os.path.exists(out_path):
        sys.exit("ERROR: Could not download gene2go. Check network access or provide manually.")
    print("  Done.")


def load_gene2go(gene2go_path, taxid):
    """Load gene2go for a specific taxid. Returns {gene_id: set(GO terms)}."""
    taxid_str = str(taxid)
    gene2go = {}
    opener = gzip.open if gene2go_path.endswith(".gz") else open
    with opener(gene2go_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            tax, gene_id, go_id = parts[0], parts[1], parts[2]
            if tax != taxid_str:
                continue
            gene2go.setdefault(gene_id, set()).add(go_id)
    return gene2go


def run_enrichment(outlier_ids, background_ids, gene2go, out_prefix):
    try:
        from goatools.obo_parser import GODag
        from goatools.go_enrichment import GOEnrichmentStudy
    except ImportError:
        sys.exit(
            "ERROR: goatools not installed.\n"
            "Install with: pip install goatools"
        )

    # Download GO OBO if needed
    obo_path = out_prefix + ".go-basic.obo"
    if not os.path.exists(obo_path):
        print("Downloading go-basic.obo...", flush=True)
        os.system(f"wget -q -O {obo_path} http://geneontology.org/ontology/go-basic.obo "
                  f"|| curl -s -o {obo_path} http://geneontology.org/ontology/go-basic.obo")

    print("Loading GO DAG...", flush=True)
    godag = GODag(obo_path)

    # Build association dict for background
    assoc = {gid: gene2go[gid] for gid in background_ids if gid in gene2go}
    study_ids = {gid for gid in outlier_ids if gid in gene2go}
    background_set = set(background_ids)

    n_outlier_with_go = len(study_ids)
    n_bg_with_go = len(assoc)
    print(f"  Outlier genes with GO annotation: {n_outlier_with_go}/{len(outlier_ids)}")
    print(f"  Background genes with GO annotation: {n_bg_with_go}/{len(background_ids)}")

    if n_outlier_with_go == 0:
        print("\nWARNING: No outlier genes have GO annotations for this taxon.")
        print("The A. palmata genome may have limited GO coverage in gene2go.")
        print("Consider using g:Profiler (https://biit.cs.ut.ee/gprofiler/) with GeneIDs.")
        return

    goeaobj = GOEnrichmentStudy(
        background_set,
        assoc,
        godag,
        propagate_counts=True,
        alpha=0.05,
        methods=["fdr_bh"],
    )

    results = goeaobj.run_study(study_ids)
    results_sig = [r for r in results if r.p_fdr_bh < 0.05 and r.study_count > 0]
    results_sig.sort(key=lambda r: r.p_fdr_bh)

    # Write all results
    all_out = out_prefix + ".go_all.tsv"
    with open(all_out, "w") as f:
        f.write("GO_id\tNS\tname\tstudy_count\tstudy_n\tpop_count\tpop_n\tp_uncorr\tp_fdr_bh\n")
        for r in sorted(results, key=lambda r: r.p_uncorrected):
            if r.study_count == 0:
                continue
            f.write(f"{r.GO}\t{r.NS}\t{r.name}\t{r.study_count}\t{r.study_n}\t"
                    f"{r.pop_count}\t{r.pop_n}\t{r.p_uncorrected:.2e}\t{r.p_fdr_bh:.2e}\n")
    print(f"\nWrote all GO results: {all_out}")

    # Write significant results
    sig_out = out_prefix + ".go_sig.tsv"
    with open(sig_out, "w") as f:
        f.write("GO_id\tNS\tname\tstudy_count\tstudy_n\tpop_count\tpop_n\tp_uncorr\tp_fdr_bh\n")
        for r in results_sig:
            f.write(f"{r.GO}\t{r.NS}\t{r.name}\t{r.study_count}\t{r.study_n}\t"
                    f"{r.pop_count}\t{r.pop_n}\t{r.p_uncorrected:.2e}\t{r.p_fdr_bh:.2e}\n")
    print(f"Wrote significant GO results (FDR<0.05): {sig_out}")

    if results_sig:
        print(f"\nTop significant GO terms (FDR < 0.05, n={len(results_sig)}):")
        print(f"  {'GO ID':<12} {'NS':<3}  {'p_fdr':>8}  {'study':>6}  {'pop':>6}  Term")
        print("  " + "-" * 70)
        for r in results_sig[:30]:
            print(f"  {r.GO:<12} {r.NS:<3}  {r.p_fdr_bh:>8.2e}  "
                  f"{r.study_count:>6}  {r.pop_count:>6}  {r.name}")
    else:
        print("\nNo GO terms significant at FDR < 0.05.")
        print("  (Common for non-model organisms with limited annotation coverage)")
        print("  Try g:Profiler or DAVID with the GeneID list for broader annotation.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outlier-ids",     required=True, help="File of outlier NCBI GeneIDs (one per line)")
    ap.add_argument("--background-gff",  required=True, help="GFF to extract background gene IDs")
    ap.add_argument("--out",             required=True, help="Output prefix")
    ap.add_argument("--taxid",           type=int, default=6131, help="NCBI taxon ID (default: 6131 = A. palmata)")
    ap.add_argument("--gene2go",         default=None, help="Path to gene2go.gz (downloaded if missing)")
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)

    # Load outlier IDs
    with open(args.outlier_ids) as f:
        outlier_ids = set(line.strip() for line in f if line.strip())
    print(f"Outlier gene IDs: {len(outlier_ids)}")

    # Load background from GFF
    print(f"Loading background gene IDs from GFF...", flush=True)
    background_ids = load_gene_ids_from_gff(args.background_gff)
    print(f"  Background: {len(background_ids):,} genes")

    # Download/load gene2go
    gene2go_path = args.gene2go or args.out + ".gene2go.gz"
    if not os.path.exists(gene2go_path):
        download_gene2go(gene2go_path)

    print(f"Loading gene2go for taxid={args.taxid}...", flush=True)
    gene2go = load_gene2go(gene2go_path, args.taxid)
    print(f"  {len(gene2go):,} A. palmata genes with GO annotations in gene2go")

    if len(gene2go) == 0:
        print(f"\nWARNING: No GO annotations found for taxid={args.taxid} in gene2go.")
        print("Possible reasons:")
        print("  1. A. palmata GeneIDs are not in gene2go (common for new assemblies)")
        print("  2. Wrong taxid (try 6131 for A. palmata)")
        print("\nFallback: use g:Profiler web tool with your gene ID list:")
        print("  https://biit.cs.ut.ee/gprofiler/gost")
        print(f"  Input: {args.outlier_ids}")
        print("  Organism: Acropora palmata (if available) or closest coral")
        # Still write background ID list for web tools
        bg_out = args.out + ".background_gene_ids.txt"
        with open(bg_out, "w") as f:
            f.write("\n".join(sorted(background_ids)))
        print(f"\nBackground gene ID list written: {bg_out}")
        return

    run_enrichment(outlier_ids, background_ids, gene2go, args.out)


if __name__ == "__main__":
    main()
