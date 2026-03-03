#!/usr/bin/env python3
"""
Population Genomics Report Generator
Reads Snakemake pipeline outputs and generates a self-contained HTML report
with embedded figures, plus saves PNG files to results/figures/.

Usage (called by Snakemake generate_report rule):
    python generate_report.py \
        --results_dir results/ \
        --metadata config/samples_aws_test.csv \
        --config config/config.yaml \
        --output results/report.html \
        --figures_dir results/figures
"""

import argparse
import base64
import gzip
import os
import sys
from io import BytesIO
from pathlib import Path

import numpy as np
import pandas as pd

# ─── Plotting setup ──────────────────────────────────────────────────────────

def setup_matplotlib():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 12,
        "axes.linewidth": 1.2,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.major.width": 1.2,
        "ytick.major.width": 1.2,
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.facecolor": "white",
    })
    return plt, gridspec

POP_COLORS = ["#2166AC", "#D6604D", "#4DAC26", "#8073AC", "#01665E", "#8C510A"]

def pop_color(pop, pops):
    idx = sorted(pops).index(pop) if pop in pops else 0
    return POP_COLORS[idx % len(POP_COLORS)]


def save_fig(fig, figures_dir, name):
    """Save figure to disk and return base64-encoded HTML <img> tag."""
    plt, _ = setup_matplotlib()
    # Save PNG
    png_path = Path(figures_dir) / f"{name}.png"
    fig.savefig(png_path, dpi=300, bbox_inches="tight", facecolor="white")
    # Encode for HTML
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight", facecolor="white")
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode()
    plt.close(fig)
    return f'<img src="data:image/png;base64,{b64}" style="max-width:100%;margin:12px 0;">'


# ─── Data loaders ─────────────────────────────────────────────────────────────

def load_metadata(path):
    df = pd.read_csv(path).set_index("sample_id")
    return df


def load_filtering_summary(results_dir):
    """Load pre-built filtering_summary.txt (Sample, Total_reads, Mapped_reads, Mapping_rate, Mean_depth)."""
    path = Path(results_dir) / "filtering" / "filtering_summary.txt"
    if not path.exists():
        return pd.DataFrame()
    df = pd.read_csv(path, sep="\t")
    df.columns = [c.lower() for c in df.columns]
    df = df.rename(columns={"sample": "sample_id"})
    return df.set_index("sample_id")


def load_fastp_stats(results_dir, samples):
    import json
    records = []
    for sample in samples:
        p = Path(results_dir) / "qc" / f"{sample}_fastp.json"
        if not p.exists():
            continue
        with open(p) as f:
            d = json.load(f)
        s = d["summary"]
        records.append({
            "sample_id": sample,
            "raw_reads":   s["before_filtering"]["total_reads"],
            "clean_reads": s["after_filtering"]["total_reads"],
            "q30_rate":    s["after_filtering"]["q30_rate"],
            "gc_content":  s["after_filtering"]["gc_content"],
        })
    if not records:
        return pd.DataFrame()
    return pd.DataFrame(records).set_index("sample_id")


def load_depth_summary(results_dir):
    path = Path(results_dir) / "qc" / "depth_summary.txt"
    if not path.exists():
        return pd.DataFrame()
    df = pd.read_csv(path, sep="\t")
    df.columns = ["sample_id", "mean_depth"]
    return df.set_index("sample_id")


def load_pca(results_dir):
    path = Path(results_dir) / "pca" / "pcangsd.cov"
    if not path.exists():
        return None, None
    cov = np.loadtxt(path)
    vals, vecs = np.linalg.eigh(cov)
    idx = np.argsort(vals)[::-1]
    vals, vecs = vals[idx], vecs[:, idx]
    pct = vals / np.sum(np.abs(vals)) * 100
    return vecs, pct


def load_admixture(results_dir, max_k=5):
    """Load PCAngsd admixture Q matrices (pcangsd_K{k}.Q)."""
    admix = {}
    for k in range(2, max_k + 1):
        p = Path(results_dir) / "admixture" / f"pcangsd_K{k}.Q"
        if p.exists():
            admix[k] = np.loadtxt(p)
    return admix


def load_fst_global(results_dir, pop_pairs):
    records = []
    for p1, p2 in pop_pairs:
        path = Path(results_dir) / "fst" / f"{p1}_vs_{p2}.fst.global"
        if not path.exists():
            continue
        with open(path) as f:
            vals = f.readline().strip().split()
        if len(vals) >= 2:
            records.append({
                "comparison": f"{p1} vs {p2}",
                "fst_unweighted": float(vals[0]),
                "fst_weighted":   float(vals[1]),
            })
    return pd.DataFrame(records) if records else pd.DataFrame()


def load_windowed_fst(results_dir, p1, p2):
    """realSFS fst stats2 output: region chr midPos Nsites fst"""
    path = Path(results_dir) / "fst" / f"{p1}_vs_{p2}.fst.windows"
    if not path.exists():
        return None
    # realSFS fst stats2 emits a 4-column header (no "fst" field name) followed
    # by 5-column data rows.  Skip the header and supply our own column names.
    df = pd.read_csv(path, sep="\t", header=None, skiprows=1,
                     names=["region", "chr", "midPos", "Nsites", "fst"])
    return df


def load_thetas(results_dir, popname):
    path = Path(results_dir) / "diversity" / f"{popname}.thetas.idx.pestPG"
    if not path.exists():
        return None
    df = pd.read_csv(path, sep="\t")
    return df


def load_sfs(results_dir, popname):
    path = Path(results_dir) / "sfs" / f"{popname}.sfs"
    if not path.exists():
        return None
    with open(path) as f:
        return np.array([float(x) for x in f.readline().strip().split()])


def load_heterozygosity(results_dir, samples):
    """realSFS 1D output: two numbers (hom, het) per sample."""
    records = []
    for sample in samples:
        path = Path(results_dir) / "heterozygosity" / f"{sample}.het"
        if not path.exists():
            continue
        with open(path) as f:
            vals = [float(x) for x in f.readline().strip().split()]
        if len(vals) >= 3:
            total = vals[0] + vals[1] + vals[2]
        elif len(vals) >= 2:
            total = vals[0] + vals[1]
        else:
            continue
        het_rate = vals[1] / total if total > 0 else 0
        records.append({"sample_id": sample,
                         "heterozygosity": het_rate,
                         "n_sites": int(total)})
    return pd.DataFrame(records).set_index("sample_id") if records else pd.DataFrame()


def load_relatedness_matrix(results_dir):
    path = Path(results_dir) / "relatedness" / "relatedness_matrix.txt"
    if not path.exists():
        return None
    return pd.read_csv(path, sep="\t", index_col=0)


def load_ld_decay(results_dir):
    path = Path(results_dir) / "ld" / "ld_decay.csv"
    if not path.exists():
        return None
    return pd.read_csv(path)


def count_snps(results_dir):
    p1 = Path(results_dir) / "angsd" / "pass1_snps.txt"
    p2 = Path(results_dir) / "angsd" / "all.mafs.gz"
    n_pass1, n_pass2 = None, None
    if p1.exists():
        n_pass1 = sum(1 for _ in open(p1))
    if p2.exists():
        with gzip.open(p2, "rt") as f:
            next(f)
            n_pass2 = sum(1 for _ in f)
    return n_pass1, n_pass2


# ─── Figures ──────────────────────────────────────────────────────────────────

def fig_mapping_rates(filt_df, metadata, figures_dir):
    plt, _ = setup_matplotlib()
    if filt_df.empty or "mapping_rate" not in filt_df.columns:
        return ""
    pops = metadata["population"].unique()
    df = filt_df.join(metadata[["population"]], how="inner").sort_values(
        ["population", "mapping_rate"])
    colors = [pop_color(p, pops) for p in df["population"]]

    fig, ax = plt.subplots(figsize=(max(10, len(df) * 0.7), 5))
    ax.bar(range(len(df)), df["mapping_rate"] * 100, color=colors, edgecolor="white")
    ax.axhline(50, color="red", linestyle="--", linewidth=1, label="50% threshold")
    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(df.index, rotation=45, ha="right", fontsize=9)
    ax.set_ylabel("Mapping Rate (%)")
    ax.set_ylim(0, 105)
    ax.set_title("Per-Sample Mapping Rate")
    # legend for populations
    from matplotlib.patches import Patch
    handles = [Patch(color=pop_color(p, pops), label=p) for p in sorted(pops)]
    ax.legend(handles=handles, title="Population")
    return save_fig(fig, figures_dir, "mapping_rates")


def fig_depth(depth_df, metadata, figures_dir):
    plt, _ = setup_matplotlib()
    if depth_df.empty:
        return ""
    pops = metadata["population"].unique()
    df = depth_df.join(metadata[["population"]], how="inner").sort_values(
        ["population", "mean_depth"])
    colors = [pop_color(p, pops) for p in df["population"]]

    fig, ax = plt.subplots(figsize=(max(10, len(df) * 0.7), 5))
    ax.bar(range(len(df)), df["mean_depth"].astype(float), color=colors, edgecolor="white")
    ax.axhline(10, color="red", linestyle="--", linewidth=1, label="10× threshold")
    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(df.index, rotation=45, ha="right", fontsize=9)
    ax.set_ylabel("Mean Depth (×)")
    ax.set_title("Per-Sample Sequencing Depth")
    from matplotlib.patches import Patch
    handles = [Patch(color=pop_color(p, pops), label=p) for p in sorted(pops)]
    ax.legend(handles=handles, title="Population")
    return save_fig(fig, figures_dir, "depth")


def fig_pca(eigenvectors, pct_var, metadata, figures_dir):
    plt, _ = setup_matplotlib()
    if eigenvectors is None:
        return ""
    samples = metadata.index.tolist()
    pops = metadata["population"].unique()

    fig, ax = plt.subplots(figsize=(8, 7))
    for pop in sorted(pops):
        mask = metadata["population"] == pop
        idx = [i for i, s in enumerate(samples) if mask.get(s, False)]
        ax.scatter(eigenvectors[idx, 0], eigenvectors[idx, 1],
                   c=pop_color(pop, pops), label=pop,
                   s=100, edgecolors="black", linewidth=0.6, alpha=0.9, zorder=3)
    for i, s in enumerate(samples):
        ax.annotate(s, (eigenvectors[i, 0], eigenvectors[i, 1]),
                    fontsize=7, alpha=0.65, xytext=(5, 5),
                    textcoords="offset points")
    ax.set_xlabel(f"PC1 ({pct_var[0]:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({pct_var[1]:.1f}% variance)")
    ax.set_title("PCA — LD-pruned SNPs")
    ax.legend(title="Population", frameon=True)
    ax.grid(True, alpha=0.25)
    return save_fig(fig, figures_dir, "pca")


def fig_admixture(admix, metadata, figures_dir):
    plt, gs_mod = setup_matplotlib()
    if not admix:
        return ""
    samples = metadata.index.tolist()
    pops = metadata["population"].unique()
    k_values = sorted(admix.keys())

    # Order samples by population then name
    order = metadata.sort_values("population").index
    order_idx = [samples.index(s) for s in order if s in samples]

    fig = plt.figure(figsize=(max(12, len(samples) * 0.6), 3 * len(k_values)))
    gs = gs_mod.GridSpec(len(k_values), 1, hspace=0.4, figure=fig)

    cluster_colors = ["#2166AC", "#D6604D", "#4DAC26", "#8073AC",
                      "#01665E", "#8C510A", "#F4A582", "#92C5DE"]

    for ki, k in enumerate(k_values):
        ax = fig.add_subplot(gs[ki])
        q = admix[k][order_idx, :]
        bottom = np.zeros(len(order_idx))
        for j in range(k):
            ax.bar(range(len(order_idx)), q[:, j], bottom=bottom,
                   color=cluster_colors[j % len(cluster_colors)],
                   width=1.0, edgecolor="none")
            bottom += q[:, j]
        ax.set_xlim(-0.5, len(order_idx) - 0.5)
        ax.set_ylim(0, 1)
        ax.set_ylabel(f"K = {k}", fontsize=11)
        ax.set_yticks([0, 0.5, 1])

        # Population separator lines
        pop_order = metadata.loc[order, "population"]
        prev = pop_order.iloc[0]
        for i, p in enumerate(pop_order):
            if p != prev:
                ax.axvline(x=i - 0.5, color="black", linewidth=1.5)
                prev = p

        if ki == len(k_values) - 1:
            ax.set_xticks(range(len(order_idx)))
            ax.set_xticklabels(list(order), rotation=90, fontsize=8)
        else:
            ax.set_xticks([])

    fig.suptitle("Admixture Analysis (PCAngsd)", fontsize=14, y=1.01)
    return save_fig(fig, figures_dir, "admixture")


def fig_fst_manhattan(p1, p2, wfst_df, figures_dir):
    plt, _ = setup_matplotlib()
    if wfst_df is None or wfst_df.empty:
        return ""
    df = wfst_df[wfst_df["fst"] >= 0].copy()
    if df.empty:
        return ""

    chroms = df["chr"].unique()
    offset = {}
    cum = 0
    for c in chroms:
        offset[c] = cum
        cum += df[df["chr"] == c]["midPos"].max()
    df["cumpos"] = df.apply(lambda r: r["midPos"] + offset[r["chr"]], axis=1)

    thresh = df["fst"].quantile(0.99)
    chrom_colors = {c: ("#2166AC" if i % 2 == 0 else "#6BAED6")
                    for i, c in enumerate(chroms)}
    point_colors = df.apply(
        lambda r: "#D6604D" if r["fst"] >= thresh else chrom_colors[r["chr"]], axis=1)

    fig, ax = plt.subplots(figsize=(16, 5))
    ax.scatter(df["cumpos"], df["fst"], c=point_colors, s=4, alpha=0.7, linewidths=0)
    ax.axhline(thresh, color="red", linestyle="--", linewidth=1,
               label=f"Top 1% (FST > {thresh:.3f})")

    centers = {c: df[df["chr"] == c]["cumpos"].median() for c in chroms}
    ax.set_xticks(list(centers.values()))
    ax.set_xticklabels([str(c).split("_")[0] for c in chroms], rotation=45, fontsize=7)
    ax.set_ylabel("FST")
    ax.set_title(f"Windowed FST — {p1.title()} vs {p2.title()} (50 kb windows)")
    ax.legend(fontsize=9)
    return save_fig(fig, figures_dir, f"fst_manhattan_{p1}_vs_{p2}")


def fig_heterozygosity(het_df, metadata, figures_dir):
    plt, _ = setup_matplotlib()
    if het_df is None or het_df.empty:
        return ""
    pops = metadata["population"].unique()
    df = het_df.join(metadata[["population"]], how="inner").sort_values(
        ["population", "heterozygosity"])
    colors = [pop_color(p, pops) for p in df["population"]]
    mean_h = df["heterozygosity"].mean()
    sd_h = df["heterozygosity"].std()

    fig, ax = plt.subplots(figsize=(max(10, len(df) * 0.7), 5))
    ax.bar(range(len(df)), df["heterozygosity"] * 100, color=colors, edgecolor="white")
    ax.axhline(mean_h * 100, color="gray", linewidth=1, linestyle="-", alpha=0.6, label="Mean")
    ax.axhline((mean_h + 2 * sd_h) * 100, color="red", linewidth=1,
               linestyle="--", alpha=0.7, label="±2 SD")
    ax.axhline((mean_h - 2 * sd_h) * 100, color="red", linewidth=1,
               linestyle="--", alpha=0.7)
    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(df.index, rotation=45, ha="right", fontsize=9)
    ax.set_ylabel("Heterozygosity (%)")
    ax.set_title("Per-Individual Heterozygosity")
    from matplotlib.patches import Patch
    handles = [Patch(color=pop_color(p, pops), label=p) for p in sorted(pops)]
    handles += [plt.Line2D([0], [0], color="gray", label="Mean"),
                plt.Line2D([0], [0], color="red", linestyle="--", label="±2 SD")]
    ax.legend(handles=handles)
    return save_fig(fig, figures_dir, "heterozygosity")


def fig_sfs(sfs_dict, figures_dir):
    plt, _ = setup_matplotlib()
    pops = [p for p, s in sfs_dict.items() if s is not None]
    if not pops:
        return ""
    fig, axes = plt.subplots(1, len(pops), figsize=(6 * len(pops), 5))
    if len(pops) == 1:
        axes = [axes]
    for i, pop in enumerate(pops):
        sfs = sfs_dict[pop]
        n = len(sfs)
        folded = sfs[: n // 2 + 1].copy()
        for j in range(1, n // 2):
            folded[j] += sfs[n - j]
        folded = folded[1:]
        axes[i].bar(range(1, len(folded) + 1), folded,
                    color=pop_color(pop, pops), edgecolor="white")
        axes[i].set_xlabel("Minor Allele Count")
        axes[i].set_ylabel("Number of Sites")
        axes[i].set_title(f"Folded SFS — {pop.title()}")
    plt.suptitle("Site Frequency Spectra", fontsize=14, y=1.02)
    return save_fig(fig, figures_dir, "sfs")


def fig_diversity(thetas_dict, figures_dir):
    plt, _ = setup_matplotlib()
    pops = [p for p, d in thetas_dict.items() if d is not None and len(d) > 0]
    if not pops:
        return ""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    metrics = [
        ("tP",     "Nucleotide Diversity (π/site)"),
        ("tW",     "Watterson's θ (per site)"),
        ("Tajima", "Tajima's D"),
    ]
    for ax, (col, title) in zip(axes, metrics):
        means, sems, labels = [], [], []
        for pop in pops:
            df = thetas_dict[pop]
            if col not in df.columns:
                continue
            if "nSites" in df.columns and col != "Tajima":
                vals = df[col] / df["nSites"]
            else:
                vals = df[col]
            means.append(vals.mean())
            sems.append(vals.sem())
            labels.append(pop.title())
        bars = ax.bar(range(len(labels)), means, yerr=sems, capsize=5,
                      color=[pop_color(p, pops) for p in pops],
                      edgecolor="white")
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels)
        ax.set_title(title)
    return save_fig(fig, figures_dir, "diversity")


def fig_kinship(rel_matrix, figures_dir):
    plt, _ = setup_matplotlib()
    if rel_matrix is None:
        return ""
    import matplotlib.colors as mcolors
    fig, ax = plt.subplots(figsize=(max(8, len(rel_matrix) * 0.6),
                                     max(7, len(rel_matrix) * 0.6)))
    mat = rel_matrix.values.astype(float)
    np.fill_diagonal(mat, 0)
    im = ax.imshow(mat, cmap="YlOrRd", vmin=0, vmax=0.5)
    ax.set_xticks(range(len(rel_matrix)))
    ax.set_yticks(range(len(rel_matrix)))
    ax.set_xticklabels(rel_matrix.columns, rotation=90, fontsize=8)
    ax.set_yticklabels(rel_matrix.index, fontsize=8)
    plt.colorbar(im, ax=ax, label="KING Kinship", shrink=0.8)
    ax.set_title("Pairwise Kinship (KING)")
    return save_fig(fig, figures_dir, "kinship")


def fig_ld_decay(ld_df, figures_dir):
    plt, _ = setup_matplotlib()
    if ld_df is None or ld_df.empty:
        return ""
    fig, ax = plt.subplots(figsize=(9, 5))
    # Expect columns: distance_kb, mean_r2 (or similar)
    xcol = next((c for c in ld_df.columns if "dist" in c.lower()), ld_df.columns[0])
    ycol = next((c for c in ld_df.columns if "r2" in c.lower() or "ld" in c.lower()),
                ld_df.columns[1])
    ax.plot(ld_df[xcol], ld_df[ycol], color="#2166AC", linewidth=2)
    ax.axhline(0.1, color="gray", linestyle="--", linewidth=1, label="r² = 0.1")
    ax.set_xlabel("Distance (kb)")
    ax.set_ylabel("Mean r²")
    ax.set_title("LD Decay")
    ax.legend()
    return save_fig(fig, figures_dir, "ld_decay")


# ─── HTML assembly ────────────────────────────────────────────────────────────

CSS = """
<style>
body { font-family: Arial, Helvetica, sans-serif; max-width: 1200px;
       margin: 0 auto; padding: 40px 40px 60px; color: #222; line-height: 1.6; }
h1   { color: #1565C0; border-bottom: 3px solid #1565C0; padding-bottom: 10px; }
h2   { color: #1565C0; border-bottom: 1px solid #ddd; padding-bottom: 6px; margin-top: 40px; }
h3   { color: #1976D2; margin-top: 24px; }
.summary { background:#E3F2FD; border-left:4px solid #1565C0; padding:18px 22px;
           margin:20px 0; border-radius:0 8px 8px 0; font-size:1.05em; }
.warn    { background:#FFF3E0; border-left:4px solid #FF9800; padding:14px 18px;
           margin:14px 0; border-radius:0 8px 8px 0; }
.ok      { background:#E8F5E9; border-left:4px solid #4CAF50; padding:14px 18px;
           margin:14px 0; border-radius:0 8px 8px 0; }
table { border-collapse:collapse; width:100%; margin:14px 0; font-size:0.93em; }
th { background:#1976D2; color:white; padding:9px 12px; text-align:left; }
td { border:1px solid #ddd; padding:8px 12px; }
tr:nth-child(even) { background:#F5F5F5; }
tr:hover { background:#E3F2FD; }
.metric { display:inline-block; background:#F5F5F5; border-radius:8px;
          padding:14px 22px; margin:8px; text-align:center; min-width:120px; }
.metric-val  { font-size:1.8em; font-weight:bold; color:#1565C0; }
.metric-lbl  { font-size:0.85em; color:#666; }
pre  { background:#F5F5F5; padding:14px; border-radius:6px; font-size:0.88em;
       white-space:pre-wrap; }
.footer { margin-top:50px; padding-top:18px; border-top:1px solid #ddd;
          color:#999; font-size:0.83em; }
</style>
"""

def metric_box(value, label):
    return f'<div class="metric"><div class="metric-val">{value}</div><div class="metric-lbl">{label}</div></div>'

def df_to_html(df, float_fmt="{:.4f}"):
    if df is None or df.empty:
        return "<p><em>No data available.</em></p>"
    return df.to_html(classes="", float_format=lambda x: float_fmt.format(x))


def generate_methods(pipeline_config, metadata):
    n = len(metadata)
    pops = metadata["population"].unique()
    return f"""Whole-genome sequencing data from {n} samples spanning {len(pops)} populations
({', '.join(p.title() for p in sorted(pops))}) were processed using a Snakemake pipeline.
Raw reads were quality-trimmed with fastp (adapter auto-detection).
Trimmed reads were mapped to the <i>Acropora palmata</i> reference genome
(GCF_964030605.1; Wellcome Sanger Institute) using BWA-MEM;
duplicates were marked with samtools markdup.
Genotype likelihoods were estimated with ANGSD (GL model {pipeline_config.get('gl_model', 1)},
minMapQ {pipeline_config.get('minmapq', 20)}, minQ {pipeline_config.get('minq', 20)},
MAF &gt; {pipeline_config.get('min_maf', 0.10)},
SNP p-value &lt; {pipeline_config.get('snp_pval', '1e-6')},
minInd {int(pipeline_config.get('min_ind_frac', 0.8)*100)}%).
Population structure was assessed by PCA and admixture analysis (K = 2–{pipeline_config.get('max_k', 5)}) using PCAngsd.
Pairwise relatedness was estimated with NgsRelate (KING statistic).
Site frequency spectra, nucleotide diversity (π), Watterson's θ, and Tajima's D
were calculated per population using realSFS and ANGSD.
Pairwise FST was estimated using the 2D-SFS approach in realSFS
with 50 kb sliding windows (10 kb step).
LD decay was estimated with ngsLD; SNPs were pruned using a greedy graph algorithm (r² threshold {pipeline_config.get('ld_r2_threshold', 0.3)})."""


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results_dir", required=True)
    parser.add_argument("--metadata",    required=True)
    parser.add_argument("--config",      default=None)
    parser.add_argument("--output",      default="results/report.html")
    parser.add_argument("--figures_dir", default="results/figures")
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    figures_dir = Path(args.figures_dir)
    figures_dir.mkdir(parents=True, exist_ok=True)

    # Load pipeline config
    pipeline_config = {}
    if args.config and Path(args.config).exists():
        import yaml
        with open(args.config) as f:
            pipeline_config = yaml.safe_load(f) or {}

    print("Loading metadata...")
    metadata = load_metadata(args.metadata)
    samples = metadata.index.tolist()
    pops = sorted(metadata["population"].unique())
    from itertools import combinations
    pop_pairs = list(combinations(pops, 2))

    print("Loading pipeline outputs...")
    filt_df  = load_filtering_summary(results_dir)
    depth_df = load_depth_summary(results_dir)
    fastp_df = load_fastp_stats(results_dir, samples)
    eigenvectors, pct_var = load_pca(results_dir)
    admix    = load_admixture(results_dir, pipeline_config.get("max_k", 5))
    fst_df   = load_fst_global(results_dir, pop_pairs)
    het_df   = load_heterozygosity(results_dir, samples)
    rel_mat  = load_relatedness_matrix(results_dir)
    ld_df    = load_ld_decay(results_dir)
    thetas   = {p: load_thetas(results_dir, p) for p in pops}
    sfs      = {p: load_sfs(results_dir, p) for p in pops}
    n_pass1, n_pass2 = count_snps(results_dir)

    # Build PCA-ordered metadata: only unrelated samples, in bamlist row order.
    # PCA eigenvector row i corresponds to sample i in unrelated_samples.txt.
    # Using full metadata for PCA plots causes IndexError when some samples
    # were excluded as clones/related.
    unrel_path = results_dir / "relatedness" / "unrelated_samples.txt"
    if unrel_path.exists():
        pca_sample_order = [l.strip() for l in open(unrel_path) if l.strip()]
        pca_meta = metadata.reindex([s for s in pca_sample_order if s in metadata.index])
    else:
        # Fallback: clip to eigenvector count if no unrelated list available
        n_pca = eigenvectors.shape[0] if eigenvectors is not None else len(samples)
        pca_meta = metadata.iloc[:n_pca]

    print("Generating figures...")
    img_mapping  = fig_mapping_rates(filt_df, metadata, figures_dir)
    img_depth    = fig_depth(depth_df, metadata, figures_dir)
    img_pca      = fig_pca(eigenvectors, pct_var, pca_meta, figures_dir)
    img_admix    = fig_admixture(admix, pca_meta, figures_dir)
    img_het      = fig_heterozygosity(het_df, metadata, figures_dir)
    img_sfs      = fig_sfs(sfs, figures_dir)
    img_div      = fig_diversity(thetas, figures_dir)
    img_kinship  = fig_kinship(rel_mat, figures_dir)
    img_ld       = fig_ld_decay(ld_df, figures_dir)
    img_fst_man  = {}
    for p1, p2 in pop_pairs:
        wfst = load_windowed_fst(results_dir, p1, p2)
        img_fst_man[(p1, p2)] = fig_fst_manhattan(p1, p2, wfst, figures_dir)

    # QC flags
    flagged = []
    if not filt_df.empty and "mapping_rate" in filt_df.columns:
        flagged = filt_df[filt_df["mapping_rate"] < 0.5].index.tolist()
    if not depth_df.empty:
        low_depth = depth_df[depth_df["mean_depth"].astype(float) < 10].index.tolist()
        flagged = list(set(flagged + low_depth))

    print("Assembling report...")
    parts = [f"<!DOCTYPE html><html><head><meta charset='utf-8'>"
             f"<title>Population Genomics Report</title>{CSS}</head><body>"]
    parts.append("<h1>Population Genomics Report — <em>Acropora palmata</em></h1>")
    parts.append(f"<p style='color:#666'>Generated from {len(samples)} samples "
                 f"({', '.join(f'{(metadata.population==p).sum()} {p.title()}' for p in pops)})</p>")

    # ── Overview metrics
    parts.append("<h2>1. Overview</h2>")
    boxes = [metric_box(len(samples), "Samples"),
             metric_box(len(pops), "Populations")]
    if n_pass1:
        boxes.append(metric_box(f"{n_pass1:,}", "Pass-1 SNPs"))
    if n_pass2:
        boxes.append(metric_box(f"{n_pass2:,}", "Pass-2 SNPs"))
    if not filt_df.empty and "mapping_rate" in filt_df.columns:
        mr = filt_df["mapping_rate"].mean()
        boxes.append(metric_box(f"{mr:.1%}", "Mean Mapping Rate"))
    if not depth_df.empty:
        md = depth_df["mean_depth"].astype(float).mean()
        boxes.append(metric_box(f"{md:.1f}×", "Mean Depth"))
    parts.append("".join(boxes))

    if flagged:
        parts.append(f'<div class="warn">⚠ {len(flagged)} sample(s) flagged for low mapping rate or depth: '
                     f'{", ".join(flagged)}</div>')
    else:
        parts.append('<div class="ok">✓ All samples passed QC thresholds</div>')

    # ── Sequencing QC
    parts.append("<h2>2. Sequencing &amp; Mapping Quality</h2>")
    parts.append("<h3>2a. Filtering Summary</h3>")
    parts.append(df_to_html(filt_df, "{:.4f}"))
    parts.append(img_mapping)
    parts.append("<h3>2b. Sequencing Depth</h3>")
    parts.append(img_depth)
    parts.append("<h3>2c. fastp Read Quality</h3>")
    parts.append(df_to_html(fastp_df, "{:.4f}"))

    # ── Population structure
    parts.append("<h2>3. Population Structure</h2>")
    parts.append("<h3>3a. PCA</h3>")
    parts.append(img_pca or "<p><em>PCA not yet available.</em></p>")
    parts.append("<h3>3b. Admixture</h3>")
    parts.append(img_admix or "<p><em>Admixture not yet available.</em></p>")

    # ── Diversity
    parts.append("<h2>4. Genetic Diversity</h2>")
    parts.append("<h3>4a. Site Frequency Spectra</h3>")
    parts.append(img_sfs or "<p><em>SFS not yet available.</em></p>")
    parts.append("<h3>4b. π, θ, Tajima's D</h3>")
    parts.append(img_div or "<p><em>Diversity stats not yet available.</em></p>")
    parts.append("<h3>4c. Individual Heterozygosity</h3>")
    parts.append(img_het or "<p><em>Heterozygosity not yet available.</em></p>")
    if het_df is not None and not het_df.empty:
        parts.append(df_to_html(het_df, "{:.6f}"))

    # ── Differentiation
    parts.append("<h2>5. Population Differentiation (FST)</h2>")
    if not fst_df.empty:
        parts.append(df_to_html(fst_df, "{:.4f}"))
    for p1, p2 in pop_pairs:
        img = img_fst_man.get((p1, p2), "")
        if img:
            parts.append(f"<h3>FST Manhattan — {p1.title()} vs {p2.title()}</h3>" + img)

    # ── Relatedness
    parts.append("<h2>6. Relatedness</h2>")
    parts.append(img_kinship or "<p><em>Relatedness not yet available.</em></p>")
    # Show close relatives if any
    rel_summary = results_dir / "relatedness" / "relatedness_summary.txt"
    if rel_summary.exists():
        df_rel = pd.read_csv(rel_summary, sep="\t")
        if len(df_rel) > 0:
            parts.append("<h3>Close Relatives / Clones</h3>")
            parts.append(df_to_html(df_rel))
        else:
            parts.append('<div class="ok">✓ No close relatives or clones detected</div>')

    # ── LD
    parts.append("<h2>7. Linkage Disequilibrium</h2>")
    parts.append(img_ld or "<p><em>LD decay not yet available.</em></p>")

    # ── Methods
    parts.append("<h2>8. Methods</h2>")
    parts.append(f"<pre>{generate_methods(pipeline_config, metadata)}</pre>")

    parts.append('<div class="footer">Generated by the coral-angsd-pipeline · '
                 'Florida Atlantic University ECOS</div>')
    parts.append("</body></html>")

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        f.write("\n".join(parts))

    print(f"Report: {out_path}")
    print(f"Figures: {figures_dir}/ ({len(list(figures_dir.glob('*.png')))} PNG files)")


if __name__ == "__main__":
    main()
