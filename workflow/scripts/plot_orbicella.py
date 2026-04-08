#!/usr/bin/env python3
"""
plot_orbicella.py — Population structure plots for RRorbicella_angsd

Generates:
  1. PCA (PC1 vs PC2, colored by species)
  2. ngsAdmix bar plots K=2,3,4
  3. Delta-K plot (ngsAdmix)
  4. LD decay curve

Usage:
    python plot_orbicella.py [--workdir /work/vollmer/orbicella_genomics]
                             [--outdir docs/figures]
"""

import argparse
import glob
import os
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec

# ── colours ──────────────────────────────────────────────────────────────────
SPECIES_COLORS = {
    "Oannularis": "#2196F3",   # blue
    "Ofaveolata":  "#FF9800",  # orange
    "Ofranksi":    "#4CAF50",  # green
}
SPECIES_LABELS = {
    "Oannularis": "O. annularis",
    "Ofaveolata":  "O. faveolata",
    "Ofranksi":    "O. franksi",
}
ADMIX_PALETTE = ["#2196F3", "#FF9800", "#4CAF50", "#E91E63",
                 "#9C27B0", "#00BCD4", "#FF5722", "#607D8B",
                 "#8BC34A", "#FFC107"]

def load_metadata(workdir):
    samples_csv = os.path.join(workdir, "config/samples.csv")
    unrelated   = os.path.join(workdir, "results/relatedness/unrelated_samples.txt")
    meta = pd.read_csv(samples_csv)
    with open(unrelated) as f:
        order = [l.strip() for l in f if l.strip()]
    meta = meta.set_index("sample_id").loc[order].reset_index()
    return meta, order

def load_pca(workdir):
    cov = np.loadtxt(os.path.join(workdir, "results/pca/pcangsd.cov"))
    vals, vecs = np.linalg.eigh(cov)
    idx  = np.argsort(vals)[::-1]
    vals = vals[idx]; vecs = vecs[:, idx]
    pve  = vals / vals.sum() * 100
    return vecs, pve

def load_q(path):
    return np.loadtxt(path)

def collect_loglikes(admix_dir, k):
    loglikes = []
    for f in glob.glob(os.path.join(admix_dir, f"ngsadmix_K{k}_rep*.log")):
        with open(f) as fh:
            for line in fh:
                m = re.search(r"best like=([-\d.]+)", line)
                if m:
                    loglikes.append(float(m.group(1)))
    return loglikes

def delta_k(admix_dir, k_range):
    means, sds, loglikes_all = {}, {}, {}
    for k in k_range:
        ll = collect_loglikes(admix_dir, k)
        if ll:
            means[k] = np.mean(ll)
            sds[k]   = np.std(ll)
            loglikes_all[k] = ll
    ks = sorted(means.keys())
    dk = {}
    for i, k in enumerate(ks):
        if i == 0 or i == len(ks) - 1:
            continue
        if sds[k] > 0:
            dk[k] = abs(means[ks[i+1]] - 2*means[k] + means[ks[i-1]]) / sds[k]
    return means, sds, dk

# ── 1. PCA ────────────────────────────────────────────────────────────────────
def plot_pca(vecs, pve, meta, outpath):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, (pc_x, pc_y) in zip(axes, [(0, 1), (0, 2)]):
        for sp, grp in meta.groupby("species"):
            idx = grp.index
            ax.scatter(vecs[idx, pc_x], vecs[idx, pc_y],
                       c=SPECIES_COLORS.get(sp, "grey"),
                       label=SPECIES_LABELS.get(sp, sp),
                       s=40, alpha=0.8, edgecolors="white", linewidths=0.3)
        ax.set_xlabel(f"PC{pc_x+1} ({pve[pc_x]:.1f}%)", fontsize=11)
        ax.set_ylabel(f"PC{pc_y+1} ({pve[pc_y]:.1f}%)", fontsize=11)
        ax.axhline(0, color="grey", lw=0.5, ls="--")
        ax.axvline(0, color="grey", lw=0.5, ls="--")
        ax.tick_params(labelsize=9)
    axes[0].legend(fontsize=9, framealpha=0.8)
    fig.suptitle("PCA — Orbicella spp. (PCAngsd, 134 unrelated samples)", fontsize=12)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved: {outpath}")

# ── 2. Admixture bar plots ────────────────────────────────────────────────────
def plot_admix_bars(admix_dir, meta, k_list, outpath, method="ngsadmix"):
    n = len(k_list)
    fig, axes = plt.subplots(n, 1, figsize=(14, 2.2 * n))
    if n == 1:
        axes = [axes]

    # sort samples: species then region
    sort_order = meta.sort_values(["species", "region"]).index.tolist()
    sorted_meta = meta.loc[sort_order].reset_index(drop=True)

    # species dividers
    sp_bounds = []
    prev = None
    for i, sp in enumerate(sorted_meta["species"]):
        if sp != prev:
            sp_bounds.append(i)
            prev = sp
    sp_bounds.append(len(sorted_meta))

    for ax, k in zip(axes, k_list):
        prefix = "ngsadmix" if method == "ngsadmix" else "pcangsd"
        qfile  = os.path.join(admix_dir, f"{prefix}_K{k}.Q")
        if not os.path.exists(qfile):
            continue
        Q = load_q(qfile)[sort_order]
        # sort columns by dominant component
        col_order = np.argsort(-Q.mean(axis=0))
        Q = Q[:, col_order]

        bottom = np.zeros(len(Q))
        for c in range(k):
            ax.bar(range(len(Q)), Q[:, c], bottom=bottom,
                   color=ADMIX_PALETTE[c % len(ADMIX_PALETTE)],
                   width=1.0, linewidth=0)
            bottom += Q[:, c]

        # species dividers + labels
        for b in sp_bounds[1:-1]:
            ax.axvline(b - 0.5, color="white", lw=1.5)
        for i in range(len(sp_bounds) - 1):
            mid = (sp_bounds[i] + sp_bounds[i+1]) / 2
            sp  = sorted_meta["species"].iloc[sp_bounds[i]]
            ax.text(mid, -0.08, SPECIES_LABELS.get(sp, sp),
                    ha="center", va="top", fontsize=8, style="italic",
                    transform=ax.get_xaxis_transform())

        ax.set_xlim(-0.5, len(Q) - 0.5)
        ax.set_ylim(0, 1)
        ax.set_yticks([0, 0.5, 1])
        ax.set_yticklabels(["0", "0.5", "1"], fontsize=8)
        ax.set_xticks([])
        ax.set_ylabel(f"K={k}", fontsize=10)

    title = "ngsAdmix" if method == "ngsadmix" else "PCAngsd"
    fig.suptitle(f"{title} admixture — Orbicella spp. (134 samples)", fontsize=12, y=1.01)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved: {outpath}")

# ── 3. Delta-K ────────────────────────────────────────────────────────────────
def plot_delta_k(admix_dir, outpath, k_range=range(1, 11)):
    means, sds, dk = delta_k(admix_dir, k_range)
    ks = sorted(means.keys())
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))

    # mean log-likelihood
    axes[0].errorbar(ks, [means[k] for k in ks], yerr=[sds[k] for k in ks],
                     fmt="o-", color="#2196F3", capsize=4, lw=1.5)
    axes[0].set_xlabel("K", fontsize=11)
    axes[0].set_ylabel("Mean log-likelihood", fontsize=11)
    axes[0].set_title("Log-likelihood per K", fontsize=11)
    axes[0].set_xticks(ks)
    axes[0].tick_params(labelsize=9)

    # delta-K
    dk_ks = sorted(dk.keys())
    axes[1].bar(dk_ks, [dk[k] for k in dk_ks], color="#FF9800", edgecolor="white")
    best_k = max(dk, key=dk.get)
    axes[1].axvline(best_k, color="red", ls="--", lw=1.5, label=f"Best K={best_k}")
    axes[1].set_xlabel("K", fontsize=11)
    axes[1].set_ylabel("ΔK", fontsize=11)
    axes[1].set_title("Evanno ΔK", fontsize=11)
    axes[1].set_xticks(dk_ks)
    axes[1].legend(fontsize=9)
    axes[1].tick_params(labelsize=9)

    fig.suptitle("ngsAdmix model selection — Orbicella spp.", fontsize=12)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved: {outpath}  (best K={best_k}, ΔK={dk[best_k]:.1f})")
    return best_k

# ── 4. LD decay ───────────────────────────────────────────────────────────────
def plot_ld_decay(workdir, outpath):
    ld = pd.read_csv(os.path.join(workdir, "results/ld/ld_decay.csv"))
    # parse bin midpoints
    mids = []
    for b in ld["bin"]:
        lo, hi = b.split("-")
        mids.append((int(lo) + int(hi)) / 2)
    ld["midpoint"] = mids

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(ld["midpoint"] / 1000, ld["mean_r2"], "o-", color="#4CAF50",
            lw=1.8, ms=5, alpha=0.9)
    ax.axhline(0.1, color="grey", ls="--", lw=1, label="r²=0.1")
    ax.set_xlabel("Distance (kb)", fontsize=11)
    ax.set_ylabel("Mean r²", fontsize=11)
    ax.set_title("LD decay — Orbicella spp. (134 samples, all chromosomes)", fontsize=11)
    ax.legend(fontsize=9)
    ax.tick_params(labelsize=9)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved: {outpath}")

# ── main ──────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workdir", default="/work/vollmer/orbicella_genomics")
    parser.add_argument("--outdir",  default="/projects/vollmer/RRorbicella_angsd/docs/figures")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    admix_dir = os.path.join(args.workdir, "results/admixture")

    print("Loading metadata...")
    meta, order = load_metadata(args.workdir)
    print(f"  {len(meta)} unrelated samples: {meta['species'].value_counts().to_dict()}")

    print("Loading PCA...")
    vecs, pve = load_pca(args.workdir)
    print(f"  PC1={pve[0]:.1f}%  PC2={pve[1]:.1f}%  PC3={pve[2]:.1f}%")

    print("Plotting PCA...")
    plot_pca(vecs, pve, meta, os.path.join(args.outdir, "pca_species.png"))

    print("Plotting delta-K...")
    best_k = plot_delta_k(admix_dir, os.path.join(args.outdir, "admixture_deltaK.png"))

    print("Plotting ngsAdmix bars K=2,3,4...")
    plot_admix_bars(admix_dir, meta, [2, 3, 4],
                    os.path.join(args.outdir, "admixture_ngsadmix_K234.png"),
                    method="ngsadmix")

    print("Plotting PCAngsd bars K=2,3,4...")
    plot_admix_bars(admix_dir, meta, [2, 3, 4],
                    os.path.join(args.outdir, "admixture_pcangsd_K234.png"),
                    method="pcangsd")

    print("Plotting LD decay...")
    plot_ld_decay(args.workdir, os.path.join(args.outdir, "ld_decay.png"))

    print("\nDone.")

if __name__ == "__main__":
    main()
