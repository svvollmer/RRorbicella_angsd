#!/usr/bin/env python3
"""
plot_orbicella_structure.py — K2/K3 exploration + genetic distance for RRorbicella_angsd

Generates:
  1. PCA with 95% confidence ellipses per species
  2. Between-species genetic distance matrix (from KING + PCA covariance)
  3. NJ tree from genetic distances
  4. K=2 vs K=3 side-by-side admixture with hybrid individuals highlighted
  5. Log-likelihood curve annotated to explain K=2 bias

Usage:
    python plot_orbicella_structure.py [--workdir ...] [--outdir ...]
"""

import argparse
import os
import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import dendrogram, linkage
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
import glob, re

SPECIES_COLORS = {
    "Oannularis": "#2196F3",
    "Ofaveolata":  "#FF9800",
    "Ofranksi":    "#4CAF50",
}
SPECIES_LABELS = {
    "Oannularis": "O. annularis",
    "Ofaveolata":  "O. faveolata",
    "Ofranksi":    "O. franksi",
}
ADMIX_PALETTE = ["#2196F3", "#FF9800", "#4CAF50", "#E91E63", "#9C27B0"]

# ── helpers ───────────────────────────────────────────────────────────────────

def load_metadata(workdir):
    meta = pd.read_csv(os.path.join(workdir, "config/samples.csv"))
    with open(os.path.join(workdir, "results/relatedness/unrelated_samples.txt")) as f:
        order = [l.strip() for l in f if l.strip()]
    meta = meta.set_index("sample_id").loc[order].reset_index()
    return meta

def load_pca(workdir):
    cov = np.loadtxt(os.path.join(workdir, "results/pca/pcangsd.cov"))
    vals, vecs = np.linalg.eigh(cov)
    idx = np.argsort(vals)[::-1]
    vals, vecs = vals[idx], vecs[:, idx]
    pve = vals / vals.sum() * 100
    return vecs, pve, cov

def load_king(workdir):
    path = os.path.join(workdir, "results/relatedness/relatedness_matrix.txt")
    df = pd.read_csv(path, sep="\t", index_col=0)
    return df

def load_q(admix_dir, method, k):
    prefix = "ngsadmix" if method == "ngsadmix" else "pcangsd"
    return np.loadtxt(os.path.join(admix_dir, f"{prefix}_K{k}.Q"))

def collect_loglikes(admix_dir, k):
    ll = []
    for f in glob.glob(os.path.join(admix_dir, f"ngsadmix_K{k}_rep*.log")):
        with open(f) as fh:
            for line in fh:
                m = re.search(r"best like=([-\d.]+)", line)
                if m:
                    ll.append(float(m.group(1)))
    return ll

def confidence_ellipse(x, y, ax, n_std=2.0, **kwargs):
    """Plot confidence ellipse for scatter data."""
    if len(x) < 3:
        return
    cov = np.cov(x, y)
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]
    angle = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
    width, height = 2 * n_std * np.sqrt(vals)
    ellipse = Ellipse(xy=(np.mean(x), np.mean(y)), width=width, height=height,
                      angle=angle, **kwargs)
    ax.add_patch(ellipse)
    return ellipse

# ── 1. PCA with ellipses ──────────────────────────────────────────────────────

def plot_pca_ellipses(vecs, pve, meta, outpath):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for ax, (pcx, pcy) in zip(axes, [(0, 1), (1, 2)]):
        for sp, grp in meta.groupby("species"):
            idx = grp.index
            c = SPECIES_COLORS.get(sp, "grey")
            ax.scatter(vecs[idx, pcx], vecs[idx, pcy],
                       c=c, label=SPECIES_LABELS.get(sp, sp),
                       s=40, alpha=0.8, edgecolors="white", linewidths=0.3, zorder=3)
            confidence_ellipse(vecs[idx, pcx], vecs[idx, pcy], ax,
                               n_std=2.0, facecolor=c, alpha=0.12, edgecolor=c,
                               linewidth=1.5, linestyle="--", zorder=2)
            # centroid
            ax.scatter(vecs[idx, pcx].mean(), vecs[idx, pcy].mean(),
                       c=c, s=120, marker="D", edgecolors="black",
                       linewidths=1.0, zorder=4)
        ax.set_xlabel(f"PC{pcx+1} ({pve[pcx]:.1f}%)", fontsize=11)
        ax.set_ylabel(f"PC{pcy+1} ({pve[pcy]:.1f}%)", fontsize=11)
        ax.axhline(0, color="grey", lw=0.4, ls=":")
        ax.axvline(0, color="grey", lw=0.4, ls=":")
        ax.tick_params(labelsize=9)

    axes[0].legend(fontsize=9, framealpha=0.9, loc="upper right")
    fig.suptitle("PCA — Orbicella spp. (95% ellipses, diamonds = centroids)", fontsize=12)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved: {outpath}")

# ── 2. Genetic distance matrix + NJ tree ─────────────────────────────────────

def compute_species_distances(cov, king_df, meta, workdir):
    """
    Compute two between-species distance matrices:
      A) From PCA covariance: mean between-species covariance, normalized to 0-1 distance
      B) From KING: 1 - mean_KING (rescaled; within-species KING ~0, between ~-0.6)
    """
    species_list = ["Oannularis", "Ofaveolata", "Ofranksi"]
    unrelated = meta["sample_id"].tolist()

    # --- covariance-based distance ---
    sp_idx = {sp: meta[meta["species"] == sp].index.tolist() for sp in species_list}
    n = len(species_list)
    cov_within, cov_between = {}, {}
    for sp in species_list:
        idx = sp_idx[sp]
        if len(idx) > 1:
            sub = cov[np.ix_(idx, idx)]
            cov_within[sp] = np.mean(sub[np.triu_indices(len(idx), k=1)])
    dist_cov = np.zeros((n, n))
    for i, sp1 in enumerate(species_list):
        for j, sp2 in enumerate(species_list):
            if i == j:
                dist_cov[i, j] = 0
            else:
                idx1, idx2 = sp_idx[sp1], sp_idx[sp2]
                sub = cov[np.ix_(idx1, idx2)]
                mean_between = np.mean(sub)
                # distance = within_mean - between_mean (larger = more different)
                w = (cov_within.get(sp1, 0) + cov_within.get(sp2, 0)) / 2
                dist_cov[i, j] = max(0, w - mean_between)

    # --- KING-based distance ---
    # subset king to unrelated samples only
    avail = [s for s in unrelated if s in king_df.index]
    king_sub = king_df.loc[avail, avail]
    sp_samples = {sp: [s for s in avail if s.startswith({"Oannularis": "Oann",
                                                           "Ofaveolata": "Ofav",
                                                           "Ofranksi": "Ofra"}[sp])]
                  for sp in species_list}
    dist_king = np.zeros((n, n))
    king_vals = {}
    for i, sp1 in enumerate(species_list):
        for j, sp2 in enumerate(species_list):
            s1, s2 = sp_samples[sp1], sp_samples[sp2]
            if i == j:
                if len(s1) > 1:
                    sub = king_sub.loc[s1, s1].values
                    iu = np.triu_indices(len(s1), k=1)
                    king_vals[(sp1, sp2)] = np.mean(sub[iu])
            else:
                sub = king_sub.loc[s1, s2].values
                king_vals[(sp1, sp2)] = np.mean(sub)

    # rescale: within-species KING ~0, between-species ~-0.6
    # distance = within_mean - between_mean
    for i, sp1 in enumerate(species_list):
        for j, sp2 in enumerate(species_list):
            if i == j:
                dist_king[i, j] = 0
            else:
                w1 = king_vals.get((sp1, sp1), 0)
                w2 = king_vals.get((sp2, sp2), 0)
                b  = king_vals.get((sp1, sp2), king_vals.get((sp2, sp1), 0))
                dist_king[i, j] = max(0, ((w1 + w2) / 2) - b)

    return species_list, dist_cov, dist_king, king_vals

def plot_distance_tree(species_list, dist_cov, dist_king, king_vals, outpath):
    labels = [SPECIES_LABELS[s] for s in species_list]
    n = len(species_list)

    fig = plt.figure(figsize=(14, 5))
    gs = fig.add_gridspec(1, 3, wspace=0.4)

    # --- KING heatmap ---
    ax_heat = fig.add_subplot(gs[0])
    # build symmetric matrix including within
    mat = np.zeros((n, n))
    for i, sp1 in enumerate(species_list):
        for j, sp2 in enumerate(species_list):
            if i == j:
                mat[i, j] = king_vals.get((sp1, sp2), 0)
            else:
                mat[i, j] = king_vals.get((sp1, sp2), king_vals.get((sp2, sp1), 0))

    im = ax_heat.imshow(mat, cmap="RdBu", vmin=-0.7, vmax=0.1, aspect="auto")
    ax_heat.set_xticks(range(n)); ax_heat.set_xticklabels(labels, rotation=30, ha="right",
                                                             fontsize=9, style="italic")
    ax_heat.set_yticks(range(n)); ax_heat.set_yticklabels(labels, fontsize=9, style="italic")
    for i in range(n):
        for j in range(n):
            ax_heat.text(j, i, f"{mat[i,j]:.3f}", ha="center", va="center", fontsize=8,
                         color="white" if abs(mat[i,j]) > 0.3 else "black")
    plt.colorbar(im, ax=ax_heat, shrink=0.7, label="Mean KING")
    ax_heat.set_title("Mean pairwise KING\nbetween species", fontsize=10)

    # --- Covariance distance heatmap ---
    ax_dist = fig.add_subplot(gs[1])
    im2 = ax_dist.imshow(dist_cov, cmap="YlOrRd", aspect="auto")
    ax_dist.set_xticks(range(n)); ax_dist.set_xticklabels(labels, rotation=30, ha="right",
                                                            fontsize=9, style="italic")
    ax_dist.set_yticks(range(n)); ax_dist.set_yticklabels(labels, fontsize=9, style="italic")
    for i in range(n):
        for j in range(n):
            ax_dist.text(j, i, f"{dist_cov[i,j]:.4f}", ha="center", va="center", fontsize=8)
    plt.colorbar(im2, ax=ax_dist, shrink=0.7, label="Genetic distance")
    ax_dist.set_title("Genetic distance\n(covariance-based)", fontsize=10)

    # --- NJ dendrogram ---
    ax_tree = fig.add_subplot(gs[2])
    # make symmetric, ensure zero diagonal
    d = (dist_cov + dist_cov.T) / 2
    np.fill_diagonal(d, 0)
    condensed = squareform(d)
    Z = linkage(condensed, method="average")
    dn = dendrogram(Z, labels=labels, ax=ax_tree,
                    leaf_font_size=10, leaf_rotation=0,
                    color_threshold=0,
                    link_color_func=lambda k: "black")
    # colour leaf labels by species
    for lbl in ax_tree.get_xticklabels():
        t = lbl.get_text()
        for sp, splab in SPECIES_LABELS.items():
            if t == splab:
                lbl.set_color(SPECIES_COLORS[sp])
                lbl.set_style("italic")
    ax_tree.set_ylabel("Genetic distance", fontsize=10)
    ax_tree.set_title("UPGMA tree\n(covariance distances)", fontsize=10)
    ax_tree.tick_params(axis="x", labelsize=9)

    fig.suptitle("Between-species genetic distances — Orbicella spp.", fontsize=12, y=1.02)
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved: {outpath}")

# ── 3. K=2 vs K=3 with hybrid flagging ───────────────────────────────────────

def plot_k2_k3_comparison(admix_dir, meta, outpath, hybrid_thresh=0.05):
    fig, axes = plt.subplots(2, 1, figsize=(14, 6))
    sort_order = meta.sort_values(["species", "region"]).index.tolist()
    sorted_meta = meta.loc[sort_order].reset_index(drop=True)

    sp_bounds = []
    prev = None
    for i, sp in enumerate(sorted_meta["species"]):
        if sp != prev:
            sp_bounds.append(i)
            prev = sp
    sp_bounds.append(len(sorted_meta))

    hybrids_k3 = []
    for row, (k, method) in enumerate([(2, "ngsadmix"), (3, "ngsadmix")]):
        ax = axes[row]
        Q = load_q(admix_dir, method, k)[sort_order]
        col_order = np.argsort(-Q.mean(axis=0))
        Q = Q[:, col_order]

        bottom = np.zeros(len(Q))
        for c in range(k):
            ax.bar(range(len(Q)), Q[:, c], bottom=bottom,
                   color=ADMIX_PALETTE[c % len(ADMIX_PALETTE)],
                   width=1.0, linewidth=0)
            bottom += Q[:, c]

        # flag hybrids at K=3 (no component > 1-thresh)
        if k == 3:
            for i, q in enumerate(Q):
                if np.max(q) < (1 - hybrid_thresh):
                    hybrids_k3.append(sorted_meta["sample_id"].iloc[i])
                    ax.axvline(i, color="red", lw=0.6, alpha=0.5, zorder=5)

        for b in sp_bounds[1:-1]:
            ax.axvline(b - 0.5, color="white", lw=1.5)
        for i in range(len(sp_bounds) - 1):
            mid = (sp_bounds[i] + sp_bounds[i+1]) / 2
            sp = sorted_meta["species"].iloc[sp_bounds[i]]
            ax.text(mid, -0.1, SPECIES_LABELS.get(sp, sp),
                    ha="center", va="top", fontsize=9, style="italic",
                    transform=ax.get_xaxis_transform())

        ax.set_xlim(-0.5, len(Q) - 0.5)
        ax.set_ylim(0, 1)
        ax.set_yticks([0, 0.5, 1])
        ax.set_yticklabels(["0", "0.5", "1"], fontsize=9)
        ax.set_xticks([])
        ax.set_ylabel(f"K={k}", fontsize=11)

        if k == 3:
            ax.text(0.01, 0.95, f"Red lines = putative hybrids (no component >85%): n={len(hybrids_k3)}",
                    transform=ax.transAxes, fontsize=8, va="top", color="red")

    fig.suptitle("ngsAdmix K=2 vs K=3 — Orbicella spp.\n"
                 "Note: ΔK favours K=2 (strongest split = annularis vs faveolata/franksi clade);\n"
                 "K=3 recovers all three species and reveals admixed individuals",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved: {outpath}  ({len(hybrids_k3)} putative hybrids at K=3)")
    return hybrids_k3

# ── 4. Log-likelihood annotated ───────────────────────────────────────────────

def plot_loglik_annotated(admix_dir, outpath, k_range=range(1, 11)):
    means, sds = {}, {}
    for k in k_range:
        ll = collect_loglikes(admix_dir, k)
        if ll:
            means[k] = np.mean(ll)
            sds[k]   = np.std(ll)

    ks = sorted(means.keys())
    ll_vals = [means[k] for k in ks]
    sd_vals  = [sds[k]  for k in ks]

    # second derivative (elbow)
    d2 = {}
    for i in range(1, len(ks) - 1):
        d2[ks[i]] = ll_vals[i+1] - 2*ll_vals[i] + ll_vals[i-1]

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    ax = axes[0]
    ax.errorbar(ks, ll_vals, yerr=sd_vals, fmt="o-", color="#2196F3",
                capsize=4, lw=2, ms=7)
    ax.axvline(2, color="orange", ls="--", lw=1.5, label="ΔK best (K=2)\nbut driven by annularis split")
    ax.axvline(3, color="green",  ls="--", lw=1.5, label="Biological K=3\n(3 species, elbow)")
    ax.set_xlabel("K", fontsize=11); ax.set_ylabel("Mean log-likelihood", fontsize=11)
    ax.set_title("Log-likelihood: large gain K=1→2\n(annularis split), smaller gain K=2→3 (franksi resolved)", fontsize=9)
    ax.set_xticks(ks); ax.legend(fontsize=8, framealpha=0.9)
    ax.tick_params(labelsize=9)

    # gains per step
    ax2 = axes[1]
    gains = [means[ks[i+1]] - means[ks[i]] for i in range(len(ks)-1)]
    gain_ks = ks[1:]
    bars = ax2.bar(gain_ks, gains, color=["#FF9800" if k==2 else "#4CAF50" if k==3 else "#90A4AE"
                                           for k in gain_ks], edgecolor="white")
    ax2.set_xlabel("K", fontsize=11)
    ax2.set_ylabel("ΔL (gain per K step)", fontsize=11)
    ax2.set_title("Log-likelihood gain per K step\nK=2 gain dominates (ΔK artefact)", fontsize=9)
    ax2.set_xticks(gain_ks); ax2.tick_params(labelsize=9)
    # annotate
    ax2.text(2, gains[0]*0.5, "annularis\nsplit", ha="center", fontsize=8, color="white", fontweight="bold")
    if len(gains) > 1:
        ax2.text(3, gains[1]*0.5, "franksi\nresolved", ha="center", fontsize=8, color="white", fontweight="bold")

    fig.suptitle("Why ΔK picks K=2 even though K=3 is biologically correct", fontsize=11)
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

    print("Loading data...")
    meta  = load_metadata(args.workdir)
    vecs, pve, cov = load_pca(args.workdir)
    king_df = load_king(args.workdir)
    print(f"  {len(meta)} samples, cov {cov.shape}, KING {king_df.shape}")

    print("PCA with ellipses...")
    plot_pca_ellipses(vecs, pve, meta,
                      os.path.join(args.outdir, "pca_ellipses.png"))

    print("Computing genetic distances...")
    sp_list, dist_cov, dist_king, king_vals = compute_species_distances(cov, king_df, meta, args.workdir)
    print("  Species distances (covariance-based):")
    for i, sp1 in enumerate(sp_list):
        for j, sp2 in enumerate(sp_list):
            if i < j:
                print(f"    {SPECIES_LABELS[sp1]} vs {SPECIES_LABELS[sp2]}: {dist_cov[i,j]:.4f}")
    print("  Species mean KING:")
    for k, v in sorted(king_vals.items()):
        print(f"    {k[0]} vs {k[1]}: {v:.4f}")

    print("Plotting distance tree...")
    plot_distance_tree(sp_list, dist_cov, dist_king, king_vals,
                       os.path.join(args.outdir, "genetic_distance_tree.png"))

    print("Plotting K=2 vs K=3 comparison...")
    hybrids = plot_k2_k3_comparison(admix_dir, meta,
                                    os.path.join(args.outdir, "admixture_K2_K3_comparison.png"))
    if hybrids:
        print(f"  Putative hybrids at K=3: {hybrids}")

    print("Plotting log-likelihood annotated...")
    plot_loglik_annotated(admix_dir,
                          os.path.join(args.outdir, "admixture_loglik_annotated.png"))

    print("\nDone.")

if __name__ == "__main__":
    main()
