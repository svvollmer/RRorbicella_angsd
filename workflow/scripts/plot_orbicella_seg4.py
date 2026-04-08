#!/usr/bin/env python3
"""
plot_orbicella_seg4.py — Diversity and FST plots for RRorbicella_angsd Segment 4

Generates:
  1. Per-chromosome diversity (pi, Watterson theta, Tajima's D) — all 3 species
  2. Sliding-window diversity (50kb) — pi and Tajima's D across genome
  3. Sliding-window FST — stacked per-chromosome Manhattan for each comparison
  4. FST summary table / heatmap

Usage:
    python plot_orbicella_seg4.py [--workdir /work/vollmer/orbicella_genomics]
                                  [--outdir docs/figures]
"""

import argparse
import glob
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── Chromosome order and labels (O. franksi jaOrbFran1.1, 15 chromosomes) ────
CHROM_ORDER = [
    "OZ116698.1","OZ116699.1","OZ116700.1","OZ116701.1","OZ116702.1",
    "OZ116703.1","OZ116704.1","OZ116705.1","OZ116706.1","OZ116707.1",
    "OZ116708.1","OZ116709.1","OZ116710.1","OZ116711.1","OZ116712.1",
]
CHROM_LABELS = {c: str(i+1) for i, c in enumerate(CHROM_ORDER)}

GROUPS = ["Oannularis_all", "Ofaveolata_all", "Ofranksi_all"]
GROUP_LABELS = {
    "Oannularis_all": "O. annularis (n=59)",
    "Ofaveolata_all": "O. faveolata (n=50)",
    "Ofranksi_all":   "O. franksi (n=24)",
}
GROUP_COLORS = {
    "Oannularis_all": "#2196F3",   # blue
    "Ofaveolata_all": "#FF9800",   # orange
    "Ofranksi_all":   "#4CAF50",   # green
}

FST_PAIRS = [
    ("Oannularis_all", "Ofaveolata_all"),
    ("Oannularis_all", "Ofranksi_all"),
    ("Ofaveolata_all", "Ofranksi_all"),
]
FST_PAIR_LABELS = {
    ("Oannularis_all","Ofaveolata_all"): "O. annularis vs O. faveolata",
    ("Oannularis_all","Ofranksi_all"):   "O. annularis vs O. franksi",
    ("Ofaveolata_all","Ofranksi_all"):   "O. faveolata vs O. franksi",
}

HIGH_Q   = 0.975
LOW_Q    = 0.025
MIN_BLOCK = 3
COL_HIGH = "#D62728"
COL_LOW  = "#2CA02C"
COL_NEU  = "#AEC7E8"


# ── Helpers ───────────────────────────────────────────────────────────────────

def load_pestpg(path):
    df = pd.read_csv(path, sep="\t", comment="#",
                     names=["win","Chr","WinCenter","tW","tP","tF","tH","tL",
                            "Tajima","fuf","fud","fayh","zeng","nSites"])
    df = df[df["Chr"].isin(CHROM_ORDER)].copy()
    df["Chr"] = pd.Categorical(df["Chr"], categories=CHROM_ORDER, ordered=True)
    return df.sort_values(["Chr","WinCenter"]).reset_index(drop=True)


def per_site(df):
    df = df.copy()
    df["tW_site"] = df["tW"] / df["nSites"]
    df["tP_site"] = df["tP"] / df["nSites"]
    return df


def load_fst_windows(path):
    df = pd.read_csv(path, sep="\t", header=0,
                     names=["region","chr","midPos","Nsites","Fst"])
    df = df[df["chr"].isin(CHROM_ORDER)].copy()
    df["chr"] = pd.Categorical(df["chr"], categories=CHROM_ORDER, ordered=True)
    return df.sort_values(["chr","midPos"]).reset_index(drop=True)


def call_fst_outliers(df):
    df = df.copy()
    high_thresh = df["Fst"].quantile(HIGH_Q)
    low_thresh  = df["Fst"].quantile(LOW_Q)
    df["outlier"] = "neutral"
    df.loc[df["Fst"] > high_thresh, "outlier"] = "high"
    df.loc[df["Fst"] < low_thresh,  "outlier"] = "low"
    df["block"] = False
    for chrom in CHROM_ORDER:
        for otype in ["high","low"]:
            idx = df.index[(df["chr"] == chrom) & (df["outlier"] == otype)].tolist()
            run = []
            for i, ix in enumerate(idx):
                if i == 0 or ix != idx[i-1] + 1:
                    if len(run) >= MIN_BLOCK:
                        df.loc[run, "block"] = True
                    run = [ix]
                else:
                    run.append(ix)
            if len(run) >= MIN_BLOCK:
                df.loc[run, "block"] = True
    return df, high_thresh, low_thresh


def genome_offsets(df, chrom_col="Chr", pos_col="WinCenter"):
    """Build cumulative genome offsets for Manhattan-style plots."""
    offsets = {}
    offset = 0
    for chrom in CHROM_ORDER:
        sub = df[df[chrom_col] == chrom]
        if len(sub):
            offsets[chrom] = offset
            offset += int(sub[pos_col].max()) + 2_000_000
    return offsets


# ── 1. Per-chromosome diversity ───────────────────────────────────────────────

def plot_chrom_diversity(workdir, outdir):
    chrom_data = {}
    for grp in GROUPS:
        f = os.path.join(workdir, f"results/diversity/{grp}.thetas.idx.pestPG")
        if not os.path.exists(f):
            print(f"  MISSING: {f}")
            continue
        df = load_pestpg(f)
        df = per_site(df)
        chrom_data[grp] = df.set_index("Chr")

    if not chrom_data:
        print("  No diversity data found.")
        return

    n_chroms = len(CHROM_ORDER)
    x = np.arange(n_chroms)
    width = 0.26

    fig, axes = plt.subplots(3, 1, figsize=(15, 12), sharex=True)
    for ax, stat, ylabel in zip(axes,
            ["tP_site", "tW_site", "Tajima"],
            ["Nucleotide diversity (π per site)",
             "Watterson's θ per site",
             "Tajima's D"]):
        for i, grp in enumerate(GROUPS):
            if grp not in chrom_data:
                continue
            vals = [chrom_data[grp].loc[c, stat] if c in chrom_data[grp].index else np.nan
                    for c in CHROM_ORDER]
            offset = (i - 1) * width
            ax.bar(x + offset, vals, width=width,
                   color=GROUP_COLORS[grp], label=GROUP_LABELS[grp],
                   edgecolor="none", alpha=0.85)
        if stat == "Tajima":
            ax.axhline(0, color="grey", lw=0.8, ls="--")
        ax.set_ylabel(ylabel, fontsize=10)
        for j in range(0, n_chroms, 2):
            ax.axvspan(j - 0.5, j + 0.5, color="#f5f5f5", zorder=0)

    axes[0].legend(fontsize=9, frameon=False)
    axes[2].set_xticks(x)
    axes[2].set_xticklabels([CHROM_LABELS[c] for c in CHROM_ORDER], fontsize=9)
    axes[2].set_xlabel("Chromosome", fontsize=10)
    fig.suptitle(
        "Per-chromosome diversity — Orbicella spp. (133 samples, genetic lineage assignment)\n"
        "Reference: O. franksi jaOrbFran1.1",
        fontsize=12, fontweight="bold")
    plt.tight_layout()
    out = os.path.join(outdir, "chrom_diversity.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  saved: {out}")


# ── 2. Sliding-window diversity ───────────────────────────────────────────────

def plot_windowed_diversity(workdir, outdir):
    # Load windowed pestPG files
    dfs = {}
    for grp in GROUPS:
        f = os.path.join(workdir, f"results/diversity/{grp}.windowed.pestPG")
        if not os.path.exists(f):
            print(f"  MISSING windowed: {f} -- skipping windowed diversity plots")
            return
        df = load_pestpg(f)
        df = per_site(df)
        df = df[df["nSites"] >= 5000]
        dfs[grp] = df

    # Build offsets from first group
    offsets = genome_offsets(list(dfs.values())[0])

    for stat, ylabel, fname in [
            ("tP_site", "Nucleotide diversity π (50 kb windows)", "windowed_pi.png"),
            ("Tajima",  "Tajima's D (50 kb windows)",             "windowed_tajima.png")]:
        fig, axes = plt.subplots(len(GROUPS), 1, figsize=(18, 9), sharey=True)
        for ax, grp in zip(axes, GROUPS):
            if grp not in dfs:
                continue
            df = dfs[grp]
            color = GROUP_COLORS[grp]
            for j, chrom in enumerate(CHROM_ORDER):
                sub = df[df["Chr"] == chrom]
                if len(sub) == 0:
                    continue
                x = sub["WinCenter"] + offsets.get(chrom, 0)
                alpha = 0.6 if j % 2 == 0 else 0.35
                ax.scatter(x, sub[stat], s=1.2, c=color, alpha=alpha, edgecolors="none")
            # chromosome dividers + labels
            for chrom in CHROM_ORDER:
                if chrom in offsets:
                    ax.axvline(offsets[chrom], color="lightgrey", lw=0.5)
            if stat == "Tajima":
                ax.axhline(0, color="grey", lw=0.8, ls="--")
            ax.set_ylabel(GROUP_LABELS[grp].split("(")[0].strip(), fontsize=9)
            ax.tick_params(labelsize=8)

        # chromosome labels on bottom axis
        df0 = list(dfs.values())[0]
        for chrom in CHROM_ORDER:
            if chrom in offsets:
                sub = df0[df0["Chr"] == chrom]
                if len(sub):
                    mid = sub["WinCenter"].mean() + offsets[chrom]
                    axes[-1].text(mid, axes[-1].get_ylim()[0] * 1.02,
                                  CHROM_LABELS[chrom], ha="center", fontsize=7)
        axes[-1].set_xlabel("Genomic position (chromosome)", fontsize=9)
        axes[-1].set_xticks([])
        fig.suptitle(
            f"{ylabel} — Orbicella spp. (133 samples)\n"
            "Reference: O. franksi jaOrbFran1.1",
            fontsize=11, fontweight="bold")
        plt.tight_layout()
        out = os.path.join(outdir, fname)
        fig.savefig(out, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"  saved: {out}")


# ── 3. FST sliding windows ────────────────────────────────────────────────────

def plot_fst_windows(workdir, outdir):
    for g1, g2 in FST_PAIRS:
        fname = f"{g1}_vs_{g2}.fst.windows"
        fpath = os.path.join(workdir, f"results/fst/{fname}")
        if not os.path.exists(fpath):
            print(f"  MISSING: {fpath} -- skipping")
            continue

        label = FST_PAIR_LABELS[(g1, g2)]
        df = load_fst_windows(fpath)
        df, high_thresh, low_thresh = call_fst_outliers(df)

        n_high = (df["outlier"] == "high").sum()
        n_low  = (df["outlier"] == "low").sum()
        print(f"  {label}: {len(df)} windows | {n_high} high-FST | {n_low} low-FST "
              f"| thresholds: >{high_thresh:.3f} / <{low_thresh:.3f}")

        # Stacked per-chromosome FST fill plot
        fig, axes = plt.subplots(len(CHROM_ORDER), 1,
                                 figsize=(16, 1.4 * len(CHROM_ORDER)),
                                 sharex=False)
        fig.subplots_adjust(hspace=0.08, left=0.07, right=0.97, top=0.95, bottom=0.03)

        for i, chrom in enumerate(CHROM_ORDER):
            ax  = axes[i]
            sub = df[df["chr"] == chrom].copy()
            if len(sub) == 0:
                ax.set_visible(False)
                continue
            pos_mb = sub["midPos"] / 1e6
            fst    = sub["Fst"]
            mask_neu  = sub["outlier"] == "neutral"
            mask_high = sub["outlier"] == "high"
            mask_low  = sub["outlier"] == "low"

            ax.fill_between(pos_mb, fst, where=mask_neu,  color=COL_NEU,  alpha=0.6,  lw=0)
            ax.fill_between(pos_mb, fst, where=mask_high, color=COL_HIGH, alpha=0.85, lw=0)
            ax.fill_between(pos_mb, fst, where=mask_low,  color=COL_LOW,  alpha=0.85, lw=0)
            ax.plot(pos_mb, fst, color="white", lw=0.2, alpha=0.4)
            ax.axhline(high_thresh,     color=COL_HIGH, lw=0.7, ls="--", alpha=0.6)
            ax.axhline(low_thresh,      color=COL_LOW,  lw=0.7, ls="--", alpha=0.6)
            ax.axhline(df["Fst"].median(), color="0.5", lw=0.5, ls=":",  alpha=0.5)
            ax.set_xlim(0, pos_mb.max() * 1.01)
            ax.set_ylim(-0.02, min(1.02, df["Fst"].quantile(0.999) * 1.5))
            ax.set_ylabel(f"Chr {CHROM_LABELS[chrom]}", fontsize=8, rotation=0,
                          labelpad=28, va="center")
            ax.tick_params(labelsize=7)
            if i < len(CHROM_ORDER) - 1:
                ax.set_xticklabels([])

        axes[-1].set_xlabel("Position (Mb)", fontsize=9)
        fig.text(0.01, 0.5, "FST", va="center", rotation="vertical", fontsize=10)

        legend_handles = [
            mpatches.Patch(color=COL_HIGH, label=f"High FST (top 2.5%, >{high_thresh:.3f})"),
            mpatches.Patch(color=COL_LOW,  label=f"Low FST (bottom 2.5%, <{low_thresh:.3f})"),
            mpatches.Patch(color=COL_NEU,  label="Neutral"),
        ]
        axes[0].legend(handles=legend_handles, fontsize=8, loc="upper right",
                       framealpha=0.8, ncol=3)
        axes[0].set_title(
            f"Sliding-window FST: {label}  (50 kb windows, 10 kb step)",
            fontsize=11, pad=6)

        pair_tag = f"{g1.replace('_all','')}_vs_{g2.replace('_all','')}"
        out = os.path.join(outdir, f"fst_windows_{pair_tag}.png")
        fig.savefig(out, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"  saved: {out}")


# ── 4. FST summary heatmap ────────────────────────────────────────────────────

def plot_fst_summary(workdir, outdir):
    sp_labels = {
        "Oannularis_all": "O. annularis",
        "Ofaveolata_all": "O. faveolata",
        "Ofranksi_all":   "O. franksi",
    }
    fst_vals = {}
    for g1, g2 in FST_PAIRS:
        f = os.path.join(workdir, f"results/fst/{g1}_vs_{g2}.fst.global")
        if not os.path.exists(f):
            continue
        vals = open(f).read().split()
        if len(vals) >= 2:
            fst_vals[(g1, g2)] = float(vals[0]) / float(vals[1])

    if not fst_vals:
        print("  No FST global files found.")
        return

    groups = GROUPS
    n = len(groups)
    gidx = {g: i for i, g in enumerate(groups)}
    mat = np.full((n, n), np.nan)
    for (g1, g2), fst in fst_vals.items():
        i, j = gidx[g1], gidx[g2]
        mat[i, j] = mat[j, i] = fst

    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.imshow(mat, cmap="YlOrRd", vmin=0, vmax=max(fst_vals.values()) * 1.1, aspect="auto")
    tick_labels = [sp_labels[g] for g in groups]
    ax.set_xticks(range(n)); ax.set_xticklabels(tick_labels, rotation=30, ha="right", fontsize=9, style="italic")
    ax.set_yticks(range(n)); ax.set_yticklabels(tick_labels, fontsize=9, style="italic")
    plt.colorbar(im, ax=ax, shrink=0.8, label="Weighted FST")
    for i in range(n):
        for j in range(n):
            if not np.isnan(mat[i, j]):
                ax.text(j, i, f"{mat[i,j]:.3f}", ha="center", va="center", fontsize=10, fontweight="bold")
            elif i != j:
                ax.text(j, i, "pending", ha="center", va="center", fontsize=8, color="grey")
    ax.set_title("Pairwise FST — Orbicella spp.\n(genetic lineage assignments, n=133)",
                 fontsize=11, fontweight="bold")
    plt.tight_layout()
    out = os.path.join(outdir, "fst_summary.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  saved: {out}")


# ── 5. Diversity summary bar chart ───────────────────────────────────────────

def plot_diversity_summary(workdir, outdir):
    rows = []
    for grp in GROUPS:
        f = os.path.join(workdir, f"results/diversity/{grp}.thetas.idx.pestPG")
        if not os.path.exists(f):
            continue
        df = load_pestpg(f)
        ns = df["nSites"].sum()
        rows.append({
            "group": grp,
            "label": GROUP_LABELS[grp],
            "pi":    df["tP"].sum() / ns,
            "tW":    df["tW"].sum() / ns,
            "tajD":  (df["Tajima"] * df["nSites"]).sum() / ns,
        })
    if not rows:
        return
    data = pd.DataFrame(rows)

    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    for ax, col, xlabel in zip(axes,
            ["pi", "tW", "tajD"],
            ["Nucleotide diversity (π)", "Watterson's θ", "Tajima's D"]):
        colors = [GROUP_COLORS[g] for g in data["group"]]
        bars = ax.barh(data["label"], data[col], color=colors, edgecolor="none", alpha=0.85)
        ax.set_xlabel(xlabel, fontsize=10)
        if col == "tajD":
            ax.axvline(0, color="grey", lw=0.8, ls="--")
        for bar, val in zip(bars, data[col]):
            ax.text(bar.get_width() + abs(data[col].max()) * 0.02, bar.get_y() + bar.get_height()/2,
                    f"{val:.4f}", va="center", fontsize=9)
        ax.invert_yaxis()
        ax.tick_params(labelsize=9)
    fig.suptitle("Genome-wide diversity — Orbicella spp. (per site)",
                 fontsize=12, fontweight="bold")
    plt.tight_layout()
    out = os.path.join(outdir, "diversity_summary.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  saved: {out}")


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workdir", default="/work/vollmer/orbicella_genomics")
    parser.add_argument("--outdir",  default="/projects/vollmer/RRorbicella_angsd/docs/figures")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print("1. Per-chromosome diversity...")
    plot_chrom_diversity(args.workdir, args.outdir)

    print("2. Sliding-window diversity (pi, Tajima's D)...")
    plot_windowed_diversity(args.workdir, args.outdir)

    print("3. FST sliding-window plots...")
    plot_fst_windows(args.workdir, args.outdir)

    print("4. FST summary heatmap...")
    plot_fst_summary(args.workdir, args.outdir)

    print("5. Diversity summary bar chart...")
    plot_diversity_summary(args.workdir, args.outdir)

    print("\nDone.")


if __name__ == "__main__":
    main()
