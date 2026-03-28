"""
plot_fst_windows.py — FST sliding-window analysis with Z-score outlier detection,
three stacked chromosome plot styles, Manhattan plot, per-chromosome panels,
and gene intersection.

Usage:
    python plot_fst_windows.py \
        --fst results/fst/lineageB_FL_vs_lineageA_FL.fst.windows \
        --gff /work/vollmer/RR_heat-tolerance/Acropora/reference/jaAcrPala1.3_genomic.gff.gz \
        --out docs/figures/fst_acer_vs_apal_FL
"""
import argparse
import gzip
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

CHROM_ORDER = [
    "NC_133882.1","NC_133883.1","NC_133884.1","NC_133885.1",
    "NC_133886.1","NC_133887.1","NC_133888.1","NC_133889.1",
    "NC_133890.1","NC_133891.1","NC_133892.1","NC_133893.1",
    "NC_133894.1","NC_133895.1",
]
CHROM_LABELS = {c: str(i+1) for i, c in enumerate(CHROM_ORDER)}

HIGH_Q = 0.975   # top 2.5%
LOW_Q  = 0.025   # bottom 2.5%
MIN_BLOCK = 3    # minimum consecutive outlier windows to call a block

COL_HIGH = "#D62728"
COL_LOW  = "#2CA02C"
COL_NEU  = "#AEC7E8"


# ── Z-score outlier calling ──────────────────────────────────────────────────

def per_chrom_zscore(df):
    """Modified Z-score computed per chromosome (relative to local background)."""
    z = np.zeros(len(df))
    for chrom in CHROM_ORDER:
        idx = df.index[df["chr"] == chrom]
        vals = df.loc[idx, "Fst"].values
        med = np.median(vals)
        mad = np.median(np.abs(vals - med))
        if mad == 0:
            mad = np.mean(np.abs(vals - med)) or 1e-9
        z[df.index.get_indexer(idx)] = 0.6745 * (vals - med) / mad
    return z


def call_outliers(df):
    """
    Classify windows using genome-wide quantile thresholds.
    Per-chromosome Z-scores are computed for display/reference.
    Also identifies outlier blocks (MIN_BLOCK+ consecutive outlier windows).
    """
    df = df.copy()
    high_thresh = df["Fst"].quantile(HIGH_Q)
    low_thresh  = df["Fst"].quantile(LOW_Q)
    df["zscore"]  = per_chrom_zscore(df)
    df["outlier"] = "neutral"
    df.loc[df["Fst"] > high_thresh, "outlier"] = "high"
    df.loc[df["Fst"] < low_thresh,  "outlier"] = "low"
    df["high_thresh"] = high_thresh
    df["low_thresh"]  = low_thresh

    # Mark outlier blocks (MIN_BLOCK+ consecutive same-type outlier windows)
    df["block"] = False
    for chrom in CHROM_ORDER:
        for otype in ["high", "low"]:
            mask = (df["chr"] == chrom) & (df["outlier"] == otype)
            idx = df.index[mask].tolist()
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


def get_outlier_spans(sub, pos_col, otype):
    """Return list of (start_mb, width_mb) spans for contiguous outlier blocks."""
    mask = (sub["outlier"] == otype) & sub["block"]
    pos = sub.loc[mask, pos_col].values
    if len(pos) == 0:
        return []
    # step size
    step = np.median(np.diff(sub[pos_col].values)) if len(sub) > 1 else 10000
    half = step / 2
    spans = []
    start = pos[0]
    prev  = pos[0]
    for p in pos[1:]:
        if p - prev > step * 1.5:
            spans.append((start - half, prev - start + step))
            start = p
        prev = p
    spans.append((start - half, prev - start + step))
    return spans


# ── GFF parsing ──────────────────────────────────────────────────────────────

def parse_gff_genes(gff_path):
    """Extract gene records from GFF3 (gzipped or plain)."""
    records = []
    opener = gzip.open if gff_path.endswith(".gz") else open
    with opener(gff_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            chrom, start, end, attrs = parts[0], int(parts[3]), int(parts[4]), parts[8]
            if chrom not in CHROM_ORDER:
                continue
            name    = re.search(r"Name=([^;]+)", attrs)
            desc    = re.search(r"description=([^;]+)", attrs)
            gene_id = re.search(r"gene=([^;]+)", attrs)
            biotype = re.search(r"gene_biotype=([^;]+)", attrs)
            records.append({
                "chr":         chrom,
                "start":       start,
                "end":         end,
                "name":        name.group(1) if name else "",
                "gene_id":     gene_id.group(1) if gene_id else "",
                "description": desc.group(1).replace("%20", " ") if desc else "",
                "biotype":     biotype.group(1) if biotype else "",
            })
    return pd.DataFrame(records)


def find_genes_in_windows(outlier_df, genes_df, win_half=25000):
    """For each outlier window return overlapping genes."""
    results = []
    for _, row in outlier_df.iterrows():
        win_start = row["midPos"] - win_half
        win_end   = row["midPos"] + win_half
        hits = genes_df[
            (genes_df["chr"] == row["chr"]) &
            (genes_df["end"]   >= win_start) &
            (genes_df["start"] <= win_end)
        ].copy()
        if len(hits):
            hits["window_chr"]    = row["chr"]
            hits["window_midPos"] = row["midPos"]
            hits["window_Fst"]    = row["Fst"]
            hits["window_zscore"] = row["zscore"]
            hits["outlier_type"]  = row["outlier"]
            results.append(hits)
    if results:
        return pd.concat(results, ignore_index=True).drop_duplicates(
            subset=["chr","start","end"])
    return pd.DataFrame()


# ── Shared axis setup helpers ─────────────────────────────────────────────────

def _stacked_fig():
    fig, axes = plt.subplots(
        len(CHROM_ORDER), 1,
        figsize=(16, 1.4 * len(CHROM_ORDER)),
        sharex=False,
    )
    fig.subplots_adjust(hspace=0.08, left=0.07, right=0.97, top=0.96, bottom=0.04)
    return fig, axes


def _label_axes(axes, ylabel):
    axes[-1].set_xlabel("Position (Mb)", fontsize=10)
    for ax in axes[:-1]:
        ax.set_xticklabels([])
    fig = axes[0].get_figure()
    fig.text(0.01, 0.5, ylabel, va="center", rotation="vertical", fontsize=10)


# ── Plot style 1: raw FST fill ────────────────────────────────────────────────

def plot_stacked_fst(df, out, high_thresh, low_thresh):
    """Fill-between coloured by outlier class; y = raw FST."""
    fig, axes = _stacked_fig()

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

        ax.fill_between(pos_mb[mask_neu],  fst[mask_neu],  color=COL_NEU,  alpha=0.6, lw=0)
        ax.fill_between(pos_mb[mask_high], fst[mask_high], color=COL_HIGH, alpha=0.85, lw=0)
        ax.fill_between(pos_mb[mask_low],  fst[mask_low],  color=COL_LOW,  alpha=0.85, lw=0)
        ax.plot(pos_mb, fst, color="white", lw=0.2, alpha=0.4)

        ax.axhline(high_thresh,        color=COL_HIGH, lw=0.7, ls="--", alpha=0.6)
        ax.axhline(low_thresh,         color=COL_LOW,  lw=0.7, ls="--", alpha=0.6)
        ax.axhline(df["Fst"].median(), color="0.5",    lw=0.5, ls=":",  alpha=0.5)

        ax.set_xlim(0, pos_mb.max() * 1.01)
        ax.set_ylim(-0.02, 1.02)
        ax.set_ylabel(f"Chr {CHROM_LABELS[chrom]}", fontsize=8, rotation=0,
                      labelpad=28, va="center")
        ax.tick_params(labelsize=7)

    _label_axes(axes, "FST (unweighted)")
    legend_handles = [
        mpatches.Patch(color=COL_HIGH, label=f"High FST (top 2.5%, >{high_thresh:.3f})"),
        mpatches.Patch(color=COL_LOW,  label=f"Low FST (bottom 2.5%, <{low_thresh:.3f})"),
        mpatches.Patch(color=COL_NEU,  label="Neutral"),
    ]
    axes[0].legend(handles=legend_handles, fontsize=8, loc="upper right",
                   framealpha=0.8, ncol=3)
    axes[0].set_title(
        "Sliding-window FST: A. cervicornis FL vs A. palmata FL  "
        "(50 kb windows, 10 kb step)  —  Style 1: raw FST",
        fontsize=11, pad=6,
    )
    fig.savefig(out + "_stacked_fst.png", dpi=180, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out}_stacked_fst.png")


# ── Plot style 2: per-chromosome Z-score fill ─────────────────────────────────

def plot_stacked_zscore(df, out):
    """Fill-between coloured by outlier class; y = per-chrom modified Z-score."""
    fig, axes = _stacked_fig()

    # Per-chrom Z-score outlier thresholds (same quantile ranks as FST)
    # We already classified outliers by FST quantile; reuse that, just plot Z on y-axis.
    # Also compute per-chrom Z-score quantile thresholds for reference lines.
    z_high_q = df.loc[df["outlier"] == "high", "zscore"].min()
    z_low_q  = df.loc[df["outlier"] == "low",  "zscore"].max()
    z_abs_max = df["zscore"].abs().quantile(0.999)
    ylim = max(z_abs_max, 4) * 1.05

    for i, chrom in enumerate(CHROM_ORDER):
        ax  = axes[i]
        sub = df[df["chr"] == chrom].copy()
        if len(sub) == 0:
            ax.set_visible(False)
            continue

        pos_mb = sub["midPos"] / 1e6
        zsc    = sub["zscore"]

        mask_neu  = sub["outlier"] == "neutral"
        mask_high = sub["outlier"] == "high"
        mask_low  = sub["outlier"] == "low"

        ax.fill_between(pos_mb[mask_neu],  zsc[mask_neu],  color=COL_NEU,  alpha=0.55, lw=0)
        ax.fill_between(pos_mb[mask_high], zsc[mask_high], color=COL_HIGH, alpha=0.85, lw=0)
        ax.fill_between(pos_mb[mask_low],  zsc[mask_low],  color=COL_LOW,  alpha=0.85, lw=0)
        ax.plot(pos_mb, zsc, color="white", lw=0.2, alpha=0.4)

        ax.axhline(0,       color="0.5",    lw=0.6, ls=":",  alpha=0.6)
        ax.axhline(z_high_q, color=COL_HIGH, lw=0.7, ls="--", alpha=0.6)
        ax.axhline(z_low_q,  color=COL_LOW,  lw=0.7, ls="--", alpha=0.6)

        ax.set_xlim(0, pos_mb.max() * 1.01)
        ax.set_ylim(-ylim, ylim)
        ax.set_ylabel(f"Chr {CHROM_LABELS[chrom]}", fontsize=8, rotation=0,
                      labelpad=28, va="center")
        ax.tick_params(labelsize=7)

    _label_axes(axes, "Modified Z-score (per-chromosome)")
    legend_handles = [
        mpatches.Patch(color=COL_HIGH, label="High FST (top 2.5% genome-wide)"),
        mpatches.Patch(color=COL_LOW,  label="Low FST (bottom 2.5% genome-wide)"),
        mpatches.Patch(color=COL_NEU,  label="Neutral"),
    ]
    axes[0].legend(handles=legend_handles, fontsize=8, loc="upper right",
                   framealpha=0.8, ncol=3)
    axes[0].set_title(
        "Sliding-window FST: A. cervicornis FL vs A. palmata FL  "
        "(50 kb windows, 10 kb step)  —  Style 2: per-chromosome Z-score",
        fontsize=11, pad=6,
    )
    fig.savefig(out + "_stacked_zscore.png", dpi=180, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out}_stacked_zscore.png")


# ── Plot style 3: line trace + outlier block spans ────────────────────────────

def plot_stacked_regions(df, out, high_thresh, low_thresh):
    """
    Thin gray FST line trace with semi-transparent spans for outlier blocks
    (MIN_BLOCK+ consecutive windows). Cleanest for showing where the regions are.
    """
    fig, axes = _stacked_fig()

    for i, chrom in enumerate(CHROM_ORDER):
        ax  = axes[i]
        sub = df[df["chr"] == chrom].copy()
        if len(sub) == 0:
            ax.set_visible(False)
            continue

        pos_mb = sub["midPos"] / 1e6
        fst    = sub["Fst"]
        step   = np.median(np.diff(pos_mb.values)) if len(sub) > 1 else 0.01

        # Gray FST line
        ax.plot(pos_mb, fst, color="0.55", lw=0.5, alpha=0.8)

        # Shade all outlier windows (single + block) lightly
        mask_high = sub["outlier"] == "high"
        mask_low  = sub["outlier"] == "low"
        ax.fill_between(pos_mb[mask_high], fst[mask_high],
                        color=COL_HIGH, alpha=0.25, lw=0)
        ax.fill_between(pos_mb[mask_low],  fst[mask_low],
                        color=COL_LOW,  alpha=0.25, lw=0)

        # Dark spans for blocks (≥MIN_BLOCK consecutive)
        for otype, col in [("high", COL_HIGH), ("low", COL_LOW)]:
            spans = get_outlier_spans(sub, "midPos", otype)
            for (s_bp, w_bp) in spans:
                s_mb = s_bp / 1e6
                w_mb = w_bp / 1e6
                ax.axvspan(s_mb, s_mb + w_mb, color=col, alpha=0.45, lw=0)

        # Threshold lines
        ax.axhline(high_thresh,        color=COL_HIGH, lw=0.8, ls="--", alpha=0.7)
        ax.axhline(low_thresh,         color=COL_LOW,  lw=0.8, ls="--", alpha=0.7)
        ax.axhline(df["Fst"].median(), color="0.4",    lw=0.5, ls=":",  alpha=0.5)

        ax.set_xlim(0, pos_mb.max() * 1.01)
        ax.set_ylim(-0.02, 1.02)
        ax.set_ylabel(f"Chr {CHROM_LABELS[chrom]}", fontsize=8, rotation=0,
                      labelpad=28, va="center")
        ax.tick_params(labelsize=7)

    _label_axes(axes, "FST (unweighted)")
    legend_handles = [
        mpatches.Patch(color=COL_HIGH, label=f"High FST block (≥{MIN_BLOCK} windows, >{high_thresh:.3f})"),
        mpatches.Patch(color=COL_LOW,  label=f"Low FST block (≥{MIN_BLOCK} windows, <{low_thresh:.3f})"),
        mpatches.Patch(facecolor=COL_HIGH, alpha=0.25, label="Single outlier window (high)"),
        mpatches.Patch(facecolor=COL_LOW,  alpha=0.25, label="Single outlier window (low)"),
    ]
    axes[0].legend(handles=legend_handles, fontsize=7.5, loc="upper right",
                   framealpha=0.8, ncol=2)
    axes[0].set_title(
        "Sliding-window FST: A. cervicornis FL vs A. palmata FL  "
        "(50 kb windows, 10 kb step)  —  Style 3: line + outlier blocks",
        fontsize=11, pad=6,
    )
    fig.savefig(out + "_stacked_regions.png", dpi=180, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out}_stacked_regions.png")


# ── Summary tables ───────────────────────────────────────────────────────────

def print_outlier_summary(df, gene_df, out, high_thresh, low_thresh):
    high_df = df[df["outlier"] == "high"].sort_values("Fst", ascending=False)
    low_df  = df[df["outlier"] == "low"].sort_values("Fst")

    print(f"\n{'='*65}")
    print(f"  Outlier summary  (top/bottom {(1-HIGH_Q)*100:.1f}% genome-wide quantiles)")
    print(f"  Genome-wide: median FST = {df['Fst'].median():.4f}  "
          f"mean = {df['Fst'].mean():.4f}")
    print(f"  High threshold: FST > {high_thresh:.4f}  "
          f"Low threshold: FST < {low_thresh:.4f}")
    print(f"{'='*65}")

    print(f"\n  High-FST windows: {len(high_df)}  (candidates: inversions / barriers)")
    print(f"  Low-FST  windows: {len(low_df)}   (candidates: introgression tracts)\n")

    print("  Per-chromosome breakdown:")
    for chrom in CHROM_ORDER:
        n_high = len(high_df[high_df["chr"] == chrom])
        n_low  = len(low_df[low_df["chr"]  == chrom])
        if n_high + n_low > 0:
            print(f"    Chr {CHROM_LABELS[chrom]:>2}: {n_high} high, {n_low} low")

    print(f"\n  Top 20 HIGH-FST windows:")
    _print_windows(high_df.head(20))
    print(f"\n  Top 20 LOW-FST windows:")
    _print_windows(low_df.head(20))

    if len(gene_df):
        print(f"\n{'='*65}")
        print(f"  Genes in HIGH-FST outlier windows ({len(gene_df[gene_df['outlier_type']=='high'])} unique genes)")
        print(f"{'='*65}")
        _print_genes(gene_df[gene_df["outlier_type"] == "high"])

        print(f"\n{'='*65}")
        print(f"  Genes in LOW-FST outlier windows ({len(gene_df[gene_df['outlier_type']=='low'])} unique genes)")
        print(f"{'='*65}")
        _print_genes(gene_df[gene_df["outlier_type"] == "low"])

        for otype in ["high", "low"]:
            sub = gene_df[gene_df["outlier_type"] == otype][[
                "window_chr","window_midPos","window_Fst","window_zscore",
                "chr","start","end","gene_id","description","biotype"
            ]].rename(columns={"window_chr":"fst_chr","window_midPos":"fst_midPos",
                                "window_Fst":"FST","window_zscore":"Z"})
            sub["Chr"]    = sub["fst_chr"].map(CHROM_LABELS).apply(lambda x: f"Chr{x}")
            sub["pos_Mb"] = (sub["fst_midPos"] / 1e6).round(2)
            fname = out + f"_genes_{otype}_fst.tsv"
            sub.to_csv(fname, sep="\t", index=False)
            print(f"\n  Gene table saved: {fname}")


def _print_windows(df):
    for _, r in df.iterrows():
        print(f"    Chr{CHROM_LABELS[r['chr']]:>2}  {r['midPos']/1e6:6.2f} Mb  "
              f"FST={r['Fst']:.4f}  Z={r['zscore']:+.1f}  N={r['Nsites']}")


def _print_genes(df):
    shown = set()
    asc   = r_outlier(df) == "low"
    for _, r in df.sort_values("window_Fst", ascending=asc).iterrows():
        key = (r["chr"], r["start"])
        if key in shown:
            continue
        shown.add(key)
        desc = r["description"][:55] if r["description"] else "uncharacterized"
        gene = r["gene_id"] if r["gene_id"] else r["name"]
        print(f"    Chr{CHROM_LABELS[r['chr']]:>2}  {r['start']/1e6:.2f}–{r['end']/1e6:.2f} Mb  "
              f"{gene:<20}  {desc}")
        if len(shown) >= 40:
            print("    ... (see TSV for full list)")
            break


def r_outlier(df):
    return df["outlier_type"].iloc[0] if len(df) else "high"


# ── Main ─────────────────────────────────────────────────────────────────────

def load_windows(path):
    df = pd.read_csv(path, sep="\t", header=0,
                     names=["region","chr","midPos","Nsites","Fst"])
    df = df[df["chr"].isin(CHROM_ORDER)].copy()
    df["chr"] = pd.Categorical(df["chr"], categories=CHROM_ORDER, ordered=True)
    df = df.sort_values(["chr","midPos"]).reset_index(drop=True)
    return df


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fst",  required=True, help="FST windows file")
    parser.add_argument("--gff",  required=True, help="GFF3 annotation (gzipped OK)")
    parser.add_argument("--out",  default="fst_windows", help="Output prefix")
    args = parser.parse_args()

    print(f"Loading FST windows: {args.fst}")
    df = load_windows(args.fst)
    df, high_thresh, low_thresh = call_outliers(df)

    n_high  = (df["outlier"] == "high").sum()
    n_low   = (df["outlier"] == "low").sum()
    n_block = df["block"].sum()
    print(f"Windows: {len(df)} total  |  {n_high} high-FST  |  {n_low} low-FST  "
          f"|  {n_block} in outlier blocks (≥{MIN_BLOCK} consecutive)")
    print(f"Thresholds: FST > {high_thresh:.4f} (high) | FST < {low_thresh:.4f} (low)")

    print(f"\nLoading GFF: {args.gff}")
    genes = parse_gff_genes(args.gff)
    print(f"Loaded {len(genes)} gene records from GFF")

    print("\nIntersecting outlier windows with genes...")
    outlier_df = df[df["outlier"] != "neutral"]
    gene_hits  = find_genes_in_windows(outlier_df, genes)

    print_outlier_summary(df, gene_hits, args.out, high_thresh, low_thresh)

    print("\nGenerating stacked chromosome plots (3 styles)...")
    plot_stacked_fst(df, args.out, high_thresh, low_thresh)
    plot_stacked_zscore(df, args.out)
    plot_stacked_regions(df, args.out, high_thresh, low_thresh)

    df[df["outlier"] == "high"].to_csv(args.out + "_high_fst.tsv", sep="\t", index=False)
    df[df["outlier"] == "low"].to_csv(args.out + "_low_fst.tsv",  sep="\t", index=False)
    print(f"Outlier window tables saved.")


if __name__ == "__main__":
    main()
