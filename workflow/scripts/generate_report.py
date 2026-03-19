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

# Consistent scheme used across all plots
REGION_ORDER   = ["FL", "PA", "BON"]
REGION_COLORS  = {"FL": "#2166AC", "PA": "#4DAC26", "BON": "#D6604D"}
SPECIES_MARKERS = {"Acervicornis": "o", "Apalmata": "^"}
# AC = solid fill, AP = lighter/hatched
SPECIES_ALPHA  = {"Acervicornis": 0.85, "Apalmata": 0.45}

def sp_abbr(species):
    if not species:
        return "?"
    return "AC" if "cerv" in species.lower() else "AP"

def group_label(species, region):
    return f"{sp_abbr(species)}_{region}"

def group_color(species, region):
    base = REGION_COLORS.get(region, "#888888")
    # Lighten Apal groups slightly so AC/AP are visually distinct at same region color
    if sp_abbr(species) == "AP":
        import matplotlib.colors as mc
        r, g, b = mc.to_rgb(base)
        # blend toward white
        r, g, b = [0.4 + 0.6 * x for x in (r, g, b)]
        return mc.to_hex((r, g, b))
    return base


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
    """Load per-sample depth from individual {sample}.depth.txt files in qc/."""
    qc_dir = Path(results_dir) / "qc"
    records = []
    for p in sorted(qc_dir.glob("*.depth.txt")):
        sample = p.stem.replace(".depth", "")
        try:
            val = float(p.read_text().strip().split()[0])
            records.append({"sample_id": sample, "mean_depth": val})
        except Exception:
            pass
    if not records:
        # fallback: depth_summary.txt
        path = qc_dir / "depth_summary.txt"
        if path.exists():
            df = pd.read_csv(path, sep="\t")
            df.columns = ["sample_id", "mean_depth"]
            return df.set_index("sample_id")
        return pd.DataFrame()
    return pd.DataFrame(records).set_index("sample_id")


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
    """Per-sample depth as violin plots grouped by AC_FL/AP_FL etc."""
    plt, _ = setup_matplotlib()
    if depth_df is None or depth_df.empty:
        return ""

    spec_col = "species" if "species" in metadata.columns else None
    reg_col  = "region"  if "region"  in metadata.columns else None
    join_cols = [c for c in [spec_col, reg_col] if c]
    df = depth_df.join(metadata[join_cols], how="inner")

    species_list = ["Acervicornis", "Apalmata"] if spec_col else [None]
    regions      = [r for r in REGION_ORDER if reg_col and r in df[reg_col].values]
    groups = [(sp, reg) for sp in species_list for reg in regions
              if len(df[(df[spec_col] == sp) & (df[reg_col] == reg)]) > 0]
    xlabels = [group_label(sp, reg) for sp, reg in groups]
    colors  = [group_color(sp, reg) for sp, reg in groups]

    fig, ax = plt.subplots(figsize=(max(8, len(groups) * 1.2), 5))

    all_depths = []
    for gi, (sp, reg) in enumerate(groups):
        mask = pd.Series([True] * len(df), index=df.index)
        if spec_col: mask &= df[spec_col] == sp
        if reg_col:  mask &= df[reg_col]  == reg
        vals = df.loc[mask, "mean_depth"].astype(float).values
        all_depths.extend(vals)
        if len(vals) < 2:
            ax.scatter([gi], vals, color=colors[gi], s=60, zorder=3)
            continue
        vp = ax.violinplot(vals, positions=[gi], widths=0.7,
                           showmedians=True, showextrema=True)
        for pc in vp["bodies"]:
            pc.set_facecolor(colors[gi])
            pc.set_alpha(0.75)
            pc.set_edgecolor("white")
        for part in ["cmedians", "cmins", "cmaxes", "cbars"]:
            vp[part].set_color("black")
            vp[part].set_linewidth(1.2)

    ax.axhline(10, color="red", linestyle="--", linewidth=1, label="10× threshold")
    mean_d = np.mean(all_depths)
    ax.axhline(mean_d, color="gray", linestyle="-", linewidth=1,
               alpha=0.7, label=f"Mean {mean_d:.1f}×")
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels(xlabels, fontsize=10)
    ax.set_ylabel("Mean depth (×)")
    ax.set_title("Sequencing depth per sample (grouped by species × region)")
    ax.legend(fontsize=9)
    return save_fig(fig, figures_dir, "depth")


REGION_COLORS = {"FL": "#2166AC", "PA": "#4DAC26", "BON": "#D6604D"}
SPECIES_MARKERS = {"Acervicornis": "o", "Apalmata": "^"}

def fig_pca(eigenvectors, pct_var, metadata, figures_dir):
    """PCA: shape = species (circle=Acer, triangle=Apal), color = region (FL/PA/BON)."""
    plt, _ = setup_matplotlib()
    if eigenvectors is None:
        return ""
    samples = metadata.index.tolist()

    region_col  = "region"  if "region"  in metadata.columns else None
    species_col = "species" if "species" in metadata.columns else None

    fig, ax = plt.subplots(figsize=(8, 7))

    # Plot by species × region combination
    species_list = sorted(metadata[species_col].unique()) if species_col else [None]
    region_list  = sorted(metadata[region_col].unique())  if region_col  else [None]

    for sp in species_list:
        for reg in region_list:
            mask = pd.Series([True] * len(metadata), index=metadata.index)
            if sp  and species_col: mask &= metadata[species_col] == sp
            if reg and region_col:  mask &= metadata[region_col]  == reg
            idx = [i for i, s in enumerate(samples) if s in metadata.index and mask.get(s, False)]
            if not idx:
                continue
            sp_short  = sp.replace("Acervicornis", "Acer").replace("Apalmata", "Apal") if sp else ""
            marker = SPECIES_MARKERS.get(sp, "o") if sp else "o"
            color  = REGION_COLORS.get(reg, "#888888") if reg else "#888888"
            label  = f"{sp_short} — {reg}" if sp_short and reg else (sp_short or reg or "")
            ax.scatter(eigenvectors[idx, 0], eigenvectors[idx, 1],
                       c=color, marker=marker, label=label,
                       s=60, edgecolors="white", linewidth=0.4, alpha=0.85, zorder=3)

    ax.set_xlabel(f"PC1 ({pct_var[0]:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({pct_var[1]:.1f}% variance)")
    ax.set_title("PCA — LD-pruned SNPs")
    ax.grid(True, alpha=0.25)

    # Legend: shapes for species, filled circles for region
    from matplotlib.lines import Line2D
    legend_handles = []
    for sp in species_list:
        m = SPECIES_MARKERS.get(sp, "o") if sp else "o"
        sp_short = sp_abbr(sp) if sp else "?"
        legend_handles.append(Line2D([0], [0], marker=m, color="w", markerfacecolor="#444",
                                     markersize=9, label=sp_short, markeredgecolor="#444"))
    for reg in sorted(region_list, key=lambda r: REGION_ORDER.index(r) if r in REGION_ORDER else 99):
        c = REGION_COLORS.get(reg, "#888") if reg else "#888"
        legend_handles.append(Line2D([0], [0], marker="o", color="w", markerfacecolor=c,
                                     markersize=9, label=reg or "unknown", markeredgecolor="white"))
    ax.legend(handles=legend_handles, frameon=True, fontsize=9,
              title="shape=species  color=region")
    return save_fig(fig, figures_dir, "pca")


def fig_admixture(admix, metadata, figures_dir):
    """
    K=2 and K=3 admixture shown as violin plots.
    x-axis: AC_FL, AC_PA, AC_BON, AP_FL, AP_PA, AP_BON
    y-axis: Q proportion for each ancestry component.
    Color = region, shade = species (dark=AC, light=AP).
    """
    plt, _ = setup_matplotlib()
    if not admix:
        return ""

    samples   = metadata.index.tolist()
    spec_col  = "species"    if "species"    in metadata.columns else None
    reg_col   = "region"     if "region"     in metadata.columns else None

    # Build ordered group list: AC_FL, AC_PA, AC_BON, AP_FL, AP_PA, AP_BON
    species_list = ["Acervicornis", "Apalmata"] if spec_col else [None]
    regions      = [r for r in REGION_ORDER
                    if reg_col and r in metadata[reg_col].values]
    groups = [(sp, reg) for sp in species_list for reg in regions
              if len([s for s in samples
                      if s in metadata.index
                      and (not spec_col or metadata.loc[s, spec_col] == sp)
                      and (not reg_col  or metadata.loc[s, reg_col]  == reg)]) > 0]
    xlabels = [group_label(sp, reg) for sp, reg in groups]
    colors  = [group_color(sp, reg) for sp, reg in groups]

    k_show = [k for k in [2, 3] if k in admix]
    if not k_show:
        k_show = sorted(admix.keys())[:2]

    fig, axes = plt.subplots(1, len(k_show), figsize=(7 * len(k_show), 5), squeeze=False)

    for ki, k in enumerate(k_show):
        ax    = axes[0][ki]
        q_all = admix[k]

        # For each component, collect per-group Q values
        # Identify which component most strongly corresponds to each ancestry
        # by ranking components by mean Q in the first group
        component_means = [np.mean([q_all[samples.index(s), j]
                                    for s in samples
                                    if s in metadata.index]) for j in range(k)]
        comp_order = np.argsort(component_means)[::-1]

        positions = np.arange(len(groups))
        width = 0.8 / k

        for ci, comp in enumerate(comp_order):
            offset = (ci - (k - 1) / 2) * width
            for gi, (sp, reg) in enumerate(groups):
                idx = [samples.index(s) for s in samples
                       if s in metadata.index
                       and (not spec_col or metadata.loc[s, spec_col] == sp)
                       and (not reg_col  or metadata.loc[s, reg_col]  == reg)]
                if not idx:
                    continue
                vals = q_all[idx, comp]
                if len(vals) < 2:
                    ax.bar(gi + offset, vals[0], width=width * 0.9,
                           color=colors[gi], alpha=0.5 + 0.5 * (ci == 0))
                    continue
                vp = ax.violinplot(vals, positions=[gi + offset],
                                   widths=width * 0.9,
                                   showmedians=True, showextrema=False)
                for pc in vp["bodies"]:
                    pc.set_facecolor(colors[gi])
                    pc.set_alpha(0.85 if ci == 0 else 0.4)
                    pc.set_edgecolor("white")
                vp["cmedians"].set_color("black")
                vp["cmedians"].set_linewidth(1.5)

        ax.set_xticks(positions)
        ax.set_xticklabels(xlabels, rotation=30, ha="right", fontsize=9)
        ax.set_ylim(0, 1)
        ax.set_ylabel("Admixture proportion (Q)", fontsize=10)
        ax.set_title(f"K = {k}", fontsize=11, fontweight="bold")
        ax.axhline(0.5, color="gray", linestyle="--", linewidth=0.8, alpha=0.5)
        # Species separator
        if spec_col and len(species_list) > 1:
            n_ac = sum(1 for sp, _ in groups if sp_abbr(sp) == "AC")
            if 0 < n_ac < len(groups):
                ax.axvline(x=n_ac - 0.5, color="black", linewidth=1.5)

    fig.suptitle("Admixture proportions — K=2 and K=3 (PCAngsd)", fontsize=13)
    plt.tight_layout()
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

    spec_col = "species" if "species" in metadata.columns else None
    reg_col  = "region"  if "region"  in metadata.columns else None
    join_cols = [c for c in [spec_col, reg_col] if c]
    df = het_df.join(metadata[join_cols], how="inner")
    if join_cols:
        sort_key = [c for c in [spec_col, reg_col] if c]
        df["_reg_order"] = df[reg_col].map({r: i for i, r in enumerate(REGION_ORDER)}) if reg_col else 0
        df = df.sort_values(([spec_col] if spec_col else []) + ["_reg_order", "heterozygosity"])

    # Assign color per sample by species × region
    def _bar_color(row):
        sp  = row.get(spec_col, "") if spec_col else ""
        reg = row.get(reg_col,  "") if reg_col  else ""
        return group_color(sp, reg)

    colors = [_bar_color(row) for _, row in df.iterrows()]
    mean_h = df["heterozygosity"].mean()
    sd_h   = df["heterozygosity"].std()

    fig, ax = plt.subplots(figsize=(14, 5))
    ax.bar(range(len(df)), df["heterozygosity"] * 100, color=colors, edgecolor="none", width=1.0)
    ax.axhline(mean_h * 100, color="gray", linewidth=1, linestyle="-", alpha=0.6, label="Mean")
    ax.axhline((mean_h + 2 * sd_h) * 100, color="red", linewidth=1,
               linestyle="--", alpha=0.7, label="±2 SD")
    ax.axhline((mean_h - 2 * sd_h) * 100, color="red", linewidth=1, linestyle="--", alpha=0.7)

    # Group boundary lines + AC_FL / AP_FL style x-axis labels
    if join_cols:
        cur_grp  = (df.iloc[0].get(spec_col,""), df.iloc[0].get(reg_col,""))
        start    = 0
        ticks, tick_labels = [], []
        rows = list(df.iterrows())
        for i, (_, row) in enumerate(rows):
            grp = (row.get(spec_col,""), row.get(reg_col,""))
            if grp != cur_grp:
                ax.axvline(x=i - 0.5, color="black", linewidth=1.2, alpha=0.4)
                ticks.append((start + i) / 2)
                tick_labels.append(group_label(*cur_grp))
                start, cur_grp = i, grp
        ticks.append((start + len(df)) / 2)
        tick_labels.append(group_label(*cur_grp))
        ax.set_xticks(ticks)
        ax.set_xticklabels(tick_labels, fontsize=9)
    else:
        ax.set_xticks([])

    ax.set_xlim(-0.5, len(df) - 0.5)
    ax.set_ylabel("Heterozygosity (%)")
    ax.set_title("Per-Individual Heterozygosity")
    from matplotlib.lines import Line2D
    handles = [Line2D([0],[0], marker="o", color="w", markerfacecolor=group_color(sp, reg),
                      markersize=9, label=group_label(sp, reg), markeredgecolor="white")
               for sp in (["Acervicornis","Apalmata"] if spec_col else [""])
               for reg in ([r for r in REGION_ORDER if reg_col and r in df[reg_col].values] if reg_col else [""])]
    handles += [Line2D([0],[0], color="gray",  label="Mean"),
                Line2D([0],[0], color="red", linestyle="--", label="±2 SD")]
    ax.legend(handles=handles, fontsize=8, ncol=4)
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


def fig_kinship(rel_matrix, metadata, figures_dir):
    """
    Kinship heatmaps split by species, then by population within species.
    One subplot per species, samples ordered by population.
    """
    plt, gs_mod = setup_matplotlib()
    if rel_matrix is None or metadata is None:
        return ""

    species_col = "species" if "species" in metadata.columns else None
    pop_col     = "population" if "population" in metadata.columns else None

    # Build ordered sample list per species
    def ordered_samples_for(species=None):
        if species and species_col:
            meta = metadata[metadata[species_col] == species]
        else:
            meta = metadata
        if pop_col:
            meta = meta.sort_values([pop_col, meta.index.name or "sample_id"])
        return [s for s in meta.index if s in rel_matrix.index]

    species_list = sorted(metadata[species_col].unique()) if species_col else [None]
    n_sp = len(species_list)

    fig, axes = plt.subplots(1, n_sp, figsize=(min(10, 5 * n_sp), 5 * n_sp))
    if n_sp == 1:
        axes = [axes]

    for ax, sp in zip(axes, species_list):
        samps = ordered_samples_for(sp)
        if not samps:
            continue
        mat = rel_matrix.loc[samps, samps].values.astype(float)
        np.fill_diagonal(mat, 0)
        im = ax.imshow(mat, cmap="YlOrRd", vmin=0, vmax=0.5, aspect="auto")
        title = sp if sp else "All"
        ax.set_title(f"Kinship — {title}", fontsize=11)

        # Draw population boundary lines
        if pop_col and species_col:
            pops = metadata.loc[samps, pop_col]
            prev = pops.iloc[0]
            ticks, tick_labels = [], []
            start = 0
            for i, p in enumerate(pops):
                if p != prev:
                    ax.axvline(x=i - 0.5, color="white", linewidth=1)
                    ax.axhline(y=i - 0.5, color="white", linewidth=1)
                    ticks.append((start + i) / 2)
                    tick_labels.append(prev)
                    start = i
                    prev = p
            ticks.append((start + len(pops)) / 2)
            tick_labels.append(prev)
            ax.set_xticks(ticks)
            ax.set_xticklabels(tick_labels, rotation=45, ha="right", fontsize=8)
            ax.set_yticks(ticks)
            ax.set_yticklabels(tick_labels, fontsize=8)
        else:
            ax.set_xticks([])
            ax.set_yticks([])

        plt.colorbar(im, ax=ax, label="KING Kinship", shrink=0.6)

    fig.suptitle("Pairwise Kinship (KING)", fontsize=13)
    plt.tight_layout()
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


# ─── Lineage loaders and figure ───────────────────────────────────────────────

def load_lineage_assignments(results_dir):
    """Load lineage_assignments.txt → dict sample_id: lineage."""
    path = Path(results_dir) / "admixture" / "lineage_assignments.txt"
    if not path.exists():
        return {}
    assignments = {}
    with open(path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                assignments[parts[0]] = parts[1]
    return assignments


def fig_lineage(admix_k2, lineage_assignments, metadata, figures_dir):
    """
    K=2 Q bar plot with samples colored and sorted by lineage assignment.
    lineageA=blue, lineageB=red, admixed=gray.
    """
    plt, _ = setup_matplotlib()
    if admix_k2 is None or not lineage_assignments:
        return ""

    unrel_path = Path(figures_dir).parent / "relatedness" / "unrelated_samples.txt"
    if not unrel_path.exists():
        return ""
    with open(unrel_path) as f:
        samples = [l.strip() for l in f if l.strip()]

    lineage_color = {"lineageA": "#2166AC", "lineageB": "#D6604D", "admixed": "#AAAAAA"}

    def sort_key(i):
        s   = samples[i]
        lin = lineage_assignments.get(s, "admixed")
        order = {"lineageA": 0, "lineageB": 1, "admixed": 2}
        pop = (metadata.loc[s, "population"]
               if s in metadata.index and "population" in metadata.columns else "")
        return (order.get(lin, 2), pop, s)

    order_idx       = sorted(range(len(samples)), key=sort_key)
    ordered_samples = [samples[i] for i in order_idx]
    q               = admix_k2[order_idx, :]

    fig, ax = plt.subplots(figsize=(max(14, len(samples) * 0.5), 4))

    for j in range(q.shape[0]):
        s   = ordered_samples[j]
        lin = lineage_assignments.get(s, "admixed")
        ax.bar(j, q[j, 0], color=lineage_color.get(lin, "#AAAAAA"), width=1.0, edgecolor="none")
        ax.bar(j, q[j, 1], bottom=q[j, 0], color="#DDDDDD",         width=1.0, edgecolor="none")

    prev_lin = lineage_assignments.get(ordered_samples[0], "admixed")
    for j in range(1, len(ordered_samples)):
        lin = lineage_assignments.get(ordered_samples[j], "admixed")
        if lin != prev_lin:
            ax.axvline(x=j - 0.5, color="black", linewidth=2)
            prev_lin = lin

    ax.set_xlim(-0.5, len(samples) - 0.5)
    ax.set_ylim(0, 1)
    ax.set_xticks(range(len(ordered_samples)))
    ax.set_xticklabels(ordered_samples, rotation=90, fontsize=7)
    ax.set_ylabel("Q (K=2)")
    ax.set_title("K=2 Admixture — Lineage Assignment")
    ax.axhline(0.8, color="black", linestyle="--", linewidth=1, alpha=0.5, label="Q=0.80 threshold")

    from matplotlib.patches import Patch
    ax.legend(handles=[
        Patch(color=lineage_color["lineageA"], label="lineageA"),
        Patch(color=lineage_color["lineageB"], label="lineageB"),
        Patch(color=lineage_color["admixed"],  label="admixed"),
    ], loc="upper right", fontsize=9)

    return save_fig(fig, figures_dir, "lineage_assignment")


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results_dir", required=True)
    parser.add_argument("--metadata",    required=True)
    parser.add_argument("--config",      default=None)
    parser.add_argument("--output",      default="results/report.html")
    parser.add_argument("--figures_dir", default="results/figures")
    parser.add_argument("--fst-comparisons", nargs="*", dest="fst_comparisons", default=None,
                        help="FST comparisons as 'g1_vs_g2' strings (space-separated)")
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

    # Parse FST comparisons — from CLI arg or fall back to pop_pairs
    fst_comparisons = []
    if args.fst_comparisons:
        for comp in args.fst_comparisons:
            parts = comp.split("_vs_", 1)
            if len(parts) == 2:
                fst_comparisons.append(tuple(parts))
    if not fst_comparisons:
        fst_comparisons = pop_pairs

    print("Loading pipeline outputs...")
    filt_df  = load_filtering_summary(results_dir)
    depth_df = load_depth_summary(results_dir)
    fastp_df = load_fastp_stats(results_dir, samples)
    eigenvectors, pct_var = load_pca(results_dir)
    admix    = load_admixture(results_dir, pipeline_config.get("max_k", 5))
    fst_df   = load_fst_global(results_dir, fst_comparisons)
    het_df   = load_heterozygosity(results_dir, samples)
    rel_mat  = load_relatedness_matrix(results_dir)
    ld_df    = load_ld_decay(results_dir)
    lineage_assignments = load_lineage_assignments(results_dir)

    # Diversity: load thetas for all FST groups + standard pops
    all_groups = list({g for pair in fst_comparisons for g in pair} | set(pops))
    thetas = {g: load_thetas(results_dir, g) for g in all_groups}
    sfs    = {p: load_sfs(results_dir, p) for p in pops}  # SFS plot uses pop groups only
    n_pass1, n_pass2 = count_snps(results_dir)

    # Build PCA-ordered metadata (eigenvector row i = sample i in unrelated_samples.txt)
    unrel_path = results_dir / "relatedness" / "unrelated_samples.txt"
    if unrel_path.exists():
        pca_sample_order = [l.strip() for l in open(unrel_path) if l.strip()]
        pca_meta = metadata.reindex([s for s in pca_sample_order if s in metadata.index])
    else:
        n_pca = eigenvectors.shape[0] if eigenvectors is not None else len(samples)
        pca_meta = metadata.iloc[:n_pca]

    print("Generating figures...")
    img_mapping  = fig_mapping_rates(filt_df, metadata, figures_dir)
    img_depth    = fig_depth(depth_df, metadata, figures_dir)
    img_pca      = fig_pca(eigenvectors, pct_var, pca_meta, figures_dir)
    img_admix    = fig_admixture(admix, pca_meta, figures_dir)
    img_lineage  = fig_lineage(admix.get(2), lineage_assignments, pca_meta, figures_dir)
    img_het      = fig_heterozygosity(het_df, metadata, figures_dir)
    img_sfs      = fig_sfs(sfs, figures_dir)
    img_div      = fig_diversity({p: thetas[p] for p in pops if thetas.get(p) is not None}, figures_dir)
    img_kinship  = fig_kinship(rel_mat, metadata, figures_dir)
    img_ld       = fig_ld_decay(ld_df, figures_dir)
    img_fst_man  = {}
    for g1, g2 in fst_comparisons:
        wfst = load_windowed_fst(results_dir, g1, g2)
        img_fst_man[(g1, g2)] = fig_fst_manhattan(g1, g2, wfst, figures_dir)

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

    # Sample counts by group
    spec_col = "species" if "species" in metadata.columns else None
    reg_col  = "region"  if "region"  in metadata.columns else None
    group_counts = []
    if spec_col and reg_col:
        for sp in ["Acervicornis", "Apalmata"]:
            for reg in REGION_ORDER:
                n = int(((metadata[spec_col] == sp) & (metadata[reg_col] == reg)).sum())
                if n > 0:
                    group_counts.append(f"{group_label(sp, reg)} n={n}")
    sample_desc = " | ".join(group_counts) if group_counts else f"{len(samples)} samples"

    boxes = [metric_box(len(samples), "Total samples"),
             metric_box(sum(1 for l in open(unrel_path) if l.strip()) if unrel_path.exists() else "—", "Unrelated")]
    if n_pass1:
        boxes.append(metric_box(f"{n_pass1:,}", "Pass-1 SNPs"))
    if n_pass2:
        boxes.append(metric_box(f"{n_pass2:,}", "Pass-2 SNPs"))
    if not depth_df.empty:
        md = depth_df["mean_depth"].astype(float).mean()
        boxes.append(metric_box(f"{md:.1f}×", "Mean depth"))
    parts.append("".join(boxes))
    parts.append(f"<p style='color:#555;font-size:0.95em'>{sample_desc}</p>")

    if flagged:
        parts.append(f'<div class="warn">⚠ {len(flagged)} sample(s) flagged for low depth (&lt;10×): '
                     f'{", ".join(flagged[:20])}{"…" if len(flagged) > 20 else ""}</div>')
    else:
        parts.append('<div class="ok">✓ All samples passed depth QC (≥10×)</div>')

    # ── Sequencing QC
    parts.append("<h2>2. Sequencing &amp; Mapping Quality</h2>")
    parts.append("<h3>2a. Per-sample sequencing depth</h3>")
    if not depth_df.empty:
        d = depth_df["mean_depth"].astype(float)
        parts.append(f"<p>n={len(d)} samples | mean={d.mean():.1f}× | median={d.median():.1f}× | "
                     f"min={d.min():.1f}× | max={d.max():.1f}× | "
                     f"samples &lt;10×: {(d < 10).sum()}</p>")
    parts.append(img_depth)
    parts.append("<h3>2b. Filtering Summary</h3>")
    parts.append(df_to_html(filt_df, "{:.4f}"))
    parts.append(img_mapping)
    parts.append("<h3>2c. fastp Read Quality</h3>")
    parts.append(df_to_html(fastp_df, "{:.4f}"))

    # ── Population structure
    parts.append("<h2>3. Population Structure</h2>")
    parts.append("<h3>3a. PCA</h3>")
    parts.append(img_pca or "<p><em>PCA not yet available.</em></p>")
    parts.append("<h3>3b. Admixture (all K)</h3>")
    parts.append(img_admix or "<p><em>Admixture not yet available.</em></p>")
    if img_lineage:
        parts.append("<h3>3c. Lineage Assignment (K=2)</h3>")
        parts.append(img_lineage)
        if lineage_assignments:
            counts = {}
            for lin in lineage_assignments.values():
                counts[lin] = counts.get(lin, 0) + 1
            parts.append("<p>" + "  |  ".join(
                f"<b>{lin}</b>: {n}" for lin, n in sorted(counts.items())
            ) + "</p>")

    # ── Diversity
    parts.append("<h2>4. Genetic Diversity</h2>")
    parts.append("<h3>4a. Site Frequency Spectra</h3>")
    parts.append(img_sfs or "<p><em>SFS not yet available.</em></p>")
    parts.append("<h3>4b. π, θ, Tajima's D (by population)</h3>")
    parts.append(img_div or "<p><em>Diversity stats not yet available.</em></p>")
    parts.append("<h3>4c. Individual Heterozygosity</h3>")
    parts.append(img_het or "<p><em>Heterozygosity not yet available.</em></p>")
    if het_df is not None and not het_df.empty:
        parts.append(df_to_html(het_df, "{:.6f}"))

    # ── Differentiation
    parts.append("<h2>5. Population Differentiation (FST)</h2>")
    if not fst_df.empty:
        parts.append(df_to_html(fst_df, "{:.4f}"))
    else:
        parts.append("<p><em>FST not yet available.</em></p>")
    for g1, g2 in fst_comparisons:
        img = img_fst_man.get((g1, g2), "")
        if img:
            parts.append(f"<h3>FST Manhattan — {g1} vs {g2}</h3>" + img)

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
