#!/usr/bin/env python3
"""
Population Genomics Report Generator
Reads Snakemake pipeline outputs and generates HTML report
Usage: python generate_report.py --results_dir results/ --metadata config/samples.csv --output report.html
"""

import argparse
import os
import json
import numpy as np
import pandas as pd
from pathlib import Path

# ============================================================================
# CONFIGURATION
# ============================================================================

class ReportConfig:
    """Report settings and thresholds for flagging"""
    def __init__(self, config_path=None):
        # Defaults
        self.min_mapping_rate = 0.50
        self.min_depth = 10
        self.depth_warn_threshold = 5
        self.het_sd_threshold = 2
        self.kinship_threshold = 0.2
        self.tstv_green = (1.8, 2.3)
        self.tstv_yellow = (1.5, 2.5)
        
        if config_path:
            self.load(config_path)
    
    def load(self, path):
        import yaml
        with open(path) as f:
            cfg = yaml.safe_load(f)
        for k, v in cfg.items():
            if hasattr(self, k):
                setattr(self, k, v)


# ============================================================================
# DATA LOADERS — read pipeline outputs
# ============================================================================

def load_metadata(path):
    """Load sample metadata CSV"""
    df = pd.read_csv(path)
    df = df.set_index("sample_id")
    return df

def load_mapping_stats(results_dir, samples):
    """Parse samtools flagstat outputs"""
    records = []
    for sample in samples:
        path = Path(results_dir) / "qc" / f"{sample}.dedup_metrics.txt"
        if not path.exists():
            continue
        with open(path) as f:
            lines = f.readlines()
        total = int(lines[0].split()[0])
        mapped = int(lines[4].split()[0])
        duplicates = int(lines[3].split()[0])
        rate = mapped / total if total > 0 else 0
        dup_rate = duplicates / total if total > 0 else 0
        records.append({
            "sample_id": sample,
            "total_reads": total,
            "mapped_reads": mapped,
            "mapping_rate": rate,
            "duplicate_rate": dup_rate,
        })
    return pd.DataFrame(records).set_index("sample_id")

def load_fastp_stats(results_dir, samples):
    """Parse fastp JSON outputs"""
    records = []
    for sample in samples:
        path = Path(results_dir) / "qc" / f"{sample}_fastp.json"
        if not path.exists():
            continue
        with open(path) as f:
            data = json.load(f)
        s = data["summary"]
        records.append({
            "sample_id": sample,
            "raw_reads": s["before_filtering"]["total_reads"],
            "clean_reads": s["after_filtering"]["total_reads"],
            "raw_bases": s["before_filtering"]["total_bases"],
            "clean_bases": s["after_filtering"]["total_bases"],
            "q30_rate": s["after_filtering"]["q30_rate"],
            "gc_content": s["after_filtering"]["gc_content"],
            "adapter_trimmed_pct": data.get("adapter_cutting", {}).get("adapter_trimmed_reads", 0) / max(s["before_filtering"]["total_reads"], 1),
        })
    return pd.DataFrame(records).set_index("sample_id")

def load_pca(results_dir):
    """Load PCAngsd covariance matrix and decompose"""
    cov_path = Path(results_dir) / "pca" / "pcangsd.cov"
    if not cov_path.exists():
        return None, None
    cov = np.loadtxt(cov_path)
    eigenvalues, eigenvectors = np.linalg.eigh(cov)
    # Sort descending
    idx = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]
    # Percent variance explained
    pct_var = eigenvalues / np.sum(np.abs(eigenvalues)) * 100
    return eigenvectors, pct_var

def load_admixture(results_dir, max_k=5):
    """Load NGSadmix Q matrices for each K"""
    admix = {}
    likelihoods = {}
    for k in range(2, max_k + 1):
        qpath = Path(results_dir) / "admixture" / f"K{k}.qopt"
        logpath = Path(results_dir) / "admixture" / f"K{k}.log"
        if qpath.exists():
            admix[k] = np.loadtxt(qpath)
        if logpath.exists():
            with open(logpath) as f:
                for line in f:
                    if "best like=" in line:
                        ll = float(line.strip().split("best like=")[1].split()[0])
                        likelihoods[k] = ll
    return admix, likelihoods

def load_fst(results_dir, pop_pairs):
    """Load global FST values"""
    records = []
    for p1, p2 in pop_pairs:
        path = Path(results_dir) / "fst" / f"{p1}_vs_{p2}.fst.global"
        if not path.exists():
            continue
        with open(path) as f:
            line = f.readline().strip()
        vals = line.split()
        if len(vals) >= 2:
            fst_unweight = float(vals[0])
            fst_weight = float(vals[1])
            records.append({
                "pop1": p1,
                "pop2": p2,
                "fst_unweighted": fst_unweight,
                "fst_weighted": fst_weight,
            })
    return pd.DataFrame(records)

def load_windowed_fst(results_dir, p1, p2):
    """Load windowed FST for Manhattan plot"""
    path = Path(results_dir) / "fst" / f"{p1}_vs_{p2}.fst.global.windows"
    if not path.exists():
        return None
    df = pd.read_csv(path, sep="\t", header=None,
                     names=["region", "chr", "midpoint", "n_sites", "fst"])
    return df

def load_thetas(results_dir, popname):
    """Load windowed diversity statistics"""
    path = Path(results_dir) / "diversity" / f"{popname}.thetas.idx.pestPG"
    if not path.exists():
        return None
    df = pd.read_csv(path, sep="\t")
    return df

def load_sfs(results_dir, popname):
    """Load site frequency spectrum"""
    path = Path(results_dir) / "sfs" / f"{popname}.sfs"
    if not path.exists():
        return None
    with open(path) as f:
        vals = [float(x) for x in f.readline().strip().split()]
    return np.array(vals)

def load_heterozygosity(results_dir, samples):
    """Load per-individual heterozygosity from realSFS output"""
    records = []
    for sample in samples:
        path = Path(results_dir) / "heterozygosity" / f"{sample}.het"
        if not path.exists():
            continue
        with open(path) as f:
            vals = [float(x) for x in f.readline().strip().split()]
        if len(vals) >= 2:
            hom = vals[0]
            het = vals[1]
            total = hom + het
            het_rate = het / total if total > 0 else 0
            records.append({
                "sample_id": sample,
                "homozygous_sites": int(hom),
                "heterozygous_sites": int(het),
                "heterozygosity": het_rate,
            })
    return pd.DataFrame(records).set_index("sample_id")

def load_ld_decay(results_dir):
    """Load LD decay data from ngsLD output"""
    path = Path(results_dir) / "ld" / "ld_decay.csv"
    if not path.exists():
        return None
    return pd.read_csv(path)


# ============================================================================
# QC FLAGGING
# ============================================================================

def flag_samples(mapping_stats, het_stats, config):
    """Generate sample flags based on QC thresholds"""
    flags = []
    
    for sample in mapping_stats.index:
        sample_flags = []
        
        # Mapping rate
        mr = mapping_stats.loc[sample, "mapping_rate"]
        if mr < config.min_mapping_rate:
            sample_flags.append(f"low_mapping_rate ({mr:.1%})")
        
        # Add depth flags when depth data available
        
        flags.append({
            "sample_id": sample,
            "flags": "; ".join(sample_flags) if sample_flags else "pass",
            "status": "flagged" if sample_flags else "pass",
        })
    
    # Heterozygosity outliers
    if het_stats is not None and len(het_stats) > 2:
        mean_het = het_stats["heterozygosity"].mean()
        sd_het = het_stats["heterozygosity"].std()
        for sample in het_stats.index:
            h = het_stats.loc[sample, "heterozygosity"]
            for f in flags:
                if f["sample_id"] == sample:
                    if h > mean_het + config.het_sd_threshold * sd_het:
                        f["flags"] += f"; excess_heterozygosity ({h:.4f})"
                        f["status"] = "flagged"
                    elif h < mean_het - config.het_sd_threshold * sd_het:
                        f["flags"] += f"; low_heterozygosity ({h:.4f})"
                        f["status"] = "flagged"
    
    return pd.DataFrame(flags).set_index("sample_id")


# ============================================================================
# EVANNO K SELECTION
# ============================================================================

def evanno_method(likelihoods):
    """Calculate delta K using Evanno method"""
    if len(likelihoods) < 3:
        return None
    
    ks = sorted(likelihoods.keys())
    lnk = [likelihoods[k] for k in ks]
    
    # First derivative
    lnk_prime = [lnk[i+1] - lnk[i] for i in range(len(lnk)-1)]
    
    # Second derivative
    lnk_dprime = [abs(lnk_prime[i+1] - lnk_prime[i]) for i in range(len(lnk_prime)-1)]
    
    # Delta K (using absolute value of second derivative)
    # Can't calculate SD from single runs, so just use |L''(K)|
    delta_k = {ks[i+1]: lnk_dprime[i] for i in range(len(lnk_dprime))}
    
    return delta_k


# ============================================================================
# PLOTTING FUNCTIONS
# ============================================================================

def setup_plotting():
    """Import and configure plotting libraries"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    
    # Publication style
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica"],
        "font.size": 11,
        "axes.linewidth": 1.2,
        "xtick.major.width": 1.2,
        "ytick.major.width": 1.2,
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
    })
    return plt, gridspec


def plot_mapping_rates(mapping_stats, metadata, config, output_path):
    """Barplot of per-sample mapping rates with flag threshold"""
    plt, _ = setup_plotting()
    fig, ax = plt.subplots(figsize=(12, 5))
    
    df = mapping_stats.join(metadata[["population"]])
    df = df.sort_values(["population", "mapping_rate"])
    
    colors = df["population"].map(lambda x: "#2196F3" if x == "florida" else "#FF9800")
    
    ax.bar(range(len(df)), df["mapping_rate"] * 100, color=colors, edgecolor="white", linewidth=0.5)
    ax.axhline(y=config.min_mapping_rate * 100, color="red", linestyle="--", linewidth=1, label=f"Flag threshold ({config.min_mapping_rate:.0%})")
    
    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(df.index, rotation=45, ha="right", fontsize=9)
    ax.set_ylabel("Mapping Rate (%)")
    ax.set_title("Per-Sample Mapping Rate")
    ax.legend()
    ax.set_ylim(0, 105)
    
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def plot_depth_distribution(mapping_stats, metadata, config, output_path):
    """Barplot of per-sample mean depth"""
    plt, _ = setup_plotting()
    fig, ax = plt.subplots(figsize=(12, 5))
    
    # Estimate depth from mapped bases / genome size
    genome_size = 334e6  # A. palmata
    df = mapping_stats.join(metadata[["population"]])
    if "mean_depth" not in df.columns:
        # Approximate from mapped reads * read length / genome size
        df["mean_depth"] = (df["mapped_reads"] * 150) / genome_size
    
    df = df.sort_values(["population", "mean_depth"])
    colors = df["population"].map(lambda x: "#2196F3" if x == "florida" else "#FF9800")
    
    ax.bar(range(len(df)), df["mean_depth"], color=colors, edgecolor="white", linewidth=0.5)
    ax.axhline(y=config.min_depth, color="red", linestyle="--", linewidth=1, label=f"Flag threshold ({config.min_depth}x)")
    
    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(df.index, rotation=45, ha="right", fontsize=9)
    ax.set_ylabel("Mean Depth (x)")
    ax.set_title("Per-Sample Sequencing Depth")
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def plot_pca(eigenvectors, pct_var, metadata, output_path, color_by="population"):
    """PCA scatter plot colored by metadata variable"""
    plt, _ = setup_plotting()
    fig, ax = plt.subplots(figsize=(8, 7))
    
    samples = metadata.index.tolist()
    groups = metadata[color_by].unique()
    
    palette = ["#2196F3", "#FF9800", "#4CAF50", "#E91E63", "#9C27B0", "#00BCD4", "#FF5722", "#607D8B"]
    color_map = {g: palette[i % len(palette)] for i, g in enumerate(groups)}
    
    for group in groups:
        mask = metadata[color_by] == group
        idx = [i for i, s in enumerate(samples) if mask[s]]
        ax.scatter(eigenvectors[idx, 0], eigenvectors[idx, 1],
                  c=color_map[group], label=group, s=80, edgecolors="black",
                  linewidth=0.5, alpha=0.85, zorder=3)
    
    ax.set_xlabel(f"PC1 ({pct_var[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({pct_var[1]:.1f}%)")
    ax.set_title("Principal Component Analysis")
    ax.legend(title=color_by.replace("_", " ").title(), frameon=True)
    ax.grid(True, alpha=0.3)
    
    # Add sample labels
    for i, sample in enumerate(samples):
        ax.annotate(sample, (eigenvectors[i, 0], eigenvectors[i, 1]),
                   fontsize=7, alpha=0.6, xytext=(5, 5),
                   textcoords="offset points")
    
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def plot_admixture(admix_dict, metadata, output_path, order_by="population"):
    """Stacked admixture barplots for all K values"""
    plt, gridspec = setup_plotting()
    
    k_values = sorted(admix_dict.keys())
    n_k = len(k_values)
    
    fig = plt.figure(figsize=(14, 3 * n_k))
    gs = gridspec.GridSpec(n_k, 1, hspace=0.3)
    
    samples = metadata.index.tolist()
    sort_idx = metadata.sort_values([order_by]).index
    sample_order = [samples.index(s) for s in sort_idx]
    
    palette = ["#2196F3", "#FF9800", "#4CAF50", "#E91E63", "#9C27B0", "#00BCD4", "#FF5722", "#607D8B"]
    
    for ki, k in enumerate(k_values):
        ax = fig.add_subplot(gs[ki])
        q = admix_dict[k][sample_order, :]
        
        bottom = np.zeros(len(samples))
        for j in range(k):
            ax.bar(range(len(samples)), q[:, j], bottom=bottom,
                  color=palette[j % len(palette)], width=1.0, edgecolor="none")
            bottom += q[:, j]
        
        ax.set_xlim(-0.5, len(samples) - 0.5)
        ax.set_ylim(0, 1)
        ax.set_ylabel(f"K={k}")
        
        if ki == n_k - 1:
            ax.set_xticks(range(len(samples)))
            ax.set_xticklabels([sort_idx[i] for i in range(len(samples))],
                             rotation=90, fontsize=7)
        else:
            ax.set_xticks([])
        
        # Add population separators
        pops = metadata.loc[sort_idx, order_by]
        prev = pops.iloc[0]
        for i, p in enumerate(pops):
            if p != prev:
                ax.axvline(x=i - 0.5, color="black", linewidth=1.5)
                prev = p
    
    fig.suptitle("Admixture Analysis", fontsize=14, y=1.02)
    plt.savefig(output_path)
    plt.close()


def plot_k_selection(likelihoods, delta_k, output_path):
    """Plot likelihood and delta K for K selection"""
    plt, _ = setup_plotting()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Log-likelihood vs K
    ks = sorted(likelihoods.keys())
    lls = [likelihoods[k] for k in ks]
    ax1.plot(ks, lls, "o-", color="#2196F3", linewidth=2, markersize=8)
    ax1.set_xlabel("K")
    ax1.set_ylabel("Log-Likelihood")
    ax1.set_title("Log-Likelihood vs K")
    ax1.set_xticks(ks)
    
    # Delta K
    if delta_k:
        dk_ks = sorted(delta_k.keys())
        dk_vals = [delta_k[k] for k in dk_ks]
        best_k = dk_ks[np.argmax(dk_vals)]
        ax2.plot(dk_ks, dk_vals, "o-", color="#FF9800", linewidth=2, markersize=8)
        ax2.set_xlabel("K")
        ax2.set_ylabel("ΔK")
        ax2.set_title(f"Evanno ΔK (Best K = {best_k})")
        ax2.set_xticks(dk_ks)
        ax2.axvline(x=best_k, color="red", linestyle="--", alpha=0.5)
    
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def plot_fst_manhattan(windowed_fst, output_path):
    """Manhattan plot of windowed FST"""
    plt, _ = setup_plotting()
    fig, ax = plt.subplots(figsize=(16, 5))
    
    if windowed_fst is None:
        return
    
    df = windowed_fst.copy()
    df = df[df["fst"] >= 0]  # Remove negative FST
    
    # Color alternating chromosomes
    chroms = df["chr"].unique()
    chrom_colors = {c: "#2196F3" if i % 2 == 0 else "#90CAF9" for i, c in enumerate(chroms)}
    
    # Calculate cumulative position
    chrom_offset = {}
    offset = 0
    for chrom in chroms:
        chrom_offset[chrom] = offset
        offset += df[df["chr"] == chrom]["midpoint"].max()
    
    df["cumpos"] = df.apply(lambda r: r["midpoint"] + chrom_offset[r["chr"]], axis=1)
    
    # Outlier threshold (top 1%)
    threshold = df["fst"].quantile(0.99)
    
    colors = df.apply(lambda r: "#E91E63" if r["fst"] >= threshold else chrom_colors[r["chr"]], axis=1)
    
    ax.scatter(df["cumpos"], df["fst"], c=colors, s=3, alpha=0.6)
    ax.axhline(y=threshold, color="red", linestyle="--", linewidth=1, alpha=0.7, label=f"Top 1% (FST > {threshold:.3f})")
    
    # Chromosome labels
    chrom_centers = {c: df[df["chr"] == c]["cumpos"].median() for c in chroms}
    ax.set_xticks(list(chrom_centers.values()))
    ax.set_xticklabels([str(c).replace("NC_", "chr") for c in chroms], rotation=45, fontsize=8)
    
    ax.set_ylabel("FST")
    ax.set_title("Windowed FST — Florida vs Panama")
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def plot_heterozygosity(het_stats, metadata, config, output_path):
    """Per-individual heterozygosity barplot"""
    plt, _ = setup_plotting()
    fig, ax = plt.subplots(figsize=(12, 5))
    
    df = het_stats.join(metadata[["population"]])
    df = df.sort_values(["population", "heterozygosity"])
    
    colors = df["population"].map(lambda x: "#2196F3" if x == "florida" else "#FF9800")
    
    mean_het = df["heterozygosity"].mean()
    sd_het = df["heterozygosity"].std()
    
    ax.bar(range(len(df)), df["heterozygosity"] * 100, color=colors, edgecolor="white", linewidth=0.5)
    ax.axhline(y=(mean_het + config.het_sd_threshold * sd_het) * 100, color="red", linestyle="--", linewidth=1, alpha=0.7, label=f"+{config.het_sd_threshold} SD")
    ax.axhline(y=(mean_het - config.het_sd_threshold * sd_het) * 100, color="red", linestyle="--", linewidth=1, alpha=0.7, label=f"-{config.het_sd_threshold} SD")
    ax.axhline(y=mean_het * 100, color="gray", linestyle="-", linewidth=1, alpha=0.5, label="Mean")
    
    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(df.index, rotation=45, ha="right", fontsize=9)
    ax.set_ylabel("Heterozygosity (%)")
    ax.set_title("Per-Individual Heterozygosity")
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def plot_sfs(sfs_dict, output_path):
    """Folded SFS per population"""
    plt, _ = setup_plotting()
    n_pops = len(sfs_dict)
    fig, axes = plt.subplots(1, n_pops, figsize=(6 * n_pops, 5))
    if n_pops == 1:
        axes = [axes]
    
    palette = ["#2196F3", "#FF9800", "#4CAF50"]
    
    for i, (popname, sfs) in enumerate(sfs_dict.items()):
        ax = axes[i]
        n = len(sfs)
        # Fold the SFS
        folded = sfs[:n//2+1].copy()
        for j in range(1, n//2):
            folded[j] += sfs[n-j]
        # Remove monomorphic
        folded = folded[1:]
        
        ax.bar(range(1, len(folded)+1), folded, color=palette[i % len(palette)], edgecolor="white")
        ax.set_xlabel("Minor Allele Count")
        ax.set_ylabel("Number of Sites")
        ax.set_title(f"Folded SFS — {popname.title()}")
    
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def plot_diversity_comparison(thetas_dict, output_path):
    """Compare pi, theta, Tajima's D between populations"""
    plt, _ = setup_plotting()
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    pops = list(thetas_dict.keys())
    palette = ["#2196F3", "#FF9800", "#4CAF50", "#E91E63"]
    
    metrics = [
        ("tP", "Nucleotide Diversity (π)", axes[0]),
        ("tW", "Watterson's θ", axes[1]),
        ("Tajima", "Tajima's D", axes[2]),
    ]
    
    for col, title, ax in metrics:
        means = []
        sems = []
        for popname in pops:
            df = thetas_dict[popname]
            if df is not None and col in df.columns:
                # Normalize by number of sites per window
                if "nSites" in df.columns:
                    vals = df[col] / df["nSites"]
                else:
                    vals = df[col]
                means.append(vals.mean())
                sems.append(vals.sem())
            else:
                means.append(0)
                sems.append(0)
        
        bars = ax.bar(range(len(pops)), means, yerr=sems, capsize=5,
                     color=[palette[i % len(palette)] for i in range(len(pops))],
                     edgecolor="white", linewidth=0.5)
        ax.set_xticks(range(len(pops)))
        ax.set_xticklabels([p.title() for p in pops])
        ax.set_title(title)
        ax.set_ylabel(title)
    
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def plot_kinship_heatmap(kinship_df, output_path):
    """Pairwise kinship heatmap"""
    plt, _ = setup_plotting()
    
    samples = sorted(set(kinship_df["sample_1"].tolist() + kinship_df["sample_2"].tolist()))
    n = len(samples)
    mat = np.zeros((n, n))
    
    sample_idx = {s: i for i, s in enumerate(samples)}
    for _, row in kinship_df.iterrows():
        i = sample_idx[row["sample_1"]]
        j = sample_idx[row["sample_2"]]
        mat[i, j] = row["kinship"]
        mat[j, i] = row["kinship"]
    
    fig, ax = plt.subplots(figsize=(10, 9))
    im = ax.imshow(mat, cmap="YlOrRd", vmin=0, vmax=0.5)
    
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(samples, rotation=90, fontsize=8)
    ax.set_yticklabels(samples, fontsize=8)
    
    plt.colorbar(im, ax=ax, label="Kinship Coefficient", shrink=0.8)
    ax.set_title("Pairwise Kinship")
    
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


# ============================================================================
# CSV EXPORT — all intermediate deliverables
# ============================================================================

def export_csvs(results_dir, output_dir, metadata, mapping_stats, fastp_stats,
                het_stats, fst_global, flags):
    """Export all deliverable CSVs"""
    os.makedirs(output_dir, exist_ok=True)
    
    # Master sample table
    master = metadata.copy()
    if mapping_stats is not None:
        master = master.join(mapping_stats, how="left")
    if fastp_stats is not None:
        master = master.join(fastp_stats, how="left")
    if het_stats is not None:
        master = master.join(het_stats, how="left")
    if flags is not None:
        master = master.join(flags, how="left")
    master.to_csv(Path(output_dir) / "sample_metadata_with_results.csv")
    
    # Individual CSVs
    if mapping_stats is not None:
        mapping_stats.to_csv(Path(output_dir) / "qc" / "mapping_stats.csv")
    if fastp_stats is not None:
        fastp_stats.to_csv(Path(output_dir) / "qc" / "sequencing_stats.csv")
    if het_stats is not None:
        het_stats.to_csv(Path(output_dir) / "individuals" / "heterozygosity.csv")
    if fst_global is not None:
        fst_global.to_csv(Path(output_dir) / "fst" / "global_fst.csv", index=False)
    if flags is not None:
        flags.to_csv(Path(output_dir) / "individuals" / "sample_flags.csv")


# ============================================================================
# AUTO-GENERATED INTERPRETATION TEXT
# ============================================================================

def generate_executive_summary(metadata, mapping_stats, het_stats, fst_global,
                                admix_likelihoods, delta_k, flags):
    """Auto-generate executive summary text"""
    n_samples = len(metadata)
    pops = metadata["population"].unique()
    n_pops = len(pops)
    pop_counts = metadata["population"].value_counts()
    
    summary = []
    summary.append(f"We analyzed {n_samples} samples across {n_pops} populations "
                   f"({', '.join(f'{pop_counts[p]} {p}' for p in pops)}).")
    
    if mapping_stats is not None:
        mean_mr = mapping_stats["mapping_rate"].mean()
        summary.append(f"Mean mapping rate was {mean_mr:.1%}.")
    
    if flags is not None:
        n_flagged = (flags["status"] == "flagged").sum()
        if n_flagged > 0:
            summary.append(f"{n_flagged} sample(s) were flagged for quality concerns.")
        else:
            summary.append("All samples passed quality checks.")
    
    if delta_k:
        best_k = max(delta_k, key=delta_k.get)
        summary.append(f"Population structure analysis supports K={best_k} ancestral populations.")
    
    if fst_global is not None and len(fst_global) > 0:
        for _, row in fst_global.iterrows():
            summary.append(f"Weighted FST between {row['pop1']} and {row['pop2']} "
                          f"is {row['fst_weighted']:.4f}.")
    
    return " ".join(summary)


def generate_methods_text(config, metadata):
    """Auto-generate publication-ready methods paragraph"""
    n = len(metadata)
    pops = metadata["population"].unique()
    
    methods = f"""Whole genome sequencing data from {n} samples were processed using 
a custom Snakemake pipeline. Raw reads were quality-trimmed using fastp with adapter 
auto-detection. Trimmed reads were mapped to the Acropora palmata reference genome 
(GCF_964030605.1, Wellcome Sanger Institute) using BWA-MEM. Duplicate reads were 
marked using samtools markdup. 

Genotype likelihoods were estimated using ANGSD (GL model {config.get('gl_model', 1)}) 
with minimum mapping quality {config.get('minmapq', 20)}, minimum base quality 
{config.get('minq', 20)}, minimum minor allele frequency {config.get('min_maf', 0.10)}, 
SNP p-value threshold {config.get('snp_pval', '1e-6')}, and requiring data in at least 
{int(config.get('min_ind_frac', 0.8) * 100)}% of individuals per site. Maximum per-site 
depth was set to the empirical mean plus two standard deviations to exclude putative 
paralogous regions.

Population structure was assessed using PCAngsd for principal component analysis and 
NGSadmix for admixture analysis (K=2 through K={config.get('max_k', 5)}). The optimal 
number of ancestral populations was determined using the Evanno delta-K method. Pairwise 
relatedness was estimated using NgsRelate.

Site frequency spectra were estimated per population using ANGSD and realSFS. Nucleotide 
diversity (pi), Watterson's theta, and Tajima's D were calculated in 50 kb sliding 
windows with 10 kb steps. Pairwise FST between populations was estimated using the 
2D-SFS approach in realSFS with 50 kb sliding windows."""
    
    return methods


# ============================================================================
# HTML REPORT GENERATION
# ============================================================================

def generate_html_report(title, sections, output_path):
    """Assemble sections into a single HTML report"""
    
    html = f"""<!DOCTYPE html>
<html>
<head>
    <title>{title}</title>
    <style>
        body {{
            font-family: Arial, Helvetica, sans-serif;
            max-width: 1200px;
            margin: 0 auto;
            padding: 40px;
            color: #333;
            line-height: 1.6;
        }}
        h1 {{
            color: #1565C0;
            border-bottom: 3px solid #1565C0;
            padding-bottom: 10px;
        }}
        h2 {{
            color: #1976D2;
            border-bottom: 1px solid #ddd;
            padding-bottom: 8px;
            margin-top: 40px;
        }}
        h3 {{
            color: #1E88E5;
        }}
        .executive-summary {{
            background: #E3F2FD;
            border-left: 4px solid #1565C0;
            padding: 20px;
            margin: 20px 0;
            border-radius: 0 8px 8px 0;
            font-size: 1.05em;
        }}
        .flag-warning {{
            background: #FFF3E0;
            border-left: 4px solid #FF9800;
            padding: 15px;
            margin: 15px 0;
            border-radius: 0 8px 8px 0;
        }}
        .flag-ok {{
            background: #E8F5E9;
            border-left: 4px solid #4CAF50;
            padding: 15px;
            margin: 15px 0;
            border-radius: 0 8px 8px 0;
        }}
        table {{
            border-collapse: collapse;
            width: 100%;
            margin: 15px 0;
        }}
        th, td {{
            border: 1px solid #ddd;
            padding: 10px;
            text-align: left;
        }}
        th {{
            background: #1976D2;
            color: white;
        }}
        tr:nth-child(even) {{
            background: #f9f9f9;
        }}
        tr:hover {{
            background: #E3F2FD;
        }}
        img {{
            max-width: 100%;
            margin: 15px 0;
            border: 1px solid #ddd;
            border-radius: 4px;
        }}
        .metric-box {{
            display: inline-block;
            background: #f5f5f5;
            border-radius: 8px;
            padding: 15px 25px;
            margin: 10px;
            text-align: center;
        }}
        .metric-value {{
            font-size: 2em;
            font-weight: bold;
            color: #1565C0;
        }}
        .metric-label {{
            font-size: 0.9em;
            color: #666;
        }}
        .traffic-green {{ color: #4CAF50; font-weight: bold; }}
        .traffic-yellow {{ color: #FF9800; font-weight: bold; }}
        .traffic-red {{ color: #E91E63; font-weight: bold; }}
        .footer {{
            margin-top: 50px;
            padding-top: 20px;
            border-top: 1px solid #ddd;
            color: #999;
            font-size: 0.85em;
        }}
    </style>
</head>
<body>
    <h1>{title}</h1>
"""
    
    for section in sections:
        html += section
    
    html += """
    <div class="footer">
        Generated by ECOS Population Genomics Pipeline<br>
        School of Environmental, Coastal, and Ocean Sustainability<br>
        Florida Atlantic University
    </div>
</body>
</html>"""
    
    with open(output_path, "w") as f:
        f.write(html)


# ============================================================================
# MAIN
# ============================================================================

def main():
    parser = argparse.ArgumentParser(description="Generate population genomics report")
    parser.add_argument("--results_dir", required=True, help="Path to pipeline results/")
    parser.add_argument("--metadata", required=True, help="Path to sample metadata CSV")
    parser.add_argument("--config", default=None, help="Path to config.yaml")
    parser.add_argument("--output", default="report.html", help="Output HTML report path")
    parser.add_argument("--output_dir", default="deliverables", help="Directory for CSV deliverables")
    args = parser.parse_args()
    
    # Load config
    report_config = ReportConfig(args.config)
    
    # Load pipeline config
    pipeline_config = {}
    config_path = Path(args.results_dir).parent / "config" / "config.yaml"
    if config_path.exists():
        import yaml
        with open(config_path) as f:
            pipeline_config = yaml.safe_load(f)
    
    # Load metadata
    metadata = load_metadata(args.metadata)
    samples = metadata.index.tolist()
    pops = metadata["population"].unique().tolist()
    
    # Create figure directory
    fig_dir = Path(args.output_dir) / "figures"
    os.makedirs(fig_dir, exist_ok=True)
    os.makedirs(Path(args.output_dir) / "qc", exist_ok=True)
    os.makedirs(Path(args.output_dir) / "fst", exist_ok=True)
    os.makedirs(Path(args.output_dir) / "individuals", exist_ok=True)
    os.makedirs(Path(args.output_dir) / "structure", exist_ok=True)
    os.makedirs(Path(args.output_dir) / "diversity", exist_ok=True)
    os.makedirs(Path(args.output_dir) / "relatedness", exist_ok=True)
    
    # Load all data
    print("Loading pipeline outputs...")
    mapping_stats = load_mapping_stats(args.results_dir, samples)
    fastp_stats = load_fastp_stats(args.results_dir, samples)
    eigenvectors, pct_var = load_pca(args.results_dir)
    admix, likelihoods = load_admixture(args.results_dir, pipeline_config.get("max_k", 5))
    het_stats = load_heterozygosity(args.results_dir, samples)
    
    from itertools import combinations
    pop_pairs = list(combinations(sorted(pops), 2))
    fst_global = load_fst(args.results_dir, pop_pairs)
    
    thetas_dict = {p: load_thetas(args.results_dir, p) for p in pops}
    sfs_dict = {p: load_sfs(args.results_dir, p) for p in pops}
    
    # QC flagging
    flags = flag_samples(mapping_stats, het_stats, report_config)
    
    # K selection
    delta_k = evanno_method(likelihoods) if likelihoods else None
    
    # Generate plots
    print("Generating figures...")
    
    if len(mapping_stats) > 0:
        plot_mapping_rates(mapping_stats, metadata, report_config, fig_dir / "mapping_rates.png")
        plot_depth_distribution(mapping_stats, metadata, report_config, fig_dir / "depth_distribution.png")
    
    if eigenvectors is not None:
        plot_pca(eigenvectors, pct_var, metadata, fig_dir / "pca_population.png", "population")
        # Additional PCA colorings for each grouping variable
        for col in metadata.columns:
            if col not in ["sra_accession", "fastq_path"] and metadata[col].nunique() > 1:
                plot_pca(eigenvectors, pct_var, metadata, fig_dir / f"pca_{col}.png", col)
    
    if admix:
        plot_admixture(admix, metadata, fig_dir / "admixture.png")
    
    if likelihoods and delta_k:
        plot_k_selection(likelihoods, delta_k, fig_dir / "k_selection.png")
    
    if het_stats is not None and len(het_stats) > 0:
        plot_heterozygosity(het_stats, metadata, report_config, fig_dir / "heterozygosity.png")
    
    for p1, p2 in pop_pairs:
        wfst = load_windowed_fst(args.results_dir, p1, p2)
        if wfst is not None:
            plot_fst_manhattan(wfst, fig_dir / f"fst_manhattan_{p1}_vs_{p2}.png")
    
    if any(v is not None for v in sfs_dict.values()):
        plot_sfs(sfs_dict, fig_dir / "sfs.png")
    
    if any(v is not None for v in thetas_dict.values()):
        plot_diversity_comparison(thetas_dict, fig_dir / "diversity_comparison.png")
    
    # Generate text
    exec_summary = generate_executive_summary(metadata, mapping_stats, het_stats,
                                               fst_global, likelihoods, delta_k, flags)
    methods_text = generate_methods_text(pipeline_config, metadata)
    
    # Build HTML sections
    print("Building report...")
    sections = []
    
    # Executive summary
    sections.append(f'<div class="executive-summary"><strong>Summary:</strong> {exec_summary}</div>')
    
    # Sample overview
    sections.append("<h2>1. Sample Overview</h2>")
    sections.append(metadata.to_html(classes=""))
    
    # Sequencing QC
    sections.append("<h2>2. Sequencing Quality</h2>")
    if fastp_stats is not None and len(fastp_stats) > 0:
        sections.append(fastp_stats.to_html(classes="", float_format="%.2f"))
    
    # Mapping
    sections.append("<h2>3. Mapping & Coverage</h2>")
    if len(mapping_stats) > 0:
        sections.append('<img src="figures/mapping_rates.png">')
        sections.append('<img src="figures/depth_distribution.png">')
        sections.append(mapping_stats.to_html(classes="", float_format="%.4f"))
    
    # Flags
    sections.append("<h2>4. Sample QC Flags</h2>")
    n_flagged = (flags["status"] == "flagged").sum() if flags is not None else 0
    if n_flagged > 0:
        sections.append(f'<div class="flag-warning">{n_flagged} sample(s) flagged. Review recommended before interpretation.</div>')
    else:
        sections.append('<div class="flag-ok">All samples passed quality checks.</div>')
    if flags is not None:
        sections.append(flags.to_html(classes=""))
    
    # Population structure
    sections.append("<h2>5. Population Structure</h2>")
    sections.append("<h3>5a. K Selection</h3>")
    if likelihoods and delta_k:
        best_k = max(delta_k, key=delta_k.get)
        sections.append(f"<p>Best K = <strong>{best_k}</strong> based on Evanno ΔK method.</p>")
        sections.append('<img src="figures/k_selection.png">')
    
    sections.append("<h3>5b. Admixture</h3>")
    if admix:
        sections.append('<img src="figures/admixture.png">')
    
    sections.append("<h3>5c. PCA</h3>")
    if eigenvectors is not None:
        sections.append('<img src="figures/pca_population.png">')
    
    # Diversity
    sections.append("<h2>6. Genetic Diversity</h2>")
    sections.append('<img src="figures/sfs.png">')
    sections.append('<img src="figures/diversity_comparison.png">')
    
    # FST
    sections.append("<h2>7. Population Differentiation</h2>")
    if fst_global is not None and len(fst_global) > 0:
        sections.append(fst_global.to_html(classes="", index=False, float_format="%.4f"))
        for p1, p2 in pop_pairs:
            manhattan_path = fig_dir / f"fst_manhattan_{p1}_vs_{p2}.png"
            if manhattan_path.exists():
                sections.append(f'<img src="figures/fst_manhattan_{p1}_vs_{p2}.png">')
    
    # Heterozygosity
    sections.append("<h2>8. Individual Statistics</h2>")
    if het_stats is not None and len(het_stats) > 0:
        sections.append('<img src="figures/heterozygosity.png">')
        sections.append(het_stats.to_html(classes="", float_format="%.6f"))
    
    # Methods
    sections.append("<h2>9. Methods</h2>")
    sections.append(f"<p>{methods_text}</p>")
    
    # Generate report
    title = "Population Genomics Report"
    generate_html_report(title, sections, args.output)
    
    # Export CSVs
    print("Exporting deliverable CSVs...")
    export_csvs(args.results_dir, args.output_dir, metadata, mapping_stats,
                fastp_stats, het_stats, fst_global, flags)
    
    # Save methods
    with open(Path(args.output_dir) / "methods" / "methods.txt", "w") as f:
        f.write(methods_text)
    
    # Save PCA coordinates
    if eigenvectors is not None:
        pca_df = pd.DataFrame(eigenvectors[:, :10],
                             index=samples,
                             columns=[f"PC{i+1}" for i in range(10)])
        pca_df = pca_df.join(metadata)
        pca_df.to_csv(Path(args.output_dir) / "structure" / "pca_coordinates.csv")
    
    # Save admixture proportions
    for k, q in admix.items():
        q_df = pd.DataFrame(q, index=samples,
                           columns=[f"Cluster_{i+1}" for i in range(k)])
        q_df = q_df.join(metadata)
        q_df.to_csv(Path(args.output_dir) / "structure" / f"admixture_K{k}.csv")
    
    print(f"Report saved to {args.output}")
    print(f"Deliverables saved to {args.output_dir}/")


if __name__ == "__main__":
    main()
