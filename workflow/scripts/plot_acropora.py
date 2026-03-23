import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import os

WORKDIR = "/work/vollmer/acropora_genomics"
OUTDIR  = f"{WORKDIR}/results/plots"
os.makedirs(OUTDIR, exist_ok=True)

samples = pd.read_csv(f"{WORKDIR}/config/samples_RR.csv")
# Ground truth sample list = 253 in lineage_assignments.txt
la = pd.read_csv(f"{WORKDIR}/results/admixture/lineage_assignments.txt", sep="\t",
                 comment="#", header=None,
                 names=["sample_id","lineage","Q_A","Q_B"])
meta = la.merge(samples[["sample_id","species","population","region"]], on="sample_id", how="left")
# Row order in Q matrices matches row order in lineage_assignments.txt
meta = meta.reset_index(drop=True)

# ── colour palette ──────────────────────────────────────────────────────────
POP_COLORS = {
    "FL_RR":"#e6194b","FL_CRF":"#f58231","FL_MOTE":"#ffe119",
    "FL_FWC":"#bfef45","FL_upper":"#3cb44b","FL_lower":"#42d4f4",
    "FL_middle":"#4363d8","FL_biscayne":"#911eb4",
    "BON":"#f032e6","PA":"#a9a9a9",
}
SP_COLORS  = {"Apalmata":"#2196F3","Acervicornis":"#FF5722"}
LIN_COLORS = {"lineageA":"#2196F3","lineageB":"#FF5722"}

# ── 1. PCA ───────────────────────────────────────────────────────────────────
cov = np.loadtxt(f"{WORKDIR}/results/pca/pcangsd.cov")
vals, vecs = np.linalg.eigh(cov)
idx = np.argsort(vals)[::-1]
vals, vecs = vals[idx], vecs[:, idx]
pve = vals / vals.sum() * 100

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
for ax, (pc1, pc2) in zip(axes, [(0,1),(0,2)]):
    for sp, grp in meta.groupby("species"):
        ax.scatter(vecs[grp.index, pc1], vecs[grp.index, pc2],
                   c=SP_COLORS[sp], label=sp.replace("A","A. "), alpha=0.7, s=25, edgecolors="none")
    ax.set_xlabel(f"PC{pc1+1} ({pve[pc1]:.1f}%)")
    ax.set_ylabel(f"PC{pc2+1} ({pve[pc2]:.1f}%)")
    ax.axhline(0, c="grey", lw=0.4); ax.axvline(0, c="grey", lw=0.4)
axes[0].legend(frameon=False, fontsize=9)
fig.suptitle("PCAngsd — 253 unrelated samples", fontweight="bold")
plt.tight_layout()
fig.savefig(f"{OUTDIR}/pca_species.png", dpi=150, bbox_inches="tight")
plt.close()

# PCA coloured by population
fig, ax = plt.subplots(figsize=(7, 6))
for pop, grp in meta.groupby("population"):
    c = POP_COLORS.get(pop, "#aaaaaa")
    ax.scatter(vecs[grp.index, 0], vecs[grp.index, 1],
               c=c, label=pop, alpha=0.75, s=25, edgecolors="none")
ax.set_xlabel(f"PC1 ({pve[0]:.1f}%)"); ax.set_ylabel(f"PC2 ({pve[1]:.1f}%)")
ax.axhline(0, c="grey", lw=0.4); ax.axvline(0, c="grey", lw=0.4)
ax.legend(fontsize=7, ncol=2, frameon=False, bbox_to_anchor=(1,1), loc="upper left")
fig.suptitle("PCAngsd — coloured by population", fontweight="bold")
plt.tight_layout()
fig.savefig(f"{OUTDIR}/pca_population.png", dpi=150, bbox_inches="tight")
plt.close()
print("PCA done")

# ── 2. Admixture bar plots K2-K5 ────────────────────────────────────────────
# Sort samples: lineageA first, then by population
meta_sorted = meta.sort_values(["lineage","population","sample_id"]).reset_index(drop=True)

fig, axes = plt.subplots(4, 1, figsize=(14, 8), sharex=True)
for i, K in enumerate(range(2, 6)):
    Q = np.loadtxt(f"{WORKDIR}/results/admixture/pcangsd_K{K}.Q")
    Q_sorted = Q[meta_sorted.index]
    ax = axes[i]
    bottom = np.zeros(len(meta_sorted))
    colors = plt.cm.tab10(np.linspace(0, 0.9, K))
    for k in range(K):
        ax.bar(range(len(meta_sorted)), Q_sorted[:, k], bottom=bottom,
               color=colors[k], width=1.0, linewidth=0)
        bottom += Q_sorted[:, k]
    ax.set_ylabel(f"K={K}", fontsize=9)
    ax.set_ylim(0, 1); ax.set_yticks([0, 0.5, 1])
    ax.tick_params(left=True, bottom=False)

# Population tick labels
pops = meta_sorted["population"].values
boundaries = [0] + list(np.where(np.array(pops[:-1]) != np.array(pops[1:]))[0] + 1) + [len(pops)]
ax = axes[-1]
ax.set_xlim(-0.5, len(meta_sorted) - 0.5)
ax.set_xticks([(boundaries[i]+boundaries[i+1])//2 for i in range(len(boundaries)-1)])
ax.set_xticklabels([pops[b] for b in boundaries[:-1]], rotation=45, ha="right", fontsize=8)

# Lineage separator
lin = meta_sorted["lineage"].values
sep = np.where(np.array(lin[:-1]) != np.array(lin[1:]))[0][0]
for ax in axes:
    ax.axvline(sep + 0.5, color="black", lw=1.5, linestyle="--")

fig.suptitle("PCAngsd admixture — K2 to K5", fontweight="bold")
plt.tight_layout()
fig.savefig(f"{OUTDIR}/admixture_K2_K5.png", dpi=150, bbox_inches="tight")
plt.close()
print("Admixture done")

# ── 3. Diversity (Watterson's theta per site) ────────────────────────────────
POP_ORDER = ["FL_RR","FL_CRF","FL_MOTE","FL_FWC","FL_upper","FL_lower",
             "FL_middle","FL_biscayne","PA","BON",
             "lineageA_FL","lineageA_BON","lineageB_FL","lineageB_PA","lineageB_BON"]
div_rows = []
for pop in POP_ORDER:
    f = f"{WORKDIR}/results/diversity/{pop}.thetas.idx.pestPG"
    if not os.path.exists(f): continue
    df = pd.read_csv(f, sep="\t", comment="#",
                     names=["win","Chr","WinCenter","tW","tP","tF","tH","tL",
                            "Tajima","fuf","fud","fayh","zeng","nSites"])
    nsites = df["nSites"].sum()
    tW = df["tW"].sum() / nsites
    tP = df["tP"].sum() / nsites
    tajD = (df["Tajima"] * df["nSites"]).sum() / nsites
    div_rows.append({"pop":pop,"tW":tW,"tP":tP,"tajD":tajD})
div = pd.DataFrame(div_rows)

fig, axes = plt.subplots(1, 3, figsize=(14, 5))
for ax, col, label in zip(axes,
    ["tW","tP","tajD"],
    ["Watterson's θ per site","Nucleotide diversity π per site","Tajima's D"]):
    colors = ["#2196F3" if "lineageA" in p or p in ["FL_RR"] else
              "#FF5722" if "lineageB" in p else
              "#aaaaaa" for p in div["pop"]]
    ax.barh(range(len(div)), div[col], color=colors)
    ax.set_yticks(range(len(div))); ax.set_yticklabels(div["pop"], fontsize=9)
    ax.set_xlabel(label, fontsize=9)
    ax.axvline(0, c="grey", lw=0.5)
    ax.invert_yaxis()

patches = [mpatches.Patch(color="#2196F3",label="Apal / lineageA"),
           mpatches.Patch(color="#FF5722",label="Acer / lineageB"),
           mpatches.Patch(color="#aaaaaa",label="Population")]
axes[0].legend(handles=patches, fontsize=8, frameon=False)
fig.suptitle("Genome-wide diversity (θ, π, Tajima's D)", fontweight="bold")
plt.tight_layout()
fig.savefig(f"{OUTDIR}/diversity_summary.png", dpi=150, bbox_inches="tight")
plt.close()
print("Diversity done")

# ── 4. FST heatmap (completed comparisons) ───────────────────────────────────
fst_data = {
    ("lineageA_FL","lineageA_BON"): (0.0077, 0.0533),
    ("lineageB_PA","lineageB_BON"): (0.0086, 0.0933),
    ("lineageB_BON","lineageA_BON"): (0.0194, 0.4556),
}
# also try to load any .fst.global files present
import glob
for f in glob.glob(f"{WORKDIR}/results/fst/*.fst.global"):
    name = os.path.basename(f).replace(".fst.global","")
    p1, p2 = name.replace("_vs_","||").split("||") if "_vs_" in name else (None,None)
    if p1:
        vals_line = open(f).read().split()
        if len(vals_line) >= 2:
            fst_data[(p1, p2)] = (float(vals_line[0]), float(vals_line[1]))

pops_fst = sorted(set([p for pair in fst_data for p in pair]))
n = len(pops_fst); pidx = {p:i for i,p in enumerate(pops_fst)}
mat_uw = np.full((n,n), np.nan); mat_w = np.full((n,n), np.nan)
for (p1,p2),(uw,w) in fst_data.items():
    i,j = pidx.get(p1), pidx.get(p2)
    if i is not None and j is not None:
        mat_uw[i,j] = mat_uw[j,i] = uw
        mat_w[i,j]  = mat_w[j,i]  = w

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
for ax, mat, title in zip(axes, [mat_uw, mat_w], ["Fst (unweighted)","Fst (weighted)"]):
    im = ax.imshow(mat, cmap="YlOrRd", vmin=0, aspect="auto")
    ax.set_xticks(range(n)); ax.set_xticklabels(pops_fst, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(n)); ax.set_yticklabels(pops_fst, fontsize=8)
    plt.colorbar(im, ax=ax, shrink=0.8)
    ax.set_title(title)
    for i in range(n):
        for j in range(n):
            if not np.isnan(mat[i,j]):
                ax.text(j, i, f"{mat[i,j]:.3f}", ha="center", va="center", fontsize=7)
fig.suptitle("Pairwise FST (completed comparisons)", fontweight="bold")
plt.tight_layout()
fig.savefig(f"{OUTDIR}/fst_heatmap.png", dpi=150, bbox_inches="tight")
plt.close()
print("FST done")


# ── 5. Per-chromosome pi, theta, Tajima's D ──────────────────────────────────
CHROM_POPS = ["lineageA_FL", "lineageB_FL", "lineageB_PA"]
CHROM_LABELS = {
    "NC_133882.1":"1","NC_133883.1":"2","NC_133884.1":"3","NC_133885.1":"4",
    "NC_133886.1":"5","NC_133887.1":"6","NC_133888.1":"7","NC_133889.1":"8",
    "NC_133890.1":"9","NC_133891.1":"10","NC_133892.1":"11","NC_133893.1":"12",
    "NC_133894.1":"13","NC_133895.1":"14",
}
CHROM_COLORS = {"lineageA_FL":"#2196F3","lineageB_FL":"#FF5722","lineageB_PA":"#FF9800"}
CHR14 = "NC_133895.1"

fig, axes = plt.subplots(3, 3, figsize=(15, 12), sharey="row")
for col, pop in enumerate(CHROM_POPS):
    f = f"{WORKDIR}/results/diversity/{pop}.thetas.idx.pestPG"
    df = pd.read_csv(f, sep="\t", comment="#",
                     names=["win","Chr","WinCenter","tW","tP","tF","tH","tL",
                            "Tajima","fuf","fud","fayh","zeng","nSites"])
    df = df.sort_values("Chr")
    df["tW_site"] = df["tW"] / df["nSites"]
    df["tP_site"] = df["tP"] / df["nSites"]
    df["label"]   = df["Chr"].map(CHROM_LABELS)
    colors = ["#d32f2f" if c == CHR14 else CHROM_COLORS[pop] for c in df["Chr"]]

    for row, (col_name, ylabel) in enumerate([
            ("tW_site", "Watterson's theta per site"),
            ("tP_site", "pi per site"),
            ("Tajima",  "Tajima's D")]):
        ax = axes[row][col]
        ax.bar(range(len(df)), df[col_name], color=colors, edgecolor="none", width=0.7)
        ax.set_xticks(range(len(df)))
        ax.set_xticklabels(df["label"], fontsize=8)
        ax.set_xlabel("Chromosome", fontsize=8)
        if col == 0:
            ax.set_ylabel(ylabel, fontsize=9)
        if row == 0:
            ax.set_title(pop.replace("lineage", "Lineage ").replace("_", " "), fontsize=10, fontweight="bold")
        if col_name == "Tajima":
            ax.axhline(0, c="grey", lw=0.8, ls="--")
        if CHR14 in list(df["Chr"]):
            idx14 = list(df["Chr"]).index(CHR14)
            val = df[col_name].iloc[idx14]
            offset = val * 0.08 if abs(val) > 0.0001 else 0.0001
            ax.annotate("14", xy=(idx14, val), xytext=(idx14, val + offset),
                        ha="center", fontsize=7, color="#d32f2f", fontweight="bold")

red_patch = mpatches.Patch(color="#d32f2f", label="Chr 14 (outlier)")
axes[0][0].legend(handles=[red_patch], fontsize=8, frameon=False)
fig.suptitle("Per-chromosome diversity -- theta, pi, Tajima's D\n(Chr 14 highlighted)", fontweight="bold")
plt.tight_layout()
fig.savefig(f"{OUTDIR}/chrom_diversity.png", dpi=150, bbox_inches="tight")
plt.close()
print("Per-chromosome diversity done")

# ── 6. Windowed pi / theta across chromosomes (if available) ─────────────────
windowed_files = [f"{WORKDIR}/results/diversity/{p}.windowed.pestPG" for p in CHROM_POPS]
if all(os.path.exists(f) for f in windowed_files):
    CHROM_ORDER = ["NC_133882.1","NC_133883.1","NC_133884.1","NC_133885.1",
                   "NC_133886.1","NC_133887.1","NC_133888.1","NC_133889.1",
                   "NC_133890.1","NC_133891.1","NC_133892.1","NC_133893.1",
                   "NC_133894.1","NC_133895.1"]

    # Build genome-wide x-offsets from first population
    df0 = pd.read_csv(windowed_files[0], sep="\t", comment="#",
                      names=["win","Chr","WinCenter","tW","tP","tF","tH","tL",
                             "Tajima","fuf","fud","fayh","zeng","nSites"])
    chrom_offsets = {}
    offset = 0
    for chrom in CHROM_ORDER:
        sub = df0[df0["Chr"] == chrom]
        if len(sub):
            chrom_offsets[chrom] = offset
            offset += int(sub["WinCenter"].max()) + 1_000_000

    for stat, ylabel, fname in [
            ("tP_site", "pi per site (50kb windows)", "windowed_pi.png"),
            ("Tajima",  "Tajima's D (50kb windows)",  "windowed_tajima.png")]:
        fig, axes2 = plt.subplots(len(CHROM_POPS), 1, figsize=(18, 10), sharey=True)
        for ax, pop in zip(axes2, CHROM_POPS):
            f = f"{WORKDIR}/results/diversity/{pop}.windowed.pestPG"
            df = pd.read_csv(f, sep="\t", comment="#",
                             names=["win","Chr","WinCenter","tW","tP","tF","tH","tL",
                                    "Tajima","fuf","fud","fayh","zeng","nSites"])
            df = df[df["nSites"] >= 5000].copy()
            if stat == "tP_site":
                df["tP_site"] = df["tP"] / df["nSites"]
            df["x"] = df.apply(
                lambda r: r["WinCenter"] + chrom_offsets.get(r["Chr"], 0), axis=1)
            c = CHROM_COLORS[pop]
            for i, chrom in enumerate(CHROM_ORDER):
                sub = df[df["Chr"] == chrom]
                color = "#d32f2f" if chrom == CHR14 else c
                alpha = 0.7 if chrom == CHR14 else (0.5 if i % 2 == 0 else 0.3)
                ax.scatter(sub["x"], sub[stat], s=1.5, c=color, alpha=alpha, edgecolors="none")
                if chrom in chrom_offsets:
                    ax.axvline(chrom_offsets[chrom], c="lightgrey", lw=0.5)
            if stat == "Tajima":
                ax.axhline(0, c="grey", lw=0.8, ls="--")
            ax.set_ylabel(pop.replace("lineage","Lin ").replace("_"," "), fontsize=9)
        # chromosome labels on bottom axis
        for chrom in CHROM_ORDER:
            if chrom in chrom_offsets:
                sub = df[df["Chr"] == chrom]
                if len(sub):
                    mid = sub["x"].mean()
                    col = "#d32f2f" if chrom == CHR14 else "black"
                    axes2[-1].text(mid, axes2[-1].get_ylim()[0],
                                   CHROM_LABELS[chrom], ha="center", fontsize=7, color=col)
        axes2[-1].set_xlabel("Genomic position", fontsize=9)
        axes2[-1].set_xticks([])
        title = ("Nucleotide diversity pi" if stat == "tP_site" else "Tajima's D")
        fig.suptitle(f"{title} -- 50kb windows across genome\n(Chr 14 in red)", fontweight="bold")
        plt.tight_layout()
        fig.savefig(f"{OUTDIR}/{fname}", dpi=150, bbox_inches="tight")
        plt.close()
        print(f"{fname} done")
else:
    print("Windowed .pestPG files not yet available -- skipping windowed plots")

print(f"\nAll plots saved to {OUTDIR}/")
