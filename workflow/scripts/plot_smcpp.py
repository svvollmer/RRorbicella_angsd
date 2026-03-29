#!/usr/bin/env python3
"""
plot_smcpp.py — Annotated SMC++ Ne(t) plot for A. palmata and A. cervicornis.

Reads CSV output from `smc++ plot --csv` and produces a publication-quality
Ne(t) figure with:
  - Log-log axes
  - Bottleneck and present-day Ne annotations
  - Key geological/ecological events noted on x-axis
  - Assumption box (mu, generation time)

Usage:
    python plot_smcpp.py \
        docs/figures/apal_Ne.csv \
        docs/figures/acer_Ne.csv \
        docs/figures/smcpp_Ne_annotated.png
"""
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# ── Parameters ────────────────────────────────────────────────────────────────
MU        = 1.8e-8   # per site per generation (Acropora; see note)
GEN_TIME  = 5        # years per generation (branching Acropora sexual maturity)

COL_APAL = "#2166ac"   # blue  — A. palmata
COL_ACER = "#d6604d"   # red   — A. cervicornis

def load_csv(path):
    df = pd.read_csv(path)
    # smc++ outputs plot_type='path' for the spline; exclude knot/step rows
    df = df[df["plot_type"] == "path"].copy()
    # t=0 appears twice (step function artifact); keep unique x
    df = df.drop_duplicates(subset="x")
    df = df.sort_values("x").reset_index(drop=True)
    return df

def main(apal_csv, acer_csv, outpng):
    apal = load_csv(apal_csv)
    acer = load_csv(acer_csv)

    # Exclude t=0 from log axis (plot from first nonzero t onward)
    apal_plot = apal[apal["x"] > 0]
    acer_plot = acer[acer["x"] > 0]

    # ── Present-day values (t closest to 0 but > 0) ──────────────────────────
    apal_present = apal.iloc[0]["y"]   # first row is t=0
    acer_present = acer.iloc[0]["y"]

    # ── Bottleneck (minimum Ne) ───────────────────────────────────────────────
    apal_bn_idx = apal_plot["y"].idxmin()
    acer_bn_idx = acer_plot["y"].idxmin()
    apal_bn_t, apal_bn_ne = apal_plot.loc[apal_bn_idx, ["x","y"]]
    acer_bn_t, acer_bn_ne = acer_plot.loc[acer_bn_idx, ["x","y"]]

    # ── Plot ──────────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(9, 5.5))

    ax.plot(apal_plot["x"], apal_plot["y"], color=COL_APAL, lw=2.2,
            label=r"$A.\ palmata$" + f"  (n=49 FL, Q≥0.99)")
    ax.plot(acer_plot["x"], acer_plot["y"], color=COL_ACER, lw=2.2,
            label=r"$A.\ cervicornis$" + f"  (n=82 FL, Q≥0.99)")

    # Present-day Ne markers
    ax.axvline(apal_plot["x"].iloc[0], color=COL_APAL, lw=0.8, ls=":", alpha=0.6)
    ax.axvline(acer_plot["x"].iloc[0], color=COL_ACER, lw=0.8, ls=":", alpha=0.6)

    # Annotate present-day Ne
    ax.annotate(f"Ne={apal_present:,.0f}",
                xy=(apal_plot["x"].iloc[0], apal_present),
                xytext=(1.6e3, apal_present * 1.25),
                color=COL_APAL, fontsize=8,
                arrowprops=dict(arrowstyle="-", color=COL_APAL, lw=0.8))
    ax.annotate(f"Ne={acer_present:,.0f}",
                xy=(acer_plot["x"].iloc[0], acer_present),
                xytext=(1.6e3, acer_present * 0.6),
                color=COL_ACER, fontsize=8,
                arrowprops=dict(arrowstyle="-", color=COL_ACER, lw=0.8))

    # Shade LGM window (~20–26 kya)
    ax.axvspan(20000, 26000, alpha=0.08, color="gray", label="LGM (~20–26 kya)")

    # Annotate LGM
    ax.text(23000, ax.get_ylim()[0] if ax.get_ylim()[0] > 0 else 3000,
            "LGM", ha="center", va="bottom", fontsize=7.5, color="gray",
            rotation=90)

    # ── Axes formatting ───────────────────────────────────────────────────────
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Years ago", fontsize=11)
    ax.set_ylabel("Effective population size (Ne)", fontsize=11)
    ax.set_title(r"$Acropora$ population size history (SMC++)", fontsize=12)
    ax.set_xlim(1e3, 6e5)
    ax.set_ylim(2e3, 8e4)

    ax.legend(loc="upper left", fontsize=9, framealpha=0.85)

    # ── Assumption box ────────────────────────────────────────────────────────
    assume_txt = (
        f"Assumptions:\n"
        f"  μ = {MU:.1e} per site/gen\n"
        f"    (Acropora; see Matz et al.)\n"
        f"  Generation time = {GEN_TIME} yr\n"
        f"  Cubic spline, folded SFS\n"
        f"  NOTE: model.final.json fit with old μ=3.4×10⁻⁸\n"
        f"  Re-run pending with corrected μ"
    )
    ax.text(0.99, 0.03, assume_txt, transform=ax.transAxes,
            fontsize=7.5, va="bottom", ha="right",
            bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="gray", alpha=0.85))

    # ── Key statistics annotation ─────────────────────────────────────────────
    stats_txt = (
        f"Bottleneck:\n"
        f"  {r'$A.p.$'} Ne={apal_bn_ne:,.0f} at ~{apal_bn_t/1000:.0f} kya\n"
        f"  {r'$A.c.$'} Ne={acer_bn_ne:,.0f} at ~{acer_bn_t/1000:.0f} kya"
    )
    ax.text(0.01, 0.03, stats_txt, transform=ax.transAxes,
            fontsize=7.5, va="bottom", ha="left",
            bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="gray", alpha=0.85))

    plt.tight_layout()
    plt.savefig(outpng, dpi=180, bbox_inches="tight")
    print(f"Saved: {outpng}")

    # Print summary stats
    print(f"\nSummary:")
    print(f"  A. palmata  present Ne: {apal_present:>8,.0f}")
    print(f"  A. palmata  bottleneck: {apal_bn_ne:>8,.0f}  at {apal_bn_t/1000:.0f} kya")
    print(f"  A. palmata  ancient Ne: {apal['y'].max():>8,.0f}  at {apal.loc[apal['y'].idxmax(),'x']/1000:.0f} kya")
    print(f"  A. cerv.    present Ne: {acer_present:>8,.0f}")
    print(f"  A. cerv.    bottleneck: {acer_bn_ne:>8,.0f}  at {acer_bn_t/1000:.0f} kya")
    print(f"  A. cerv.    ancient Ne: {acer['y'].max():>8,.0f}  at {acer.loc[acer['y'].idxmax(),'x']/1000:.0f} kya")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3])
