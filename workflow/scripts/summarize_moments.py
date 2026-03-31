#!/usr/bin/env python3
"""
summarize_moments.py — Parse moments 2-pop inference JSONs and convert
dimensionless parameters to physical units.

moments parameters are dimensionless, scaled by ancestral Ne (N_anc):
  nu_i  = Ne_i / N_anc
  T     = t / (2 * N_anc)   [generations]
  m_ij  = migration_rate * 2 * N_anc

N_anc is estimated from Watterson's theta by default:
  theta_W = n_seg / (L * a1)   where a1 = sum(1/i for i=1..n_total-1)
  N_anc   = theta_W / (4 * mu)

WARNING: for interspecific comparisons the joint 2D SFS only retains sites
segregating in both populations — fixed differences are excluded. This
dramatically underestimates theta and therefore N_anc. Use --n-anc or
--n-anc-table to supply externally calibrated values (e.g. from SMC++).

Usage:
    python summarize_moments.py results/demography/inference/ \\
        --L 217341919 --mu 1.8e-8 --gen 5

    # Override N_anc globally:
    python summarize_moments.py ... --n-anc 15000

    # Per-comparison N_anc from JSON {"comparison_name": N_anc, ...}:
    python summarize_moments.py ... --n-anc-table n_anc.json
"""
import json, sys, os, math, argparse

def harmonic(n):
    return sum(1/i for i in range(1, n))

def summarize(inference_dir, L, mu, gen_time, n_anc_global=None, n_anc_table=None):
    files = sorted(f for f in os.listdir(inference_dir) if f.endswith(".moments.json"))

    n_anc_source = "theta"
    if n_anc_global:
        n_anc_source = "supplied (global)"
    elif n_anc_table:
        n_anc_source = "supplied (per-comparison)"

    print(f"Assumptions: mu={mu:.1e}/site/gen, gen_time={gen_time}yr, L={L:,} callable sites")
    print(f"N_anc source: {n_anc_source}\n")
    print(f"{'Comparison':<24} {'Best model':<8} {'dAIC_2nd':>9}  "
          f"{'N_anc':>7} {'src':>5}  {'Ne1':>7}  {'Ne2':>7}  "
          f"{'T_split(yr)':>12}  {'m12(per gen)':>13}  {'m21(per gen)':>13}")
    print("-" * 125)

    rows = []
    for fname in files:
        with open(os.path.join(inference_dir, fname)) as f:
            d = json.load(f)

        comp   = d["comparison"]
        ns     = d["ns"]            # [n1, n2] haploid projections
        n_seg  = d["n_segregating"]
        best   = d["best_model"]
        models = d["models"]

        # N_anc: use supplied value if available, otherwise Watterson's theta
        if n_anc_global:
            N_anc = float(n_anc_global)
            src = "ext"
        elif n_anc_table and comp in n_anc_table:
            N_anc = float(n_anc_table[comp])
            src = "ext"
        else:
            n_total = sum(ns)
            a1 = harmonic(n_total)
            theta_w = n_seg / (L * a1)
            N_anc = theta_w / (4 * mu)
            src = "θ"

        # Skip incomplete runs
        if best is None or not models:
            print(f"{comp:<24} {'(no results)'}")
            continue

        # Best model params
        p = models[best]["params"]
        best_aic = models[best]["aic"]

        # Second-best AIC for dAIC calculation
        aics = sorted((v["aic"], k) for k, v in models.items())
        second_aic = aics[1][0] if aics[0][1] == best else aics[0][0]
        dAIC_2nd = second_aic - best_aic

        # Convert to physical units
        nu1 = p.get("nu1", p.get("nu", None))
        nu2 = p.get("nu2", None)

        # T: use T, T1 (isolation-with-migration period), or T1+T2
        if "T" in p:
            T_dim = p["T"]
        elif "T1" in p and "T2" in p:
            T_dim = p["T1"] + p["T2"]   # total time since split
        else:
            T_dim = None

        Ne1 = int(nu1 * N_anc) if nu1 else None
        Ne2 = int(nu2 * N_anc) if nu2 else None
        T_yr = T_dim * 2 * N_anc * gen_time if T_dim else None

        # Migration rates (per generation = m / (2*N_anc))
        m12_pg = p["m12"] / (2 * N_anc) if "m12" in p else (p["m"] / (2*N_anc) if "m" in p else None)
        m21_pg = p["m21"] / (2 * N_anc) if "m21" in p else m12_pg

        row = dict(comp=comp, best=best, dAIC_2nd=dAIC_2nd,
                   N_anc=int(N_anc), src=src, Ne1=Ne1, Ne2=Ne2, T_yr=T_yr,
                   m12=m12_pg, m21=m21_pg, params=p, model_name=best)
        rows.append(row)

        m12_str = f"{m12_pg:.2e}" if m12_pg else "—"
        m21_str = f"{m21_pg:.2e}" if m21_pg and "m21" in p else "—"
        T_str   = f"{T_yr:>12,.0f}" if T_yr else f"{'—':>12}"
        print(f"{comp:<24} {best:<8} {dAIC_2nd:>9.1f}  "
              f"{int(N_anc):>7,} {src:>5}  "
              f"{Ne1 if Ne1 else 0:>7,}  "
              f"{Ne2 if Ne2 else 0:>7,}  "
              f"{T_str}  "
              f"{m12_str:>13}  {m21_str:>13}")

    print()
    print("Parameter key:")
    print("  N_anc   = ancestral Ne estimated from Watterson's theta")
    print("  Ne1/Ne2 = current Ne of pop1/pop2 (nu * N_anc)")
    print("  T_split = time since split (or total divergence time for AM model)")
    print("  m12     = migration FROM pop2 INTO pop1 per generation")
    print("  m21     = migration FROM pop1 INTO pop2 per generation")
    print()
    print("Model codes:")
    print("  SI   = strict isolation (no migration)")
    print("  IM   = isolation with symmetric migration")
    print("  IM_a = isolation with asymmetric migration  ← most common winner")
    print("  AM   = ancestral migration (ancient gene flow, then isolated)")
    print("  SC   = secondary contact (isolated, then reconnected)")
    print()
    print("Detailed parameters for best models:")
    print("-" * 60)
    for r in rows:
        print(f"\n{r['comp']}  [{r['best']}]  N_anc={r['N_anc']:,} ({r['src']})")
        pop_ids = None
        # try to recover pop_ids from JSON for direction labels
        json_path = os.path.join(inference_dir, r['comp'] + ".moments.json")
        if os.path.exists(json_path):
            with open(json_path) as jf:
                jd = json.load(jf)
                pop_ids = jd.get("pop_ids")
        for k, v in r["params"].items():
            unit = ""
            if k.startswith("nu"):   unit = f"  → Ne={v * r['N_anc']:,.0f}"
            if k == "T":             unit = f"  → {v * 2 * r['N_anc'] * gen_time:,.0f} yr"
            if k in ("T1","T2"):     unit = f"  → {v * 2 * r['N_anc'] * gen_time:,.0f} yr"
            if k.startswith("m"):    unit = f"  → {v/(2*r['N_anc']):.2e}/gen"
            direction = ""
            if pop_ids:
                if k == "m12":   direction = f"  ({pop_ids[1]} → {pop_ids[0]})"
                elif k == "m21": direction = f"  ({pop_ids[0]} → {pop_ids[1]})"
            print(f"  {k:6s} = {v:12.6f}{unit}{direction}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("inference_dir")
    ap.add_argument("--L",           type=float, default=217341919)
    ap.add_argument("--mu",          type=float, default=1.8e-8)
    ap.add_argument("--gen",         type=float, default=5)
    ap.add_argument("--n-anc",       type=float, default=None,
                    help="Override N_anc globally (e.g. from SMC++)")
    ap.add_argument("--n-anc-table", type=str,   default=None,
                    help="JSON file mapping comparison name to N_anc")
    args = ap.parse_args()

    n_anc_table = None
    if args.n_anc_table:
        with open(args.n_anc_table) as f:
            n_anc_table = json.load(f)

    summarize(args.inference_dir, args.L, args.mu, args.gen,
              n_anc_global=args.n_anc, n_anc_table=n_anc_table)
