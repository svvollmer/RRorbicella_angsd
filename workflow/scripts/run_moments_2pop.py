#!/usr/bin/env python3
"""
run_moments_2pop.py — Two-population demographic inference with moments.

Fits 5 models to a folded 2D SFS (output of realSFS dadi):
  SI   — strict isolation (split, no migration ever)
  IM   — isolation with migration, symmetric
  IM_a — isolation with migration, asymmetric
  SC   — secondary contact: split → isolation → secondary contact
  AM   — ancient migration: split → migration → isolation

Each model is fit with --restarts random starts; best log-likelihood kept.
AIC used for model comparison. All times in units of 2*Ne generations;
all population sizes relative to ancestral Ne.

Usage:
    python run_moments_2pop.py \\
        --sfs data.dadi.txt \\
        --projections projections.json \\
        --pop-ids pop1 pop2 \\
        --out results/demography/inference/comparison_name \\
        --restarts 20 \\
        --models SI IM IM_a SC AM
"""

import argparse
import json
import os
import sys
import time
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

try:
    import moments
except ImportError:
    sys.exit(
        "ERROR: moments not installed.\n"
        "Install from source: pip install git+https://github.com/MomentsLD/moments.git"
    )


# ============================================================================
# Model definitions
# All take (params, ns) and return a folded 2D Spectrum.
# ns = [n1, n2] in haploids (projected sample sizes).
# ============================================================================

def model_SI(params, ns):
    """
    Strict Isolation: split, no migration.
    params: nu1, nu2, T
      nu1, nu2 — relative sizes after split (relative to ancestral Ne)
      T        — time since split (2Ne generations)
    """
    nu1, nu2, T = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T)
    return fs.fold()


def model_IM(params, ns):
    """
    Isolation with Migration, symmetric.
    params: nu1, nu2, T, m
      m — symmetric migration rate (2*Ne*m per generation)
    """
    nu1, nu2, T, m = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T, m=np.array([[0, m], [m, 0]]))
    return fs.fold()


def model_IM_a(params, ns):
    """
    Isolation with Migration, asymmetric.
    params: nu1, nu2, T, m12, m21
      m12 — migration pop1 → pop2
      m21 — migration pop2 → pop1
    """
    nu1, nu2, T, m12, m21 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T, m=np.array([[0, m12], [m21, 0]]))
    return fs.fold()


def model_SC(params, ns):
    """
    Secondary Contact: split → isolation → secondary contact.
    Relevant for recent hybridization after a period of isolation.
    params: nu1, nu2, T1, T2, m
      T1 — duration of isolation
      T2 — duration of secondary contact (ongoing migration)
      m  — symmetric migration rate during contact
    """
    nu1, nu2, T1, T2, m = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T1)
    fs.integrate([nu1, nu2], T2, m=np.array([[0, m], [m, 0]]))
    return fs.fold()


def model_AM(params, ns):
    """
    Ancient Migration: split with early gene flow → isolation.
    Relevant for ancient hybridization that has since stopped.
    params: nu1, nu2, T1, T2, m
      T1 — duration of ancient migration phase
      T2 — duration of isolation (more recent, no gene flow)
      m  — symmetric migration rate in ancient phase
    """
    nu1, nu2, T1, T2, m = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T1, m=np.array([[0, m], [m, 0]]))
    fs.integrate([nu1, nu2], T2)
    return fs.fold()


# name → (func, param_names, n_params)
MODELS = {
    "SI":  (model_SI,  ["nu1", "nu2", "T"],              3),
    "IM":  (model_IM,  ["nu1", "nu2", "T", "m"],         4),
    "IM_a":(model_IM_a,["nu1", "nu2", "T", "m12", "m21"],5),
    "SC":  (model_SC,  ["nu1", "nu2", "T1", "T2", "m"],  5),
    "AM":  (model_AM,  ["nu1", "nu2", "T1", "T2", "m"],  5),
}

# Log-uniform sampling ranges for random starts.
# Migration lower bound is 1e-5: critical for inter-species comparisons where true
# m may be near-zero (rare introgression). A bound of 0.01 would prevent the
# optimizer from finding solutions in the plausible range and bias SI vs IM AIC
# comparisons — the IM model would always report m=0.01 at the lower bound,
# making it look like migration is supported even when it isn't.
PARAM_RANGES = {
    "nu1": (0.01,  10.0),
    "nu2": (0.01,  10.0),
    "T":   (0.01,   5.0),
    "T1":  (0.01,   4.0),
    "T2":  (0.001,  2.0),
    "m":   (1e-5,  10.0),   # allow near-zero gene flow
    "m12": (1e-5,  10.0),   # allow near-zero gene flow
    "m21": (1e-5,  10.0),   # allow near-zero gene flow
}


def random_start(param_names, rng):
    p = []
    for name in param_names:
        lo, hi = PARAM_RANGES[name]
        p.append(np.exp(rng.uniform(np.log(lo), np.log(hi))))
    return p


def _single_restart(p0, fs_data, model_name, ns):
    """Top-level function for ProcessPoolExecutor (must be picklable)."""
    func, param_names, k = MODELS[model_name]
    lb = [1e-6 if n in ("m", "m12", "m21") else 1e-4 for n in param_names]
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            popt = moments.Inference.optimize_log(
                p0, fs_data, func,
                lower_bound=lb,
                upper_bound=[100.0] * k,
                verbose=0,
                maxiter=200,
            )
        ll = moments.Inference.ll_multinom(func(popt, ns), fs_data)
        if np.isfinite(ll):
            return ll, list(popt)
    except Exception:
        pass
    return -np.inf, None


def fit_model(fs_data, model_name, restarts, ns, rng, verbose=False, threads=1):
    func, param_names, k = MODELS[model_name]
    p0_list = [random_start(param_names, rng) for _ in range(restarts)]

    best_ll = -np.inf
    best_params = None

    if threads > 1:
        with ProcessPoolExecutor(max_workers=threads) as pool:
            futures = {pool.submit(_single_restart, p0, fs_data, model_name, ns): i
                       for i, p0 in enumerate(p0_list)}
            for fut in as_completed(futures):
                i = futures[fut]
                ll, popt = fut.result()
                if popt is not None and ll > best_ll:
                    best_ll = ll
                    best_params = popt
                    if verbose:
                        print(f"  [{model_name}] restart {i+1}/{restarts}: ll={ll:.4f}", flush=True)
                elif verbose and popt is None:
                    print(f"  [{model_name}] restart {i+1}/{restarts}: failed", flush=True)
    else:
        for i, p0 in enumerate(p0_list):
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    # lower_bound per-parameter: migration params allow 1e-6 to
                    # detect near-zero inter-species gene flow; size/time params 1e-4
                    lb = [1e-6 if n in ("m","m12","m21") else 1e-4
                          for n in param_names]
                    popt = moments.Inference.optimize_log(
                        p0, fs_data, func,
                        lower_bound=lb,
                        upper_bound=[100.0] * k,
                        verbose=0,
                        maxiter=200,
                    )
                ll = moments.Inference.ll_multinom(func(popt, ns), fs_data)
                if np.isfinite(ll) and ll > best_ll:
                    best_ll = ll
                    best_params = list(popt)
                    if verbose:
                        print(f"  [{model_name}] restart {i+1}/{restarts}: ll={ll:.4f}", flush=True)
            except Exception as e:
                if verbose:
                    print(f"  [{model_name}] restart {i+1}/{restarts}: failed ({e})", flush=True)

    return best_ll, best_params


def aic(ll, k):
    return 2 * k - 2 * ll



def load_realsfs_2dsfs(path, n_ind1, n_ind2, ns, pop_ids):
    """
    Load a realSFS 2D SFS flat file (output of: realSFS saf1.idx saf2.idx -fold 1).
    Dimensions: (2*n_ind1+1) x (2*n_ind2+1). Projects to ns and folds.
    Some realSFS versions output two rows (iteration log + ML estimate); if so,
    the last row is the ML estimate.
    """
    data = np.fromstring(open(path).read(), sep=" ")
    dim1 = 2 * n_ind1 + 1
    dim2 = 2 * n_ind2 + 1
    expected = dim1 * dim2
    if data.size == 2 * expected:
        # Two-row output: take last row (ML estimate)
        data = data[expected:]
    elif data.size != expected:
        raise ValueError(
            f"Expected {dim1}x{dim2}={expected} values for n_ind=({n_ind1},{n_ind2}), "
            f"got {data.size}. Check --n-ind values."
        )
    fs = moments.Spectrum(data.reshape((dim1, dim2)), pop_ids=pop_ids)
    return fs.project(ns).fold()


def main():
    parser = argparse.ArgumentParser(description="Two-population demographic inference with moments")
    parser.add_argument("--sfs",         required=True,              help="2D SFS file (dadi or realSFS format)")
    parser.add_argument("--sfs-format",  default="realsfs",
                        choices=["realsfs", "dadi"],            help="SFS file format (default: realsfs)")
    parser.add_argument("--n-ind",       nargs=2, type=int,     help="Actual n_individuals per pop (required for realsfs format)")
    parser.add_argument("--projections", required=True,              help="JSON {group: n_haploids}")
    parser.add_argument("--pop-ids",     nargs=2,                    help="Population IDs")
    parser.add_argument("--out",         required=True,              help="Output prefix")
    parser.add_argument("--restarts",    type=int,   default=20,     help="Random restarts per model")
    parser.add_argument("--models",      nargs="+",  default=list(MODELS.keys()),
                        choices=list(MODELS.keys()))
    parser.add_argument("--seed",        type=int,   default=42)
    parser.add_argument("--threads",     type=int,   default=1,      help="Parallel workers for restarts")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)

    # Load projections
    with open(args.projections) as f:
        proj = json.load(f)
    groups = list(proj.keys())
    ns = [proj[g] for g in groups]
    pop_ids = args.pop_ids or groups

    print(f"Loading SFS: {args.sfs}")
    print(f"Populations: {pop_ids[0]} (n={ns[0]}) vs {pop_ids[1]} (n={ns[1]})")

    if args.sfs_format == "realsfs":
        if not args.n_ind:
            parser.error("--n-ind N1 N2 required when --sfs-format realsfs")
        fs_data = load_realsfs_2dsfs(args.sfs, args.n_ind[0], args.n_ind[1], ns, pop_ids)
    else:
        try:
            fs_raw = moments.Spectrum.from_file(args.sfs, pop_ids=pop_ids)
        except TypeError:
            fs_raw = moments.Spectrum.from_file(args.sfs)
            if hasattr(fs_raw, "pop_ids"):
                fs_raw.pop_ids = pop_ids
        fs_data = fs_raw.project(ns).fold()
    print(f"SFS projected to {ns}, S={fs_data.S():.0f} segregating sites")
    print()

    results = {}
    for model_name in args.models:
        _, param_names, k = MODELS[model_name]
        print(f"Fitting {model_name} ({k} params, {args.restarts} restarts)...", flush=True)
        t0 = time.time()
        best_ll, best_params = fit_model(fs_data, model_name, args.restarts, ns, rng, args.verbose, args.threads)
        elapsed = time.time() - t0

        if best_params is None:
            print(f"  {model_name}: ALL RESTARTS FAILED")
            results[model_name] = {"ll": None, "aic": None, "k": k, "params": None, "elapsed_s": elapsed}
        else:
            aic_val = aic(best_ll, k)
            results[model_name] = {
                "ll":       float(best_ll),
                "aic":      float(aic_val),
                "k":        k,
                "params":   {n: float(v) for n, v in zip(param_names, best_params)},
                "elapsed_s": round(elapsed, 1),
            }
            print(f"  {model_name}: ll={best_ll:.4f}  AIC={aic_val:.2f}  ({elapsed:.0f}s)")

    # Rank by AIC
    ranked = sorted(
        [(m, r) for m, r in results.items() if r.get("aic") is not None],
        key=lambda x: x[1]["aic"]
    )
    best_model = ranked[0][0] if ranked else None
    best_aic   = ranked[0][1]["aic"] if ranked else None

    # Write JSON
    output = {
        "comparison":      os.path.basename(args.out),
        "pop_ids":         pop_ids,
        "ns":              ns,
        "n_segregating":   float(fs_data.S()),
        "restarts":        args.restarts,
        "best_model":      best_model,
        "best_aic":        best_aic,
        "models":          results,
    }
    with open(args.out + ".moments.json", "w") as f:
        json.dump(output, f, indent=2)

    # Write human-readable summary
    with open(args.out + ".bestfit.txt", "w") as f:
        f.write(f"Two-population demographic inference — {os.path.basename(args.out)}\n")
        f.write(f"Populations: {pop_ids[0]} (n={ns[0]}) vs {pop_ids[1]} (n={ns[1]})\n")
        f.write(f"Segregating sites: {fs_data.S():.0f}  |  Restarts: {args.restarts}\n\n")
        f.write(f"{'Model':<8}  {'k':>3}  {'log-L':>12}  {'AIC':>12}  {'dAIC':>8}\n")
        f.write("-" * 50 + "\n")
        for name, res in ranked:
            delta = res["aic"] - best_aic
            f.write(f"{name:<8}  {res['k']:>3}  {res['ll']:>12.4f}  {res['aic']:>12.4f}  {delta:>8.2f}\n")
        f.write("\n")
        if best_model and results[best_model]["params"]:
            f.write(f"Best model: {best_model}\n")
            f.write("Parameters:\n")
            for pname, pval in results[best_model]["params"].items():
                f.write(f"  {pname}: {pval:.6f}\n")

    # Print summary
    print()
    print("=" * 50)
    print(f"Model comparison — {os.path.basename(args.out)}")
    print("=" * 50)
    print(f"{'Model':<8}  {'k':>3}  {'log-L':>12}  {'AIC':>12}  {'dAIC':>8}")
    print("-" * 50)
    for name, res in ranked:
        delta = res["aic"] - best_aic
        marker = " <-- BEST" if name == best_model else ""
        print(f"{name:<8}  {res['k']:>3}  {res['ll']:>12.4f}  {res['aic']:>12.4f}  {delta:>8.2f}{marker}")
    print("=" * 50)
    print(f"\nResults: {args.out}.moments.json")
    print(f"Summary: {args.out}.bestfit.txt")


if __name__ == "__main__":
    main()
