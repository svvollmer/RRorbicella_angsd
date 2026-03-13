"""
================================================================================
common.smk — Shared configuration, sample loading, and helper functions
Included by all segment Snakefiles via `include: "common.smk"`
================================================================================
"""

import os
import re
import pandas as pd
from itertools import combinations

configfile: "config/config.yaml"

# ---------------------------------------------------------------------------
# Shell PATH setup for local runs (conda env).
# Leave local_conda_env unset on HPC/AWS — Singularity handles PATH.
# ---------------------------------------------------------------------------
_conda_env = config.get("local_conda_env", "")
if _conda_env and os.path.isdir(_conda_env):
    shell.executable("/bin/bash")
    shell.prefix(f"export PATH={_conda_env}/bin:/usr/bin:/bin:/usr/sbin:/sbin:$PATH; ")

# ---------------------------------------------------------------------------
# Sample table
# ---------------------------------------------------------------------------
samples = pd.read_csv(config["samples_csv"])
SAMPLES = samples["sample_id"].tolist()
POPS    = samples[config.get("primary_grouping", "population")].unique().tolist()
POP_PAIRS = list(combinations(sorted(POPS), 2))

# Prevent {popname} wildcard from matching "all" or other non-population strings
_pop_pattern = "|".join(re.escape(p) for p in POPS)
wildcard_constraints:
    popname  = _pop_pattern,
    popname1 = _pop_pattern,
    popname2 = _pop_pattern,

# ---------------------------------------------------------------------------
# Approved sample count (available after Gate 1; falls back to full CSV count)
# ---------------------------------------------------------------------------
N_SAMPLES      = len(SAMPLES)
_approved_path = "results/qc/samples_approved.txt"
N_APPROVED = (
    sum(1 for l in open(_approved_path)
        if l.strip() and not l.startswith('#') and not l.startswith('sample'))
    if os.path.exists(_approved_path) else N_SAMPLES
)

MIN_DEPTH_TOTAL   = N_APPROVED * 2        # 2 reads × N samples minimum
MIN_IND_DISCOVERY = int(N_APPROVED * 0.6) # relaxed — pass 1 only

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def get_sra(wc):
    return samples.loc[samples["sample_id"] == wc.sample, "sra_accession"].values[0]

def get_pop_samples(wc):
    """CRAMs for a given population wildcard (used by make_bamlist_pop)."""
    s = samples.loc[
        samples[config.get("primary_grouping", "population")] == wc.popname,
        "sample_id"
    ].tolist()
    return expand("results/bams/{sample}.filtered.cram", sample=s)

def load_filter_params():
    """
    Read MAF and minInd from filter_params.yaml (written by filter_select.py).
    Returns (min_maf, min_ind) as (float, int).
    """
    import yaml
    params_file = "results/angsd/filter_params.yaml"
    if not os.path.exists(params_file):
        return None, None
    with open(params_file) as f:
        p = yaml.safe_load(f)
    return float(p["min_maf"]), int(p["min_ind"])
