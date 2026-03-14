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
samples = pd.read_csv(config["samples_csv"], keep_default_na=False)

# Required columns
SAMPLES   = samples["sample_id"].tolist()
POPS      = samples[config.get("primary_grouping", "population")].unique().tolist()
POP_PAIRS = list(combinations(sorted(POPS), 2))

# Optional: species — enables multi-species analyses (demography, cross-species FST)
SPECIES = samples["species"].unique().tolist() if "species" in samples.columns else None

# Optional: lineage — enables lineage-stratified FST and demography lineage models
HAS_LINEAGE = (
    "lineage" in samples.columns
    and samples["lineage"].replace("", pd.NA).notna().any()
)

# Optional: admix_fraction — enables clean/admixed sample splits for demography
HAS_ADMIX = (
    "admix_fraction" in samples.columns
    and samples["admix_fraction"].replace("", pd.NA).notna().any()
)

# Optional: phenotypes — enables WGAS
HAS_HEAT    = "heat_tolerance_score"     in samples.columns and samples["heat_tolerance_score"].replace("", pd.NA).notna().any()
HAS_DISEASE = "disease_resistance_score" in samples.columns and samples["disease_resistance_score"].replace("", pd.NA).notna().any()

# Optional: spatial — enables IBD / seascape analyses
HAS_SPATIAL = (
    "lat" in samples.columns and "lon" in samples.columns
    and samples["lat"].replace("", pd.NA).notna().any()
)

# Log which optional analyses are enabled
import sys as _sys
_enabled  = [k for k, v in [("lineage-stratified FST", HAS_LINEAGE),
                              ("demography admix-split", HAS_ADMIX),
                              ("WGAS heat",    HAS_HEAT),
                              ("WGAS disease", HAS_DISEASE),
                              ("spatial/IBD",  HAS_SPATIAL)] if v]
_disabled = [k for k, v in [("lineage-stratified FST", HAS_LINEAGE),
                              ("demography admix-split", HAS_ADMIX),
                              ("WGAS heat",    HAS_HEAT),
                              ("WGAS disease", HAS_DISEASE),
                              ("spatial/IBD",  HAS_SPATIAL)] if not v]
if _enabled:
    print(f"[common.smk] Optional analyses ENABLED:  {', '.join(_enabled)}", file=_sys.stderr)
if _disabled:
    print(f"[common.smk] Optional analyses DISABLED: {', '.join(_disabled)} (add columns to samples.csv to enable)", file=_sys.stderr)

# ---------------------------------------------------------------------------
# BAM/CRAM format — set bam_ext: "bam" or "cram" in config.yaml
# IDX_EXT: index extension (.bam.bai or .cram.crai)
# ---------------------------------------------------------------------------
BAM_EXT = config.get("bam_ext", "cram")
IDX_EXT = f"{BAM_EXT}.bai" if BAM_EXT == "bam" else f"{BAM_EXT}.crai"

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
    """BAMs/CRAMs for a given population wildcard (used by make_bamlist_pop)."""
    s = samples.loc[
        samples[config.get("primary_grouping", "population")] == wc.popname,
        "sample_id"
    ].tolist()
    return expand(f"results/bams/{{sample}}.filtered.{BAM_EXT}", sample=s)

def get_species_samples(species):
    """Sample IDs for a given species (requires species column)."""
    if not SPECIES:
        raise ValueError("species column not present in samples.csv")
    return samples.loc[samples["species"] == species, "sample_id"].tolist()

def get_clean_samples(pop, species=None, admix_threshold=None):
    """
    Sample IDs passing admixture filter (admix_fraction < threshold).
    Falls back to all samples if admix_fraction column absent.
    """
    threshold = admix_threshold or config.get("demography", {}).get("admixture_threshold", 0.05)
    mask = samples["population"] == pop
    if species:
        mask &= samples["species"] == species
    if HAS_ADMIX:
        admix = pd.to_numeric(samples["admix_fraction"], errors="coerce")
        mask &= admix.fillna(1.0) < threshold
    return samples.loc[mask, "sample_id"].tolist()

def get_lineage_samples(lineage, pop=None):
    """Sample IDs for a given lineage label (requires lineage column)."""
    if not HAS_LINEAGE:
        raise ValueError("lineage column not present or empty in samples.csv")
    mask = samples["lineage"] == lineage
    if pop:
        mask &= samples["population"] == pop
    return samples.loc[mask, "sample_id"].tolist()

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
