# Results

Multi-species *Acropora* population genomics — *A. cervicornis* and *A. palmata* across
Florida, Panama, and Bonaire, aligned to the *A. palmata* reference genome (Vollmer Lab).

> **Status (2026-03-17):** 290-sample Discovery HPC run — Segment 2 complete (SNP discovery,
> relatedness/clonality). Segment 3 (PCA, admixture) and Segment 4 (diversity, FST) pending.

---

## 290-Sample Production Run

**Samples:** 290 — *A. cervicornis* (178) and *A. palmata* (112) from Florida, Panama, and Bonaire
(`config/samples_RR.csv`). 10 samples failed QC and were excluded before SNP calling.

**Reference:** *Acropora palmata* genome assembly (14 chromosomes, scaffold size ≥ 10 Mb).

**HPC:** Northeastern Discovery, SLURM `short` partition, Singularity containers.

---

### Sequencing QC

| Metric | Value |
|--------|-------|
| Samples submitted | 300 |
| Samples passing QC | 290 |
| Chromosomes analyzed | 14 |
| Sites scanned (pass 1) | — |
| SNPs discovered (pass 1) | — |
| SNPs retained (pass 2, MAF > 0.05, ≥ 80% individuals) | 2,034,805 |

*Pass 1 site count to be added from filter_report.txt.*

---

### Relatedness and Clonality

Pairwise KING kinship coefficients estimated with ngsRelate from high-MAF SNPs.
Clone threshold: KING ≥ 0.45. Close-relative threshold: KING ≥ 0.20.

Two complementary runs were performed and compared:

| Run | Sites | Strategy | Clone pairs | Close relatives |
|-----|-------|----------|-------------|-----------------|
| All-vs-all | 100,000 (MAF > 0.30, random subsample) | All 290 × 290 pairs | 74 | — |
| Grouped (per species × region) | ~783,000 (MAF > 0.30, all sites) | Within-group only (6 parallel jobs) | 74 | 120 |

**Clone results are identical across both runs**, confirming the 100K-subsampled all-vs-all
approach was sufficient for clone detection. The grouped run with ~8× more sites added
resolution for close-relative detection (120 pairs vs. fewer in the subsampled run).

One sample (RR_FL_Apal_015) was excluded as a technical artifact before analysis — it appeared
as a close relative (KING 0.20–0.45) to >60 FL *A. palmata* samples within the grouped run,
a pattern biologically impossible and indicative of a mapping or contamination artifact.
After exclusion, all cross-population FL Apal "relatives" disappeared.

#### Clonality by group

| Group | N | Clone pairs | Excluded | N genets | Genet diversity |
|-------|---|-------------|----------|----------|-----------------|
| *A. cervicornis* BON | 25 | 0 | 0 | 25 | 1.000 |
| *A. cervicornis* FL | 104 | 28 | 13 | 91 | 0.875 |
| *A. cervicornis* PA | 49 | 32 | 10 | 39 | 0.796 |
| *A. palmata* BON | 25 | 0 | 0 | 25 | 1.000 |
| *A. palmata* FL | 88 | 11 | 11 | 77 | 0.875 |
| *A. palmata* PA | 9 | 3 | 2 | 7 | 0.778 |
| **Total** | **300** | **74** | **36** | **264** | **0.880** |

Zero cross-species or cross-population clone pairs detected — clonality is entirely within-group
as expected biologically.

**Key patterns:**
- Bonaire shows 100% genet diversity in both species — consistent with wild reef sampling of
  distinct colonies.
- Panama has the lowest genet diversity (79–78%). The Ac_PA_HS cluster is a single genotype
  sampled 8 times (HS1–8, KING 0.485–0.492); Ac_PA_CK41/42/44 form a 3-sample clone cluster.
  Consistent with nursery/aquaculture sourcing.
- Florida shows intermediate clonality (12.5% both species), reflecting restoration nursery
  collections (CRF TV genotype distributed to multiple sites; MOTE nursery propagation).
  The CRF_TV genotype cluster spans 7 samples: AC_CRF_TV_68/71/78, Ac_FL_K1/K2/U14, RR_FL_Acer_K3
  (KING 0.484–0.492).
- **253 unrelated samples retained** for downstream population structure and diversity analyses.

#### Close relatives (grouped run, 783K sites)

120 within-group pairs with KING 0.20–0.45, flagged but retained in analyses.

Notable clusters:
- **BON_Acer_AC8/11/23/24** — a 4-member kin cluster in Bonaire wild reef *A. cervicornis*
  (KING 0.21–0.30). Likely a family group from a single wild colony neighborhood; biologically
  interesting given Bonaire's 100% genet diversity overall.
- **RR_FL_Apal_009 × RR_FL_Apal_015 KING=0.384** — the highest non-clone KING value (excluding
  the artifact sample). Warrants re-checking the 015 exclusion decision; 009 shows elevated
  relatedness to several other FL Apal samples (016, 050, 179, 180) suggesting a kin network.
- **RR_FL_Apal_158 × RR_FL_Apal_200 KING=0.313**, **011 × 158 KING=0.299** — a Florida Apal
  kin cluster (001, 011, 158, 200 all mutually elevated).
- **AC_MOTE_45 × AC_MOTE_8 KING=0.250**, **AC_FWC_MR_20 × AC_MOTE_10 KING=0.238** — first-degree
  relatives within MOTE and cross-collection, potentially from shared broodstock.

---

### Population Structure

*PCA and admixture (K=2–5) pending — Segment 3.*

Expected: K=2 as primary structure based on prior 96-sample *A. cervicornis* run (Science paper).
Two-lineage structure (lineage A / lineage B) does not split cleanly by geography — both lineages
are present in FL and PA, consistent with the Vollmer Lab hybridization findings.

---

### Genetic Diversity

*Per-population π, θ_W, and Tajima's D pending — Segment 4.*

Prior 96-sample *A. cervicornis* results (reference only):

| Metric | Florida | Panama |
|--------|---------|--------|
| Tajima's D (mean) | +0.74 | +0.53 |
| θ_π | 0.00351 | 0.00354 |

Positive Tajima's D likely driven by two-lineage admixture structure inflating
intermediate-frequency alleles. Full 290-sample run will separate by species and lineage.

---

### FST and Population Differentiation

*Genome-wide FST and per-comparison SFS pending — Segment 4.*

Planned comparisons (from `config/config.yaml`):
- FL vs. PA (*A. cervicornis*, species-wide)
- FL vs. PA (*A. palmata*, species-wide)
- FL vs. BON and PA vs. BON (both species)
- Within-lineage geography: lineageA_FL vs. lineageA_PA; lineageB_FL vs. lineageB_PA
- Cross-lineage: lineageA_all vs. lineageB_all

Prior 5-sample pilot (reference only):

| Population Pair | FST (unweighted) | FST (weighted) |
|----------------|-----------------|----------------|
| *A. cervicornis* FL vs. PA | 0.127 | 0.162 |

---

### Demographic Inference

*Moments/dadi model fitting pending — Segment 6.*

Planned: 2-pop and 4-pop models for FL/PA/BON × species comparisons.
Scripts (`run_moments_2pop.py`, `run_moments_4pop.py`) pending implementation.

---

## Previous Runs

### 96-Sample AWS Run (*A. cervicornis* FL + PA only)

**Status:** Complete through diversity/PCA/admixture. FST killed (no `-maxIter` cap).
**Key findings:**
- K=2 confirmed as primary structure; two lineages present in both FL and PA
- Tajima's D positive both populations (FL +0.74, PA +0.53); Chr 14 outlier (FL D=+1.10)
- θ_π nearly identical FL (0.00351) vs PA (0.00354)
- 78 unrelated samples retained after clone gate (15 excluded)

### 5-Sample Pilot (*A. cervicornis* FL + PA)

**Status:** Complete.

| Metric | Value |
|--------|-------|
| SNPs retained (MAF > 0.10, 80% individuals) | 2,062,496 |
| FST FL vs. PA (unweighted / weighted) | 0.127 / 0.162 |
| Tajima's D Florida | +0.41 |
| Tajima's D Panama | +0.29 |

---

## Reproducibility

All results reproducible from BAMs in `/projects/vollmer/RR_heat-tolerance/Acropora/2_mapping.bwa/`
and `config/samples_RR.csv` using:

```bash
cd /work/vollmer/acropora_genomics
bash run.sh 2a   # SNP discovery + relatedness (stops at clone gate)
# → python3 /projects/vollmer/coral-angsd-pipeline/workflow/scripts/clone_approve.py
bash run.sh 2b   # subset beagle to unrelated samples
bash run.sh 2c   # grouped ngsRelate validation (optional, comparison run)
bash run.sh 3    # PCA, admixture, LD
bash run.sh 4    # SAF, SFS, diversity, FST
```

Reference genome: *Acropora palmata* assembly (`reference.fna`).
SNP coordinates are in *A. palmata* reference space — liftover to a native *A. cervicornis*
assembly pending.
