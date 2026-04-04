# Results

Three-species *Orbicella* population genomics — *O. annularis*, *O. faveolata*, and *O. franksi*
from the Reef Renewal Florida restoration program, aligned to the *O. franksi* reference genome
(jaOrbFran1.1 / GCA_964199315.1).

> **Status (2026-04-04):** Segments 1–3 complete. 134 unrelated samples. PCA + admixture done.
> Lineage assignment gate pending → Segment 4 (diversity/FST) next.

---

## Samples

| Species | N | Status |
|---------|---|--------|
| *O. annularis* | 63 | Complete |
| *O. faveolata* | 57 | Complete (incl. Ofra_26 reclassified) |
| *O. franksi* | 21 | Complete (Ofra_26 excluded as mislabel) |
| **Total** | **140** | |

**Reference:** *O. franksi* jaOrbFran1.1 (15 chromosomes)

---

## Sequencing QC

| Metric | Value |
|--------|-------|
| Samples passing QC | 140/140 |
| Chromosomes analyzed | 15 |
| SNPs scanned (pass 1) | 28,063,923 |
| SNPs retained (pass 2, MAF > 0.05, ≥ 80% individuals, n=112) | ~7,233,371 |

SNP filter params: `docs/outputs/filter_params.yaml`

---

## Relatedness and Clonality

Pairwise all-vs-all ngsRelate. Clone threshold: KING ≥ 0.45. No close relatives detected.

| Clone pair | KING | Excluded | Reason |
|-----------|------|----------|--------|
| Oann_72 / Oann_73 | 0.4949 | Oann_72 | lower depth (22.9x vs 39.9x) |
| Oann_75 / Oann_77 | 0.4959 | Oann_77 | lower depth (28.9x vs 39.1x) |
| Oann_84 / Oann_85 | 0.4949 | Oann_85 | lower depth (40.2x vs 47.1x) |
| Oann_95 / Oann_96 | 0.4955 | Oann_95 | lower depth (29.8x vs 31.0x) |
| Ofav_106 / Ofra_26 | 0.4929 | **Ofra_26** | **MISLABEL** — confirmed *O. faveolata* |
| Ofav_92 / Ofav_93 | 0.4853 | Ofav_92 | lower depth (29.6x vs 32.8x) |

**134 unrelated samples retained** (57 *O. annularis*, 56 *O. faveolata*, 21 *O. franksi*).

### Ofra_26 mislabel

Ofra_26 (labeled *O. franksi*, FL_RR Upper Keys / East Turtle Shoal) is a confirmed *O. faveolata* mislabel.
Evidence from KING matrix: mean KING vs Ofav=−0.038 (within-species), vs Ofra=−0.521 (between-species),
KING vs Ofav_106=+0.493 (clone — same individual). Decision: exclude Ofra_26, retain Ofav_106.

Full audit: `docs/outputs/relatedness/clone_exclusions.txt`

---

## Population Structure

**Status:** Complete — Segment 3 (134 unrelated samples).

### PCA (PCAngsd)

Covariance matrix: `docs/outputs/pca/pcangsd.cov`

*Plots pending lineage assignment.*

### Admixture

ngsAdmix K=1–10 (20 replicates each) + PCAngsd K=2–10.
Q matrices: `docs/outputs/admixture/`

**Expectation:** K=3 should cleanly separate the three species. Admixture between
*O. annularis* and *O. faveolata* has been documented; K=3 or K=4 may reveal hybrids.

*Delta-K and structure plots pending lineage assignment gate.*

### LD

LD decay: `docs/outputs/ld/ld_decay.csv`
LD prune window: `docs/outputs/ld/ld_prune_window.txt`

---

## Genetic Diversity

*Pending Segment 4.*

---

## FST and Population Differentiation

*Pending Segment 4.*

**Planned comparisons:**

| Comparison | Expected FST | Notes |
|------------|-------------|-------|
| *O. annularis* vs *O. faveolata* | moderate | closest pair; documented hybridization |
| *O. annularis* vs *O. franksi* | higher | more divergent morphology |
| *O. faveolata* vs *O. franksi* | higher | more divergent morphology |

---

## Demographic Inference (moments)

*Pending Segment 6.*

---

## Population Size History (SMC++)

*Pending Segment 7.*

Generation time assumed ~10 yr for massive corals. μ = 1.8×10⁻⁸ (Matz et al.).
**Expectation:** All three species bottlenecked at LGM (~18–20 kya).

---

## Reproducibility

```bash
cd /work/vollmer/orbicella_genomics
bash run.sh 2a   # SNP discovery + relatedness (stops at clone gate)
# → python3 workflow/scripts/clone_approve.py
bash run.sh 2b   # subset BEAGLE to unrelated samples
bash run.sh 3    # LD, PCA, admixture
bash run.sh 4    # SAF, SFS, diversity, FST
```

Pipeline code: `/projects/vollmer/RRorbicella_angsd/`
Reference: jaOrbFran1.1 (`/projects/vollmer/RR_heat-tolerance/Orbicella/reference/`)
