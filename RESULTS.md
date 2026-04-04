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

**Status:** Complete — Segment 3 (134 unrelated samples, 2026-04-04).

Figures: `docs/figures/`

### LD Decay

![LD decay](figures/ld_decay.png)

LD decays from r²≈0.25 at <1 kb to r²≈0.10 by ~8–10 kb, reaching r²≈0.05 by 35 kb and
plateauing near background (~0.05) beyond 50 kb. This is notably **faster than *Acropora***
(which crosses r²=0.1 at ~50–100 kb), consistent with larger historical effective population
sizes in massive corals (*Orbicella* Ne likely 2–5× higher than *Acropora* based on diversity).
The LD pruning window was set accordingly to capture independent SNPs across the genome.

### PCA (PCAngsd)

![PCA by species](figures/pca_species.png)

Covariance matrix: `docs/outputs/pca/pcangsd.cov`

- **PC1 (22.8%)** separates *O. annularis* (negative scores) from *O. faveolata* + *O. franksi*
  (positive scores). The large variance explained reflects deep species-level divergence.
- **PC2 (4.9%)** separates *O. franksi* (high positive) from *O. faveolata* (near zero), with
  *O. annularis* spread along PC2 — suggesting within-species geographic or admixture structure.
- **PC3 (1.8%)** separates within *O. annularis*.
- Notably, *O. faveolata* and *O. franksi* are closer to each other on PC1 than either is to
  *O. annularis* — consistent with phylogenetics (*O. faveolata* and *O. franksi* are sister species).
- Several *O. faveolata* and *O. franksi* individuals plot closer to the opposing species cluster,
  suggesting admixture or past introgression — consistent with documented natural hybridization.

### Admixture

![ngsAdmix delta-K](figures/admixture_deltaK.png)

#### Model selection (ngsAdmix, K=1–10, 20 reps)

| K | Mean log-likelihood | ΔK |
|---|--------------------|----|
| 1 | −4,286,604 | — |
| **2** | **−3,553,251** | **26,440,299 ← best** |
| 3 | −3,415,338 | 132,893 |
| 4 | −3,371,496 | — |
| 5–10 | gradual improvement | low, noisy |

**Best K=2** by Evanno ΔK — overwhelmingly dominant (ΔK at K=2 is ~200× larger than K=3).
This is expected: the deepest signal is *O. annularis* vs the *faveolata*/*franksi* clade.

#### K=2, 3, 4 bar plots

![ngsAdmix K2-4](figures/admixture_ngsadmix_K234.png)
![PCAngsd K2-4](figures/admixture_pcangsd_K234.png)

**K=2:** *O. annularis* forms one cluster (blue); *O. faveolata* + *O. franksi* share the other
(orange). A few *O. annularis* individuals show minor orange ancestry, and vice versa —
consistent with documented *O. annularis × O. faveolata* hybridization.

**K=3:** The three species cleanly separate into distinct clusters, with *O. franksi* resolved
as a third component (green). Several individuals in *O. faveolata* and *O. franksi* show
substantial admixture proportions from the other species — likely hybrids or recent backcrosses.
*O. annularis* remains largely pure at K=3.

**K=4:** Sub-structure emerges within *O. franksi* (splitting into two components: green + pink),
and admixed individuals in *O. faveolata* are more finely resolved. Whether this reflects true
intraspecific population structure within *O. franksi* or overfitting should be assessed after
the lineage assignment gate.

**Interpretation:** The predominant signal is species identity. Admixture between *O. faveolata*
and *O. franksi* is real and visible at K=3 — these two sister species are known to hybridize
in Florida reefs. The Reef Renewal collection samples appear to include some admixed individuals,
which has implications for restoration genetics (admixed colonies may perform differently under
stress). Full hybrid quantification pending the lineage assignment gate.

---

## Genetic Diversity

*Pending Segment 4.*

---

## FST and Population Differentiation

*Pending Segment 4.*

**Planned comparisons:**

| Comparison | Expected FST | Notes |
|------------|-------------|-------|
| *O. annularis* vs *O. faveolata* | moderate–high | documented hybridization may lower FST |
| *O. annularis* vs *O. franksi* | high | more divergent morphology |
| *O. faveolata* vs *O. franksi* | moderate | sister species; hybridization expected |

---

## Demographic Inference (moments)

*Pending Segment 6.*

---

## Population Size History (SMC++)

*Pending Segment 7.*

Generation time assumed ~10 yr for massive corals. μ = 1.8×10⁻⁸ (Matz et al.).
**Expectation:** All three species bottlenecked at LGM (~18–20 kya). *O. faveolata*
post-1980s collapse may appear as terminal Ne decline if resolution sufficient.

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

Plots: `python workflow/scripts/plot_orbicella.py`

Pipeline code: `/projects/vollmer/RRorbicella_angsd/`
Reference: jaOrbFran1.1 (`/projects/vollmer/RR_heat-tolerance/Orbicella/reference/`)
