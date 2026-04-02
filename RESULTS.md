# Results

Three-species *Orbicella* population genomics — *O. annularis*, *O. faveolata*, and *O. franksi*
from the Reef Renewal Florida restoration program, aligned to the *O. franksi* reference genome
(jaOrbFran1.1 / GCA_964199315.1).

> **Status (2026-04-02):** Segment 1 running. Results pending.

---

## Samples

| Species | N (target) | Input | Status |
|---------|-----------|-------|--------|
| *O. annularis* | 63 | 46 new FASTQ + 17 old BAM | Seg 1 running (46 new) |
| *O. faveolata* | 55 | 46 new FASTQ + 9 old BAM | Seg 1 running (46 new) |
| *O. franksi* | 22 | old BAM only | Pending ingest |

**Reference:** *O. franksi* jaOrbFran1.1 (14 chromosomes; all scaffolds ≥ 10 Mb retained)

---

## Sequencing QC

*Pending Segment 1 completion.*

---

## Relatedness and Clonality

*Pending Segment 2.*

---

## Population Structure

*Pending Segment 3.*

**Expectation:** K=3 should cleanly separate the three species. Admixture between
*O. annularis* and *O. faveolata* has been documented; K=3 or K=4 may reveal hybrid
individuals. PCAngsd and NGSAdmix will both be run for cross-validation.

---

## Genetic Diversity

*Pending Segment 4.*

---

## FST and Population Differentiation

*Pending Segment 4.*

**Planned comparisons (all three inter-species pairs):**

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

**Expectation:** All three species bottlenecked at LGM (~18–20 kya). *O. faveolata*
post-1980s collapse may appear as a terminal Ne decline if resolution is sufficient.
Generation time assumed ~10 yr for massive corals.

---

## Reproducibility

```bash
cd /work/vollmer/orbicella_genomics
bash run.sh 1    # BAM ingestion + QC
bash run.sh 2a   # SNP discovery + relatedness (stops at clone gate)
bash run.sh 2b   # subset BEAGLE to unrelated samples
bash run.sh 3    # PCA + admixture
bash run.sh 4    # SAF + SFS + diversity + FST
```

Pipeline code: `/projects/vollmer/RRorbicella_angsd/`
Reference: jaOrbFran1.1 (`/projects/vollmer/RR_heat-tolerance/Orbicella/reference/`)
