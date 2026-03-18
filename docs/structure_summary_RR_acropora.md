# Population Structure Summary — RR Acropora Dataset
**Run:** acropora_genomics (Discovery HPC, 2026-03-17)
**Samples:** 290 total → 253 unrelated after clone filtering
**Reference:** *A. palmata* jaAcrPala1.3
**SNPs:** 29,127 LD-pruned sites used for PCA/admixture (MAF > 0.05)

---

## PCA Eigenvalues

| PC | Eigenvalue | % Variance | Cumulative % |
|----|-----------|-----------|-------------|
| PC1 | 148.99 | 45.1% | 45.1% |
| PC2 | 10.23 | 3.1% | 48.2% |
| PC3 | 5.91 | 1.8% | 50.0% |
| PC4 | 3.96 | 1.2% | 51.2% |
| PC5 | 2.57 | 0.8% | 52.0% |
| PC6 | 1.37 | 0.4% | 52.4% |
| PC7 | 1.32 | 0.4% | 52.8% |
| PC8 | 1.16 | 0.3% | 53.2% |

PC1 captures the interspecific split (*A. palmata* vs *A. cervicornis*). PCs 2–5 carry within-species geographic/source structure; PC6+ is noise floor (~1.0–1.4).

---

## K=2 — Species Split

Cleanly separates *A. palmata* (lineageA) from *A. cervicornis* (lineageB). Zero admixed individuals.

| Group | n | Q_Apal (mean) | Q_Apal (min–max) | Assignment |
|-------|---|--------------|-----------------|-----------|
| Apal (FL_RR, BON, PA) | 105 | 0.991 | 0.976–1.000 | lineageA |
| Acer CRF (FL_CRF, FL_MOTE, FL_FWC) | 36 | 0.002 | 0.000–0.012 | lineageB |
| Acer wild FL+PA (Ac_FL, Ac_PA) | 75 | 0.016 | 0.000–0.049 | lineageB |
| Acer RR+BON (RR_FL_Acer, BON_Acer) | 35 | 0.017 | 0.000–0.034 | lineageB |
| Acer PA outgroup (AC_PA_CK14_8, AC_PA_Tetas_10) | 2 | 0.038 | 0.031–0.044 | lineageB |

**Interpretation:** K=2 is the relevant split for hybridization questions (Apal × Acer). No hybrid individuals detected in this dataset.

---

## K=3 — Within-Acer Geographic/Source Structure

Apal remains monolithic (Q3 ≈ 0.99). Within *A. cervicornis*, two clusters emerge.

| Group | n | Q1 (CRF) | Q2 (RR/BON) | Q3 (Apal) |
|-------|---|---------|------------|---------|
| Apal | 105 | 0.005 | 0.005 | **0.990** |
| Acer CRF (nursery/restoration) | 36 | **0.948** | 0.047 | 0.005 |
| Acer wild FL+PA | 75 | 0.537 | 0.456 | 0.007 |
| Acer RR+BON | 35 | 0.312 | **0.686** | 0.002 |

**Interpretation:** K=3 reveals structure within *A. cervicornis*: CRF restoration stock (FL_CRF, FL_MOTE) clusters separately from Reef Renewal + Bonaire samples. Wild FL+PA samples are intermediate (~50/50), suggesting they draw from both gene pools or represent the ancestral source.

---

## K=4 and K=5

K=4 resolves CRF restoration as a distinct cluster (Q4 ≈ 0.93) and further separates RR+BON from wild FL+PA. K=5 begins fragmenting Apal with no clear biological signal — K=4 is the interpretable stopping point.

---

## Lineage Assignments (used downstream)

Written to `results/admixture/lineage_assignments.txt` (K=2, Q threshold = 0.80):
- **lineageA** (*A. palmata*): 105 samples
- **lineageB** (*A. cervicornis*): 148 samples
- **admixed**: 0 samples

---

## Demography Comparisons (Segment 6)

Groups restricted geographically to avoid Wahlund effect in SFS (e.g. `apal_fl` = FL_RR Apal only, not all Apal).

### K2 path — Hybridization / Apal–Acer species divergence
Models: SI, IM, IM_a, SC (secondary contact), AM (ancient migration)

| Comparison | Pop 1 | n | Pop 2 | n |
|-----------|-------|---|-------|---|
| `apal_vs_acer_fl` | Apal FL (FL_RR) | 76 | Acer FL (all FL pops) | 85 |
| `apal_vs_acer_bon` | Apal BON | 24 | Acer BON | 24 |

### K3 path — Within-Acer geographic structure
Models: SI, IM, IM_a, SC, AM

| Comparison | Pop 1 | n | Pop 2 | n |
|-----------|-------|---|-------|---|
| `acer_fl_vs_pa` | Acer FL (all) | 85 | Acer PA | 39 |
| `acer_fl_vs_bon` | Acer FL (all) | 85 | Acer BON | 24 |
| `acer_pa_vs_bon` | Acer PA | 39 | Acer BON | 24 |
