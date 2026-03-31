# Population Genomics Results — *Acropora* RR Dataset

300 samples (290 unrelated) spanning *A. palmata* and *A. cervicornis* across Florida, Panama, and Bonaire.
See [structure_summary_RR_acropora.md](structure_summary_RR_acropora.md) for detailed tables.

---

## 1. PCA

PC1 (45.1%) separates species; PC2 (3.1%) captures within-*A. cervicornis* geographic structure.

![PCA](figures/pca.png)

---

## 2. Admixture (K=2–5)

Bar plots sorted by species (Apal first) then region (FL → PA → BON). Violin plots show
per-group ancestry distributions (species × region) for K=2 and K=3.

![PCAngsd admixture K2–K5](figures/admixture_K2_K5.png)
![PCAngsd violin K2 and K3](figures/admixture_violin_pcangsd.png)
![NGSAdmix admixture K2–K5](figures/ngsadmix_admixture_K2_K5.png)
![NGSAdmix violin K2 and K3](figures/admixture_violin_ngsadmix.png)
![PCAngsd vs NGSAdmix violin comparison](figures/admixture_violin_compare.png)

---

## 3. Lineage Assignment (K=2)

K=2 cleanly separates species: lineageA = *A. palmata* (105), lineageB = *A. cervicornis* (148), 0 admixed.

![Lineage Assignment](figures/lineage_assignment.png)

---

## 4. Individual Heterozygosity

Per-sample heterozygosity grouped by species × population. Mean ± 2 SD shown.

![Heterozygosity](figures/heterozygosity.png)

---

## 5. Pairwise Kinship (KING)

Kinship heatmaps split by species, samples ordered by population within each species.

![Kinship](figures/kinship.png)

---

## 6. LD Decay

LD decay estimated with ngsLD on 29,127 LD-pruned SNPs.

![LD Decay](figures/ld_decay.png)

---

## 7. Site Frequency Spectra

Folded 1D SFS per population. *(Pending Segment 4 completion)*

![SFS](figures/sfs.png)

---

## 8. Genetic Diversity (π, θ, Tajima's D)

*(Pending Segment 4 completion)*

![Diversity](figures/diversity.png)

---

## 9. FST

*(Pending Segment 4 completion — FST Manhattan plots will appear here)*

---

## Demography

### Seg 6 — Two-population moments inference (COMPLETE)

Five comparisons fit with 5 models (SI, IM, IM_a, AM, SC), 20 restarts each.
μ = 1.8×10⁻⁸/site/gen, generation time = 5 yr, L = 217,341,919 callable sites.

**Migration convention:** m[i,j] = FROM pop j INTO pop i (moments forward-time diffusion).
In all tables below, m12 = FROM pop2 INTO pop1; m21 = FROM pop1 INTO pop2.

#### Model selection

| Comparison | Best model | ΔAIC vs 2nd | Interpretation |
|---|---|---|---|
| apal_vs_acer_fl | **IM_a** | 378 | Ongoing asymmetric gene flow; apal→acer dominant |
| apal_vs_acer_bon | **AM** | 463 | Ancient gene flow, now ceased; interspecific Bonaire contact is historical |
| acer_fl_vs_bon | **IM_a** | 1156 | Strong ongoing FL→BON dispersal |
| acer_fl_vs_pa | **IM_a** | 84 | Ongoing FL→PA dispersal |
| acer_pa_vs_bon | **IM_a** | 106 | Ongoing PA→BON dispersal |

IM_a wins 4/5 comparisons; AM wins the only interspecific Bonaire comparison.
Secondary contact (SC) never wins, arguing against recent secondary contact as the primary model.

#### Parameter summary (physical units)

N_anc estimated from Watterson's θ on the projected 2D SFS. Ne values relative to that N_anc.

| Comparison | N_anc | Ne(pop1) | Ne(pop2) | T_split (yr) | m12 (/gen) | m21 (/gen) | Ratio m21/m12 |
|---|---|---|---|---|---|---|---|
| apal_vs_acer_fl | 190 | 305 (apal_fl) | 183 (acer_fl) | ~8,250 | 2.0×10⁻⁴ (acer→apal) | 1.4×10⁻³ (apal→acer) | **~7×** |
| apal_vs_acer_bon | 1,053 | 76 (apal_bon) | 83 (acer_bon) | ~3,975 (T1+T2) | 3.9×10⁻³ (symm.) | — | — |
| acer_fl_vs_bon | 927 | 332 (acer_fl) | 121 (acer_bon) | ~1,090 | 1.2×10⁻³ (BON→FL) | 1.8×10⁻² (FL→BON) | **~15×** |
| acer_fl_vs_pa | 94 | 15 (acer_fl) | 13 (acer_pa) | ~345 | 4.7×10⁻² (PA→FL) | 2.3×10⁻¹ (FL→PA) | **~5×** |
| acer_pa_vs_bon | 1,302 | 128,691 (acer_pa) | 288 (acer_bon) | ~590 | 4.1×10⁻⁴ (BON→PA) | 2.2×10⁻³ (PA→BON) | **~5×** |

*Note: N_anc from projected SFS Watterson's θ — treat absolute Ne as order-of-magnitude; use SMC++ Ne(t) for publication Ne values.*

#### Key biological findings

1. **Interspecific hybridization (apal × acer):** IM_a best model for FL (apal→acer ~7× acer→apal), AM best for Bonaire (ancient migration, now ceased). Supports ongoing hybridization in Florida but not Bonaire. apal→acer direction is consistent with *A. prolifera* backcrossing into *A. cervicornis*.

2. **Florida as intraspecific source:** In all three intraspecific acer comparisons, Florida disperses *to* PA and BON (m21 >> m12), not the reverse. FL→BON ~15×, FL→PA ~5×, PA→BON ~5×. Consistent with Florida as a regional source population under prevailing Caribbean circulation.

3. **Ancient vs. ongoing interspecific contact:** The Bonaire apal–acer comparison (AM) shows gene flow ceased, while Florida (IM_a) shows it continues — geographic differentiation in hybrid zone dynamics.

### Seg 7 — SMC++ effective population size (COMPLETE)

Florida-only high-purity samples (Q ≥ 0.99): apal n=49, acer n=82.
μ = 1.8×10⁻⁸, generation time = 5 yr.

| Time (ya) | *A. palmata* Ne | *A. cervicornis* Ne |
|---|---|---|
| Present | ~28,700 | ~12,300 |
| LGM bottleneck (~17–18 kya) | ~8,200 | ~7,000 |
| Ancient peak (~170–213 kya) | ~73,300 | ~95,100 |

Both species bottleneck at LGM (~17–18 kya); apal recovered 4× post-LGM vs acer 1.7×.
Ancient Ne: acer > apal — acer was historically the larger/more diverse species.

### SMC++ split (species divergence time)

`smc++ split` was attempted using the apal and acer Florida Ne(t) models with all 14 chromosomes
(3 distinguished pairs each). It failed with `OverflowError: math range error` — the apal/acer
divergence (~2–4 Mya) is beyond SMC++'s estimable time window, which saturates at ~500 kya.
Species divergence time is covered by the moments IM_a T parameter for `apal_vs_acer_fl` (~8,250 yr
scaled; interpret cautiously given small N_anc from projected SFS). Literature places apal/acer
divergence at ~1–3 Mya; SMC++ split is not a suitable tool for this pair.

---

## FST Outlier Genes and GO Enrichment (acer vs apal, Florida)

FST outlier windows (top 1%, sliding window) were intersected with the *A. palmata* jaAcrPala1.3
gene annotation. Significant GO terms tested by Fisher's exact test (goatools, FDR-BH < 0.05).

**Files:** `docs/selection/acer_vs_apal_FL.genes.tsv`, `acer_vs_apal_FL_GO.go_sig.tsv`

### Outlier summary

- 324 FST outlier windows (FST > ~0.94 weighted)
- 705 genes overlapping outlier windows
- 293 significant GO terms (FDR-BH < 0.05)

### Top FST outlier genes (by FST)

| Gene | Description | FST |
|---|---|---|
| LOC141897654 | CD59 glycoprotein-like | 0.958 |
| LOC141863290 | E3 ubiquitin-protein ligase RNF113A-like | 0.953 |
| LOC141873931 | CCR4-NOT transcription complex subunit 6-like | 0.951 |
| LOC141863732 | Helicase ARIP4-like | 0.949 |
| LOC141873566 | Ras-related and estrogen-regulated growth inhibitor-like | 0.949 |
| LOC141874147 | Adenosine receptor A2b-like | 0.947 |
| LOC141877289 | Activating signal cointegrator 1-like | 0.948 |

### Key enriched GO categories (FDR-BH < 0.05, most specific terms)

| GO term | Category | Name | Study/Pop | FDR |
|---|---|---|---|---|
| GO:0007165 | BP | Signal transduction | 63/479 | 3.8e-09 |
| GO:0140110 | MF | Transcription regulator activity | 28/479 | 3.1e-07 |
| GO:0016491 | MF | Oxidoreductase activity | 28/479 | 1.1e-06 |
| GO:0032502 | BP | Developmental process | 24/479 | 1.1e-06 |
| GO:0016301 | MF | Kinase activity | 24/479 | 2.9e-05 |
| GO:0030547 | MF | Signaling receptor inhibitor activity | 5/479 | 3.9e-06 |
| GO:0022836 | MF | Gated channel activity | 17/479 | 1.7e-05 |
| GO:0005216 | MF | Ion channel activity | 21/479 | 2.9e-05 |
| GO:0048856 | BP | Anatomical structure development | 15/479 | 8.8e-05 |
| GO:0005615 | CC | Extracellular space | 18/479 | 2.3e-04 |

### Biological interpretation

FST outliers are highly enriched for genes involved in signal transduction, receptor signaling,
ion channels, and transcription regulation — consistent with species-level differentiation in
cell communication, developmental pathways, and stress response. CD59 glycoprotein (top outlier)
is a complement inhibitor relevant to immune defense. Adenosine receptor A2b is involved in
inflammation and stress signaling. CCR4-NOT complex regulates mRNA stability and transcription.
These categories align with known reproductive isolation and adaptive divergence loci in coral
comparative genomics.
