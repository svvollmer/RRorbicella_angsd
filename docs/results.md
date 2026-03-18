# Population Genomics Results — *Acropora* RR Dataset

300 samples (290 unrelated) spanning *A. palmata* and *A. cervicornis* across Florida, Panama, and Bonaire.
See [structure_summary_RR_acropora.md](structure_summary_RR_acropora.md) for detailed tables.

---

## 1. PCA

PC1 (45.1%) separates species; PC2 (3.1%) captures within-*A. cervicornis* geographic structure.

![PCA](figures/pca.png)

---

## 2. Admixture (K=2–5)

Faceted by geographic region (FL, PA, BON). Horizontal bars sorted by species then population within each panel.

![Admixture](figures/admixture.png)

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

Five two-population comparisons queued for moments demographic inference (Segment 6):

| Comparison | Groups | Question |
|---|---|---|
| apal_vs_acer_fl | apal_fl vs acer_fl | FL hybridization — IM, SC, or AM? |
| apal_vs_acer_bon | apal_bon vs acer_bon | BON hybridization replicate |
| acer_fl_vs_pa | acer_fl vs acer_pa | FL–Panama divergence |
| acer_fl_vs_bon | acer_fl vs acer_bon | FL–Bonaire divergence |
| acer_pa_vs_bon | acer_pa vs acer_bon | Panama–Bonaire divergence |
