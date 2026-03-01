# Results

Results from the *Acropora palmata* population genomics analysis (Florida vs. Panama).

> **Status**: 10-sample AWS test run in progress. Full 96-sample results pending.
> Figures will be added here after the production run completes.

---

## 5-Sample Pilot Results

Pilot run using 5 samples (Ac_FL_M10, Ac_FL_M5, Ac_PA_CK140, Ac_PA_HS1, Ac_PA_Tet1).

### Sequencing QC

| Metric | Value |
|--------|-------|
| Sites scanned (pass 1) | 839,000,000 |
| SNPs discovered (pass 1) | 6,514,518 |
| SNPs retained (pass 2, MAF > 0.10, 80% individuals) | 2,062,496 |

### Population Differentiation

| Metric | Florida | Panama |
|--------|---------|--------|
| Tajima's D | +0.41 | +0.29 |

| Population Pair | FST (unweighted) | FST (weighted) |
|----------------|-----------------|----------------|
| Florida vs. Panama | 0.127 | 0.162 |

### Relatedness

All pairwise kinship coefficients below threshold — no clones or first-degree relatives detected in the pilot sample.

---

## 10-Sample AWS Test Run

*Results pending — run in progress.*

Samples: 5 Florida (Ac_FL_K1, Ac_FL_K2, Ac_FL_KW15, Ac_FL_U31, Ac_FL_U39) + 5 Panama (Ac_PA_CK141, Ac_PA_CK142, Ac_PA_Tet5, Ac_PA_Tet6, Ac_PA_HS3)

---

## Full 96-Sample Production Run

*Pending 10-sample test validation.*

### Population Structure

*PCA and admixture plots will be added here.*

### Genetic Diversity

*Per-population π, θ, and Tajima's D table will be added here.*

### FST and Outlier Loci

*Genome-wide FST Manhattan plot and outlier gene table will be added here.*

### Relatedness

*Pairwise kinship heatmap will be added here.*

---

## Reproducibility

All results are fully reproducible from the raw SRA accessions in `config/samples.csv` using:

```bash
snakemake --snakefile workflow/Snakefile --profile profiles/aws \
  --config samples_csv=config/samples.csv
```

Reference genome: *Acropora palmata* assembly (see `config/config.yaml` for path).
