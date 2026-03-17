# Results

Results from the *Acropora cervicornis* population genomics analysis (Florida vs. Panama), aligned to the *A. palmata* reference genome.

> **Status**: 10-sample AWS test run in progress. Full 96-sample results pending.
> Figures will be added here after the production run completes.

---

## 5-Sample Pilot Results

Pilot run using 5 *A. cervicornis* samples (Ac_FL_M10, Ac_FL_M5, Ac_PA_CK140, Ac_PA_HS1, Ac_PA_Tet1).

### Sequencing QC

| Metric | Value |
|--------|-------|
| Sites scanned (pass 1) | 839,000,000 |
| SNPs discovered (pass 1) | 6,514,518 |
| SNPs retained (pass 2, MAF > 0.10, 80% individuals) | 2,062,496 |

~0.8% of scanned sites are variable — typical for a within-species coral resequencing dataset. The ~3-fold reduction from pass 1 to pass 2 reflects tighter quality filters (MAF > 0.10, 80% individual coverage), retaining the most informative and well-genotyped SNPs for population inference.

### Population Differentiation

| Population Pair | FST (unweighted) | FST (weighted) |
|----------------|-----------------|----------------|
| Florida vs. Panama | 0.127 | 0.162 |

FST of ~0.13–0.16 indicates strong genetic differentiation between Florida and Panama populations. For context, FST values in the range of 0.05–0.15 are typically interpreted as moderate-to-high differentiation in marine invertebrates with planktonic larval dispersal. Values above 0.10 in a coral species with ~5–10 day larval competency windows are consistent with effective isolation across the ~2,000 km Florida–Panama distance, likely maintained by the Caribbean current system limiting northward larval transport from Panama to Florida. The weighted FST (0.162) accounts for variation in sample size across sites and is the more conservative estimate.

### Genetic Diversity

| Metric | Florida | Panama |
|--------|---------|--------|
| Tajima's D | +0.41 | +0.29 |

Both populations show positive Tajima's D, indicating an excess of intermediate-frequency alleles relative to the neutral expectation. Positive values can reflect:
- **Recent population bottleneck** followed by recovery — rare alleles are lost, intermediate alleles remain
- **Balancing selection** maintaining variation at certain loci
- **Population structure** within the sampled region creating apparent excess heterozygosity

Positive Tajima's D has been reported in other Caribbean coral populations and may reflect historical bottlenecks from bleaching events or disease outbreaks. Florida's slightly higher D (+0.41 vs +0.29) could indicate stronger or more recent demographic contraction, consistent with documented population declines in the Florida Keys reef tract. These patterns should be interpreted cautiously with only 5 pilot samples — the full 96-sample run will provide more robust estimates.

### Relatedness

All pairwise kinship coefficients below threshold — no clones or first-degree relatives detected in the pilot sample.

Absence of clones is expected given that samples were collected from visually distinct colonies across reef sites. The low relatedness confirms that the sampling design successfully captured unrelated individuals, which is important for population structure analyses that assume independence among samples. The full 96-sample analysis will include a formal relatedness screen to flag any cryptic duplicates or close kin before downstream analyses.

---

## 10-Sample AWS Test Run

*Results pending — run in progress.*

Samples: 5 Florida (Ac_FL_K1, Ac_FL_K2, Ac_FL_KW15, Ac_FL_U31, Ac_FL_U39) + 5 Panama (Ac_PA_CK141, Ac_PA_CK142, Ac_PA_Tet5, Ac_PA_Tet6, Ac_PA_HS3)

---

## Full 290-Sample Production Run

*Segment 2 (SNP discovery, relatedness) complete. Segments 3–6 in progress.*

Dataset: 290 samples — *A. cervicornis* and *A. palmata* across Florida, Panama, and Bonaire.
SNPs: 2,034,805 passing filters (MAF > 0.05, ≥ 80% individuals genotyped, 14 chromosomes).

### Relatedness and Clonality

Pairwise kinship estimated with ngsRelate (KING coefficient) from 100,000 randomly subsampled
high-MAF (> 0.30) SNPs. Clone threshold: KING ≥ 0.45. Close-relative threshold: KING ≥ 0.20.

One sample (RR_FL_Apal_015) was excluded as a technical artifact — it appeared as a close
relative to >60 other samples across populations and species, a pattern inconsistent with
biology and indicative of a mapping or contamination artifact.

| Group | N | Clone pairs | Excluded | N genets | Genet diversity |
|-------|---|-------------|----------|----------|-----------------|
| *Acervicornis* BON | 25 | 0 | 0 | 25 | 1.000 |
| *Acervicornis* FL | 104 | 28 | 13 | 91 | 0.875 |
| *Acervicornis* PA | 49 | 32 | 10 | 39 | 0.796 |
| *Apalmata* BON | 25 | 0 | 0 | 25 | 1.000 |
| *Apalmata* FL | 88 | 11 | 11 | 77 | 0.875 |
| *Apalmata* PA | 9 | 3 | 2 | 7 | 0.778 |
| **Total** | **300** | **74** | **36** | **264** | **0.880** |

Zero cross-species or cross-population clone pairs were detected, confirming clonality is
entirely within-group as expected biologically.

**Key patterns:**
- Bonaire collections show 100% genet diversity in both species, consistent with wild reef sampling.
- Panama has the lowest genet diversity (79–78%), consistent with nursery/aquaculture sourcing
  (Ac_PA_HS1–8 are a single genotype sampled 8 times; Ac_PA_CK41/42/44 are a 3-sample clone cluster).
- Florida shows intermediate clonality (12.5% both species), reflecting reef restoration nursery
  collections (CRF, MOTE) where genotypes are propagated clonally.
- 253 unrelated samples retained for downstream population structure and diversity analyses.

### Population Structure

*PCA and admixture results pending (Segment 3).*

### Genetic Diversity

*Per-population π, θ, and Tajima's D pending (Segment 4).*

### FST and Outlier Loci

*Genome-wide FST and outlier annotation pending (Segment 4).*

---

## Reproducibility

All results are fully reproducible from the raw SRA accessions in `config/samples.csv` using:

```bash
snakemake --snakefile workflow/Snakefile --profile profiles/aws \
  --config samples_csv=config/samples.csv
```

Samples: *Acropora cervicornis* WGS data (Vollmer Lab, *Science*).
Reference genome: *Acropora palmata* assembly used as mapping reference (see `config/config.yaml`). SNP coordinates are in *A. palmata* reference space pending liftover to a native *A. cervicornis* assembly.
