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

## Full 96-Sample Production Run

*Pending 10-sample test validation.*

### Population Structure

*PCA and admixture plots will be added here.*

With 96 samples across Florida and Panama, we expect PC1 to cleanly separate the two regions (consistent with FST ~0.15). Secondary PCs may reveal structure within Florida (Keys vs. offshore reefs) or within Panama (Caribbean vs. Pacific-facing sites if sampled). Admixture at K=2 should recover the FL/PA split; K=3–5 will test for finer substructure or admixed individuals near the range edges.

### Genetic Diversity

*Per-population π, θ, and Tajima's D table will be added here.*

We anticipate nucleotide diversity (π) in the range of 1–5 × 10⁻³ based on published *A. palmata* resequencing studies. Lower π in Florida relative to Panama would be consistent with the documented population collapse in Florida since the 1980s (>90% cover loss). Tajima's D patterns from the pilot (FL > PA) will be tested with greater power in the full dataset.

### FST and Outlier Loci

*Genome-wide FST Manhattan plot and outlier gene table will be added here.*

Genome-wide FST outliers (top 1% of windowed FST) will be annotated against the *A. palmata* gene models. Candidate loci of interest include those involved in thermal tolerance, symbiont recognition (lectin pathway), calcification, and immunity — pathways previously implicated in coral local adaptation. Outlier windows on large chromosomes (≥ 10 Mb) will be prioritized given the scaffold-level assembly quality.

### Relatedness

*Pairwise kinship heatmap will be added here.*

Pairs with KING kinship > 0.20 will be flagged as related (≥ second-degree); pairs > 0.45 classified as clonal. Any flagged pairs will be removed from downstream PCA and FST analyses to avoid inflating structure signals.

---

## Reproducibility

All results are fully reproducible from the raw SRA accessions in `config/samples.csv` using:

```bash
snakemake --snakefile workflow/Snakefile --profile profiles/aws \
  --config samples_csv=config/samples.csv
```

Samples: *Acropora cervicornis* WGS data (Vollmer Lab, *Science*).
Reference genome: *Acropora palmata* assembly used as mapping reference (see `config/config.yaml`). SNP coordinates are in *A. palmata* reference space pending liftover to a native *A. cervicornis* assembly.
