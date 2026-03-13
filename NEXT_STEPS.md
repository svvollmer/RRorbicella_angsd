# Next Steps

Current state as of 2026-03-13. Segmented pipeline built and ready.
96-sample production run completed through diversity/PCA/admixture.
FST and report not yet complete — re-run from Segment 2 with corrected parameters.

---

## To resume the 96-sample run

All 96 CRAMs are in S3. Segment 1 is done. Start from Segment 2:

```bash
# 1. Launch c6i.16xlarge on-demand, run aws_setup.sh
bash scripts/aws_setup.sh

# 2. Run Segment 2 (pass 1 ~2h, then SNP filter gate, pass 2 ~30 min, relatedness)
bash scripts/run_segment.sh 2

# 3. SNP Filter Gate — runs automatically after pass 1, then:
python workflow/scripts/filter_select.py
# → shows SNP count table at MAF × minInd grid
# → recommend: MAF 0.05, minInd 0.90 (Science paper match)
# → writes results/angsd/filter_params.yaml
# → re-run: bash scripts/run_segment.sh 2

# 4. Clone Gate — after ngsRelate:
python workflow/scripts/clone_approve.py
# → re-run: bash scripts/run_segment.sh 2  (subsets beagle)

# 5. Segment 3: LD, PCA, admixture (c6i.8xlarge)
bash scripts/run_segment.sh 3

# 6. Lineage Gate — after admixture:
python workflow/scripts/lineage_assign.py
# → shows log-likelihoods for K=2..5, select K
# → shows Q table, confirm/adjust threshold, writes lineage_assignments.txt

# 7. Segment 4: SAF, SFS, diversity, FST (c6i.16xlarge)
bash scripts/run_segment.sh 4

# 8. Segment 5: report (run locally)
bash scripts/run_segment.sh 5 --profile local
```

---

## Parameter corrections applied in segmented pipeline

| Parameter | Old | New | Location |
|-----------|-----|-----|----------|
| min_maf | 0.10 | user-selected via filter gate | filter_params.yaml |
| min_ind_frac | 0.80 | user-selected via filter gate | filter_params.yaml |
| -setMinDepthInd 2 | missing from pass2 | added | Snakefile.snps angsd_gl_snps |
| realSFS -maxIter 200 | missing | added | Snakefile.diversity realsfs_2d |
| individual_saf -anc | missing (Bug 21) | added | Snakefile.diversity individual_saf |

---

## Segmented pipeline files

| File | Purpose |
|------|---------|
| `workflow/common.smk` | Shared config/samples/helpers |
| `workflow/Snakefile.align` | Segment 1: download, trim, map, QC |
| `workflow/Snakefile.snps` | Segment 2: pass1, filter gate, pass2, relatedness |
| `workflow/Snakefile.structure` | Segment 3: LD, PCA, admixture |
| `workflow/Snakefile.diversity` | Segment 4: SAF, SFS, diversity, lineage FST |
| `workflow/Snakefile.report` | Segment 5: annotation, HTML report |
| `scripts/run_segment.sh` | Launcher: S3 sync + snakemake + auto-stop |
| `workflow/scripts/filter_select.py` | SNP Filter Gate: select MAF + minInd |
| `workflow/scripts/lineage_assign.py` | Gate 3: select K, assign lineages |

---

## Remaining analysis work

### nonrepeat_sites.txt + -rf interaction (Bug 9 workaround still active)
SAF rules use `pass1_snps.txt` instead of full nonrepeat sites for SAF.
For publication, debug ANGSD 0.940 behavior with large binary sites index
+ BED-format -rf, then switch back to `nonrepeat_sites.txt` for unbiased
SFS/diversity. Monomorphic sites are required for accurate π/θ/D.

### HTML report
`generate_report.py` updated to handle multiple FST comparisons and lineage
assignment figure. Needs FST results to fully render.

### Chr 14 outlier investigation
NC_133895.1 has highest Tajima's D in both populations (FL D=+1.10, PA D=+0.74)
and elevated per-site π. Check gene content and repeat structure — candidate
sex chromosome or major introgression region.

### Within-lineage diversity
Now built into Segment 4 via `fst_comparisons` config — thetas are computed
for all FST groups including lineageA_florida, lineageA_panama etc.

---

## FST comparisons (config/config.yaml)

```yaml
fst_comparisons:
  - [florida, panama]                       # overall geographic
  - [lineageA_florida, lineageA_panama]     # geography within pure lineage
  - [lineageB_florida, lineageB_panama]     # geography within hybrid lineage
  - [lineageA_all, lineageB_all]            # lineage divergence pooled
```

Approximate group sizes from 96-sample run (may shift with corrected params):

| Group | N |
|-------|---|
| FL lineageA | ~9 |
| FL lineageB | ~30 |
| FL admixed | ~1 |
| PA lineageA | ~22 |
| PA lineageB | ~9 |
| PA admixed | ~4 |

FL lineageA (n≈9) is the smallest group — flag low N in paper.
