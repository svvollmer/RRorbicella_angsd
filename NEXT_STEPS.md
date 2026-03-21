# Next Steps

Current state as of 2026-03-21. 290-sample Discovery HPC run.
Segments 2–3 complete. Segment 4 (diversity/FST) mostly complete. Moments 2D SFS running.

---

## Current pipeline status

| Segment | Status | Notes |
|---------|--------|-------|
| Seg 2: SNP discovery + relatedness | ✅ Complete | 2,034,805 SNPs; 253 unrelated samples |
| Seg 3: PCA + admixture | ✅ Complete | PCAngsd + NGSAdmix K=1–10, 20 reps |
| Seg 4: Diversity + FST | 🔄 Mostly complete | 3/6 FST comparisons done; FL pairs pending |
| Seg 6: Demography (prototype) | 🔄 2D SFS running | 5 jobs running; watcher auto-submitting moments |

---

## Immediate next steps

### 1. Moments 2D SFS — check tonight / tomorrow morning

```bash
ssh s.vollmer@login.discovery.neu.edu
# Check 2D SFS files
ls -lh /work/vollmer/acropora_genomics/results/demography/sfs/*.2dsfs
# Check watcher log
tail -20 /work/vollmer/acropora_genomics/logs/moments_watcher.log
# Check queue — non-FL jobs at 24h limit, FL jobs with -maxIter 100
squeue -u s.vollmer
```

Non-FL jobs (`acer_pa_vs_bon`, `apal_vs_acer_bon`) hit 24h wall at ~19:30 on 2026-03-21.
If they timed out with 0-byte output, relaunch:
```bash
cd /work/vollmer/acropora_genomics
bash run_2dsfs.sh   # submits acer_pa_vs_bon and apal_vs_acer_bon
```

FL jobs (`acer_fl_vs_acer_pa`, `acer_fl_vs_acer_bon`, `apal_fl_vs_apal_bon`) have `-maxIter 100`
and should finish faster — check output size.

### 2. When moments 2D SFS lands

The watcher auto-submits moments jobs. Verify fits ran:
```bash
ls -lh /work/vollmer/acropora_genomics/results/demography/moments/
tail -20 /work/vollmer/acropora_genomics/logs/moments_watcher.log
```

If watcher missed a job (restart), submit manually:
```bash
# Example for acer_pa_vs_bon (39 PA, 24 BON samples, project to 40/38 haploids)
/home/s.vollmer/.conda/envs/snakemake2/bin/python \
  /work/vollmer/acropora_genomics/scripts/run_moments_2pop.py \
  --sfs /work/vollmer/acropora_genomics/results/demography/sfs/acer_pa_vs_bon.2dsfs \
  --sfs-format realsfs --n-ind 39 24 \
  --pop1 acer_pa --pop2 acer_bon \
  --outdir /work/vollmer/acropora_genomics/results/demography/moments/acer_pa_vs_bon
```

### 3. Remaining FST (Segment 4)

Three FST comparisons still pending (lineageB_FL vs FL/PA/BON pairs for Acer):
```bash
# Check seg4 log
tail -50 /work/vollmer/acropora_genomics/logs/run_seg4.log
# Restart seg4 if needed
bash run.sh 4
```

The inter-species FL comparison (lineageB_FL vs lineageA_FL) has been **skipped** — scientifically
redundant given 0 admixed individuals; repeatedly times out at 24h.

### 4. Update figures and RESULTS.md when Seg4 + moments complete

Once FST and moments results are in hand:
- Pull FST heatmap and diversity plots from HPC
- Update RESULTS.md diversity and FST sections (currently showing "pending")
- Update moments section when fits complete
- Commit and push

---

## Pending analysis decisions

### Admixture interpretation for paper
- **Primary figure**: NGSAdmix K=2 (best by Evanno with K=1 anchor; delta-K = 116M)
- **Secondary figure**: NGSAdmix K=3 (geography within *A. cervicornis*)
- PCAngsd and NGSAdmix fully agree at K=2 (0 admixed in both)
- PCAngsd shows more spread in violin plots — expected (continuous PCA model vs discrete cluster model)

### Inter-lineage FL FST
Skipped from pipeline. The BON species FST (weighted = 0.456) already establishes the species-level
divergence. The FL equivalent is biologically uninformative beyond that.

### Moments demographic inference — prototype → pipeline
Once prototype fits validate on the 5 pairwise comparisons:
1. Add `moments_2pop` rule to `Snakefile.demography`
2. Add Segment 6 to `scripts/run_segment.sh`
3. Run full 290-sample demography on Discovery or AWS

Comparisons planned, with scientific role:

| Comparison | Type | Restarts | Key question |
|------------|------|----------|-------------|
| `apal_vs_acer_bon` | **Inter-species** | 50 | Is gene flow between species rare? SI vs IM AIC |
| `acer_pa_vs_bon` | Intra-species | 20 | Acer PA↔BON geographic gene flow |
| `acer_fl_vs_acer_pa` | Intra-species | 20 | Acer FL↔PA gene flow + split time |
| `acer_fl_vs_acer_bon` | Intra-species | 20 | Acer FL↔BON gene flow + split time |
| `apal_fl_vs_apal_bon` | Intra-species | 20 | Apal FL↔BON gene flow + split time |

**Inter-species gene flow (apal_vs_acer_bon)** is the primary test for the paper narrative.
K=2 admixture shows 0 admixed individuals → expectation is SI model best or very low m in IM.
Migration lower bound set to 1e-5 (not 0.01) so the optimizer can find near-zero solutions.
If SI wins: "isolation is strict — introgression rare or absent."
If IM wins with low m: "rare gene flow quantified — consistent with near-complete reproductive isolation."
BON used (not FL) because BON has 100% genet diversity, no clonal contamination, smallest N → cleanest 2D SFS.

---

## Longer-term analysis work

### Chr 14 outlier investigation
NC_133895.1 has highest Tajima's D in both FL and PA populations in the 96-sample run
(FL D=+1.10, PA D=+0.74) and elevated per-site π. Check gene content and repeat structure —
candidate sex chromosome or major introgression hotspot. Investigate with 290-sample diversity
results once Seg4 FST is complete.

### 4-pop moments model
After 2-pop prototype validates, fit 4-pop model:
`run_moments_4pop.py` (not yet written) for lineageA_FL + lineageB_FL + Acer_PA + Acer_BON.

### HTML report (Snakefile.report)
`generate_report.py` updated to handle multiple FST comparisons and lineage assignment.
Needs full Seg4 FST results to render completely. Run after FST complete:
```bash
bash run.sh 5
```

### Liftover to native A. cervicornis assembly
All SNP coordinates are in *A. palmata* reference space. Liftover to a native *A. cervicornis*
assembly pending — not yet prioritized.

---

## How to check status (quick reference)

```bash
# Queue
squeue -u s.vollmer

# Seg4 FST progress
tail -30 /work/vollmer/acropora_genomics/logs/run_seg4.log

# Moments watcher
tail -20 /work/vollmer/acropora_genomics/logs/moments_watcher.log

# 2D SFS file sizes
ls -lh /work/vollmer/acropora_genomics/results/demography/sfs/*.2dsfs

# NGSAdmix (complete — K=1-10 .Q files all present)
ls /work/vollmer/acropora_genomics/results/admixture/ngsadmix_K*.Q

# Restart seg4 if needed
cd /work/vollmer/acropora_genomics && bash run.sh 4
```

---

## Pipeline files reference

| File | Purpose |
|------|---------|
| `workflow/common.smk` | Shared config/samples/helpers |
| `workflow/Snakefile.align` | Segment 1: download, trim, map, QC |
| `workflow/Snakefile.snps` | Segment 2: pass1, filter gate, pass2, relatedness |
| `workflow/Snakefile.structure` | Segment 3: LD, PCA, NGSAdmix K=1–10 |
| `workflow/Snakefile.diversity` | Segment 4: SAF, SFS, diversity, lineage FST |
| `workflow/Snakefile.report` | Segment 5: annotation, HTML report |
| `workflow/Snakefile.demography` | Segment 6: Gate 4+5, SAF, SFS, dadi/moments |
| `scripts/run_segment.sh` | Launcher: S3 sync + snakemake + auto-stop |
| `workflow/scripts/filter_select.py` | SNP Filter Gate: select MAF + minInd |
| `workflow/scripts/clone_approve.py` | Clone Gate: approve exclusions |
| `workflow/scripts/lineage_assign.py` | Gate 3: select K, assign lineages |
| `workflow/scripts/run_moments_2pop.py` | Moments 2-pop model fitting (prototype) |
| `workflow/scripts/plot_ngsadmix.py` | NGSAdmix delta-K + bar plots |
| `workflow/scripts/plot_admix_compare.py` | PCAngsd vs NGSAdmix comparison panels |
| `workflow/scripts/plot_admix_violin.py` | Admixture violin plots by species×region |

HPC scripts (not yet in pipeline, `/work/vollmer/acropora_genomics/scripts/`):
- `run_2dsfs.sh` — realSFS 2D SFS for non-FL pairs (24h, no maxIter cap)
- `run_2dsfs_fl.sh` — realSFS 2D SFS for FL pairs (-maxIter 100)
- `watch_and_submit_moments.sh` — auto-submits moments fits when 2D SFS files land
