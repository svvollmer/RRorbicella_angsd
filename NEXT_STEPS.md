# Next Steps

Current state as of 2026-03-13. Segmented pipeline built and ready.
96-sample production run completed through diversity/PCA/admixture.
FST and report not yet complete — re-run from Segment 2 with corrected parameters.

---

## Running on FAU KoKo HPC

### One-time setup

```bash
# 1. Edit the SLURM profile with your account name
nano profiles/slurm/config.yaml
# → change FILL_IN_YOUR_ACCOUNT to your KoKo allocation (e.g. ecos)
# → add your scratch path to singularity-args if not under /blue or /scratch

# 2. Stage data from S3 to KoKo scratch (run from KoKo login node)
#    Requires AWS CLI and credentials configured on KoKo
SCRATCH="/blue/yourgroup/yourusername/coral-angsd"   # adjust path
mkdir -p $SCRATCH

# Reference
aws s3 sync s3://coral-angsd-728009587639/reference/ $SCRATCH/reference/

# BAMs (large — ~500 GB, use a batch job or screen session)
aws s3 sync s3://coral-angsd-728009587639/results-production/results/bams/ \
    $SCRATCH/results/bams/ --exclude '*.fastq*' --exclude '*.raw.bam'

# Pass 1 outputs (to skip pass 1 re-run)
aws s3 sync s3://coral-angsd-728009587639/results-production/results/angsd/ \
    $SCRATCH/results/angsd/ \
    --include 'pass1.mafs.gz' --include 'pass1_snps.txt*' \
    --include 'nonrepeat_sites.txt*' --include 'depth_thresholds.txt'

# QC + filtering
aws s3 sync s3://coral-angsd-728009587639/results-production/results/qc/ \
    $SCRATCH/results/qc/
aws s3 sync s3://coral-angsd-728009587639/results-production/results/filtering/ \
    $SCRATCH/results/filtering/

# 3. Clone pipeline and link scratch
git clone https://github.com/svvollmer/coral-angsd-pipeline.git $SCRATCH/pipeline
cd $SCRATCH/pipeline
ln -s $SCRATCH/results results          # point results/ at scratch
ln -s $SCRATCH/reference reference      # point reference/ at scratch
```

### Running on KoKo

```bash
# From the pipeline directory on KoKo:

# SNP filter gate (interactive — run on login node, needs pass1.mafs.gz)
python workflow/scripts/filter_select.py

# Then submit each segment via SLURM profile:
bash scripts/run_segment.sh 2 --profile slurm   # pass 2 + relatedness
# clone gate → re-run
bash scripts/run_segment.sh 3 --profile slurm   # LD/PCA/admixture
# lineage gate → run_segment.sh 4
bash scripts/run_segment.sh 4 --profile slurm   # SAF/SFS/FST
bash scripts/run_segment.sh 5 --profile local   # report (login node is fine)
```

### KoKo notes
- Edit `config/config.yaml`: comment out `local_conda_env` line (Singularity handles tools)
- Segment 4 (SAF) may need `himem` partition — if jobs OOM, add
  `slurm_partition: "himem"` to the SAF rule's `resources:` block
- Gate scripts (`filter_select.py`, `clone_approve.py`, `lineage_assign.py`)
  are interactive — run on login node, not in a job

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
| realSFS -maxIter 200 | FST 2D only | replaced with -nSites 20000000 (2026-03-14) | Snakefile.diversity realsfs_2d |
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
| `workflow/Snakefile.demography` | Segment 6: demographic inference (Gate 4+5, SAF, SFS, dadi/moments) |
| `scripts/run_segment.sh` | Launcher: S3 sync + snakemake + auto-stop |
| `workflow/scripts/filter_select.py` | SNP Filter Gate: select MAF + minInd |
| `workflow/scripts/lineage_assign.py` | Gate 3: select K, assign lineages |

---

## Remaining analysis work

### nonrepeat_sites.txt (Bug 9 — RESOLVED)
SAF rules correctly use `nonrepeat_sites.txt` (monomorphic + variable sites).
The `-rf chromosomes.bed` interaction bug was fixed by removing `-rf` and relying
on the binary sites index alone. Do NOT revert to `pass1_snps.txt` — SNP-only SAF
inflates π and biases Tajima's D. Confirmed by Gemini pipeline review 2026-03-14.

### Segment 6: demographic inference (Snakefile.demography)
Skeleton complete (Gate 4, Gate 5, SAF, 1D SFS, dadi-format SFS). Pending:
- `workflow/scripts/run_moments_2pop.py` — 2-pop dadi/moments model fitting
- `workflow/scripts/run_moments_4pop.py` — 4-pop moments model fitting
- Add Segment 6 to `scripts/run_segment.sh` launcher
- Requires BAMs synced from HPC/S3 to AWS, run on c6i.16xlarge (~4-6 hrs Phase 1)

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
