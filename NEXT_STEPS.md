# Next Steps

Current state and planned work as of 2026-03-11. Pipeline has completed a full 96-sample
production run through diversity/PCA/admixture. FST not yet complete. Major architectural
redesign planned before next production run.

---

## Immediate fixes (apply before next run)

These are small targeted changes to the existing pipeline. Do these first.

### 1. Parameter corrections
```yaml
# config/config.yaml
min_maf: 0.05          # was 0.10 — match Science paper, standard for population genomics
min_ind_frac: 0.90     # was 0.80 — match Science paper, reduces missing data bias
```

```python
# workflow/Snakefile line ~71
MIN_IND_FINAL = int(N_APPROVED * 0.9)   # was 0.8
```

### 2. Add -setMinDepthInd 2 to pass2
Pass 1 and the SAF rules both use `-setMinDepthInd 2` but `angsd_gl_snps` (pass 2) does not.
Without it, a sample with 1 read at a site counts toward `minInd` despite a highly uncertain GL.

```bash
# workflow/Snakefile — angsd_gl_snps shell block, add:
-setMinDepthInd 2 \
```

### 3. Cap realSFS 2D with -maxIter 200
Without this, the 2D SFS EM ran for 20+ hours without converging (default has no iteration cap).
200 iterations is sufficient for FST estimation.

```bash
# workflow/Snakefile — realsfs_2d shell block, change to:
realSFS {input.fl} {input.pa} -P 8 -fold 1 -maxIter 200 > {output}
```

### 4. Fix IAM role for auto-stop and SNS
The `coral-pipeline-ec2` IAM role is still missing reliable `ec2:StopInstances` and
`sns:Publish` permissions in some execution contexts. Verify the policy is attached to
the correct role ARN and test with a short run before the next production job.

---

## Pipeline architecture redesign

The current single-instance design is functional but not optimal for cost, flexibility,
or client use. Planned redesign segments the pipeline into stages with different instance
types and clean stopping points.

### Segment design

```
Segment 1 — Alignment          (embarrassingly parallel, spot-friendly)
Segment 2 — SNP Discovery      (high memory, serial — ANGSD pass1/pass2 + ngsRelate)
Segment 3 — Population Structure (moderate compute — PCA, admixture, LD)
Segment 4 — Diversity & FST    (high memory — SAF, realSFS, thetas, windowed FST)
Segment 5 — Report             (trivial — can run locally or on smallest instance)
```

Each segment:
- Has its own Snakefile (`workflow/Snakefile.align`, `workflow/Snakefile.genotype`, etc.)
- Reads inputs from S3, writes outputs to S3
- Launched by `scripts/run_segment.sh <N> --mode aws|hpc`
- Stops instance (or SLURM job) cleanly on completion or at a human gate

### Two deployment modes

**HPC mode** (lab use — data already on FAU/NU servers):
- Segment 1 runs on SLURM (free allocation, data local)
- CRAMs synced to S3 after alignment
- Segments 2–5 run on AWS (memory requirements exceed typical HPC allocation)

**AWS mode** (client use — fully self-contained):
- Segment 1: spot fleet of c6i.2xlarge instances (~$4 for 96 samples vs ~$25 on single large)
  - Samples distributed via SQS queue; each worker pulls, aligns, pushes CRAM, terminates
- Segments 2–5: on-demand instances appropriate to memory needs

### Estimated AWS cost (96 samples, full pipeline)
| Segment | Instance | Time | Cost |
|---|---|---|---|
| Alignment (16× spot) | c6i.2xlarge spot | ~3 hrs | ~$4 |
| SNP discovery + genotyping | c6i.16xlarge on-demand | ~6 hrs | ~$27 |
| Population structure | c6i.8xlarge on-demand | ~3 hrs | ~$7 |
| Diversity + FST | c6i.16xlarge on-demand | ~4 hrs | ~$18 |
| S3 storage (1 month) | | | ~$5 |
| **Total** | | **~16 hrs** | **~$61** |

### Human gates (redesigned)
Gates currently use local file checks that break on instance restart.
Redesigned gates check S3 for approval files and send SNS email with presigned report URL.

```
Gate 1 (QC):    pipeline stops → SNS email with QC report link
                user runs: python scripts/gates/qc_approve.py --s3
                writes samples_approved.txt to S3
                user starts instance → pipeline resumes

Gate 2 (Clone): pipeline stops → SNS email with clone report link
                user runs: python scripts/gates/clone_approve.py --s3
                writes unrelated_samples.txt to S3
                user starts instance → pipeline resumes

Gate 3 (Lineage): NEW — after PCA/admixture
                pipeline stops → SNS email with admixture Q plots
                user runs: python scripts/gates/lineage_assign.py --q-threshold 0.80
                writes lineage_assignments.txt to S3 (sample → lineageA/lineageB/admixed)
                user starts instance → FST runs with lineage groupings
```

---

## Analysis redesign — lineage and population FST

The K=2 admixture result and *A. cervicornis × A. palmata* hybridization context means
FL vs PA FST alone is not the right comparison. Planned FST design:

### Lineage assignment
After admixture (Segment 3 Gate 3), assign each unrelated sample:
- Q > 0.80 component 1 → Lineage A (pure / predominantly *A. cervicornis*)
- Q > 0.80 component 2 → Lineage B (hybrid / *A. cervicornis × A. palmata*)
- 0.20 ≤ Q ≤ 0.80 → Admixed (excluded from lineage FST; reported separately)

### Pairwise comparisons
Define in config rather than hardcoded:

```yaml
# config/config.yaml
fst_comparisons:
  - [florida, panama]              # overall geographic (Science paper comparison)
  - [lineageA_florida, lineageA_panama]   # geography within pure lineage
  - [lineageB_florida, lineageB_panama]   # geography within hybrid lineage
  - [lineageA_all, lineageB_all]          # lineage divergence pooled

lineage_q_threshold: 0.80
```

Each comparison:
- Gets its own SAF (per group bamlist)
- Gets its own 1D realSFS + thetas (diversity within group)
- Gets its own 2D realSFS (with -maxIter 200) + FST stats + windowed FST

### Approximate group sizes from 96-sample run
| Group | N | Notes |
|---|---|---|
| FL lineage A | ~9 | FL_B/M series |
| FL lineage B | ~30 | FL_U series |
| FL admixed | ~1 | FL_U75 |
| PA lineage A | ~22 | SR, Tet, HS9/10 |
| PA lineage B | ~9 | CK43–49, CK410 |
| PA admixed | ~4 | CK140, CK1410, CK142, CK144 |

FL lineage A (n=9) is the smallest group — workable but flag low N in paper.

---

## Input flexibility — local FASTQs and S3

The pipeline currently assumes SRA accessions. Needs to accept raw FASTQs already
on HPC servers or uploaded to S3 by clients.

### samples.csv extension
```csv
sample_id,population,input_type,fastq_r1,fastq_r2,sra_accession
Ac_FL_M10,florida,sra,,,SRR24007640
Ac_PA_HS1,panama,local,/lab/data/Ac_PA_HS1_R1.fastq.gz,/lab/data/Ac_PA_HS1_R2.fastq.gz,
Client_01,florida,s3,s3://client-bucket/Client_01_R1.fastq.gz,s3://client-bucket/Client_01_R2.fastq.gz,
```

### Three ingestion rules → one downstream path
```
download_sra   (input_type == sra)   → results/fastq/{sample}_1.fastq.gz
stage_local    (input_type == local) → results/fastq/{sample}_1.fastq.gz
stage_s3       (input_type == s3)    → results/fastq/{sample}_1.fastq.gz
                                           ↓
                                    fastp → BWA → CRAM (unchanged)
```

---

## Other improvements flagged

- **ANGSD nonrepeat_sites.txt + -rf interaction** (Bug 9 workaround): SAF rules currently
  use pass1_snps.txt instead of the full nonrepeat sites file. For production, debug the
  ANGSD 0.940 behavior with large binary sites index + BED-format -rf and switch back to
  nonrepeat_sites.txt for SAF. Monomorphic sites are needed for unbiased SFS/diversity.

- **Report**: HTML report never generated for 96-sample run. Needs FST to complete first.
  Once FST redesign is in place, update `generate_report.py` to handle multiple FST
  comparisons and add the lineage assignment figure.

- **Chr 14 outlier investigation**: NC_133895.1 has the highest Tajima's D in both
  populations (FL D=+1.10, PA D=+0.74) and elevated per-site π. Check gene content
  and repeat structure — candidate sex chromosome or major introgression region.

- **Within-lineage diversity**: Rerun thetas/Tajima's D with lineage-specific bamlists
  to get unconfounded diversity estimates (current values are inflated by the two-lineage
  admixture structure).
