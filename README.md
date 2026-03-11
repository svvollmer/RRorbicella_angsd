# Coral ANGSD Pipeline

Population genomics pipeline for *Acropora cervicornis* (staghorn coral) using a two-pass ANGSD workflow. Designed to run on a local workstation, SLURM HPC, or AWS EC2 with Singularity containers.

**Vollmer Lab — Florida Atlantic University**

---

## Study System

96 *A. cervicornis* whole-genome resequencing samples from the Vollmer Lab's *Science* paper, spanning two regions:
- **Florida** (FL): Florida Keys reef system
- **Panama** (PA): Caribbean Panama

**Reference genome:** The *Acropora palmata* genome is used as the mapping reference, as it is higher quality than currently available *A. cervicornis* assemblies. *A. palmata* and *A. cervicornis* are sister species (hybridize to form *A. prolifera*) with high synteny, making cross-species mapping appropriate for population genomic inference.

**Future liftover:** SNP coordinates are currently in *A. palmata* reference space. When a chromosome-scale *A. cervicornis* assembly becomes available, a liftover will be needed to translate positions into native *A. cervicornis* coordinates. This dataset may also serve as a reference panel for future *A. cervicornis* population genomic work.

Research questions: population structure, gene flow, local adaptation, and relatedness across the Florida–Panama range.

---

## Pipeline Strategy

### Why two-pass ANGSD?

ANGSD works on genotype likelihoods rather than hard genotype calls, making it well-suited for low-to-moderate coverage coral resequencing data. A single-pass approach requires choosing SNP filters before knowing the data; a two-pass approach lets the data inform the filters:

1. **Pass 1 — Discovery**: Scan the whole genome with relaxed filters (60% of individuals, low MAF) to identify candidate SNPs. Fast because it only outputs a SNP position list.
2. **Pass 2 — Genotype Likelihoods**: Revisit only those SNP positions and compute full genotype likelihood distributions (Beagle format) across all samples. Sites are retained if ≥80% of individuals have reads and the estimated allele frequency from the likelihoods exceeds 0.10. No genotypes are ever called.

This avoids computing full GL arrays for every site in the genome, which would be prohibitively slow.

### Disk strategy

All intermediate files are marked `temp()` in Snakemake and deleted automatically once no longer needed:

```
FASTQ (raw)      → deleted after fastp
FASTQ (trimmed)  → deleted after BWA
BAM (raw)        → deleted after dedup
BAM (dedup)      → deleted after CRAM filter
CRAM (filtered)  → kept permanently (~50% smaller than BAM)
```

Peak disk usage during a 96-sample run is bounded by how many samples are being processed simultaneously, not by total sample count.

---

## Pipeline Steps

### Stage 1 — Read Processing

| Step | Tool | What it does |
|------|------|--------------|
| `download_sra` | fasterq-dump | Downloads paired-end reads from NCBI SRA. On AWS EC2 (us-east-1), NCBI routes to their free S3 bucket so transfer costs nothing. Raw FASTQs deleted after trimming. |
| `fastp` | fastp | Adapter trimming, poly-G tail removal, per-cycle quality trimming. Produces per-sample HTML/JSON QC report. |
| `bwa_map` | BWA-MEM + samtools | Aligns trimmed reads to the *A. palmata* reference genome with read groups for downstream compatibility. |
| `mark_duplicates` | samtools markdup | Marks PCR duplicates in-stream without removing them (ANGSD handles duplicates internally). |
| `filter_mapped` | samtools view | Retains only properly paired, uniquely mapped reads (MAPQ ≥ 20) on chromosomes ≥ 10 Mb. Stored as CRAM with embedded reference path for portability. |

### Stage 2A — SNP Discovery and Population Structure

| Step | Tool | What it does |
|------|------|--------------|
| `angsd_discover_snps` | ANGSD | Pass 1: scans all CRAMs genome-wide to identify candidate SNP positions. Data quality thresholds require ≥60% of individuals to have reads at a site and depth within mean ± 2SD — ensuring sufficient data for reliable likelihood estimation without calling genotypes. |
| `angsd_gl_snps` | ANGSD | Pass 2: computes genotype likelihood distributions (Beagle format) at pass-1 positions. Sites must have data in ≥80% of individuals and estimated allele frequency >0.10 from the likelihoods. ANGSD never calls genotypes — all downstream analyses (PCA, FST, diversity) work directly from the likelihood distributions. |
| `compute_depth_thresholds` | Python | Calculates data-driven depth bounds (mean + 2SD of observed per-sample coverage) passed to both ANGSD passes. Sites exceeding the upper bound are excluded as likely multi-mapping artifacts from collapsed repeats. |
| `ngsld_estimate` | ngsLD | Estimates pairwise LD (r²) between all SNP pairs within 50 kb windows. |
| `ld_decay` | Python | Calculates the LD decay curve by distance bin to determine the pruning window size. |
| `ld_prune` | Python | Greedy graph-based LD pruning: removes one SNP from each pair with r² > 0.3. Pure Python — no external dependency. |
| `pcangsd` | PCAngsd | PCA and admixture proportions (K = 2–5) directly from genotype likelihoods on LD-pruned SNPs. |
| `ngsrelate` | ngsRelate | Estimates pairwise kinship (KING) from genotype likelihoods. Flags related pairs and putative clones. |

### Stage 2B — Diversity and Differentiation

| Step | Tool | What it does |
|------|------|--------------|
| `angsd_saf` | ANGSD | Computes site allele frequency likelihoods (SAF) per population using all sites — not just SNPs — for unbiased diversity estimates. |
| `realSFS` | ANGSD | Optimizes 1D and 2D site frequency spectra from SAF likelihoods via expectation-maximization. |
| `angsd_thetas` | ANGSD | Computes per-site θ and summarizes: nucleotide diversity (π), Watterson's θ, and Tajima's D per population. |
| `fst_estimate` | realSFS | Genome-wide weighted FST between each population pair from the joint 2D SFS. |
| `fst_windows` | realSFS | Sliding-window FST across the genome for Manhattan plots and outlier detection. |
| `heterozygosity` | ANGSD | Per-individual heterozygosity from single-sample SAF likelihoods. |

### Stage 3 — QC and Reporting

| Step | Tool | What it does |
|------|------|--------------|
| `multiqc` | MultiQC | Aggregates fastp QC reports across all samples into a single interactive report. |
| `depth_summary` | samtools coverage | Per-sample mean depth and breadth of coverage from CRAMs. |
| `filtering_summary` | samtools flagstat | Read counts at each filtering stage per sample. |
| `fst_outlier_annotation` | Python | Annotates FST outlier windows with nearby genes from the GFF annotation. |
| `generate_report` | Python/matplotlib | Self-contained HTML report with embedded figures + separate 300 DPI PNGs saved to `results/figures/`. |

---

## Running the Pipeline

### Requirements

- [Snakemake](https://snakemake.readthedocs.io) ≥ 9.0
- [Singularity/Apptainer](https://apptainer.org) ≥ 1.3 (for containerized rules)
- `fasterq-dump` (SRA toolkit) installed on the host — the container version has TLS issues on EC2
- `ngsLD` compiled from source on the host — not in bioconda

### Local workstation

```bash
snakemake --snakefile workflow/Snakefile --profile profiles/local
```

Set `local_conda_env` in `config/config.yaml` to your conda environment path.

### SLURM HPC

```bash
snakemake --snakefile workflow/Snakefile --profile profiles/slurm
```

### AWS EC2

Use `scripts/run_pipeline.sh` — it handles pre-flight checks, auto-retry on transient failures
(up to 20 retries, 5 min apart), S3 sync, SNS notifications, and auto-stop:

```bash
cd ~/coral-angsd-pipeline
nohup bash scripts/run_pipeline.sh >> logs/run.log 2>&1 &
```

Set up a fresh instance first with:

```bash
bash scripts/aws_setup.sh
```

#### EC2 cost options

The pipeline currently runs on an on-demand EC2 instance. Two options exist to reduce cost
significantly for future runs:

**Option 1 — EC2 Spot Instance (~70% cheaper, minimal changes)**

Request a spot instance at launch instead of on-demand. Current spot pricing (us-east-1):

| Instance | On-demand | Spot | Savings |
|----------|-----------|------|---------|
| c6i.16xlarge | $2.72/hr | ~$0.83–0.95/hr | ~70% |
| c6i.8xlarge  | $1.36/hr | ~$0.41–0.58/hr | ~65% |

With `rerun-incomplete: true` in the AWS profile and results syncing to S3, a spot interruption
means launching a fresh instance and restarting — Snakemake picks up from completed jobs.
To make this fully automatic, pair with an EC2 Auto Scaling Group (min/max/desired = 1, spot)
so AWS replaces an interrupted instance without manual intervention.

To request spot at launch via CLI:

```bash
aws ec2 run-instances \
  --instance-market-options '{"MarketType":"spot"}' \
  --instance-type c6i.16xlarge \
  ...
```

**Option 2 — AWS Batch with Spot (more setup, per-job retry)**

Snakemake's `--executor awsbatch` submits each rule as an individual Batch job on spot instances.
Batch handles interruptions and retries at the job level, so a spot interruption only restarts
the affected job (not the whole pipeline).

Blockers to resolve before using Batch:
1. `download_sra` has no container — fasterq-dump has TLS issues in containers on EC2.
   Needs a custom Docker image with SRA toolkit + AWS CLI baked in, pushed to ECR.
2. All I/O must go through S3 — Snakemake's S3 storage plugin already partially set up
   (`profiles/aws/config.yaml`), but every rule's inputs/outputs need to be S3-backed.
3. Compute environment setup: Batch job queue, spot compute environment, IAM roles.

Both options require that the S3 bucket (`coral-angsd-728009587639`) remains the central
store for results between runs — which is already the case.

---

## Configuration

Key parameters in `config/config.yaml`:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `samples_csv` | `config/samples_test.csv` | Sample sheet with `sample_id`, `sra_accession`, `population` columns |
| `reference` | `reference/apalmata_genome.fasta` | Reference genome (soft-masked for repeat detection) |
| `min_maf` | 0.10 | Minor allele frequency threshold for pass-2 SNPs |
| `min_ind_frac` | 0.80 | Fraction of individuals required at each site (pass 2); pass 1 uses 0.60 |
| `ld_max_dist` | 50,000 bp | Max distance for pairwise LD estimation (50 kb) |
| `ld_r2_threshold` | 0.3 | r² threshold for LD pruning |
| `max_k` | 5 | Maximum K for PCAngsd admixture analysis (runs K=2 through max_k) |
| `min_scaffold_size` | 10,000,000 | Only analyze scaffolds ≥ 10 Mb (chromosomes only) |

**Depth filtering** (ANGSD best practice): `-setMaxDepth` and `-setMaxDepthInd` are not hardcoded. The `compute_depth_thresholds` rule calculates mean + 2SD of observed per-sample coverage after alignment and passes those values to both ANGSD passes. This prevents inflated read counts in repetitive or multi-mapping regions from driving false SNP calls without relying on an arbitrary ceiling.

---

## Repository Structure

```
workflow/
  Snakefile               # Main pipeline
  scripts/
    generate_report.py    # HTML report generator
config/
  config.yaml             # Main configuration
  samples.csv             # 96-sample production sheet
  samples_test.csv        # 5-sample local test
  samples_aws_test.csv    # 10-sample AWS test
profiles/
  local/config.yaml       # Local workstation settings
  slurm/config.yaml       # FAU/Northeastern HPC settings
  aws/config.yaml         # AWS EC2 settings
reference/                # Reference genome (not tracked in git)
results/                  # All outputs (not tracked in git)
```

---

## Samples

See `config/samples.csv` for the full 96-sample manifest. Each row:

| Column | Example | Description |
|--------|---------|-------------|
| `sample_id` | `Ac_FL_M10` | Unique identifier |
| `sra_accession` | `SRR24007640` | NCBI SRA accession |
| `population` | `florida` | Population label |
| `region` | `Florida` | Collection region |
