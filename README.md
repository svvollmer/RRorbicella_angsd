# Coral ANGSD Pipeline

Population genomics pipeline for *Acropora palmata* (elkhorn coral) using a two-pass ANGSD workflow. Designed to run on a local workstation, SLURM HPC, or AWS EC2 with Singularity containers.

**Vollmer Lab — Florida Atlantic University**

---

## Study System

96 *A. palmata* colonies sampled across two regions:
- **Florida** (FL): Florida Keys reef system
- **Panama** (PA): Caribbean Panama

Research questions: population structure, gene flow, local adaptation, and relatedness across the Florida–Panama range.

---

## Pipeline Strategy

### Why two-pass ANGSD?

ANGSD works on genotype likelihoods rather than hard genotype calls, making it well-suited for low-to-moderate coverage coral resequencing data. A single-pass approach requires choosing SNP filters before knowing the data; a two-pass approach lets the data inform the filters:

1. **Pass 1 — Discovery**: Scan the whole genome with relaxed filters (60% of individuals, low MAF) to identify candidate SNPs. Fast because it only outputs a SNP position list.
2. **Pass 2 — Genotype Likelihoods**: Revisit only those SNP positions with strict filters (80% of individuals, MAF > 0.10), computing full genotype likelihood arrays (Beagle format) for downstream analysis.

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
| `angsd_discover_snps` | ANGSD | Pass 1: scans all CRAMs genome-wide with relaxed filters (60% individuals, `-doMaf 1`) to identify candidate SNP positions. Max depth set to mean + 2SD of observed coverage. |
| `angsd_gl_snps` | ANGSD | Pass 2: computes genotype likelihoods at pass-1 SNPs with strict filters (80% individuals, MAF > 0.10, max depth = mean + 2SD). Outputs Beagle-format GL matrix. |
| `compute_depth_thresholds` | Python | Calculates per-individual and total max depth thresholds as mean + 2SD of observed sample coverage. ANGSD best practice — avoids hard-coded ceilings that fail to adapt to actual library depth. |
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

```bash
snakemake --snakefile workflow/Snakefile \
  --cores 16 --jobs 4 --keep-going --latency-wait 120 \
  --rerun-incomplete --rerun-triggers mtime \
  --use-singularity \
  --singularity-args '--bind /home/ubuntu --bind /usr/local/bin --bind /opt' \
  --config samples_csv=config/samples_aws_test.csv \
  > snakemake.log 2>&1 &
```

Comment out `local_conda_env` in `config/config.yaml` on AWS — containers provide the correct PATH.

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
