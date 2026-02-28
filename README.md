# Coral ANGSD Pipeline

Two-pass ANGSD population genomics workflow for *Acropora palmata* (96 samples: Florida + Panama).

## What it does

```
Stage 1: Data prep (per sample)
  SRA download → fastp QC → BWA-MEM mapping → samtools dedup → quality filter
  → QC metrics (flagstat, depth, multiqc)

Stage 2A: Population structure (SNP-based)
  ANGSD pass1 (SNP discovery) → ANGSD pass2 (genotype likelihoods)
  → ngsLD (pairwise LD within 500 kb) → LD pruning (r² > 0.3)
  → PCAngsd (PCA + admixture K=2–5) → ngsRelate (kinship)

Stage 2B: Diversity & divergence (all-sites)
  ANGSD SAF per population → realSFS (1D + 2D SFS)
  → thetaStat (π, θ, Tajima's D) → FST global + windowed
  → Individual heterozygosity per sample

Stage 3: Annotation & report
  GFF → gene BED → FST outliers × genes → HTML report
```

## Dependencies

All tools run via Singularity containers (no local installs needed beyond Snakemake + Singularity).

- Snakemake ≥ 7.0
- Singularity or Apptainer

## Configuration

| File | Purpose |
|------|---------|
| `config/config.yaml` | Main pipeline settings (reference path, filters, MAF, etc.) |
| `config/samples.csv` | Full 96-sample dataset |
| `config/samples_test.csv` | 5-sample subset (Florida + Panama) for local testing |

## Usage

### Local workstation (test run, 5 samples with existing BAMs)

```bash
# Edit config/config.yaml: set samples_csv to config/samples_test.csv
snakemake --snakefile workflow/Snakefile --profile profiles/local
```

### Local workstation (full 96-sample run)

```bash
# Edit config/config.yaml: set samples_csv to config/samples.csv
snakemake --snakefile workflow/Snakefile --profile profiles/local
```

### FAU HPC or Northeastern HPC (SLURM)

```bash
# Edit profiles/slurm/config.yaml: set slurm_partition for your cluster
# FAU KoKo: general, himem, gpu
# Northeastern Discovery: short, long, himem

module load singularity  # or apptainer
snakemake --snakefile workflow/Snakefile --profile profiles/slurm
```

### AWS EC2 (single large instance)

```bash
# Launch e.g. c6i.16xlarge (64 vCPU, 128 GB RAM) with Amazon Linux 2
# Install Snakemake + Singularity/Apptainer
snakemake --snakefile workflow/Snakefile --profile profiles/aws
```

## Reference genome

The reference and annotation files are too large for git. See `resources/README.txt` for setup.

Required path (set in `config/config.yaml`):
```
reference/GCF_021335395.2_Apul_2.0_genomic.fna
```

## Dry run (check DAG without executing)

```bash
snakemake --snakefile workflow/Snakefile --profile profiles/local --dry-run
```

## Directory structure

```
coral-angsd-pipeline/
├── workflow/
│   ├── Snakefile
│   └── scripts/
├── config/
│   ├── config.yaml
│   ├── samples.csv           ← 96-sample full run
│   └── samples_test.csv      ← 5-sample test run
├── profiles/
│   ├── local/config.yaml
│   ├── slurm/config.yaml
│   └── aws/config.yaml
├── resources/README.txt      ← reference genome setup instructions
└── archive/                  ← old Snakefile versions
```

## Key design notes

- **Two-pass ANGSD**: Pass 1 discovers SNPs (relaxed filters), Pass 2 computes genotype likelihoods only at discovered SNPs (strict filters). ~10× faster than single-pass.
- **LD pruning before PCA/Admixture**: ngsLD computes r² within 500 kb; SNPs with r² > 0.3 are removed before PCAngsd runs.
- **All-sites diversity**: SAF/SFS analysis uses all sites at discovered positions (not SNP-filtered), giving unbiased π and Tajima's D.
- **Singularity containers**: Each rule uses a Biocontainers image. Run with `--use-singularity` (included in all profiles).
