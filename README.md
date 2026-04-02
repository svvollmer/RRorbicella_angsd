# RRorbicella_angsd

Population genomics pipeline for three *Orbicella* species (massive brain corals) from the Reef Renewal Florida restoration program. Built on the Vollmer Lab coral ANGSD pipeline framework.

**Vollmer Lab — Northeastern University**

---

## Study System

**Species:**
- *Orbicella annularis* (lobed star coral) — 63 samples
- *Orbicella faveolata* (mountainous star coral) — 55 samples
- *Orbicella franksi* (boulder star coral) — 22 samples

**Populations:** Florida Keys (Reef Renewal restoration sites)

**Reference genome:** *Orbicella franksi* chromosome-scale assembly
[jaOrbFran1.1 / GCA_964199315.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_964199315.1/)
Gene annotation: BRAKER-predicted GTF (`braker.gtf`)

**Research questions:** Population structure and inter-species divergence across the *Orbicella* species complex; identification of genomic barriers and introgression tracts between species pairs; demographic history and comparative Ne(t) trajectories; genomic basis of heat tolerance in Reef Renewal restoration stock.

---

## Analytical Design

This is a **three-species comparison** pipeline. All major downstream analyses — FST, demographic inference, and SMC++ — are run for all three inter-species pairs:

| Comparison | Type |
|------------|------|
| *O. annularis* vs *O. faveolata* | inter-species |
| *O. annularis* vs *O. franksi* | inter-species |
| *O. faveolata* vs *O. franksi* | inter-species |

Within-species geographic comparisons will be added if multi-population sampling is available.

---

## Pipeline Strategy

### Two-pass ANGSD

1. **Pass 1 — Discovery:** Scan the whole genome to find candidate SNP positions (relaxed filters, 60% individuals).
2. **Pass 2 — Genotype Likelihoods:** Compute full genotype likelihood distributions (Beagle format) at pass-1 positions (strict filters, 80% individuals, MAF > 0.05). No genotypes are ever hard-called.

### Input types

All 140 samples are pre-existing data — no SRA downloads needed.

| Input type | N | Species | Description |
|------------|---|---------|-------------|
| `local` | 92 | Oann (46), Ofav (46) | Raw FASTQs available; re-aligned from scratch through pipeline for consistency. FASTQs at `/work/vollmer/sequencing_archive/250625_ReefRenewal-FL_Orbicella_WGS/` |
| `local_bam` | 48 | Oann (17), Ofav (9), Ofrank (22) | Raw FASTQs not available; using pre-existing BAMs from the heat tolerance project pipeline (`/projects/vollmer/RR_heat-tolerance/Orbicella/2_mapping.bwa/`) |

Note: the pre-existing BAMs were aligned from the same raw reads using a different pipeline (Picard MarkDuplicates, custom clipping). Re-alignment from FASTQ is preferred where possible for consistency — same samtools markdup, same MAPQ/chromosome filters, same read group tags.

---

## Pipeline Segments

| Segment | Snakefile | Description |
|---------|-----------|-------------|
| 1 | `Snakefile.align` | BAM ingestion or FASTQ → align → QC; Gate 1 (depth/mapping QC) |
| 2a | `Snakefile.snps` | Pass 1 SNP discovery + GL + relatedness; stops at clone gate |
| 2b | `Snakefile.snps` | Subset BEAGLE to unrelated samples |
| 3 | `Snakefile.structure` | LD pruning, PCA, NGSAdmix K=1–10 |
| 4 | `Snakefile.diversity` | SAF, SFS, θ, Tajima's D, pairwise FST (all 3 species pairs) |
| 6 | `Snakefile.demography` | Moments 2-pop models for all 3 inter-species pairs |
| 7 | `Snakefile.smcpp` | SMC++ Ne(t) for all 3 species |

---

## Running the Pipeline

All segments launch via the `run.sh` wrapper on the HPC:

```bash
# On Discovery HPC, from /work/vollmer/orbicella_genomics/
bash run.sh 1        # BAM ingestion + QC
# → python workflow/scripts/qc_approve.py  (Gate 1)
bash run.sh 2a       # SNP discovery + relatedness
# → python workflow/scripts/clone_approve.py  (clone gate)
bash run.sh 2b       # Subset BEAGLE to unrelated samples
bash run.sh 3        # PCA + admixture
# → python workflow/scripts/lineage_assign.py  (lineage gate)
bash run.sh 4        # Diversity + FST
```

**HPC setup:**
- Cluster: Northeastern Discovery (`s.vollmer@login.discovery.neu.edu`)
- Pipeline code: `/projects/vollmer/RRorbicella_angsd/`
- Working directory: `/work/vollmer/orbicella_genomics/`
- Snakemake env: `/home/s.vollmer/.conda/envs/snakemake2/bin/snakemake`
- Profile: `profiles/discovery`

---

## Configuration

Key parameters in `config/config.yaml`:

| Parameter | Value | Notes |
|-----------|-------|-------|
| `reference` | jaOrbFran1.1 | O. franksi genome |
| `primary_grouping` | `species` | FST and diversity grouped by species |
| `min_maf` | 0.05 | Minor allele frequency floor |
| `min_ind_frac` | 0.80 | 80% individuals required per site |
| `fst_comparisons` | all 3 pairs | Oann/Ofav, Oann/Ofrank, Ofav/Ofrank |
| `max_k` | 10 | NGSAdmix K=1–10, 20 reps each |
| `clone_threshold` | 0.45 | KING kinship |

---

## Repository Structure

```
workflow/
  Snakefile.align       # Segment 1: ingestion, align, QC
  Snakefile.snps        # Segment 2: SNP discovery, GL, relatedness
  Snakefile.structure   # Segment 3: LD, PCA, NGSAdmix
  Snakefile.diversity   # Segment 4: SAF, SFS, diversity, FST
  Snakefile.demography  # Segment 6: moments demographic inference
  Snakefile.smcpp       # Segment 7: SMC++ Ne(t)
  common.smk            # Shared config, sample loading, helpers
  scripts/              # Python helper scripts
config/
  config.yaml           # Main configuration
  samples.csv           # 140-sample manifest (all input types)
profiles/
  discovery/            # Northeastern Discovery HPC SLURM profile
docs/
  figures/              # Output plots
```
