# Changelog

Development history of the RRorbicella_angsd pipeline.

---

## v0.2 — Segment 1 launch, samples_new.csv, HPC config (2026-04-02)

**HPC working directory:** `/work/vollmer/orbicella_genomics/`
**Pipeline code:** `/projects/vollmer/RRorbicella_angsd/`

### What changed

- Launched Segment 1 on Discovery HPC for 92 new FL samples (Oann + Ofav)
  - New FASTQs at `/work/vollmer/sequencing_archive/250625_ReefRenewal-FL_Orbicella_WGS/`
  - 86/92 fastp trimming already complete from prior partial run; 6 remaining submitted
  - 92 bwa_map jobs running across 50 SLURM nodes at launch
- HPC config uses `config/samples_new.csv` (92 `local` FASTQ samples, Oann + Ofav only)
- Cleared 5,413 orphaned `samtools sort` temp files from previous crashed run
- Ofranksi (22 samples) not in new sequencing batch; will ingest via `local_bam` from
  `/projects/vollmer/RR_heat-tolerance/Orbicella/2_mapping.bwa/` in a follow-up run

---

## v0.1 — Pipeline initialization (2026-03-xx)

**Forked from:** `coral-angsd-pipeline` (Acropora multi-species framework)

### What was set up

- `config/config.yaml` configured for jaOrbFran1.1 reference, three-species design
  - `primary_grouping: species` (FST/diversity grouped by species, not geography)
  - `fst_comparisons`: all 3 pairwise inter-species combinations
  - `max_k: 10`, `admixture_replicates: 20`
- `config/samples.csv`: 140-sample manifest
  - 63 *O. annularis* (17 `local_bam` + 46 `merge_bam`)
  - 55 *O. faveolata* (9 `local_bam` + 46 `merge_bam`)
  - 22 *O. franksi* (`local_bam` only)
- `profiles/discovery/`: Northeastern Discovery HPC SLURM profile
- `run.sh` configured for `/projects/vollmer/RRorbicella_angsd/` pipeline path
  and `/work/vollmer/orbicella_genomics/` working directory
- Reference genome symlinked in working dir as `reference.fna` (BWA index pre-built)
