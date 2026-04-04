# Changelog

Development history of the RRorbicella_angsd pipeline.

---

## v0.2 — Segment 1 launch, samples_new.csv, HPC config (2026-04-02)

**HPC working directory:** `/work/vollmer/orbicella_genomics/`
**Pipeline code:** `/projects/vollmer/RRorbicella_angsd/`

### What changed

- Launched Segment 1 on Discovery HPC for 92 FASTQ-available samples (46 Oann + 46 Ofav)
  - FASTQs at `/work/vollmer/sequencing_archive/250625_ReefRenewal-FL_Orbicella_WGS/`
  - 86/92 fastp trimming already complete from prior partial run; 6 remaining submitted
  - 92 bwa_map jobs running across 50 SLURM nodes at launch
- HPC config uses `config/samples_new.csv` (92 `local` FASTQ samples, Oann + Ofav only)
- Cleared 5,413 orphaned `samtools sort` temp files from previous crashed run
- 48 BAM-only samples (17 Oann + 9 Ofav + 22 Ofranksi) pending `local_bam` ingest;
  no raw FASTQs preserved — must use pre-existing BAMs from
  `/projects/vollmer/RR_heat-tolerance/Orbicella/2_mapping.bwa/`

---

## v0.1 — Pipeline initialization (2026-03-xx)

**Forked from:** `coral-angsd-pipeline` (Acropora multi-species framework)

### What was set up

- `config/config.yaml` configured for jaOrbFran1.1 reference, three-species design
  - `primary_grouping: species` (FST/diversity grouped by species, not geography)
  - `fst_comparisons`: all 3 pairwise inter-species combinations
  - `max_k: 10`, `admixture_replicates: 20`
- `config/samples.csv`: 140-sample manifest
  - 63 *O. annularis* (46 `local` FASTQ-available + 17 `local_bam` BAM-only)
  - 55 *O. faveolata* (46 `local` FASTQ-available + 9 `local_bam` BAM-only)
  - 22 *O. franksi* (`local_bam` BAM-only; no raw FASTQs available)
- `profiles/discovery/`: Northeastern Discovery HPC SLURM profile
- `run.sh` configured for `/projects/vollmer/RRorbicella_angsd/` pipeline path
  and `/work/vollmer/orbicella_genomics/` working directory
- Reference genome symlinked in working dir as `reference.fna` (BWA index pre-built)

---

## v0.3 — Segments 2a/2b/3 complete (2026-04-03 → 2026-04-04)

**HPC working directory:** `/work/vollmer/orbicella_genomics/`
**Samples:** 140 total (63 Oannularis, 57 Ofaveolata, 21 Ofranksi)

### Segment 2a — SNP discovery + relatedness

- Two-pass ANGSD SNP calling across 15 chromosomes
- SNP filters: min_maf=0.05, min_ind=112 (80% of 140) → ~7.2M SNPs passing
- ngsRelate all-vs-all → 6 clone pairs detected, 6 excluded, **134 unrelated retained**
- **Ofra_26 confirmed mislabel** — labeled O. franksi but clusters with O. faveolata:
  mean KING vs Ofav=−0.038, vs Ofra=−0.521; KING vs Ofav_106=+0.493 (clone-level)
  Decision: exclude Ofra_26, retain Ofav_106

### Bugs fixed

- ngsRelate group keys with commas/spaces (e.g. NFWF, Middle Keys) broke sbatch --wrap:
  sanitize group names in NGSRELATE_GROUPS dict (replace ,  → _)
- ngsLD binary CXXABI mismatch on compute nodes: pass conda env LD_LIBRARY_PATH in subprocess
- ngsLD binary missing: symlinked from snakemake2 conda env to bin/ngsLD

### Segment 2b — subset BEAGLE + grouped ngsRelate

- BEAGLE subset to 134 unrelated samples
- Grouped ngsRelate per (species × region) — 14 groups, no additional clones found

### Segment 3 — LD, PCA, admixture (240 steps)

- LD estimation (ngsLD) across 15 chromosomes
- LD pruning → LD-pruned BEAGLE for admixture
- PCA (PCAngsd covariance matrix)
- ngsAdmix K=1–10, 20 replicates each
- PCAngsd admixture K=2–10
- All outputs committed to docs/outputs/
