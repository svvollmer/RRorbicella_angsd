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

---

## v0.4 — Lineage assignment, Seg4 (diversity/FST), Seg6 launch (2026-04-04 → 2026-04-06)

**HPC working directory:** 
**Samples:** 133 for seg4+ (Ofra_6 excluded as hybrid)

### Lineage assignment gate

- K=3 Q matrix analyzed to assign all 134 unrelated samples to genetic lineages
- **10 field label mislabels identified** (7.5% error rate, 7.9% including Ofra_26):
  - Oann_106 → Ofaveolata; Ofav_94 → Oannularis
  - Ofav_75/76/83/84/104/105 → Ofranksi (6 samples)
  - Ofra_1/4 → Ofaveolata
- **Ofra_6:** putative franksi×faveolata hybrid (8.4% faveolata ancestry) — excluded
- Lineage assignments written to results/admixture/lineage_assignments.txt
  (Oann=59, Ofav=50, Ofra=24; total=133)
-  updated with  column (genetic assignment, mislabels corrected)
- All downstream analyses use genetic lineage assignments, not field labels

### Key biology added to RESULTS.md

- Spawning phenology note: Ofav+Ofra spawn closer in time → actual hybridization
  Oann+Ofra spawn hours apart despite lab cross-compatibility (Levitan et al.)
  → explains UPGMA topology being reference genome bias artefact, not gene flow
- Reference genome bias caveat added throughout: Oann-Ofra appear closest due to
  SNP ascertainment in Ofranksi sequence space

### Segment 4 — SAF, SFS, diversity, FST

- Per-sample SAFs + individual heterozygosity (133 samples)
- Group SAFs for Oannularis_all, Ofaveolata_all, Ofranksi_all
- Genome-wide diversity results:
  - π: Oann=0.01241, Ofav=0.01228, Ofra=0.01272 (all similar)
  - Tajima's D: Oann=−0.59, Ofav=−0.20, Ofra=−0.52
- FST: Oann-Ofra=0.180; Ofav-Ofra=0.109; Oann-Ofav=pending
- Sliding-window FST (50kb windows, 10kb step): Oann-Ofra and Ofav-Ofra done
- Windowed diversity: 47,208 windows per species (thetaStat per-chrom approach)

### Bugs fixed (seg4)

- All singularity exec calls require  — without this, SAF/thetas
  files silently fail to open (ANGSD Problems opening file error)
- thetaStat  fails on genome-wide merged .thetas.idx — fix: run
  thetaStat per-chromosome SAF → windowed pestPG, then concatenate
- short partition allows 2-day wall time; Snakemake caps at runtime=1440min (24h)
  — realsfs 2D for Oann-Ofav needed manual sbatch with --time=2-00:00:00 and 16 CPUs

### Segment 6 — Moments demographic inference (RUNNING as of 2026-04-06)

- 3 pairwise comparisons: oann_vs_ofav, oann_vs_ofra, ofav_vs_ofra
- Models: SI, IM, IM_a (asymmetric), SC (secondary contact), AM (ancient migration)
- 20 random restarts each; AIC model selection
- SAFs regenerated without MAF filter (monomorphic sites required for unbiased SFS)
- Snakefile.demography updated: get_demog_bamlist extended to support lineage key
  in population_groups config (enables filtering by genetic lineage, not field label)

### Plots added to docs/figures/

- chrom_diversity.png — π/θ/TajimaD per chromosome, all 3 species
- diversity_summary.png — genome-wide summary bars
- fst_summary.png — pairwise FST heatmap (Oann-Ofav cell pending)
- fst_windows_Oannularis_vs_Ofranksi.png
- fst_windows_Ofaveolata_vs_Ofranksi.png
- windowed_pi.png — 50kb sliding window nucleotide diversity
- windowed_tajima.png — 50kb sliding window Tajima's D

### Pending

- Oann vs Ofav FST (SLURM job, ~17h elapsed, short partition 2-day limit)
- Seg6 moments results (running)
- Seg7 SMC++ (Snakefile.smcpp not yet written for Orbicella; gen_time=10yr, μ=1.8e-8)
