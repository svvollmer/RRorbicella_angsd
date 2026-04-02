# Next Steps

Current state as of 2026-04-02. Segment 1 running on Discovery HPC.

---

## Pipeline status

| Segment | Status | Notes |
|---------|--------|-------|
| Seg 1: BAM ingestion + QC | 🔄 Running | 92 Oann+Ofav new FASTQs; Ofranksi pending |
| Seg 2a: SNP discovery + relatedness | ⬜ Pending | |
| Seg 2b: Subset BEAGLE | ⬜ Pending | |
| Seg 3: PCA + admixture | ⬜ Pending | |
| Seg 4: Diversity + FST | ⬜ Pending | Three-species pairwise FST |
| Seg 6: Demography (moments) | ⬜ Pending | All 3 inter-species pairs |
| Seg 7: SMC++ | ⬜ Pending | Ne(t) for all 3 species |

---

## Immediate next steps

### 1. Monitor Segment 1

```bash
ssh s.vollmer@login.discovery.neu.edu
squeue -u s.vollmer
tail -f /work/vollmer/orbicella_genomics/logs/segment1_launch.log
```

Expected runtime: fastp ~30–60 min, bwa_map ~2–4h per sample (running in parallel).
Segment 1 is done when Snakemake exits successfully and writes `results/qc/samples_pass.txt`.

### 2. Gate 1 — QC review

```bash
python /projects/vollmer/RRorbicella_angsd/workflow/scripts/qc_approve.py
# → writes results/qc/samples_approved.txt
```

Flag samples with mean depth < 5× or mapping rate < 50%. *Orbicella* mapping to the
*O. franksi* reference should be high for all three species (they are closely related and
share synteny). Flag any inter-species outliers — may reflect divergent structural variants.

### 3. Add Ofranksi to the run

The 22 *O. franksi* samples (`local_bam` type) are not in the current run. After Gate 1,
either:
- Append Ofranksi rows to `samples_new.csv` with `input_type: local_bam` and
  `bam_path` pointing to the existing BAMs, then re-run Segment 1 for those samples only
- Or run a separate Segment 1 pass for Ofranksi using the `local_bam` ingest path

Ofranksi BAMs: `/projects/vollmer/RR_heat-tolerance/Orbicella/2_mapping.bwa/Ofrank_*.bwa.dedup.clip.bam`

### 4. Segment 2a — SNP discovery

```bash
bash run.sh 2a
```

Stops at clone gate. After completion:
```bash
python /projects/vollmer/RRorbicella_angsd/workflow/scripts/clone_approve.py
bash run.sh 2b
```

### 5. Segment 3 — PCA + admixture

```bash
bash run.sh 3
```

Then lineage gate:
```bash
python /projects/vollmer/RRorbicella_angsd/workflow/scripts/lineage_assign.py
```

**Key question:** Does K=3 cleanly separate the three *Orbicella* species, or is there
admixture between species pairs? The species complex is known to hybridize (*O. annularis*
× *O. faveolata* hybrids documented). Assign lineage labels post-PCA.

### 6. Segment 4 — Diversity + three-species FST

```bash
bash run.sh 4
```

FST computed for all three pairs (configured in `config.yaml`):
- *O. annularis* vs *O. faveolata*
- *O. annularis* vs *O. franksi*
- *O. faveolata* vs *O. franksi*

Also produces per-species θ_π, θ_W, and Tajima's D.

---

## Downstream analyses (post-Seg 4)

### Demographic inference (Segment 6 — moments)

2-population IM / AM / SC models for all three species pairs. Same approach as
Acropora pipeline: run SI/IM/IM_a/AM/SC as separate SLURM jobs per comparison,
compare AIC, report asymmetric migration rates and divergence times.

Key parameters:
- μ for *Orbicella*: TBD (use *A. millepora* 1.8×10⁻⁸ as first pass, or find an
  *Orbicella*-specific estimate from the literature)
- Generation time: ~10 yr (massive corals; slower than branching *Acropora*)

### Population size history (Segment 7 — SMC++)

Ne(t) curves for all three species using Florida high-purity (Q ≥ 0.95) samples.
Expected signals:
- LGM bottleneck ~18–20 kya (shared across species)
- Post-LGM recovery may differ: *O. faveolata* experienced catastrophic decline since
  1980s (white band, thermal bleaching); modern Ne may be severely depressed
- *O. annularis* colonies form large clonal patches — watch for reduced N_e from clonality
  reducing the effective sample count

### FST outlier annotation

After Seg 4 FST is complete, run sliding-window outlier analysis for each species pair
(top/bottom 2.5% 50-kb windows) and annotate against the *O. franksi* BRAKER GTF.
Candidate reproductive isolation loci and introgression tracts.

### Phenotype association (future)

Reef Renewal heat tolerance scores are available for restoration stock samples.
WGAS (genotype–phenotype association via ANGSD score test) will be enabled once
`heat_tolerance_score` column is added to samples.csv.

---

## Quick reference — HPC

```bash
# Check queue
squeue -u s.vollmer

# Segment 1 log
tail -f /work/vollmer/orbicella_genomics/logs/segment1_launch.log

# BAM outputs (when Seg 1 done)
ls /work/vollmer/orbicella_genomics/results/bams/*.filtered.* | wc -l

# Restart a segment
cd /work/vollmer/orbicella_genomics
bash run.sh 1   # or 2a, 2b, 3, 4
```

---

## Key file locations

| Item | Path |
|------|------|
| Pipeline code | `/projects/vollmer/RRorbicella_angsd/` |
| Working directory | `/work/vollmer/orbicella_genomics/` |
| Reference genome | `/projects/vollmer/RR_heat-tolerance/Orbicella/reference/GCA_964199315.1_jaOrbFran1.1_genomic.fna` |
| Old Orbicella BAMs | `/projects/vollmer/RR_heat-tolerance/Orbicella/2_mapping.bwa/` |
| New FL FASTQs | `/work/vollmer/sequencing_archive/250625_ReefRenewal-FL_Orbicella_WGS/01.RawData/` |
| Snakemake env | `/home/s.vollmer/.conda/envs/snakemake2/bin/snakemake` |
