# Changelog

Development history of the coral ANGSD pipeline from local prototype to AWS production run.

---

## v0.4 — 290-sample Discovery HPC production run (2026-03-16 → ongoing)

**Environment**: Northeastern Discovery HPC (`s.vollmer@login.discovery.neu.edu`)
**Working dir**: `/work/vollmer/acropora_genomics/`
**Pipeline**: `/projects/vollmer/coral-angsd-pipeline/`
**Samples**: 290 — *A. cervicornis* and *A. palmata* across Florida, Panama, and Bonaire (`config/samples_RR.csv`)
**Status**: Segment 2 complete (SNP discovery, GL, relatedness). Segment 3+ pending.

### What was completed

- Full two-pass ANGSD SNP discovery + GL on 14 chromosomes (290 samples, *A. palmata* reference)
- 2,034,805 SNPs passing filters (MAF > 0.05, ≥ 80% individuals genotyped)
- ngsRelate (all-vs-all, 100K MAF>0.3 sites) → 74 clone pairs across 36 excluded samples
- `clone_approve.py --yes --exclude RR_FL_Apal_015` → 253 unrelated samples retained
  - RR_FL_Apal_015 excluded as artifact (KING ≥ 0.20 with 60+ samples across species — impossible biologically)
  - Zero cross-species or cross-population clone pairs confirmed
- Subset beagle to 253 unrelated samples (`all.unrelated.beagle.gz`) — ready for Segment 3
- Grouped parallel ngsRelate launched (2026-03-17): 6 SLURM jobs (Acervicornis/Apalmata × FL/PA/BON),
  full MAF>0.3 sites (~783K), for comparison with 100K-subsampled all-vs-all run

### Clonality results (all-vs-all run)

| Group | N | Clone pairs | Excluded | N genets | Genet diversity |
|-------|---|-------------|----------|----------|-----------------|
| *Acervicornis* BON | 25 | 0 | 0 | 25 | 1.000 |
| *Acervicornis* FL | 104 | 28 | 13 | 91 | 0.875 |
| *Acervicornis* PA | 49 | 32 | 10 | 39 | 0.796 |
| *Apalmata* BON | 25 | 0 | 0 | 25 | 1.000 |
| *Apalmata* FL | 88 | 11 | 11 | 77 | 0.875 |
| *Apalmata* PA | 9 | 3 | 2 | 7 | 0.778 |
| **Total** | **300** | **74** | **36** | **264** | **0.880** |

### New rules added (Snakefile.snps)

- **`ngsrelate_group`** — per-(species×region) ngsRelate, 6 parallel SLURM jobs, full MAF>0.3 sites
  - Beagle columns subset to group members; column indices computed in Python params to avoid BusyBox shell limitations
  - Outputs: `results/relatedness/grouped/{group}.res` + `{group}.samples.txt`
- **`merge_ngsrelate`** — combines 6 `.res` files into `results/relatedness/grouped/relatedness_summary.txt`,
  remapping 0-based within-group indices to sample names
- **`filter_clones_grouped`** — produces `results/relatedness/grouped/clones_report.txt` for comparison

### New run.sh segment

- **`bash run.sh 2c`** — launches grouped ngsRelate; logs to `logs/segment2c_*.log`

### Bugs fixed

- **Bug 23**: ngsRelate all-vs-all timeout — original 60 min limit too short for 290 samples; increased to 480 min
- **Bug 24**: `ngsRelate` subsampling with `shuf` fails in BusyBox container (not installed) — replaced with awk counter
- **Bug 25**: `head -n N` in awk pipeline causes SIGPIPE with bash strict mode (`set -euo pipefail`) — replaced with awk counter `n<N && rand()<0.2`
- **Bug 26**: `mktemp --suffix` not supported in BusyBox — removed suffix argument
- **Bug 27**: `ngsrelate_group` awk column-subset program: `"\t"` and `"\n"` in Snakemake shell block interpreted as literal tab/newline by Python → embedded in awk single-quoted string → `Unexpected end of string`; fixed with `"\\t"` and `"\\n"`
- **Bug 28**: `clone_approve.py` on HPC login node: default `python` is Python 2, which chokes on UTF-8 box-drawing characters — must use `python3` explicitly
- **Bug 29**: `ngsrelate_group` segfault (exit 139) on groups containing QC-failed samples — `params.group_samples` used all CSV members (300) while `_group_cols` only added columns for bamlist members (290, QC-passed), so `N_SAMPLES > actual beagle columns` → ngsRelate reads past array bounds; fixed by adding `_group_samples_in_bamlist()` helper that filters group members to bamlist-present samples, and using it for both `group_samples` param and `{output.samples}`

### Discovery HPC setup notes

- Snakemake env: `/home/s.vollmer/.conda/envs/snakemake2/bin/snakemake` (Python 3.11, Snakemake 8.30, numpy 1.26.4)
  - Old `snakemake` env broken: numpy 2.4.3 / glibc 2.17 mismatch
- `run.sh` requires `module load singularity/3.10.3` before snakemake
- `slurm_extra` quoting: use `"'--mail-type=END,FAIL ...'"` (nested quotes) for Snakemake 8

---

## v0.3 — AWS 96-sample production run (2026-03-09 → 2026-03-11)

**Instance**: c6i.16xlarge on-demand (`i-0194f3b3de13fcd6a`)
**Status**: Complete through diversity/PCA/admixture. FST killed intentionally (see below).

### What was validated
- Full 96-sample *A. cervicornis* dataset from the Vollmer Lab *Science* paper
- All 96 CRAMs generated and stored in S3
- QC gate: 93/96 samples approved; excluded Ac_FL_M1 (0.1x, 10% mapping), Ac_FL_U72 (0x), Ac_PA_CK48 (0x, 2.4% mapping)
- Clone/relatedness gate: 15 clones excluded, 78 unrelated samples retained
  - Notable: PA-HS was an 8-way clone cluster (HS1–HS8 all KING ≥ 0.47); kept HS1 (69x)
  - FL-K1/K2/U14 were a 3-way clone; kept U14 (63x)
- ANGSD pass2: 1,201,444 SNPs at MAF ≥ 0.10, minInd 80%
- PCA/admixture K2–K5, LD decay, individual heterozygosity, 1D SFS, diversity all complete
- realSFS 2D (for FL vs PA FST) killed after 20+ hours — ran without `-maxIter` cap (see bugs)

### Key biological findings
- **K=2 confirmed** as primary structure, consistent with the Science paper
- Two-component structure does NOT split cleanly by geography (FL vs PA)
  - Lineage A: FL_B/M series + PA_SR, PA_Tet, PA_HS
  - Lineage B: FL_U series + PA_CK43–49, PA_HS1
  - Four PA_CK samples (CK140, CK1410, CK142, CK144) are ~50/50 admixed
- Interpretation: K=2 reflects *A. cervicornis × A. palmata* hybridization present at
  both locations, mixed with geographic population structure — consistent with Science paper
- Tajima's D positive in both populations (FL mean +0.74, PA mean +0.53), likely driven
  by the two-lineage admixture structure inflating intermediate-frequency alleles
- θ_π nearly identical between FL (0.00351) and PA (0.00354); θ_W slightly higher in PA
- Chr 14 (NC_133895.1) outlier: highest D in both populations (FL 1.10, PA 0.74),
  elevated per-site π — candidate sex chromosome or introgression hotspot
- LD decays from r²=0.54 at <1kb to 0.19 by 50–100kb

### Bugs fixed (production run)
- **Bug 16**: fasterq-dump not found in Snakemake subshells on AWS — `export PATH=/usr/local/bin:$PATH` in `download_sra`
- **Bug 17**: Snakemake MILP scheduler crash (PuLP/CBC IndexError) on any job failure — fixed with `scheduler: greedy` in AWS profile
- **Bug 18**: fasterq-dump exit 3 (SRA network timeout) on large accessions — added 3-attempt retry loop with backoff
- **Bug 19**: fasterq-dump streaming fails on EC2 (NCBI internal endpoint unreachable) — switched to `aws s3 cp s3://sra-pub-run-odp/sra/ACC/ACC --no-sign-request` then local fasterq-dump
- **Bug 20**: `ec2:StopInstances` missing from IAM role — added to S3PipelineAccess inline policy
- **Bug 21**: `individual_saf` missing `-anc` flag → ANGSD exits 0 with no output — added `-anc {input.ref}`
- **Bug 22**: ANGSD treats index files with equal mtime as stale — fixed with `touch` on .bin/.idx/.fai after S3 sync

### Known issues / intentional decisions
- realSFS 2D has no `-maxIter` cap → ran for 20+ hours without converging; killed intentionally
- Auto-stop via IAM still failing in some contexts despite Bug 20 fix
- FST and final HTML report not generated in this run

---

## v0.2 — AWS 10-sample test (2026-02-28 → 2026-03-02)

**Instance**: c6i.4xlarge (`i-0001a82fe42f1d26f`)
**Status**: Complete — all pipeline steps including report

### What was validated
- Full end-to-end pipeline on AWS including SRA download, alignment, ANGSD, PCA, FST, report
- Singularity containers working for all rules except `download_sra` (TLS issue) and ngsLD (not in bioconda)
- S3 as central data store; pipeline survives instance restart
- QC gate and clone gate both functional
- Report HTML generated successfully

### Bugs fixed
- **Bug 7**: BWA OOM (exit 137) on c6i.4xlarge — reduced `samtools sort -@ 16 -m 1G` to `-@ 4 -m 512M`
- **Bug 8**: `popname` wildcard matched "all.bamlist" causing dependency cycle — added `wildcard_constraints`
- **Bug 9**: ANGSD `-rf` + large `-sites` file → chromosome duplication (2–3× output) — workaround: use pass1_snps.txt (smaller) for SAF rules; production TODO to debug nonrepeat_sites.txt
- **Bug 10**: ANGSD `-setMaxDepth` unknown without `-doCounts 1` → empty SAF files — added `-doCounts 1` to SAF rules
- **Bug 11**: `-rf chromosomes.bed` + `-sites` causes tripled positions in SAF → empty 2D SFS → FST fails — removed `-rf` from `angsd_saf`; binary sites index alone restricts to chromosomes correctly
- **Bug 12**: ld_prune OOM for 10-sample test (1.1M SNPs × O(N²)) — added positional thinning fallback for LD files > 2 GB
- **Bug 13**: pruned_keep_snps.txt uses `chr:pos` format; ANGSD beagle uses `chr_pos` — fixed with `.replace(':', '_')`
- **Bug 14**: `generate_report.py` IndexError — PCA eigenvectors have rows only for unrelated samples; fixed by loading `unrelated_samples.txt` to build ordered subset
- **Bug 15**: `generate_report.py` FST midPos parsed as string — fixed with `skiprows=1` in `load_windowed_fst`

### Features added
- QC approval gate (`qc_report` → `qc_approve.py` → `samples_approved.txt`)
- Clone/relatedness filter (`filter_clones` → `clone_approve.py` → `unrelated_samples.txt`)
- Subset beagle to unrelated samples before PCA/admixture/ngsLD
- `scripts/run_pipeline.sh` orchestration script with pre-flight, retry, S3 sync, auto-stop

---

## v0.1 — Local 5-sample test (2026-02-01 → 2026-02-15)

**Environment**: Local workstation, conda env `angsd-pipeline`
**Status**: Complete

### What was validated
- Two-pass ANGSD workflow end-to-end
- ngsRelate, ngsLD, PCAngsd, realSFS all running
- 5-sample results: 2,062,496 SNPs; FST FL vs PA 0.127 unweighted / 0.162 weighted; Tajima's D FL +0.41, PA +0.29

### Bugs fixed
- **Bug 1**: ngsLD `--max_kb_dist` unit mismatch — config was in bp, ngsLD expects kb → 500,000 bp passed as 500,000 km → 168 GB file; fixed with `// 1000`
- **Bug 2**: ngsRelate invalid `-z` flag (not a site-count flag, it's a file-of-IDs flag) → "Number of genomic sites must be provided" silent exit; fixed with `-G` (beagle input) + `-L $N_SITES`
- **Bug 3**: Snakemake subshells don't inherit conda PATH → angsd/samtools not found; fixed with `shell.prefix()` toggled by `local_conda_env` config key
- **Bug 4**: `zcat`/`awk` not installed system-wide; replaced with `python3 -c "import gzip..."` one-liners
- **Bug 5**: `\t`/`\n` in shell-embedded Python strings — Snakemake expands them before bash sees them → SyntaxError; fixed with `\\t` and `\\n`
- **Bug 6**: `process_relatedness` used wrong column names `ida`/`idb`; ngsRelate output uses `a`/`b`



---

## v0.5 — Segments 3–6 production run (2026-03-17 → ongoing)

**Status**: Seg3 complete; Seg4 (diversity/FST) mostly complete; Seg6 (demography) 2D SFS running
**NGSAdmix**: K=1–10, 20 reps — COMPLETE. Best K=2 by Evanno delta-K.
**Moments**: 5 pairwise 2D SFS jobs running (external scripts), watcher auto-submitting fits when ready.

### Segments completed

#### Segment 3 — Population structure (Snakefile.structure)

- PCAngsd PCA: PC1 separates species (45.1%); PC2 captures geography within *A. cervicornis* (3.1%)
- PCAngsd admixture K=2–5: all complete
- **Lineage gate (K=2) passed**: lineageA = *A. palmata* (105), lineageB = *A. cervicornis* (148), 0 admixed
- LD pruning: `all.unrelated.ldpruned.beagle.gz` (29,127 SNPs after pruning)
- **NGSAdmix K=1–10, 20 reps** — completed 2026-03-21
  - Best K=2 by Evanno delta-K (K=1 anchor): delta-K = 116,672,186 — dominant species split signal
  - K=3 secondary (delta-K=2.06): geographic structure within *A. cervicornis* (FL vs PA/BON)
  - K=2 SD near-zero (0.021): all 20 reps converge identically — high confidence
  - 0 admixed individuals at K=2 (identical to PCAngsd; both methods agree completely)
- Admixture figures: bar plots K=2–5, violin plots K=2/K=3 by species×region, PCAngsd vs NGSAdmix comparison panels

#### Segment 4 — Diversity and FST (Snakefile.diversity)

- Individual SAF (per-sample heterozygosity): all 253 samples complete
- Merge SAF, 1D SFS, diversity (π, θ_W, Tajima's D): all 15 population groups complete
- FST global (realSFS fst): completed comparisons:

| Comparison | FST (unweighted) | FST (weighted) |
|------------|-----------------|----------------|
| lineageA FL vs lineageA BON (*Apal*) | 0.0077 | 0.0533 |
| lineageB PA vs lineageB BON (*Acer*) | 0.0086 | 0.0933 |
| lineageB BON vs lineageA BON (species) | 0.0194 | 0.4556 |
| lineageB FL vs lineageB BON (*Acer*) | pending | pending |
| lineageB FL vs lineageB PA (*Acer*) | pending | pending |
| lineageA FL vs lineageA BON (repeat) | pending | pending |

- lineageB_FL vs lineageA_FL (inter-species): **deliberately skipped** — timed out at 24h
  repeatedly; scientifically redundant given 0 admixed individuals (split is obvious, FST≈0.45
  same as BON species comparison)

### Segment 6 — Demographic inference (prototype, external scripts)

Five pairwise 2D SFS comparisons running via manual SLURM scripts:
- `acer_pa_vs_bon` (Acer PA n=39, BON n=24)
- `apal_vs_acer_bon` (Apal BON n=24, Acer BON n=24)
- `acer_fl_vs_acer_pa` (FL n=85, PA n=39) — with `-maxIter 100` to cap EM
- `acer_fl_vs_acer_bon` (FL n=85, BON n=24) — with `-maxIter 100`
- `apal_fl_vs_apal_bon` (FL n=76, BON n=24) — with `-maxIter 100`

`watch_and_submit_moments.sh` polls every 30 min and auto-submits moments fitting when each
2D SFS file becomes non-zero. Projection JSONs auto-created for FL comparisons (project to 40
haploids to reduce compute while retaining segregating sites).

Inter-species FL comparison (lineageB_FL vs lineageA_FL) skipped from moments — too large,
repeatedly times out at 24h wall time, and scientifically not needed (deep species split).

### Bugs fixed (this session)

- **Bug 30**: `merge_saf` used `-o` flag (old realSFS syntax) and `runtime=30`; correct flag is
  `-outnames`; runtime increased to 480 min; stale cached job caused crash 2026-03-18
- **Bug 31**: `fst_index` race condition — submitted before upstream `realsfs_2d` wrote output;
  resolved by pipeline restart; no code change needed
- **Bug 32**: `angsd_saf_demography`, `realsfs_demography`, `realsfs_dadi` runtime too short (480 min);
  increased to 2880 min in Snakefile.demography
- **Bug 33**: `run_moments_2pop.py` used `moments.Spectrum.from_file()` expecting dadi histogram format
  but realSFS dadi writes a per-site data table; fixed: added `load_realsfs_2dsfs()` +
  `--sfs-format realsfs|dadi` and `--n-ind N1 N2` CLI args; default is now `realsfs`
- **Bug 34**: NGSAdmix shell command used `-CPU {threads}` — NGSAdmix flag is `-P`, not `-CPU`;
  caused `Unknown arg:-CPU` exit with no output; fixed in `Snakefile.structure`
- **Bug 35**: NGSAdmix output declared as `.Q` but NGSAdmix writes `.qopt`; Snakemake waited for
  `.Q`, never found it, marked job failed; fixed: output changed to `.qopt` in `ngsadmix_rep`
  and input updated in `ngsadmix_best`
- **Bug 36**: `/work` not included in Singularity bind mounts; NGSAdmix inside container could not
  read input beagle from `/work/vollmer/...`; fixed: added `--bind /work` to `singularity-args`
  in `profiles/discovery/config.yaml`
- **Bug 37**: `latency-wait: 120` too short for NFS `/work` filesystem — NGSAdmix wrote outputs
  successfully but Snakemake couldn't detect them in time, marking jobs failed; increased to 300
- **Bug 38**: Cancelling SLURM jobs caused Snakemake to delete their partial output files; on
  relaunch Snakemake resubmitted everything; resolved by clearing stale locks
  (`rm -f .snakemake/locks/*`) before relaunch
- **Bug 39**: Sample ordering in `plot_admix_compare.py` — `sort_idx` captured *after*
  `reset_index(drop=True)`, giving `[0,1,...,252]` (no reorder); mixed species in bar plots;
  fixed: capture `sort_idx = meta_sorted.index.values` *before* `reset_index`
- **Bug 40**: `lineage` column conflict in admixture compare script — `samples_RR.csv` already has
  a `lineage` column; merging with `lineage_assignments.txt` created `lineage_x`/`lineage_y`;
  fixed: drop `lineage` from samples before merge
- **Bug 41**: Admixture violin plots entirely empty — `lineage_assignments.txt` read with default
  CSV separator (comma) but file is tab-delimited; entire tab-delimited row became dict key,
  all `lineage_map.get()` lookups returned None; fixed: use `meta['species']` column directly
  (Acervicornis/Apalmata) for group assignment — no separate file needed
- **Bug 42**: `realsfs_2d` for lineageB_FL vs lineageA_FL timed out at 24h wall time (twice);
  decision: skip this comparison — scientifically redundant given 0 admixed individuals
- **Bug 43**: NGSAdmix K=3 showed multiple likelihood modes (2 reps at −3,923,784 vs 8 at
  −3,835,871 with only 10 reps); resolved by increasing to K=1–10 with 20 reps each for
  publication quality; K=2 dominates after including K=1 as Evanno anchor
- **Bug 44**: SLURM email spam from pipeline failures — `slurm_extra` had `--mail-type=END,FAIL`;
  changed to `--mail-type=NONE` in `profiles/discovery/config.yaml`

### Profile and config changes (profiles/discovery/config.yaml)

```yaml
# Added /work bind mount (Bug 36)
singularity-args: "--bind /projects --bind /home --bind /scratch --bind /tmp --bind /work"
# Increased NFS latency tolerance (Bug 37)
latency-wait: 300   # was 120
# Disabled SLURM email notifications (Bug 44)
slurm_extra: "'--mail-type=NONE'"   # was '--mail-type=END,FAIL --mail-user=...'
```

### config/config.yaml changes

```yaml
max_k: 10                  # was 5
admixture_replicates: 20   # was 10
```

### New / modified scripts (all in /work/vollmer/acropora_genomics/scripts/)

- **`plot_ngsadmix.py`** — updated for K=1–10, 20 reps; Evanno delta-K with K=1 anchor;
  generates `ngsadmix_deltaK.png`, `ngsadmix_admixture_K2_K5.png`, `ngsadmix_K2_K5_compact.png`
- **`plot_admix_compare.py`** — PCAngsd vs NGSAdmix comparison (Bug 39/40 fixed); generates
  `admix_compare_K2.png`, `admix_compare_K2_K5.png`
- **`plot_admix_violin.py`** — admixture violin plots by species×region for PCAngsd and NGSAdmix
  K=2/K=3; generates `admixture_violin_pcangsd.png`, `admixture_violin_ngsadmix.png`,
  `admixture_violin_compare.png` (2×2 panel)
- **`run_moments_2pop.py`** — `load_realsfs_2dsfs()` for flat realSFS 2D format; `--sfs-format`,
  `--n-ind` CLI args; projects to smaller n before fitting (default 40 haploids for FL comps)
- **`run_2dsfs.sh`** — external SLURM: pairwise 2D SFS via realSFS (non-FL, no maxIter cap,
  24h limit)
- **`run_2dsfs_fl.sh`** — external SLURM: FL pairwise 2D SFS with `-maxIter 100` to cap EM
- **`watch_and_submit_moments.sh`** — polls every 30 min, auto-submits moments fitting when
  each 2D SFS becomes non-zero; running as nohup, logs to `logs/moments_watcher.log`

---

## Orbicella seg2a — clone gate (2026-04-03)

**Working dir**: `/work/vollmer/orbicella_genomics/`
**Pipeline**: `/projects/vollmer/RRorbicella_angsd/`
**Samples entering gate**: 140 (63 *O. annularis*, 57 *O. faveolata*, 21 *O. franksi*)

### SNP filters applied
- min_maf: 0.05
- min_ind: 112 (80% of 140)
- Pass 1 SNPs: 28,063,923
- Pass 2 SNPs (after filter): ~7.2M

### Clone pairs detected (6 pairs, 6 excluded, 134 retained)

| Pair | KING | Excluded | Reason |
|------|------|----------|--------|
| Oann_72 / Oann_73 | 0.4949 | Oann_72 | lower depth (22.9x vs 39.9x) |
| Oann_75 / Oann_77 | 0.4959 | Oann_77 | lower depth (28.9x vs 39.1x) |
| Oann_84 / Oann_85 | 0.4949 | Oann_85 | lower depth (40.2x vs 47.1x) |
| Oann_95 / Oann_96 | 0.4955 | Oann_95 | lower depth (29.8x vs 31.0x) |
| Ofav_106 / Ofra_26 | 0.4929 | **Ofra_26** | **MISLABEL** — see below |
| Ofav_92 / Ofav_93 | 0.4853 | Ofav_92 | lower depth (29.6x vs 32.8x) |

### Sample mislabel: Ofra_26

Ofra_26 (labeled *O. franksi*, FL_RR, Upper Keys / East Turtle Shoal) is a confirmed *O. faveolata* mislabel.

Evidence from KING relatedness matrix:
- Mean KING vs *O. faveolata* (n=57): −0.038 (within-species background)
- Mean KING vs *O. franksi* (n=20, excl. self): −0.521 (between-species background)
- KING vs Ofav_106: +0.493 (clone/identical)

Ofra_26 clusters unambiguously with *O. faveolata*, not *O. franksi*. It is the same individual as Ofav_106 (a true *O. faveolata* from FL_RR). Decision: exclude Ofra_26 (mislabel), retain Ofav_106 (confirmed Ofav).

**Action**: Flag in sample metadata. Notify data provider if possible.
