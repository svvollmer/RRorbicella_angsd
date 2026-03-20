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

## v0.5 — Segments 3-6 production run (2026-03-17 to ongoing)

**Status**: Seg3 (structure) complete; Seg4 (diversity/FST) running; Seg6 (demography) running

### Segments completed / in progress

- **Seg3 (Snakefile.structure)**: PCAngsd PCA + admixture K2-K5 complete
  - K=2 gate passed: lineageA = *A. palmata* (105 samples), lineageB = *A. cervicornis* (148 samples), 0 admixed
  - LD pruning complete: all.unrelated.ldpruned.beagle.gz
  - NGSAdmix K2-K5 (10 replicates each) launched 2026-03-20, in progress
- **Seg4 (Snakefile.diversity)**: running — diversity/FST across all 15 population groups
  - Completed: individual SAF (per-sample het), merge_saf (most groups), 1D SFS (all groups),
    diversity theta/pi/Tajimas D (15 groups), 3 of 6 FST global comparisons
  - Pending: merge_saf FL_RR (retrying), fst_index FL pairs, heterozygosity report
- **Seg6 (Snakefile.demography)**: running — 5 pairwise comparisons (K2/K3 paths)
  - Completed: bamlists, projections, all SAF files, 3 of 5 1D SFS, 2 of 5 dadi SFS
  - Pending: acer_fl + apal_fl 1D SFS (long-running SLURM jobs)

### Bugs fixed

- **Bug 30**: merge_saf used -o flag (old syntax) and runtime=30 min; should be -outnames and runtime=480;
  fixed in Snakefile.diversity; stale cached job caused crash on 2026-03-18
- **Bug 31**: fst_index race condition — submitted before upstream realsfs_2d wrote output;
  realSFS reported file looks empty; resolved by pipeline restart; no code change needed
- **Bug 32**: angsd_saf_demography, realsfs_demography, realsfs_dadi runtime too short (480 min);
  increased to 2880 min (48 h) in Snakefile.demography
- **Bug 33**: run_moments_2pop.py used moments.Spectrum.from_file() which expects dadi SFS histogram format,
  but realSFS dadi outputs a per-site data table — format mismatch (ValueError on load);
  fixed by adding load_realsfs_2dsfs() and --sfs-format realsfs|dadi + --n-ind N1 N2 CLI args;
  default is now realsfs format

### New / modified scripts

- **run_moments_2pop.py**: added load_realsfs_2dsfs() for realSFS flat 2D SFS; added --sfs-format
  and --n-ind args; backward-compatible via --sfs-format dadi
- **Prototype SLURM scripts** (manual, not yet in pipeline):
  - run_2dsfs.sh: pairwise 2D SFS via realSFS with singularity bind mount
  - run_moments_test.sh: prototype moments fitting on completed dadi comparisons

### FST global results (completed comparisons)

| Comparison | Fst (unweighted) | Fst (weighted) |
|------------|-----------------|----------------|
| lineageA FL vs lineageA BON (*Apal*) | 0.0077 | 0.0533 |
| lineageB PA vs lineageB BON (*Acer*) | 0.0086 | 0.0933 |
| lineageB BON vs lineageA BON (species) | 0.0194 | 0.4556 |
