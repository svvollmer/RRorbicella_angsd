# Changelog

Development history of the coral ANGSD pipeline from local prototype to AWS production run.

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
