#!/usr/bin/env bash
# =============================================================================
# run_kmer_compare.sh — Reference-free k-mer sharing between Orbicella species
#
# Tests whether low Oann-Ofra segregating sites (33k vs 311k for other pairs)
# reflects genuine sequence divergence or reference genome bias.
#
# Approach:
#   1. Stream reads from 5 genetically confirmed samples per species
#      directly into jellyfish count (no intermediate FASTQ files)
#   2. Dump k-mers with min count >= MINCOV to exclude sequencing errors
#   3. Compare pairwise k-mer sharing (Jaccard + containment) with Python
#
# Expected result if reference bias is the cause:
#   Oann-Ofra Jaccard ≈ Ofav-Ofra Jaccard (similar sharing despite fewer SNPs)
# Expected result if genuine divergence:
#   Oann-Ofra Jaccard << Ofav-Ofra Jaccard
#
# Samples used (5 per species, genetically confirmed, not mislabeled):
#   Oann: Oann_101, Oann_102, Oann_103, Oann_105, Oann_107
#   Ofav: Ofav_100, Ofav_102, Ofav_103, Ofav_106, Ofav_107
#   Ofra: Ofra_10,  Ofra_11,  Ofra_12,  Ofra_13,  Ofra_14
# =============================================================================
#SBATCH --job-name=kmer_orbicella
#SBATCH --partition=short
#SBATCH --time=08:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --output=/work/vollmer/orbicella_genomics/logs/kmer_compare_%j.log

set -euo pipefail

module load jellyfish/2.2.10
module load samtools/1.19.2

BAMDIR=/projects/vollmer/RR_heat-tolerance/Orbicella/2_mapping.bwa
OUTDIR=/work/vollmer/orbicella_genomics/results/kmer
SCRIPTS=/projects/vollmer/RRorbicella_angsd/workflow/scripts
PYTHON=/home/s.vollmer/.conda/envs/snakemake2/bin/python3

mkdir -p "$OUTDIR"

K=21
MINCOV=5   # min k-mer count: filters sequencing errors; real k-mers appear ~150x (5 samples x 30x)
THREADS=16
HASHSIZE=8G  # jellyfish hash; coral genome ~475Mb x 5 samples, canonical k-mers

# ---------------------------------------------------------------------------
# Sample lists (5 per species, genetically confirmed)
# ---------------------------------------------------------------------------
OANN_SAMPLES=(Oann_101 Oann_102 Oann_103 Oann_105 Oann_107)
OFAV_SAMPLES=(Ofav_100 Ofav_102 Ofav_103 Ofav_106 Ofav_107)
OFRA_SAMPLES=(Ofra_10  Ofra_11  Ofra_12  Ofra_13  Ofra_14)

# ---------------------------------------------------------------------------
# Step 1: k-mer counting per species (stream BAM -> jellyfish, no temp files)
# ---------------------------------------------------------------------------
echo "[$(date)] Counting k-mers for Oannularis..."
(for s in "${OANN_SAMPLES[@]}"; do
    samtools fastq "$BAMDIR/${s}.bwa.dedup.clip.bam" 2>/dev/null
done) | jellyfish count \
    -m $K -s $HASHSIZE -t $THREADS \
    -C \
    -o "$OUTDIR/oann.jf" \
    /dev/stdin
echo "[$(date)] Oann done. Hash: $(du -sh $OUTDIR/oann.jf)"

echo "[$(date)] Counting k-mers for Ofaveolata..."
(for s in "${OFAV_SAMPLES[@]}"; do
    samtools fastq "$BAMDIR/${s}.bwa.dedup.clip.bam" 2>/dev/null
done) | jellyfish count \
    -m $K -s $HASHSIZE -t $THREADS \
    -C \
    -o "$OUTDIR/ofav.jf" \
    /dev/stdin
echo "[$(date)] Ofav done. Hash: $(du -sh $OUTDIR/ofav.jf)"

echo "[$(date)] Counting k-mers for Ofranksi..."
(for s in "${OFRA_SAMPLES[@]}"; do
    samtools fastq "$BAMDIR/${s}.bwa.dedup.clip.bam" 2>/dev/null
done) | jellyfish count \
    -m $K -s $HASHSIZE -t $THREADS \
    -C \
    -o "$OUTDIR/ofra.jf" \
    /dev/stdin
echo "[$(date)] Ofra done. Hash: $(du -sh $OUTDIR/ofra.jf)"

# ---------------------------------------------------------------------------
# Step 2: Dump k-mers above min coverage threshold
# ---------------------------------------------------------------------------
echo "[$(date)] Dumping k-mers (min count = $MINCOV)..."
for sp in oann ofav ofra; do
    jellyfish dump -c -L $MINCOV "$OUTDIR/${sp}.jf" > "$OUTDIR/${sp}_kmers.txt"
    echo "  $sp: $(wc -l < $OUTDIR/${sp}_kmers.txt) k-mers"
done

# ---------------------------------------------------------------------------
# Step 3: Pairwise k-mer comparison
# ---------------------------------------------------------------------------
echo "[$(date)] Comparing k-mer sets..."
$PYTHON "$SCRIPTS/kmer_compare.py" \
    --oann "$OUTDIR/oann_kmers.txt" \
    --ofav "$OUTDIR/ofav_kmers.txt" \
    --ofra "$OUTDIR/ofra_kmers.txt" \
    --out  "$OUTDIR/kmer_summary.txt"

cat "$OUTDIR/kmer_summary.txt"
echo "[$(date)] Done. Results: $OUTDIR/kmer_summary.txt"
