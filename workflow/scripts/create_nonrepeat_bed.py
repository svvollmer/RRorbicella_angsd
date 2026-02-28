#!/usr/bin/env python3
"""
Extract non-repeat regions from soft-masked reference genome.
Soft-masked = lowercase bases are repeats, uppercase are non-repeats.

Much faster than the bedtools approach - processes entire sequences at once.
"""

import sys
import argparse
from pathlib import Path

def parse_fasta(fasta_file):
    """Generator that yields (chrom, sequence) from FASTA."""
    chrom = None
    seq = []
    
    with open(fasta_file) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                if chrom is not None:
                    yield chrom, ''.join(seq)
                chrom = line[1:].split()[0]  # Get chrom name, drop description
                seq = []
            else:
                seq.append(line)
        
        if chrom is not None:
            yield chrom, ''.join(seq)

def find_uppercase_regions(sequence, min_length=100):
    """
    Find runs of uppercase (non-repeat) bases.
    
    Returns list of (start, end) tuples (0-based, half-open).
    min_length: minimum region size to report (bp)
    """
    regions = []
    in_upper = False
    start = 0
    
    for i, base in enumerate(sequence):
        if base.isupper():
            if not in_upper:
                start = i
                in_upper = True
        else:
            if in_upper:
                length = i - start
                if length >= min_length:
                    regions.append((start, i))
                in_upper = False
    
    # Handle case where sequence ends in uppercase
    if in_upper:
        length = len(sequence) - start
        if length >= min_length:
            regions.append((start, len(sequence)))
    
    return regions

def main():
    parser = argparse.ArgumentParser(description='Extract non-repeat regions from soft-masked FASTA')
    parser.add_argument('fasta', help='Input soft-masked FASTA file')
    parser.add_argument('chromosomes_bed', help='BED file of chromosomes to process')
    parser.add_argument('output', help='Output BED file')
    parser.add_argument('--min-length', type=int, default=100,
                       help='Minimum region length (bp) [default: 100]')
    parser.add_argument('--merge-gap', type=int, default=50,
                       help='Merge regions within this distance (bp) [default: 50]')
    
    args = parser.parse_args()
    
    # Load chromosome list
    chromosomes = set()
    with open(args.chromosomes_bed) as f:
        for line in f:
            if line.strip():
                chromosomes.add(line.split()[0])
    
    print(f"Processing {len(chromosomes)} chromosomes...", file=sys.stderr)
    
    # Process each chromosome
    total_regions = 0
    total_bp = 0
    
    with open(args.output, 'w') as out:
        for chrom, sequence in parse_fasta(args.fasta):
            if chrom not in chromosomes:
                continue
            
            print(f"  {chrom}: {len(sequence):,} bp", file=sys.stderr)
            
            # Find uppercase regions
            regions = find_uppercase_regions(sequence, args.min_length)
            
            # Optional: merge nearby regions
            if args.merge_gap > 0 and len(regions) > 1:
                merged = [regions[0]]
                for start, end in regions[1:]:
                    prev_start, prev_end = merged[-1]
                    if start - prev_end <= args.merge_gap:
                        # Merge with previous
                        merged[-1] = (prev_start, end)
                    else:
                        merged.append((start, end))
                regions = merged
            
            # Write BED
            for start, end in regions:
                out.write(f"{chrom}\t{start}\t{end}\n")
                total_bp += (end - start)
            
            total_regions += len(regions)
            print(f"    -> {len(regions):,} non-repeat regions", file=sys.stderr)
    
    print(f"\nTotal: {total_regions:,} regions, {total_bp:,} bp", file=sys.stderr)
    print(f"Output: {args.output}", file=sys.stderr)

if __name__ == '__main__':
    main()
