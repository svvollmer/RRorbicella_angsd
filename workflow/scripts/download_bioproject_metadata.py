#!/usr/bin/env python3
"""
Download sample metadata from NCBI BioProject and create samples.csv

Usage:
    python3 download_bioproject_metadata.py PRJNA950067 --output config/samples.csv
    
Requires: biopython
Install: pip install biopython
"""

import argparse
import sys
import csv
from Bio import Entrez

# Set email for NCBI (required)
Entrez.email = "your.email@fau.edu"

def get_bioproject_samples(bioproject_id):
    """Get all SRA run accessions from a BioProject"""
    
    print(f"Fetching samples from BioProject {bioproject_id}...")
    
    # Search for SRA experiments linked to BioProject
    search_handle = Entrez.esearch(
        db="sra",
        term=f"{bioproject_id}[BioProject]",
        retmax=1000,
        usehistory="y"
    )
    search_results = Entrez.read(search_handle)
    search_handle.close()
    
    count = int(search_results["Count"])
    print(f"Found {count} SRA experiments")
    
    if count == 0:
        print(f"No samples found for {bioproject_id}")
        return []
    
    # Fetch full records
    webenv = search_results["WebEnv"]
    query_key = search_results["QueryKey"]
    
    fetch_handle = Entrez.efetch(
        db="sra",
        query_key=query_key,
        WebEnv=webenv,
        rettype="runinfo",
        retmode="text"
    )
    
    # Parse runinfo CSV
    runinfo = fetch_handle.read().decode('utf-8')
    fetch_handle.close()
    
    samples = []
    lines = runinfo.strip().split('\n')
    header = lines[0].split(',')
    
    # Find column indices
    run_idx = header.index('Run')
    sample_idx = header.index('SampleName')
    organism_idx = header.index('ScientificName')
    spots_idx = header.index('spots')
    bases_idx = header.index('bases')
    
    for line in lines[1:]:
        fields = line.split(',')
        if len(fields) < len(header):
            continue
            
        run = fields[run_idx]
        sample_name = fields[sample_idx]
        organism = fields[organism_idx]
        spots = fields[spots_idx]
        bases = fields[bases_idx]
        
        # Parse population from sample name (assumes naming like "Ac_FL_M5")
        parts = sample_name.split('_')
        if len(parts) >= 2:
            if 'FL' in parts[1]:
                population = 'florida'
                region = 'Florida'
            elif 'PA' in parts[1]:
                population = 'panama'
                region = 'Panama'
            else:
                population = 'unknown'
                region = 'Unknown'
        else:
            population = 'unknown'
            region = 'Unknown'
        
        samples.append({
            'sample_id': sample_name,
            'sra_accession': run,
            'population': population,
            'region': region,
            'organism': organism,
            'read_count': spots,
            'bases': bases
        })
    
    return samples

def write_samples_csv(samples, output_file):
    """Write samples to CSV file"""
    
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=[
            'sample_id', 'sra_accession', 'population', 'region'
        ])
        writer.writeheader()
        
        for sample in samples:
            writer.writerow({
                'sample_id': sample['sample_id'],
                'sra_accession': sample['sra_accession'],
                'population': sample['population'],
                'region': sample['region']
            })
    
    print(f"\nWrote {len(samples)} samples to {output_file}")
    
    # Print summary
    from collections import Counter
    pop_counts = Counter(s['population'] for s in samples)
    print(f"\nPopulation breakdown:")
    for pop, count in pop_counts.items():
        print(f"  {pop}: {count} samples")
    
    # Print total data size
    total_bases = sum(int(s['bases']) for s in samples)
    total_gb = total_bases / 1e9
    print(f"\nTotal data: {total_gb:.1f} GB bases")
    print(f"Estimated download: ~{total_gb * 0.3:.1f} GB compressed")

def main():
    parser = argparse.ArgumentParser(description='Download BioProject metadata from NCBI')
    parser.add_argument('bioproject', help='BioProject accession (e.g., PRJNA950067)')
    parser.add_argument('--output', default='config/samples.csv', help='Output CSV file')
    parser.add_argument('--email', default='svollmer@fau.edu', help='Your email for NCBI')
    
    args = parser.parse_args()
    
    # Set email
    Entrez.email = args.email
    
    # Download metadata
    samples = get_bioproject_samples(args.bioproject)
    
    if not samples:
        print("No samples found!")
        sys.exit(1)
    
    # Write CSV
    write_samples_csv(samples, args.output)
    
    print(f"\nReady to run pipeline with {len(samples)} samples!")
    print(f"Command: snakemake --cores 16 --jobs 4")

if __name__ == '__main__':
    main()
