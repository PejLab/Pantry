#!/usr/bin/env python3
# Assemble edit level matrix
# Based on sharedsamples_sites_matrix_FastQTL_v8.pl from https://github.com/vargasliqin/GTEx_edQTL

import argparse
import gzip
from pathlib import Path

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='Create a site-by-sample RNA editing ratio matrix from per-sample editing files'
    )
    parser.add_argument('--path_to_edit_files', required=True,
                        help='Directory containing {sample}.rnaeditlevel.tsv.gz files')
    parser.add_argument('--samples_file', required=True,
                        help='File containing sample IDs, one per line')
    parser.add_argument('--output_file', required=True,
                        help='Output TSV matrix with one row per retained edit site and one column per sample')
    parser.add_argument('--min_coverage', type=int, required=True,
                        help='Minimum read coverage required for a site to count as observed in a sample')
    parser.add_argument('--min_samples', type=int, required=True,
                        help='Minimum number of samples in which a site must meet --min_coverage to be retained')
    args = parser.parse_args()

    # Ensure path ends with slash
    path_edit_files = str(Path(args.path_to_edit_files).resolve()) + '/'

    print(f"Path to Edit Files: {path_edit_files}")
    print(f"Samples File: {args.samples_file}")
    print(f"Output File: {args.output_file}")
    print(f"Min Coverage: {args.min_coverage}")
    print(f"Min Samples: {args.min_samples}")

    # Read samples from file
    with open(args.samples_file, 'r') as f:
        samples = [line.strip() for line in f]

    sitehash = {}  # Dictionary to store sample-site ratios
    totalhash = {}  # Dictionary to store site counts
    lvlhash = {}    # Dictionary to store site levels

    # Process each sample
    for sample in samples:
        file_path = f"{path_edit_files}{sample}.rnaeditlevel.tsv.gz"

        print(f"Analyzing: {sample}")
        
        with gzip.open(file_path, 'rt') as f:
            for line in f:
                fields = line.strip().split()
                if fields[0] == '#chrom':
                    continue

                chr_name, pos, cov, edit, lvl = fields
                site = f"{chr_name}_{pos}"
                ratio = f"{edit}/{cov}"

                if int(cov) >= args.min_coverage:
                    if sample not in sitehash:
                        sitehash[sample] = {}
                    sitehash[sample][site] = ratio

                    if site in totalhash:
                        totalhash[site] += 1
                        lvlhash[site] = f"{lvlhash[site]},{ratio}"
                    else:
                        totalhash[site] = 1
                        lvlhash[site] = ratio

    print(f"Output file: {args.output_file}")

    # Check if any sites meet the minimum sample threshold
    sites_to_write = [site for site in totalhash if totalhash[site] >= args.min_samples]
    if not sites_to_write:
        raise ValueError(f"No sites found with at least {args.min_samples} samples meeting the minimum coverage threshold of {args.min_coverage}")

    # Write output matrix
    with open(args.output_file, 'w') as outfile:
        # Write header
        outfile.write("site")
        for sample in samples:
            outfile.write(f"\t{sample}")
        outfile.write("\n")

        # Write data rows
        for site in sites_to_write:
            outfile.write(site)
            for sample in samples:
                if sample in sitehash and site in sitehash[sample]:
                    outfile.write(f"\t{sitehash[sample][site]}")
                else:
                    outfile.write("\t0/0")
            outfile.write("\n")

if __name__ == "__main__":
    main()
