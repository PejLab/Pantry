#!/usr/bin/env python3
"""
Prepare Pantry phenotypes for Pheast

This script converts sample IDs to individual IDs in Pantry phenotype files.
For multi-tissue datasets like GTEx, RNA-seq samples often have IDs that differ
from individual IDs, and genetic analyses require the individual IDs to match genotypes.

Usage:
    python prepare_phenotypes.py \
        --indir ../phenotyping/output/ \
        --outdir input/phenotypes/ \
        --map input/sample_map.txt \
        --individuals input/individuals.txt
"""

import os
import argparse
import sys
import pandas as pd
import shutil
import glob
import subprocess


def load_sample_map(map_file):
    """Load sample to individual mapping from file."""
    if map_file is None:
        return None
    map_df = pd.read_csv(map_file, sep=r'\s+', header=None) 
    if map_df.shape[1] < 2:
        sys.exit(f"Error: Map file {map_file} must have at least 2 columns (sample_id, individual_id) separated by whitespace")
    sample_to_individual = dict(zip(map_df.iloc[:, 0], map_df.iloc[:, 1]))
    return sample_to_individual


def load_individuals(individuals_file):
    """Load list of individuals to include from file."""
    if individuals_file is None:
        return None
    with open(individuals_file, 'r') as f:
        individuals = set(line.strip() for line in f if line.strip())
    return individuals


def subset_bed_samples(df, ids):
    """Subset BED file to only include samples in ids."""
    n_before = df.shape[1] - 4
    df = df[list(df.columns[:4]) + [col for col in df.columns[4:] if col in ids]]
    n_after = df.shape[1] - 4
    if n_before != n_after:
        print(f"  Subset from {n_before} to {n_after} samples")
    return df


parser = argparse.ArgumentParser(description='Prepare Pantry phenotypes for Pheast.')
parser.add_argument('--indir', required=True, help='Directory containing input phenotype BED files. All files matching *.bed.gz will be processed. Any associated *.phenotype_groups.txt files will be copied to the output directory.')
parser.add_argument('--outdir', required=True, help='Directory to save output phenotype BED files. Must not exist.')
parser.add_argument('--map', help='Tab-delimited file with no header mapping sample IDs (first column) to individual IDs (second column). Sample IDs not in the mapping will be excluded. Sample IDs not in the phenotype data will be ignored.')
parser.add_argument('--individuals', help='File with list of individuals to include. If not provided, all individuals will be included. Listed individuals not in the phenotype data will be ignored.')
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=False)

sample_map = load_sample_map(args.map)
individuals = load_individuals(args.individuals)

# Find phenotype files in input directory
bed_files = glob.glob(os.path.join(args.indir, "*.bed.gz"))

if not bed_files:
    sys.exit(f"Error: No .bed.gz files found in {args.indir}")

print(f"Found {len(bed_files)} phenotype files")

# Process each phenotype file
for input_file in bed_files:
    basename = os.path.basename(input_file)
    output_file = os.path.join(args.outdir, basename)
    
    print(f"Processing {basename}...")
    df = pd.read_csv(input_file, sep='\t', dtype=str)

    if sample_map is not None:
        # Filter to samples present in mapping
        print(f"  Converting sample IDs to individual IDs")
        df = subset_bed_samples(df, sample_map.keys())
        ids = [sample_map[id] for id in df.columns[4:]]
        # Check for duplicates
        unique_ids = set(ids)
        if len(ids) != len(unique_ids):
            duplicates = [id for id in unique_ids if ids.count(id) > 1]
            raise AssertionError("Duplicate individual IDs found:", duplicates)
        df.columns = list(df.columns[:4]) + ids

    ## Subset to genotyped individuals
    if individuals is not None:
        print(f"  Subsetting to specified individuals")
        df = subset_bed_samples(df, individuals)

    # Write to uncompressed temporary file then compress with bgzip
    temp_file = output_file.replace('.gz', '')
    df.to_csv(temp_file, sep='\t', index=False)
    
    # Compress with bgzip
    print(f"  Compressing with bgzip")
    subprocess.run(['bgzip', '-f', temp_file], check=True)
    
    # Check for and copy phenotype groups file if it exists
    groups_file = input_file.replace('.bed.gz', '.phenotype_groups.txt')
    if os.path.exists(groups_file):
        output_groups_file = os.path.join(args.outdir, os.path.basename(groups_file))
        shutil.copy2(groups_file, output_groups_file)
        print(f"  Copied phenotype groups file: {output_groups_file}")

    # Create .tbi index file if input had one
    tbi_file = input_file + '.tbi'
    if os.path.exists(tbi_file):
        print(f"  Creating tabix index for {output_file}")
        os.system(f"tabix -p bed {output_file}")

print("Sample-to-individual mapping complete")
