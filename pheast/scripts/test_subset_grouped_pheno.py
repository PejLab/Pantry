import argparse
from pathlib import Path
import pandas as pd

parser = argparse.ArgumentParser(description="Subset grouped phenotypes to N groups")
parser.add_argument("--in-bed", type=Path, help="Input phenotype BED file")
parser.add_argument("--in-groups", type=Path, help="Input phenotype groups file")
parser.add_argument("--out-bed", type=Path, help="Output phenotype BED file")
parser.add_argument("--out-groups", type=Path, help="Output phenotype groups file")
parser.add_argument("--chrom", type=str, help="Subset to phenotypes on this chromosome")
parser.add_argument("--n-groups", type=int, help="Number of groups to keep")
args = parser.parse_args()

df = pd.read_csv(args.in_bed, sep='\t', dtype={'#chr': str, 'start': int, 'end': int})
groups = pd.read_csv(args.in_groups, sep='\t', names=['phenotype_id', 'gene_id'])

if args.chrom:
    df = df[df['#chr'] == args.chrom]
    groups = groups[groups['phenotype_id'].isin(df['phenotype_id'])]

if args.n_groups:
    keep = groups['gene_id'].unique()[:args.n_groups]
    groups = groups[groups['gene_id'].isin(keep)]
    df = df[df['phenotype_id'].isin(groups['phenotype_id'])]

df.to_csv(args.out_bed, sep='\t', index=False)
groups.to_csv(args.out_groups, sep='\t', index=False)
