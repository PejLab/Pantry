"""Get gene feature bins from GTF file

Extract exons from annotation, add noncoding regions, and split into bins for
computing latent features.
"""

import argparse
from pathlib import Path
from gtfparse import read_gtf
import numpy as np
import pandas as pd

def add_noncoding_regions(anno: pd.DataFrame) -> pd.DataFrame:
    """Add noncoding regions to the annotations for one gene
    
    Inserts an intron between each pair of consecutive exons, and adds
    fixed-length upstream and downstream regions.
    """
    strand = anno.iloc[0].strand
    assert strand in {'+', '-'}
    assert anno.start.unique().shape[0] == anno.shape[0]
    features, starts, ends = [], [], []
    features.append('upstream' if strand == '+' else 'downstream')
    starts.append(anno.iloc[0].start - 1000)
    ends.append(anno.iloc[0].start - 1)
    for i in range(len(anno) - 1):
        features.append('intron')
        starts.append(anno.iloc[i].end + 1)
        ends.append(anno.iloc[i + 1].start - 1)
    features.append('downstream' if strand == '+' else 'upstream')
    starts.append(anno.iloc[-1].end + 1)
    ends.append(anno.iloc[-1].end + 1000)
    nc = pd.DataFrame({
        'seqname': anno.iloc[0].seqname,
        'start': starts,
        'end': ends,
        'strand': strand,
        'gene_id': anno.iloc[0].gene_id,
        'feature': features,
    })
    anno = pd.concat([anno, nc])
    # anno = anno.sort_values(by='start')  # Sort all together for efficiency
    return anno

def split_region(region: pd.DataFrame, n_bins: int) -> pd.DataFrame:
    """Split a region into a fixed number of equal-sized bins"""
    assert region.shape[0] == 1
    region = list(region.itertuples(index=False))[0]
    starts = np.linspace(region.start, region.end + 1, n_bins + 1, dtype=int)
    bins = pd.DataFrame({
        'seqname': region.seqname,
        'start': starts[:-1],
        'end': starts[1:] - 1,
        'strand': region.strand,
        'gene_id': region.gene_id,
        'feature': region.feature,
    })
    assert bins.iloc[0].start == region.start
    assert bins.iloc[-1].end == region.end
    assert bins.start.unique().shape[0] == bins.shape[0]
    return bins

def split_regions(anno: pd.DataFrame, n_bins: int) -> pd.DataFrame:
    """Split each region into a fixed number of equal-sized bins"""
    # Filter out regions that are too small to split:
    n_regions = anno.shape[0]
    anno = anno[anno.end - anno.start >= n_bins]
    n_removed = n_regions - anno.shape[0]
    if n_removed > 0:
        print(f'Removed {n_removed} features that were too small to split')
    bins = anno.groupby(['gene_id', 'start']).apply(lambda x: split_region(x, n_bins))
    bins = bins.reset_index(drop=True)
    return bins

parser = argparse.ArgumentParser(description='Get gene feature bins from GTF file')
parser.add_argument('-g', '--gtf', type=Path, required=True, help='GTF file')
parser.add_argument('-c', '--chromosomes', type=Path, required=True, help='Chromosome lengths file, e.g. chrNameLength.txt from STAR index, to sort chromosomes.')
parser.add_argument('-n', '--n-bins', type=int, default=10, help='Number of bins to split each feature (exon, intron, etc.) into')
parser.add_argument('-o', '--output', type=Path, required=True, help='Output file (BED)')
args = parser.parse_args()

anno = read_gtf(args.gtf)
anno = anno.loc[anno['feature'] == 'exon', :]
anno = anno.loc[anno['source'] == 'ensembl_havana', :]
# anno = anno.iloc[:1000, :]

# For each gene_id, choose the transcript_id with the most exons:
tx_ids = anno.groupby('gene_id').apply(lambda x: x['transcript_id'].value_counts().index[0])
# Get exons in the selected transcript_ids:
anno = anno.loc[anno['transcript_id'].isin(tx_ids), :]

anno = anno[['seqname', 'start', 'end', 'strand', 'gene_id', 'feature']]
anno = anno.sort_values(by=['gene_id', 'start'])

# For each gene, add introns, upstream, and downstream regions:
anno = anno.groupby('gene_id').apply(add_noncoding_regions)
anno = anno.reset_index(drop=True)
anno = anno.sort_values(by=['gene_id', 'start'])
bins = split_regions(anno, n_bins=args.n_bins)

# Sort bins using chromosome order from the chromosome lengths file:
chrom = pd.read_csv(
    args.chromosomes,
    sep='\t',
    header=None,
    names=['chrom', 'length'],
    dtype={'chrom': str, 'length': int}
)
chrom = list(chrom['chrom'])
bins['seqname'] = pd.Categorical(bins['seqname'], categories=chrom, ordered=True)
bins = bins.sort_values(by=['seqname', 'start'])

# Prepare BED format:
bins['start'] = bins['start'] - 1
bins['name'] = bins['gene_id'] + '_' + bins['start'].astype(str)
bins = bins[['seqname', 'start', 'end', 'name', 'feature', 'strand']]
bins.to_csv(args.output, sep='\t', index=False, header=False)
