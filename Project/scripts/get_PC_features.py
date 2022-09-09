"""Run PCA on feature bin coverage and create BED file"""

import argparse
from pathlib import Path
from gtfparse import read_gtf
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

def norm_counts_over_gene(counts: pd.DataFrame) -> pd.DataFrame:
    """For one gene, scale coverage per sample to have the same total per sample
    
    Set the total to the median of the totals per sample.
    """
    total = counts.sum(axis=0).median()
    counts = counts.div(counts.sum(axis=0), axis=1) * total
    return counts

def pca(counts: pd.DataFrame, var_expl: float = 0.95, n_pcs: int = None) -> pd.DataFrame:
    """Perform PCA on normalized coverage for one gene
    
    Returns enough PCs to explain var_expl variance, or supply n_pcs to return a fixed number of PCs.
    """
    n_comp = var_expl if n_pcs is None else min(args.n_pcs, *counts.shape)
    samples = counts.columns
    x = counts.values.T
    # Center and scale the data:
    x = (x - x.mean(axis=0)) / x.std(axis=0)
    pca = PCA(n_components=n_comp)
    mat = pca.fit_transform(x).T
    df = pd.DataFrame(
        mat,
        index=[f'PC{i + 1}' for i in range(mat.shape[0])],
        columns=samples,
    )
    df.index.set_names('PC', inplace=True)
    return df

def load_tss(ref_anno: Path) -> pd.DataFrame:
    """Load TSS annotations from GTF file
    
    Returns TSS as the first four columns of the BED format, meaning the
    coordinates are 0-based and chromEnd is just chromStart + 1.
    """
    anno = read_gtf(ref_anno)
    anno = anno.loc[anno['feature'] == 'gene', :]
    anno['chromEnd'] = np.where(anno['strand'] == '+', anno['start'], anno['end'])
    anno['chromStart'] = anno['chromEnd'] - 1  # BED coordinates are 0-based
    anno['#chrom'] = anno['seqname']
    anno = anno.sort_values(['#chrom', 'chromStart'])
    return anno[['#chrom', 'chromStart', 'chromEnd', 'gene_id']]

parser = argparse.ArgumentParser(description='Run PCA on feature bin coverage and create BED file')
parser.add_argument('-i', '--input', type=Path, required=True, help='File containing paths to all input BED files. Base file name before first "." is sample ID')
parser.add_argument('-g', '--gtf', type=Path, required=True, help='Reference annotation GTF file to add TSS to BED output')
parser_n_comp = parser.add_mutually_exclusive_group()
parser_n_comp.add_argument('-v', '--var_expl', type=float, default=0.95, help='Return enough PCs to explain this fraction of variance')
parser_n_comp.add_argument('-n', '--n-pcs', type=int, help='Number of PCs to keep per gene')
parser.add_argument('-o', '--output', type=Path, required=True, help='Output file (BED)')
args = parser.parse_args()

infiles = pd.read_csv(args.input, header=None, names=['path'])
for i, fname in enumerate(infiles['path']):
    d = pd.read_csv(fname, sep='\t', header=None, usecols=[3, 6], names=['region', 'count'], index_col='region')
    if i == 0:
        df = pd.DataFrame(index=d.index)
        assert df.index.is_unique
    else:
        assert d.index.equals(df.index)
    sample = Path(fname).name.split('.')[0]
    df[sample] = d['count']

# Split region label to get multiindex:
df.index = df.index.str.split('_', 1, expand=True)
df.index.set_names(['gene_id', 'start'], inplace=True)

# Normalize coverage to ignore total expression
df = df.groupby('gene_id').apply(norm_counts_over_gene)
# Remove rows with NaN, for which median coverage was 0:
df = df.dropna(axis=0)
# Normalize coverage to mean coverage within each region
# Get region lengths to normalize counts:
regions = pd.read_csv(infiles['path'][0], sep='\t', header=None, usecols=[1, 2, 3], names=['start', 'end', 'region'], index_col='region')
# Split index on underscore into two columns:
regions.index = regions.index.str.split('_', 1, expand=True)
regions.index.set_names(['gene_id', 'start'], inplace=True)
regions['length'] = regions['end'] - regions['start']
df = df.div(regions['length'], axis=0)

# Variance stabilization:
df = np.sqrt(df)
# Run PCA on each gene:
df = df.loc[df.std(axis=1) > 0, :]
df = df.groupby('gene_id').apply(lambda x: pca(x, args.var_expl, args.n_pcs))

samples = list(df.columns)
df = df.reset_index()
anno = load_tss(args.gtf)
df = anno.merge(df, on='gene_id', how='inner')
df['name'] = df['gene_id'] + ':' + df['PC'].astype(str)
df = df[['#chrom', 'chromStart', 'chromEnd', 'name'] + samples]
# Rename columns for tensorQTL:
df.columns = ['#chr', 'start', 'end', 'phenotype_id'] + samples
df.to_csv(args.output, sep='\t', index=False, float_format='%g')
