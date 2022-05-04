"""Extract coverage for a set of regions from a bedgraph file"""

import argparse
from pathlib import Path
import pandas as pd
from pyBedGraph import BedGraph

def mean_coverage(bedgraph: BedGraph, regions: pd.DataFrame) -> pd.DataFrame:
    """Compute mean coverage for each region on the same chromosome"""
    seqnames = list(regions['seqname'].unique())
    assert len(seqnames) == 1
    regions = regions.copy()
    bedgraph.load_chrom_data(seqnames[0])
    # There's an error (cython?) if the dtype is int64
    means = bedgraph.stats(
        stat='mean',
        start_list=regions['start'].astype('int32').values - 1,  # I think GTF is 1-based and bedgraph is 0-based
        end_list=regions['end'].astype('int32').values,
        chrom_name=seqnames[0],
    )
    # print(means[:10])
    regions['region'] = regions.index
    regions['mean_coverage'] = means
    return regions[['region', 'mean_coverage']]
    # return means

parser = argparse.ArgumentParser(description='Extract coverage for a set of regions from a bedgraph file')
parser.add_argument('-b', '--bedgraph', type=Path, required=True, help='Bedgraph file')
parser.add_argument('-r', '--regions', type=Path, required=True, help='Regions file')
parser.add_argument('-c', '--chr-lengths', type=Path, required=True, help='Chromosome lengths file')
parser.add_argument('-o', '--output', type=Path, required=True, help='Output file')
args = parser.parse_args()

regions = pd.read_csv(args.regions, sep='\t', dtype={'seqname': str})
# regions = pd.read_csv(args.regions, sep='\t', dtype={'seqname': str}, index_col=['seqname', 'start'])
# regions = regions.iloc[:1000, :]
bg = BedGraph(args.chr_lengths, args.bedgraph, ignore_missing_bp=False)
# regions['mean_coverage'] = regions.groupby(regions['seqname']).apply(lambda x: mean_coverage(bg, x))
regions.index = regions['seqname'] + '_' + regions['start'].astype(str)
# df = regions.copy()
# print(regions)
df = regions.groupby(regions['seqname']).apply(lambda x: mean_coverage(bg, x))
df = df.reset_index(level=0, drop=True)
# print(df)
 # Join according to region label, e.g. X_100600386
regions['mean_coverage'] = df['mean_coverage']  # Fix to avoid duplicate indices by replacing seqname in region labels with gene_id
# print(regions)
# print(regions.loc['1_169827125', :])
# print(df.loc['1_169827125', :])
# assert regions['start'].equals(df['start'])
# assert regions.index.identical(df.index)
# print(len(means))
# print(means[:5])
regions = regions[['seqname', 'start', 'mean_coverage']]
regions.to_csv(args.output, sep='\t', index=False, float_format='%g')
