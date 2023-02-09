"""Assemble data into an RNA phenotype BED file"""

import argparse
from pathlib import Path
from gtfparse import read_gtf
import numpy as np
import pandas as pd

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
    anno = anno[['#chrom', 'chromStart', 'chromEnd', 'gene_id']]
    # Rename columns for tensorQTL:
    anno.columns = ['#chr', 'start', 'end', 'gene_id']
    return anno

def load_exons(ref_anno: Path) -> pd.DataFrame:
    """Load exon annotations
    
    Returns exons with start and end oriented on the gene's strand.
    """
    anno = read_gtf(ref_anno)
    anno = anno.loc[anno['feature'] == 'exon', :]
    anno['chrom'] = anno['seqname']
    anno['exonStart'] = np.where(anno['strand'] == '+', anno['start'], anno['end'])
    anno['exonEnd'] = np.where(anno['strand'] == '+', anno['end'], anno['start'])
    return anno[['gene_id', 'chrom', 'exonStart', 'exonEnd']]

def transcript_to_gene_map(ref_anno: Path) -> pd.DataFrame:
    """Load transcript IDs and corresponding gene IDs from GTF file"""
    anno = read_gtf(ref_anno)
    anno = anno.loc[anno['feature'] == 'transcript', :]
    return anno[['gene_id', 'transcript_id']]

def map_introns_to_genes(introns: list, exons: pd.DataFrame) -> pd.DataFrame:
    """Map de novo splice junctions to genes based on known exon boundaries"""
    exons['exonStart'] = exons['exonStart'].astype(str)
    exons['exonEnd'] = exons['exonEnd'].astype(str)
    df = pd.DataFrame({'intron': introns})
    df[['chrom', 'chr_start', 'chr_end', 'clu', 'cluster', 'strand']] = df['intron'].str.split(r':|_', expand=True)
    df['start'] = np.where(df['strand'] == '+', df['chr_start'], df['chr_end'])
    df['end'] = np.where(df['strand'] == '+', df['chr_end'], df['chr_start'])
    start_matches = df.merge(
        exons[['chrom', 'exonEnd', 'gene_id']].rename(columns={'exonEnd': 'start'}),
        on=['chrom', 'start'],
        how='inner'
    )
    end_matches = df.merge(
        exons[['chrom', 'exonStart', 'gene_id']].rename(columns={'exonStart': 'end'}),
        on=['chrom', 'end'],
        how='inner'
    )
    df = pd.concat([start_matches, end_matches])
    clust_genes = df.groupby('cluster', group_keys=False).agg({'gene_id': pd.Series.unique})
    clust_genes = clust_genes.reset_index().explode('gene_id')
    return clust_genes

def load_kallisto(sample_ids: list, kallisto_dir: Path, units: str) -> pd.DataFrame:
    """Assemble kallisto est_counts or tpm outputs into a table"""
    counts = []
    for i, sample in enumerate(sample_ids):
        fname = kallisto_dir / sample / 'abundance.tsv'
        d = pd.read_csv(fname, sep='\t', index_col='target_id')
        d = d[units] # est_counts or tpm
        d.name = sample
        # Store in lists and concat all at once to avoid 'PerformanceWarning: DataFrame is highly fragmented' warning
        counts.append(d)
    return pd.concat(counts, axis=1)

def assemble_alt_TSS_polyA(sample_ids: list, group1_dir: Path, group2_dir: Path, units: str, ref_anno: Path, bed: Path, min_frac: float = 0.05, max_frac: float = 0.95):
    """Assemble txrevise-based kallisto outputs into BED file
    
    By default, txrevise produces two sets of annotations per annotation type.
    Relative TSS/polyA site usage should be calculated within each set, and
    can then be combined to produce a single BED file.
    """
    df1 = load_kallisto(sample_ids, group1_dir, units)
    df2 = load_kallisto(sample_ids, group2_dir, units)
    # This assumes target_id is e.g. {gene_id}.grp_1.downstream.{transcript_id}
    df1['gene_id'] = df1.index.str.split('.').str[0]
    df2['gene_id'] = df2.index.str.split('.').str[0]
    # Calculate proportion of each transcript in each gene_id:
    gene_ids = df1['gene_id'] # groupby/apply removes gene_id, so save it to add back
    df1 = df1.groupby('gene_id', group_keys=False).apply(lambda x: x / x.sum(axis=0))
    # Remove sites with mean relative usage < `min_frac` or > `max_frac`:
    df1 = df1[(df1.mean(axis=1) >= min_frac) & (df1.mean(axis=1) <= max_frac)]
    df1 = df1.join(gene_ids, how='left')
    gene_ids = df2['gene_id']
    df2 = df2.groupby('gene_id', group_keys=False).apply(lambda x: x / x.sum(axis=0))
    # Remove sites with mean relative usage < `min_frac` or > `max_frac`:
    df2 = df2[(df2.mean(axis=1) >= min_frac) & (df2.mean(axis=1) <= max_frac)]
    df2 = df2.join(gene_ids, how='left')
    df = pd.concat([df1, df2], axis=0)
    df.index = df.index.rename('phenotype_id')
    df = df.reset_index()
    # Use gene's TSS for all of its isoforms:
    anno = load_tss(ref_anno)
    df = anno.merge(df.reset_index(), on='gene_id', how='inner')
    df = df[['#chr', 'start', 'end', 'phenotype_id'] + sample_ids]
    df.to_csv(bed, sep='\t', index=False, float_format='%g')

def assemble_expression(sample_ids: list, kallisto_dir: Path, units: str, ref_anno: Path, bed_iso: Path, bed_gene: Path, min_count: int = 10, min_frac: float = 0.05, max_frac: float = 0.95):
    """Assemble kallisto est_counts or tpm outputs into isoform- and gene-level BED files
    
    Isoform values are normalized to relative abundance in each gene. Isoforms
    with fewer than `min_count` reads in every sample are excluded
    (`est_counts` read counts are always used for this filtering).
    """
    df_iso = load_kallisto(sample_ids, kallisto_dir, units)

    # Record isoforms with any value >= `min_count` to keep:
    if units == 'est_counts':
        df_counts = df_iso
    else:
        df_counts = load_kallisto(sample_ids, kallisto_dir, 'est_counts')
    iso_enough_counts = df_counts[(df_counts >= min_count).any(axis=1)].index

    df_iso.index = df_iso.index.rename('transcript_id')
    gene_map = transcript_to_gene_map(ref_anno).set_index('transcript_id')
    df_iso = df_iso.join(gene_map, how='inner')
    # Also get gene-level expression:
    df_gene = df_iso.reset_index().drop('transcript_id', axis=1).groupby('gene_id', group_keys=True).sum()
    # Calculate proportion of each transcript in each gene_id:
    df_iso = df_iso.groupby('gene_id', group_keys=False).apply(lambda x: x / x.sum(axis=0))
    # Remove isoforms with fewer than `min_count` reads in every sample:
    df_iso = df_iso[df_iso.index.isin(iso_enough_counts)]
    # Remove isoforms with mean relative abundance < `min_frac` or > `max_frac`:
    df_iso = df_iso[(df_iso.mean(axis=1) >= min_frac) & (df_iso.mean(axis=1) <= max_frac)]
    # groupby/apply removes gene_id, so add them back
    df_iso = df_iso.join(gene_map, how='left')

    anno = load_tss(ref_anno)
    df_iso = anno.merge(df_iso.reset_index(), on='gene_id', how='inner')
    df_iso['phenotype_id'] = df_iso['gene_id'] + ':' + df_iso['transcript_id']
    df_iso = df_iso[['#chr', 'start', 'end', 'phenotype_id'] + sample_ids]
    df_iso.to_csv(bed_iso, sep='\t', index=False, float_format='%g')
    df_gene = anno.merge(df_gene.reset_index(), on='gene_id', how='inner')
    df_gene = df_gene.rename(columns={'gene_id': 'phenotype_id'})
    df_gene = df_gene[['#chr', 'start', 'end', 'phenotype_id'] + sample_ids]
    df_gene.to_csv(bed_gene, sep='\t', index=False, float_format='%g')

def assemble_splicing(counts: Path, ref_anno: Path, bed: Path, min_frac: float = 0.05, max_frac: float = 0.95):
    """Convert leafcutter output into splicing BED file"""
    df = pd.read_csv(counts, sep=' ')
    sample_ids = list(df.columns)
    cluster = df.index.str.extract(r'clu_(\d+)_', expand=False)
    df = df.groupby(cluster, group_keys=False).apply(lambda g: g / g.sum(axis=0))
    # Remove junctions with mean relative usage < `min_frac` or > `max_frac`:
    df = df[(df.mean(axis=1) >= min_frac) & (df.mean(axis=1) <= max_frac)]
    df['cluster'] = df.index.str.extract(r'clu_(\d+)_', expand=False)
    df.index = df.index.rename('intron')
    df = df.reset_index()
    exons = load_exons(ref_anno)
    genes = map_introns_to_genes(df['intron'], exons)
    df = df.merge(genes, on='cluster', how='left')
    anno = load_tss(ref_anno)
    df = anno.merge(df, on='gene_id', how='inner')
    df['phenotype_id'] = df['gene_id'] + ':' + df['intron']
    df = df[['#chr', 'start', 'end', 'phenotype_id'] + sample_ids]
    df.to_csv(bed, sep='\t', index=False, float_format='%g')

def load_featureCounts(sample_ids: list, counts_dir: Path, feature: str, min_count: int = 10) -> pd.DataFrame:
    """Assemble featureCounts outputs into a table"""
    counts = []
    for i, sample in enumerate(sample_ids):
        fname_ex = counts_dir / f'{sample}.{feature}.counts.txt'
        d = pd.read_csv(fname_ex, sep='\t', index_col='Geneid', skiprows=1)
        d = d.iloc[:, 5]
        d.name = sample
        d[d < min_count] = np.nan
        # Store in lists and concat all at once to avoid 'PerformanceWarning: DataFrame is highly fragmented' warning
        counts.append(d)
    return pd.concat(counts, axis=1)

def assemble_stability(sample_ids: list, stab_dir: Path, ref_anno: Path, bed: Path):
    """Assemble exon to intron read ratios into mRNA stability BED file"""
    exon = load_featureCounts(sample_ids, stab_dir, 'constit_exons')
    intron = load_featureCounts(sample_ids, stab_dir, 'introns')
    genes = exon.index[np.isin(exon.index, intron.index)]
    # genes = set(exon.index).intersection(intron.index)
    assert exon.loc[genes, :].index.equals(intron.loc[genes, :].index)
    assert exon.columns.equals(intron.columns)
    df = exon.loc[genes, :] / intron.loc[genes, :]
    df = df[df.isnull().mean(axis=1) <= 0.5]
    anno = load_tss(ref_anno)
    anno = anno.rename(columns={'gene_id': 'phenotype_id'})
    df.index = df.index.rename('phenotype_id')
    df = anno.merge(df.reset_index(), on='phenotype_id', how='inner')
    df.to_csv(bed, sep='\t', index=False, float_format='%g')

parser = argparse.ArgumentParser(description='Assemble data into an RNA phenotype BED file')
parser.add_argument('--type', choices=['alt_TSS_polyA', 'expression', 'splicing', 'stability'], required=True, help='Type of data to assemble')
parser.add_argument('--input', type=Path, required=False, help='Input file, for phenotype groups in which all data is already in a single file')
parser.add_argument('--input-dir', type=Path, required=False, help='Directory containing input files, for phenotype groups with per-sample input files')
parser.add_argument('--input-dir2', type=Path, required=False, help='Second input directory, for phenotype groups with per-sample input files in two directories')
parser.add_argument('--samples', type=Path, required=False, help='File listing sample IDs, for phenotype groups with per-sample input files')
parser.add_argument('--ref_anno', type=Path, required=True, help='Reference annotation file')
parser.add_argument('--output', type=Path, required=True, help='Output file ("*.bed")')
parser.add_argument('--output2', type=Path, required=False, help='Second output file ("*.bed") for phenotype groups with two output files, e.g. log2 and tpm expression')
args = parser.parse_args()

if args.samples is not None:
    samples = pd.read_csv(args.samples, sep='\t', header=None)[0].tolist()

if args.type == 'alt_TSS_polyA':
    assemble_alt_TSS_polyA(samples, args.input_dir, args.input_dir2, 'tpm', args.ref_anno, args.output, min_frac=0.05, max_frac=0.95)
elif args.type == 'expression':
    assemble_expression(samples, args.input_dir, 'tpm', args.ref_anno, args.output, args.output2, min_count=10, min_frac=0.05, max_frac=0.95)
elif args.type == 'splicing':
    assemble_splicing(args.input, args.ref_anno, args.output, min_frac=0.05, max_frac=0.95)
elif args.type == 'stability':
    assemble_stability(samples, args.input_dir, args.ref_anno, args.output)
