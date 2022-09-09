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
    anno.columns = ['#chr', 'start', 'end', 'phenotype_id']
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

def load_retros(retro_anno: Path) -> pd.DataFrame:
    """Load retroelement annotations
    
    For retroelements with multiple entries (corresponding to multiple repeats),
    returns the median of their starts as its 'TSS'.
    """
    anno = read_gtf(retro_anno)
    anno['chromEnd'] = np.where(anno['strand'] == '+', anno['start'], anno['end'])
    anno = anno[['gene_id', 'seqname', 'chromEnd']]
    anno = anno.groupby(['gene_id', 'seqname']).median().astype(int).reset_index()
    # assert len(anno['gene_id'].unique()) == len(anno['gene_id'])
    anno['chromStart'] = anno['chromEnd'] - 1  # BED coordinates are 0-based
    anno = anno.rename(columns={'seqname': '#chrom', 'gene_id': 'name'})
    anno = anno.sort_values(['#chrom', 'chromStart'])
    return anno[['#chrom', 'chromStart', 'chromEnd', 'name']]

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

def assemble_ATS_APA(sample_ids: list, kallisto_dir: Path, units: str, ref_anno: Path, bed: Path):
    """Assemble txrevise-based kallisto outputs into isoform-level BED file"""
    df = []
    for i, sample in enumerate(sample_ids):
        fname = kallisto_dir / sample / 'abundance.tsv'
        d = pd.read_csv(fname, sep='\t', index_col='target_id')
        d = d[units] # est_counts or tpm
        d.name = sample
        df.append(d)
    df = pd.concat(df, axis=1)
    df.index = df.index.rename('phenotype_id')
    # This assumes target_id is e.g. {gene_id}.grp_1.downstream.{transcript_id}
    df['gene_id'] = df.index.str.split('.').str[0]
    df = df.reset_index()
    # Use gene's TSS for all of its isoforms:
    anno = load_tss(ref_anno)
    anno = anno.rename(columns={'phenotype_id': 'gene_id'})
    df = anno.merge(df.reset_index(), on='gene_id', how='inner')
    df = df[['#chr', 'start', 'end', 'phenotype_id'] + sample_ids]
    df.to_csv(bed, sep='\t', index=False, float_format='%g')

def assemble_expression(sample_ids: list, kallisto_dir: Path, units: str, ref_anno: Path, bed_iso: Path, bed_gene: Path):
    """Assemble kallisto est_counts or tpm outputs into isoform- and gene-level BED files"""
    df = []
    for i, sample in enumerate(sample_ids):
        fname = kallisto_dir / sample / 'abundance.tsv'
        d = pd.read_csv(fname, sep='\t', index_col='target_id')
        d = d[units] # est_counts or tpm
        d.name = sample
        df.append(d)
    df = pd.concat(df, axis=1)
    # This assumes that Ensembl cDNA was used for kallisto:
    # Strip off the version number from the transcript_id:
    df['transcript_id'] = df.index.str.split('.').str[0]
    # Add gene_id and use its TSS for all of its isoforms:
    gene_map = transcript_to_gene_map(ref_anno)
    df = df.merge(gene_map, on='transcript_id', how='inner')
    # Also get gene-level expression:
    df_gene = df.drop('transcript_id', axis=1).groupby('gene_id').sum()
    anno = load_tss(ref_anno)
    anno = anno.rename(columns={'phenotype_id': 'gene_id'})
    df = anno.merge(df.reset_index(), on='gene_id', how='inner')
    df['phenotype_id'] = df['gene_id'] + ':' + df['transcript_id']
    df = df[['#chr', 'start', 'end', 'phenotype_id'] + sample_ids]
    df.to_csv(bed_iso, sep='\t', index=False, float_format='%g')
    df_gene = anno.merge(df_gene.reset_index(), on='gene_id', how='inner')
    df_gene = df_gene.rename(columns={'gene_id': 'phenotype_id'})
    df_gene = df_gene[['#chr', 'start', 'end', 'phenotype_id'] + sample_ids]
    df_gene.to_csv(bed_gene, sep='\t', index=False, float_format='%g')

def assemble_retroelements(sample_ids: list, retro_dir: Path, retro_anno: Path, bed: Path):
    """Assemble telescope output files into retroelement expression BED file"""
    for i, sample in enumerate(sample_ids):
        fname = retro_dir / f'{sample}-telescope_report.tsv'
        d = pd.read_csv(fname, sep='\t', index_col='transcript', skiprows=1)
        d = d.rename(columns={'final_count': sample})[[sample]]
        if i == 0:
            df = d
        else:
            df = df.join(d, how='outer')
    df = df.fillna(0).astype(int)
    df = df.loc[np.sum(df > 0, axis=1) > 0, :]
    anno = load_retros(retro_anno)
    # Rename columns for tensorQTL:
    anno.columns = ['#chr', 'start', 'end', 'phenotype_id']
    df.index = df.index.rename('phenotype_id')
    df = anno.merge(df.reset_index(), on='phenotype_id', how='inner')
    df.to_csv(bed, sep='\t', index=False, float_format='%g')

def assemble_splicing(counts: Path, ref_anno: Path, bed: Path):
    """Convert leafcutter output into splicing BED file"""
    df = pd.read_csv(counts, sep=' ')
    sample_ids = list(df.columns)
    df.index = df.index.rename('intron')
    df = df.reset_index()
    # df = df.rename(columns={'chrom': 'name'})
    df['cluster'] = df['intron'].str.extract(r'clu_(\d+)_', expand=False)
    for col in sample_ids:
    #     df[col] = pd.eval(df[col])  # convert fraction string to float
        df[col] = df.groupby('cluster', group_keys=False).apply(lambda g: g[col] / g[col].sum())
    exons = load_exons(ref_anno)
    genes = map_introns_to_genes(df['intron'], exons)
    df = df.merge(genes, on='cluster', how='left')
    anno = load_tss(ref_anno)
    anno = anno.rename(columns={'phenotype_id': 'gene_id'})
    df = anno.merge(df, on='gene_id', how='inner')
    df['phenotype_id'] = df['gene_id'] + ':' + df['intron']
    df = df[['#chr', 'start', 'end', 'phenotype_id'] + sample_ids]
    df.to_csv(bed, sep='\t', index=False, float_format='%g')

def assemble_stability(sample_ids: list, stab_dir: Path, ref_anno: Path, bed: Path):
    """Assemble exon to intron read ratios into mRNA stability BED file"""
    exon = []
    intron = []
    for i, sample in enumerate(sample_ids):
        fname_ex = stab_dir / f'{sample}.constit_exons.counts.txt'
        d_ex = pd.read_csv(fname_ex, sep='\t', index_col='Geneid', skiprows=1)
        fname_in = stab_dir / f'{sample}.introns.counts.txt'
        d_in = pd.read_csv(fname_in, sep='\t', index_col='Geneid', skiprows=1)
        # if i == 0:
        #     exon = pd.DataFrame(index=d_ex.index)
        #     intron = pd.DataFrame(index=d_in.index)
        # else:
        #     assert d_ex.index.equals(exon.index)
        #     assert d_in.index.equals(intron.index)
        # exon[sample] = d_ex.iloc[:, 5].copy()
        # intron[sample] = d_in.iloc[:, 5].copy()
        # Avoid 'PerformanceWarning: DataFrame is highly fragmented' warning:
        d_ex = d_ex.iloc[:, 5]
        d_in = d_in.iloc[:, 5]
        d_ex.name = sample
        d_in.name = sample
        d_ex[d_ex < 10] = np.nan
        d_in[d_in < 10] = np.nan
        exon.append(d_ex)
        intron.append(d_in)
    exon = pd.concat(exon, axis=1)
    intron = pd.concat(intron, axis=1)
    genes = exon.index[np.isin(exon.index, intron.index)]
    # genes = set(exon.index).intersection(intron.index)
    assert exon.loc[genes, :].index.equals(intron.loc[genes, :].index)
    assert exon.columns.equals(intron.columns)
    df = exon.loc[genes, :] / intron.loc[genes, :]
    df = df[df.isnull().mean(axis=1) <= 0.5]
    anno = load_tss(ref_anno)
    df.index = df.index.rename('phenotype_id')
    df = anno.merge(df.reset_index(), on='phenotype_id', how='inner')
    df.to_csv(bed, sep='\t', index=False, float_format='%g')

parser = argparse.ArgumentParser(description='Assemble data into an RNA phenotype BED file')
parser.add_argument('--type', choices=['ATS_APA', 'expression', 'retroelements', 'splicing', 'stability'], required=True, help='Type of data to assemble')
parser.add_argument('--input', type=Path, required=False, help='Input file, for phenotype groups in which all data is already in a single file')
parser.add_argument('--input-dir', type=Path, required=False, help='Directory containing input files, for phenotype groups with per-sample input files')
parser.add_argument('--samples', type=Path, required=False, help='File listing sample IDs, for phenotype groups with per-sample input files')
parser.add_argument('--ref_anno', type=Path, required=True, help='Reference annotation file')
parser.add_argument('--output', type=Path, required=True, help='Output file ("*.bed")')
parser.add_argument('--output2', type=Path, required=False, help='Second output file ("*.bed") for phenotype groups with two output files, e.g. log2 and tpm expression')
args = parser.parse_args()

if args.samples is not None:
    samples = pd.read_csv(args.samples, sep='\t', header=None)[0].tolist()

if args.type == 'ATS_APA':
    assemble_ATS_APA(samples, args.input_dir, 'tpm', args.ref_anno, args.output)
elif args.type == 'expression':
    assemble_expression(samples, args.input_dir, 'tpm', args.ref_anno, args.output, args.output2)
elif args.type == 'retroelements':
    assemble_retroelements(samples, args.input_dir, args.ref_anno, args.output)
elif args.type == 'splicing':
    assemble_splicing(args.input, args.ref_anno, args.output)
elif args.type == 'stability':
    assemble_stability(samples, args.input_dir, args.ref_anno, args.output)
