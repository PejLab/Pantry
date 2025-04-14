"""Run plink and gcta to calculate heritability per phenotype."""

import argparse
from pathlib import Path
import subprocess
import pandas as pd

def read_hsq_stats(hsq_file):
    """Read the heritability statistics from the output file of gcta."""
    hsq, se, pval = 'NA', 'NA', 'NA'
    if not Path(hsq_file).exists():
        return hsq, se, pval
    with open(hsq_file, 'r') as f:
        lines = f.read().splitlines()
    lines = [line.split('\t') for line in lines]
    for line in lines:
        if line[0] == 'V(G)/Vp':
            hsq = line[1]
            se = line[2]
        elif line[0] == 'Pval':
            pval = line[1]
    return hsq, se, pval

parser = argparse.ArgumentParser(description='Run plink and gcta to calculate heritability per phenotype.')
parser.add_argument('--bed', '-i', type=Path, help='The input BED file.')
parser.add_argument('--geno', type=str, help='Prefix of plink (bed, bim, fam) genotype files.')
parser.add_argument('--covar', type=Path, help='Covariates file (TSV, plink format, no header, first two columns are family and individual ID).')
parser.add_argument('--chrom', type=str, help='If provided, only include phenotypes on this chromosome.')
parser.add_argument('--batch', type=int, help='If provided, only include phenotypes in this batch. Currently this must be an integer 0-9, and phenotypes for which the last digit of the TSS position matches will be included.')
parser.add_argument('--grm-dir', type=Path, help='The directory that contains or will contain GRM files.')
parser.add_argument('--tmp-dir', type=Path, help='The directory that intermediate files will be written to.')
parser.add_argument('--output', '-o', type=Path, help='The output file (TSV).')
args = parser.parse_args()

df = pd.read_csv(args.bed, sep='\t', index_col='phenotype_id', dtype={'#chr': str, 'start': int, 'end': int})
if args.chrom is not None:
    df = df[df['#chr'] == args.chrom]
if args.batch is not None:
    df = df[df['end'].astype('string').str[-1] == str(args.batch)]
pos = df[['#chr', 'end']].copy()
pos['pos_min'] = pos['end'] - 250000
pos['pos_min'] = pos['pos_min'].mask(pos['pos_min'] < 1, 1)
pos['pos_max'] = pos['end'] + 250000
df = df.drop(['#chr', 'start', 'end'], axis=1)
df.index.name = None
df = df.T
df.index.name = 'sample'
phenos = list(df.columns)
df = df.reset_index()
df.insert(0, 'family', 0)  # Family must match that in the genotype files.

stats = {}
for pheno in phenos:
    prefix = args.tmp_dir / pheno.replace(':', '_')
    ph = df[['family', 'sample', pheno]]
    ph.to_csv(f'{prefix}.pheno', sep='\t', index=False, header=False)
    chrom, pos_min, pos_max = pos.loc[pheno, ['#chr', 'pos_min', 'pos_max']]
    grm = args.grm_dir / f'{chrom}_{pos_min}_{pos_max}'
    if not Path(f'{grm}.grm.bin').exists():
        cmd = [
            'plink2',
            '--bfile', args.geno,
            '--chr', chrom,
            '--from-bp', str(pos_min),
            '--to-bp', str(pos_max),
            '--pheno', f'{prefix}.pheno',
            '--keep', f'{prefix}.pheno',
            '--mind', '0.1', # Samples with higher missingness rate are removed. Added to avoid "Sample(s) present with no genotype data"
            '--make-grm-bin',
            '--out', grm,
            '--allow-extra-chr',
            '--allow-no-sex',
        ]
        subprocess.run(cmd, stdout=subprocess.DEVNULL)  # plink already creates a log
    cmd = [
        'gcta64',
        '--grm', grm,
        '--pheno', f'{prefix}.pheno',
        '--qcovar', args.covar,
        '--reml',
        '--reml-no-constrain',
        '--out', f'{prefix}.covar_h2',
    ]
    subprocess.run(cmd, stdout=subprocess.DEVNULL)  # gcta already creates a log
    stats[pheno] = read_hsq_stats(f'{prefix}.covar_h2.hsq')

with open(args.output, 'w') as f:
    f.write(f'phenotype_id\thsq\tse\tpval\n')
    for pheno in phenos:
        hsq, se, pval = stats[pheno]
        f.write(f'{pheno}\t{hsq}\t{se}\t{pval}\n')
