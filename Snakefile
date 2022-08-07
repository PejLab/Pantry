from pathlib import Path
from gtfparse import read_gtf
import numpy as np
import pandas as pd
import turnap

turnap.validate_config(config)

fastq_map = Path(config['fastq_map'])
fastq_dir = Path(config['fastq_dir'])
samples_file = Path(config['samples_file'])
read_length = config['read_length']
paired_end = config['paired_end']
project_dir = Path(config['project_dir'])
ref_genome = Path(config['ref_genome'])
ref_anno = Path(config['ref_anno'])
retro_anno = Path(config['retro_anno'])
ref_dir = project_dir / 'reference'
code_dir = Path(config['code_dir'])
threads = config['threads']
phenotypes = config['phenotypes']
genome_size = config['genome_size'] # TODO: compute from fasta file
# Currently only used for heritability and qtl:
chroms = [str(x) for x in config['chroms']]
geno_prefix = config['geno_prefix']
covar_file = Path(config['covar_file'])

# Get samples, preserving order for outputs in case it's meaningful
# with open(fastq_map, 'r') as f:
#     tmp = [line.split('\t')[-1] for line in f.read().splitlines()]
#     samples = []
#     for s in tmp:
#         if s not in samples:
#             samples.append(s)
samples = pd.read_csv(samples_file, sep='\t', header=None)[0].tolist()

outputs = []
for pheno, params in phenotypes.items():
    for f in params['files']:
        outputs.append(project_dir / f)
    # Used only for QTL mapping for now:
    params['grouped'] = any(['phenotype_groups' in f for f in params['files']])

# These steps are short and will not be submitted as cluster jobs:
localrules:
    index_bed,

include: 'steps/align.smk'
for phenotype in phenotypes:
    include: f'steps/{phenotype}.smk'
include: 'steps/heritability.smk'
include: 'steps/qtl.smk'

rule all:
    input:
        outputs,
        # expand(project_dir / 'heritability' / '{pheno}_hsq.tsv', pheno=phenotypes.keys()),
        # expand(project_dir / 'qtl' / '{pheno}.cis_independent_qtl.txt.gz', pheno=phenotypes.keys()),

rule index_bed:
    input:
        project_dir / '{pheno}.bed.gz',
    output:
        project_dir / '{pheno}.bed.gz.tbi'
    shell:
        'tabix {input}'
