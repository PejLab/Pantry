import yaml
from pathlib import Path

config = yaml.safe_load(open('TURNAP/config.yml'))
fastq_map = Path(config['fastq_map'])
fastq_dir = Path(config['fastq_dir'])
read_length = config['read_length']
paired_end = config['paired_end']
project_dir = Path(config['project_dir'])
ref_genome = Path(config['ref_genome'])
ref_anno = Path(config['ref_anno'])
ref_dir = project_dir / 'reference'
threads = config['threads']
phenotypes = config['phenotypes']

# Get samples, preserving order for outputs in case it's meaningful
with open(fastq_map, 'r') as f:
    tmp = [line.split('\t')[-1] for line in f.read().splitlines()]
    samples = []
    for s in tmp:
        if s not in samples:
            samples.append(s)

outputs = []
for pheno, params in phenotypes.items():
    for f in params['files']:
        outputs.append(project_dir / f)

include: 'steps/align.smk'
for phenotype in phenotypes:
    include: f'steps/{phenotype}.smk'

rule all:
    input:
        outputs
