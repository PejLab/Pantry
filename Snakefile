import yaml
from pathlib import Path
import pandas as pd

config = yaml.safe_load(open('TURNAP/config.yml'))
fastq_map = Path(config['fastq_map'])
fastq_dir = Path(config['fastq_dir'])
read_length = config['read_length']
paired_end = bool(config['paired_end'])
project_dir = Path(config['project_dir'])
ref_genome = Path(config['ref_genome'])
ref_anno = Path(config['ref_anno'])
ref_dir = project_dir / 'reference'
threads = config['threads']
phenotypes = config['phenotypes']

with open(fastq_map, 'r') as f:
    samples = [line.split('\t')[-1] for line in f.read().splitlines()]

include: 'steps/align.smk'
for step in phenotypes:
    include: f'steps/{step}.smk'

rule all:
    input:
        'test/expression.log2.bed',
        'test/expression.tpm.bed',
