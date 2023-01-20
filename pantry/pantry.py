"""Generate molecular phenotypes from RNA-Seq"""

from pathlib import Path
import pandas as pd

def load_fastq_map(map_file: Path, fastq_dir: Path, paired_end) -> dict:
    """Load the FASTQ file paths for each sample.

    If paired_end is True, each dictionary value is a tuple of two lists of
    paths, otherwise each value is a list of paths.
    """
    paths = {}
    with open(map_file, 'r') as f:
        for line in f.read().splitlines():
            if paired_end:
                fastq1, fastq2, sample_id = line.split('\t')
                if sample_id not in paths:
                    paths[sample_id] = ([], [])
                paths[sample_id][0].append(str(fastq_dir / fastq1))
                paths[sample_id][1].append(str(fastq_dir / fastq2))
            else:
                fastq, sample_id = line.split('\t')
                if sample_id not in paths:
                    paths[sample_id] = []
                paths[sample_id].append(str(fastq_dir / fastq))
    return paths

def validate_config(config: dict):
    """Validate the configuration"""
    if len(config.keys()) == 0:
        raise Exception('No config file provided.')
    fields = [
        'paired_end',
        'read_length',
        'fastq_dir',
        'fastq_map',
        'ref_genome',
        'ref_anno',
        'samples_file',
        'phenotypes',
        'genome_size',
    ]
    for field in fields:
        if field not in config.keys():
            raise Exception(f'{field} not in config file.')

def process_config(config: dict):
    """Prepare user config for use in the pipeline.
    
    Updates some values in place to expand paths and parse numbers.
    """
    config['genome_size'] = int(float(config['genome_size'])) # float() handles scientific notation
    paths = [
        'fastq_dir',
        'fastq_map',
        'ref_genome',
        'ref_anno',
        'samples_file',
        'intermediate_dir',
    ]
    for path in paths:
        if path in config.keys():
            config[path] = Path(config[path]).expanduser()

    samples = pd.read_csv(config['samples_file'], sep='\t', header=None)[0].tolist()
    config['samples'] = samples

    if 'intermediate_dir' not in config:
        config['intermediate_dir'] = Path('intermediate')

    config['fastq_map'] = load_fastq_map(config['fastq_map'], config['fastq_dir'], config['paired_end'])
