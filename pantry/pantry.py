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

def new_fastq_map(samples: list, fastq_dir: Path, paired_end: bool) -> dict:
    """Get the list(s) of FASTQ file paths to be created.
    
    This is used when FASTQ files are produced from input BAM files.
    If paired_end is True, each dictionary value is a tuple of two lists of
    paths, otherwise each value is a list of paths.
    """
    paths = {}
    for sample_id in samples:
        if paired_end:
            paths[sample_id] = ([], [])
            paths[sample_id][0].append(fastq_dir / f'{sample_id}.PE1.fastq.gz')
            paths[sample_id][1].append(fastq_dir / f'{sample_id}.PE2.fastq.gz')
        else:
            paths[sample_id] = [fastq_dir / f'{sample_id}.SE.fastq.gz']
    return paths

def load_bam_map(map_file: Path, bam_dir: Path) -> dict:
    """Load the BAM file path for each sample."""
    paths = {}
    with open(map_file, 'r') as f:
        for line in f.read().splitlines():
            bam, sample_id = line.split('\t')
            assert sample_id not in paths, "Only one BAM file allowed per sample."
            paths[sample_id] = bam_dir / bam
    return paths

def new_bam_map(samples: list, bam_dir: Path) -> dict:
    """Get the list of BAM file paths to be created.
    
    This is used when BAM files are produced from input FASTQ files.
    """
    paths = {}
    for sample_id in samples:
        paths[sample_id] = bam_dir / f'{sample_id}.bam'
    return paths

def validate_config(config: dict):
    """Validate the configuration"""
    if len(config.keys()) == 0:
        raise Exception('No config file provided.')
    fields = [
        'paired_end',
        'read_length',
        'ref_genome',
        'ref_anno',
        'ref_cdna',
        'samples_file',
        'phenotypes',
        'genome_size',
    ]
    for field in fields:
        if field not in config.keys():
            raise Exception(f'{field} not in config file.')
    if 'fastq_dir' in config.keys():
        if 'fastq_map' not in config.keys():
            raise Exception('Provide fastq_map to map FASTQ files to samples.')
    if 'bam_dir' in config.keys():
        if 'bam_map' not in config.keys():
            raise Exception('Provide bam_map to map BAM files to samples.')
    if 'fastq_dir' not in config.keys() and 'bam_dir' not in config.keys():
        raise Exception('Provide fastq_dir and/or bam_dir in config file.')

def process_config(config: dict):
    """Prepare user config for use in the pipeline.
    
    Updates some values in place to expand paths and parse numbers.
    """
    config['genome_size'] = int(float(config['genome_size'])) # float() handles scientific notation
    paths = [
        'fastq_dir',
        'fastq_map',
        'bam_dir',
        'bam_map',
        'ref_genome',
        'ref_anno',
        'ref_cdna',
        'retro_anno',
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

    if 'fastq_dir' in config:
        config['fastq_map'] = load_fastq_map(config['fastq_map'], config['fastq_dir'], config['paired_end'])
    else:
        config['fastq_dir'] = config['intermediate_dir'] / 'fastq'
        config['fastq_map'] = new_fastq_map(samples, config['fastq_dir'], config['paired_end'])
    if 'bam_dir' in config:
        config['bam_map'] = load_bam_map(config['bam_map'], config['bam_dir'])
    else:
        config['bam_dir'] = config['intermediate_dir'] / 'bam'
        config['bam_map'] = new_bam_map(samples, config['bam_dir'])
