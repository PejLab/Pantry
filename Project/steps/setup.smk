"""This script is loaded near the beginning of the main phenotyping Snakefile.

It validates and processes the configuration and sets global variables.
"""

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
        'modality_groups',
        'genome_size',
    ]
    for field in fields:
        if field not in config.keys():
            raise Exception(f'{field} not in config file.')
    if not config['paired_end']:
        for field in ['fragment_length_mean', 'fragment_length_mean']:
            if field not in config.keys():
                raise Exception(f'{field} estimate must be specified in config for single-end reads.')

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
        'intermediate_ref_dir',
    ]
    for path in paths:
        if path in config.keys():
            config[path] = Path(config[path]).expanduser()

    samples = pd.read_csv(config['samples_file'], sep='\t', header=None, dtype=str)[0].tolist()
    config['samples'] = samples

    if 'intermediate_dir' not in config:
        config['intermediate_dir'] = Path('intermediate')
    if 'intermediate_ref_dir' not in config:
        config['intermediate_ref_dir'] = config['intermediate_dir'] / 'reference'

    config['fastq_map'] = load_fastq_map(config['fastq_map'], config['fastq_dir'], config['paired_end'])

validate_config(config)
process_config(config)

interm_dir = config['intermediate_dir']
ref_dir = config['intermediate_ref_dir']
output_dir = Path('output')

samples_file = config['samples_file']
samples = config['samples']

paired_end = config['paired_end']
read_length = config['read_length']
fastq_dir = config['fastq_dir']
fastq_map = config['fastq_map']
if not paired_end:
    fragment_length_mean = config['fragment_length_mean']
    fragment_length_sd = config['fragment_length_sd']

ref_genome = config['ref_genome']
ref_anno = config['ref_anno']
ref_cdna = ref_dir / 'cDNA.fa.gz'

modality_groups = config['modality_groups']
genome_size = config['genome_size'] # TODO: compute from fasta file

outputs = []
for modality_group, params in modality_groups.items():
    for f in params['files']:
        outputs.append(output_dir / f)

