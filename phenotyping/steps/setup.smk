"""This script is loaded near the beginning of the main phenotyping Snakefile.

It validates and processes the configuration and sets global variables.
"""

def load_fastq_map(map_file: Path, fastq_dir: Path, paired_end=None) -> dict:
    """Load the FASTQ file paths for each sample.
    
    Returns a dictionary where each value is a tuple (paths, is_paired).
    paths is either a list (single-end) or tuple of two lists (paired-end).
    is_paired is a boolean indicating if the sample is paired-end.
    """
    paths = {}
    with open(map_file, 'r') as f:
        for line in f.read().splitlines():
            fields = line.split('\t')
            sample_id = fields[-1]
            is_paired = len(fields) == 3
            
            if sample_id not in paths:
                paths[sample_id] = (([], []), is_paired) if is_paired else ([], is_paired)
            elif paths[sample_id][1] != is_paired:
                raise ValueError(f"Sample {sample_id} has mixed paired-end and single-end reads")
            
            if is_paired:
                paths[sample_id][0][0].append(str(fastq_dir / fields[0]))
                paths[sample_id][0][1].append(str(fastq_dir / fields[1]))
            else:
                paths[sample_id][0].append(str(fastq_dir / fields[0]))
    
    # Validate against config if specified
    if paired_end is not None:
        for sample_id, (_, is_paired) in paths.items():
            if is_paired != paired_end:
                raise ValueError(f"Sample {sample_id} has {'paired' if is_paired else 'single'}-end reads, but config specifies {'paired' if paired_end else 'single'}-end")
    
    return paths

def validate_config(config: dict):
    """Check that the config has all required fields with valid values"""
    if len(config.keys()) == 0:
        raise Exception('No config file provided, or it is empty.')
    required = [
        'fastq_dir',
        'fastq_map',
        'ref_genome',
        'ref_anno',
        'samples_file',
        'modality_groups',
        'genome_size',
    ]
    for field in required:
        if field not in config:
            raise KeyError(f'Required config field missing: {field}')
    
    # paired_end is now optional
    if 'paired_end' in config:
        if not isinstance(config['paired_end'], bool):
            raise ValueError('paired_end must be True or False')
    
    if not isinstance(config['read_length'], int):
        raise ValueError('read_length must be an integer')
    
    # fragment length parameters only required for single-end samples
    if 'paired_end' in config and not config['paired_end']:
        if 'fragment_length_mean' not in config:
            raise KeyError('fragment_length_mean required for single-end data')
        if 'fragment_length_sd' not in config:
            raise KeyError('fragment_length_sd required for single-end data')

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
if paired_end is not None and not paired_end:
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

