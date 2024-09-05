"""This script is loaded near the beginning of the main Pheast Snakefile.

It validates and processes the configuration and sets global variables.
"""

def validate_config(config: dict):
    """Validate the configuration"""
    if len(config.keys()) == 0:
        raise Exception('No config file provided.')
    fields = [
        'phenotype_dir',
        'geno_prefix',
        'samples_file',
        'modalities',
        'analyses',
    ]
    for field in fields:
        if field not in config.keys():
            raise Exception(f'{field} not in config file.')

def process_config(config: dict):
    """Prepare user config for use in the pipeline.
    
    Updates some values in place to expand paths and parse numbers.
    """
    paths = [
        'phenotype_dir',
        'samples_file',
        'twas_snps',
        'intermediate_dir',
    ]
    for path in paths:
        if path in config.keys():
            config[path] = Path(config[path]).expanduser()

    config['samples'] = pd.read_csv(config['samples_file'], sep='\t', header=None)[0].tolist()

    if 'intermediate_dir' not in config:
        config['intermediate_dir'] = Path('intermediate')

validate_config(config)
process_config(config)

pheno_dir = config['phenotype_dir']
interm_dir = config['intermediate_dir']
output_dir = Path('output')
geno_prefix = config['geno_prefix']

# samples_file = Path(config['samples_file'])
# samples = pd.read_csv(samples_file, sep='\t', header=None)[0].tolist()
samples = config['samples']
geno_samples = pd.read_csv(geno_prefix + '.fam', sep=r'\s+', header=None)[1].tolist()
missing_samples = [s for s in samples if s not in geno_samples]
if len(missing_samples) > 0:
    raise Exception(f'Samples missing from genotypes: {missing_samples}')
    # print(f'WARNING: Samples missing from genotypes: {missing_samples}')

modalities = config['modalities']
for modality, params in modalities.items():
    chroms = pd.read_csv(pheno_dir / f'{modality}.bed.gz', sep='\t', usecols=[0], dtype=str)
    params['chroms'] = chroms.iloc[:, 0].unique().tolist()

analyses = config['analyses']
outputs = []
for analysis, params in analyses.items():
    for f in params['files']:
        outputs.append(expand(output_dir / analysis / f, modality=modalities.keys()))

