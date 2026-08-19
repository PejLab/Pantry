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

    # 'cross_modality' pseudomodality should only be used for cross-modality cis-QTL mapping
    if 'cross_modality' in config['modalities']:
        if any([x != 'qtl' for x in config['analyses']]) or '{modality}.trans_qtl.txt.gz' in config['analyses']['qtl']['files']:
            raise ValueError('The "cross_modality" modality should only be used with analysis "qtl" for cross-modality cis-QTL mapping.')

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
            config[path] = str(Path(config[path]).expanduser())

    config['samples'] = pd.read_csv(config['samples_file'], sep='\t', header=None, dtype=str)[0].tolist()

    if 'intermediate_dir' not in config:
        config['intermediate_dir'] = 'intermediate'

    if 'qtl_min_allele_count' not in config:
        config['qtl_min_allele_count'] = 10
    config['qtl_min_allele_count'] = int(config['qtl_min_allele_count'])

    if 'qtl_maf_threshold' in config and config['qtl_maf_threshold'] is not None:
        config['qtl_maf_threshold'] = float(config['qtl_maf_threshold'])

    if 'twas_models_per_job' not in config:
        config['twas_models_per_job'] = 200
    config['twas_models_per_job'] = int(config['twas_models_per_job'])
    if config['twas_models_per_job'] < 1:
        raise ValueError('twas_models_per_job must be at least 1')

    if 'twas_threads' not in config:
        config['twas_threads'] = 4
    config['twas_threads'] = int(config['twas_threads'])
    if config['twas_threads'] < 1:
        raise ValueError('twas_threads must be at least 1')

validate_config(config)
process_config(config)

pheno_dir = Path(config['phenotype_dir'])
interm_dir = Path(config['intermediate_dir'])
output_dir = Path('output')
geno_prefix = config['geno_prefix']

samples = config['samples']
geno_samples = pd.read_csv(geno_prefix + '.fam', sep=r'\s+', header=None, dtype=str)[1].tolist()
missing_samples = [s for s in samples if s not in geno_samples]
if len(missing_samples) > 0:
    raise Exception(f'Samples missing from genotypes: {missing_samples}')

modalities = config['modalities']
for modality, params in modalities.items():
    chroms = pd.read_csv(pheno_dir / f'{modality}.bed.gz', sep='\t', usecols=[0], dtype=str)
    params['chroms'] = chroms.iloc[:, 0].unique().tolist()

analyses = config['analyses']
outputs = []
for analysis, params in analyses.items():
    if analysis == 'qtl':
        outputs.append(output_dir / analysis / 'genotype_qc.tsv')
    for f in params['files']:
        outputs.append(expand(output_dir / analysis / f, modality=modalities.keys()))

qtl_min_allele_count = config['qtl_min_allele_count']
qtl_maf_threshold = config.get('qtl_maf_threshold')
