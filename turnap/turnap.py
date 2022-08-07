"""Generate molecular phenotypes from RNA-Seq"""

def validate_config(config: dict):
    """Validate the configuration"""
    if len(config.keys()) == 0:
        raise Exception('No config file provided.')
    fields = [
        'fastq_map',
        'fastq_dir',
        'samples_file',
        'read_length',
        'paired_end',
        'project_dir',
        'ref_genome',
        'ref_anno',
        'retro_anno',
        'code_dir',
        'threads',
        'phenotypes',
        'genome_size',
        'chroms',
        'geno_prefix',
        'covar_file',
    ]
    for field in fields:
        if field not in config.keys():
            raise Exception(f'{field} not in config file.')

def main():
    """Parse command line arguments and run the given subcommand"""
    pass


if __name__ == '__main__':
    main()
