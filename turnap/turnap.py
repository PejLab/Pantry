"""Generate molecular phenotypes from RNA-Seq"""

def validate_config(config: dict):
    """Validate the configuration"""
    if len(config.keys()) == 0:
        raise Exception('No config file provided.')
    fields = [
        'project_dir',
        'threads',
        'fastq_map',
        'fastq_dir',
        'paired_end',
        'read_length',
        'ref_genome',
        'ref_anno',
        'retro_anno',
        'samples_file',
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
