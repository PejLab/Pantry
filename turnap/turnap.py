"""Generate molecular phenotypes from RNA-Seq"""

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
        'retro_anno',
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

def main():
    """Parse command line arguments and run the given subcommand"""
    pass


if __name__ == '__main__':
    main()
