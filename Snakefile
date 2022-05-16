import yaml
from pathlib import Path
from gtfparse import read_gtf

config = yaml.safe_load(open('TURNAP/config.yml'))
fastq_map = Path(config['fastq_map'])
fastq_dir = Path(config['fastq_dir'])
read_length = config['read_length']
paired_end = config['paired_end']
project_dir = Path(config['project_dir'])
ref_genome = Path(config['ref_genome'])
ref_anno = Path(config['ref_anno'])
retro_anno = Path(config['retro_anno'])
ref_dir = project_dir / 'reference'
threads = config['threads']
phenotypes = config['phenotypes']

# Get samples, preserving order for outputs in case it's meaningful
with open(fastq_map, 'r') as f:
    tmp = [line.split('\t')[-1] for line in f.read().splitlines()]
    samples = []
    for s in tmp:
        if s not in samples:
            samples.append(s)

outputs = []
for pheno, params in phenotypes.items():
    for f in params['files']:
        outputs.append(project_dir / f)

include: 'steps/align.smk'
for phenotype in phenotypes:
    include: f'steps/{phenotype}.smk'

rule all:
    input:
        outputs
        # expand(ref_dir / 'rsem_index_ATS_APA' / '{type}.transcripts.fa',
        #        type=['grp_1.upstream', 'grp_2.upstream', 'grp_1.downstream', 'grp_2.downstream']),
        # expand(
        #     project_dir / 'ATS_APA' / '{sample_id}..{type}.isoforms.results.gz',
        #     sample_id=samples,
        #     type=['grp_1.upstream', 'grp_2.upstream', 'grp_1.downstream', 'grp_2.downstream']
        # )

def load_tss(ref_anno: Path) -> pd.DataFrame:
    """Load TSS annotations from GTF file
    
    Returns TSS as the first four columns of the BED format, meaning the
    coordinates are 0-based and chromEnd is just chromStart + 1.
    """
    anno = read_gtf(ref_anno)
    anno = anno.loc[anno['feature'] == 'gene', :]
    anno['chromEnd'] = np.where(anno['strand'] == '+', anno['start'], anno['end'])
    anno['chromStart'] = anno['chromEnd'] - 1  # BED coordinates are 0-based
    anno['#chrom'] = anno['seqname']
    anno = anno.sort_values(['#chrom', 'chromStart'])
    return anno[['#chrom', 'chromStart', 'chromEnd', 'gene_id']]
