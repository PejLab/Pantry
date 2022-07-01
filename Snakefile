import yaml
from pathlib import Path
from gtfparse import read_gtf

config = yaml.safe_load(open('config.yml'))
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

# TODO: Get chromosome list automatically from genotypes or BED file.
# Currently only used for heritability and qtl.
# CHROMS = [f'chr{i}' for i in range(1, 23)] # + ['chrX']
# GENO = 'data/genotype/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup'
# COVAR = 'data/genotype/GEUVADIS.445_samples.covariates.txt'
CHROMS = [str(x) for x in range(1, 21)] + ['X']
GENO = 'data_Brain/geno'
COVAR = 'data_Brain/covar.txt'

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
    # Used only for QTL mapping for now:
    params['grouped'] = any(['phenotype_groups' in f for f in params['files']])

# These steps are short and will not be submitted as cluster jobs:
localrules:
    index_bed,

include: 'steps/align.smk'
for phenotype in phenotypes:
    include: f'steps/{phenotype}.smk'
include: 'steps/heritability.smk'
include: 'steps/qtl.smk'

rule all:
    input:
        outputs,
        # expand(project_dir / 'heritability' / '{pheno}_hsq.tsv', pheno=phenotypes.keys()),
        # expand(project_dir / 'qtl' / '{pheno}.cis_independent_qtl.txt.gz', pheno=phenotypes.keys()),

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
    anno = anno[['#chrom', 'chromStart', 'chromEnd', 'gene_id']]
    # Rename columns for tensorQTL:
    anno.columns = ['#chr', 'start', 'end', 'phenotype_id']
    return anno

rule index_bed:
    input:
        project_dir / '{pheno}.bed.gz',
    output:
        project_dir / '{pheno}.bed.gz.tbi'
    shell:
        'tabix {input}'
