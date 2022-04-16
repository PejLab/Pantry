"""Generate molecular phenotypes from RNA-Seq"""

import argparse
import os
from pathlib import Path
import subprocess
from typing import List
import numpy as np
import pandas as pd
from gtfparse import read_gtf
from . import leafcutter_cluster_regtools_py3


def main():
    """Parse command line arguments and run the given subcommand"""
    parser = argparse.ArgumentParser(description='Generate molecular phenotypes from RNA-Seq')
    subparsers = parser.add_subparsers(title='Commands', required=True)

    parser_prep = subparsers.add_parser('prepare')
    parser_prep.add_argument('--out-dir', type=Path, required=True, help='Name of directory in which to write output')
    parser_prep.add_argument('--ref-seq', type=Path, required=True, metavar='PATH', help='Path of the genome reference fasta file')
    parser_prep.add_argument('--ref-anno', type=Path, required=True, metavar='PATH', help='Path of the genome reference annotation GTF file')
    parser_prep.add_argument('--ref-dir', type=Path, metavar='PATH', help='To compute reference files, such as for STAR and RSEM, only once and reuse for multiple datasets, specify the directory path, which will be created if it does not yet exist. By default, new reference files are written in {out-dir}/reference/.')
    parser_prep.add_argument('--read-length', type=int, required=True, metavar='N', help='RNA-Seq read length (used for preparing alignment index)')
    parser_prep.add_argument('--rerun-all', action='store_true', help='Unless this flag is provided, each step is only run if its output files do not yet exist')
    parser_prep.add_argument('-t', '--threads', type=int, default=1, metavar='N', help='Number of threads to use when possible (default 1)')
    parser_prep.add_argument('-n', '--dry-run', action='store_true', help='Include this flag to do a trial run, printing steps that would run but not running them')
    parser_prep.set_defaults(func=prepare)

    parser_extract = subparsers.add_parser('extract')
    parser_extract.add_argument('id', help='Sample ID that will be used in output files')
    paired_group = parser_extract.add_mutually_exclusive_group(required=True)
    paired_group.add_argument('--single-end', action='store_false', dest='paired_end', help='Reads are single-end. Provide one or more FASTQ files whose reads will be combined into one BAM file.')
    paired_group.add_argument('--paired-end', action='store_true', help='Reads are paired-end. Provide one or more pairs of FASTQ files whose reads will be combined into one BAM file.')
    # fastq_group = parser_extract.add_mutually_exclusive_group(required=True)
    # # fastq_group.add_argument('--fastq', nargs='+', type=Path, metavar='PATH', help='FASTQ file path(s) for the sample, separated by spaces. For paired-end reads, provide paths in comma-separated pairs, e.g. "--fastq fileA_1.fastq.gz,fileA_2.fastq.gz fileB_1.fastq.gz,fileB_2.fastq.gz"')
    # fastq_group.add_argument('--fastq', nargs='+', type=Path, metavar='PATH', help='For single-end, FASTQ file path(s) for the sample, separated by spaces"')
    # fastq_paired_group = fastq_group.add_argument_group()
    # fastq_paired_group.add_argument('--fastq1', nargs='+', type=Path, metavar='PATH', help='For paired-end, first-in-each-pair FASTQ file path(s) for the sample, separated by spaces')
    # fastq_paired_group.add_argument('--fastq2', nargs='+', type=Path, metavar='PATH', help='For paired-end, second-in-each-pair FASTQ file path(s) for the sample, separated by spaces')
    # fastq_group.add_argument('--fastq-map', type=Path, metavar='PATH', help='As an alternative to the --fastq file list, provide the path of the file containing FASTQ file paths and the name of the sample they map to, separated by tabs, with no header. Only the files mapping to the current sample ID will be used. For paired-end reads, include a row for each pair and have two FASTQ columns plus the sample name.')
    parser_extract.add_argument('--fastq-map', type=Path, required=True, metavar='PATH', help='The path of a file containing FASTQ file paths and the name of the sample they map to, separated by tabs, with no header. Only the files mapping to the current sample ID will be used. For paired-end reads, include a row for each pair and have two FASTQ columns plus the sample name.')
    parser_extract.add_argument('--fastq-dir', type=Path, default=Path('.'), metavar='PATH', help='For convenience, FASTQ file paths given by --fastq or --fastq_map can be relative to this directory. Defaults to the current directory.')
    parser_extract.add_argument('--out-dir', type=Path, required=True, metavar='PATH', help='Name of directory in which to write output')
    parser_extract.add_argument('--ref-dir', type=Path, metavar='PATH', help='Path of reference directory if different from the default {out-dir}/reference/')
    parser_extract.add_argument('--rerun-all', action='store_true', help='Unless this flag is provided, each step is only run if its output files do not yet exist')
    parser_extract.add_argument('-t', '--threads', type=int, default=1, metavar='N', help='Number of threads to use when possible (default 1)')
    parser_extract.add_argument('-n', '--dry-run', action='store_true', help='Include this flag to do a trial run, printing steps that would run but not running them')
    parser_extract.set_defaults(func=extract)

    parser_combine = subparsers.add_parser('combine')
    parser_combine.add_argument('--out-dir', type=Path, metavar='PATH', help='Name of main directory containing all outputs')
    parser_combine.add_argument('--ref-dir', type=Path, metavar='PATH', help='Path of reference directory if different from the default {out-dir}/reference/')
    parser_combine.add_argument('--ref-anno', type=Path, metavar='PATH', help='Path of the genome reference annotation GTF file')
    parser_combine.set_defaults(func=combine)

    args = parser.parse_args()
    args.func(args)


def prepare(args):
    """Prepare reference files and output directories"""
    message('TURNAP: Preparing reference files and output directories')
    ref_seq = args.ref_seq
    assert ref_seq.exists()
    ref_anno = args.ref_anno
    assert ref_anno.exists()
    out_dir = args.out_dir
    ref_dir = args.ref_dir if args.ref_dir is not None else out_dir / 'reference'
    if not args.dry_run:
        out_dir.mkdir(parents=True, exist_ok=True)
        ref_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / 'bam').mkdir(exist_ok=True)
        (out_dir / 'expression').mkdir(exist_ok=True)
        (out_dir / 'splicing').mkdir(exist_ok=True)
        (out_dir / 'stability').mkdir(exist_ok=True)
    if not (ref_dir / 'star_index' / 'SAindex').exists() or args.rerun_all:
        message('Preparing STAR index')
        if not args.dry_run:
            prepare_star_index(ref_seq, ref_anno, ref_dir, args.read_length, args.threads)
    exon_gtf = ref_dir / 'constit_exons.gtf'
    intron_gtf = ref_dir / 'introns.gtf'
    if not (exon_gtf.exists() and intron_gtf.exists()) or args.rerun_all:
        message('Extracting exons and introns from GTF')
        script = Path(__file__).parent.parent / 'src' / 'exons_introns_from_gtf.sh'
        cmd = ['sh', script, str(ref_anno), str(exon_gtf), str(intron_gtf)]
        if not args.dry_run:
            subprocess.run(cmd)
    collapsed_gtf = ref_dir / 'collapsed.gtf'
    gene_bed = ref_dir / 'genes.bed'
    if not gene_bed.exists() or args.rerun_all:
        message('Extracting genes from GTF')
        if not args.dry_run:
            script = Path(__file__).parent.parent / 'src' / 'collapse_annotation.py'
            subprocess.run(['python3', script, ref_anno, collapsed_gtf], check=True)
            cmd = f'gtf2bed < {collapsed_gtf} | awk \'{{{{ if ($8 == "gene") {{{{ print }}}}}}}}\' > {gene_bed}'
            subprocess.run(cmd, shell=True)
    if not (ref_dir / 'rsem_index' / 'rsem_ref.chrlist').exists() or args.rerun_all:
        message('Preparing RSEM reference')
        if not args.dry_run:
            prepare_rsem_reference(ref_seq, ref_anno, ref_dir, args.threads)


def extract(args):
    """Extract phenotypes for one sample"""
    message(f'TURNAP: Extracting phenotypes for sample {args.id}')
    out_dir = args.out_dir
    ref_dir = args.ref_dir if args.ref_dir is not None else out_dir / 'reference'
    bam_dir = out_dir / 'bam'
    bam = bam_dir / f'{args.id}.Aligned.sortedByCoord.out.bam'
    bam_trans = bam_dir / f'{args.id}.Aligned.toTranscriptome.out.bam'
    expr_dir = out_dir / 'expression'
    splice_dir = out_dir / 'splicing'
    stab_dir = out_dir / 'stability'

    if not bam.exists() or args.rerun_all:
        fastqs = get_fastq_paths(args)
        message('Running STAR')
        prefix = str(bam_dir / args.id) + '.'
        if not args.dry_run:
            run_star(fastqs, ref_dir, prefix, args.threads)
            subprocess.run(['samtools', 'index', str(bam)])

    if not (expr_dir / f'{args.id}.genes.results').exists() or args.rerun_all:
        message('Running RSEM')
        if not args.dry_run:
            run_rsem(bam_trans, args.paired_end, ref_dir, str(expr_dir / args.id), args.threads)

    junc = splice_dir / f'{args.id}.junc'
    if not junc.exists() or args.rerun_all:
        message('Counting splice junctions')
        if not args.dry_run:
            run_junctions(bam, junc)

    exons = stab_dir / f'{args.id}.constit_exons.counts.txt'
    if not exons.exists() or os.stat(exons).st_size == 0 or args.rerun_all:
        message('Counting constitutive exon reads for RNA stability')
        exon_gtf = ref_dir / 'constit_exons.gtf'
        if not args.dry_run:
            # run_htseq(bam, 'exon', exon_gtf, exons, args.threads)
            run_featureCounts(bam, args.paired_end, 'exon', exon_gtf, exons, args.threads)

    introns = stab_dir / f'{args.id}.introns.counts.txt'
    if not introns.exists() or os.stat(introns).st_size == 0 or args.rerun_all:
        message('Counting intron reads for RNA stability')
        intron_gtf = ref_dir / 'introns.gtf'
        if not args.dry_run:
            # run_htseq(bam, 'intron', intron_gtf, introns, args.threads)
            run_featureCounts(bam, args.paired_end, 'intron', intron_gtf, introns, args.threads)


def combine(args):
    """Combine phenotypes for multiple samples"""
    message('TURNAP: Assembling phenotype tables')
    out_dir = args.out_dir
    ref_dir = args.ref_dir if args.ref_dir is not None else out_dir / 'reference'
    ref_anno = args.ref_anno
    expr_dir = out_dir / 'expression'
    splice_dir = out_dir / 'splicing'
    stab_dir = out_dir / 'stability'
    samples = [x for x in expr_dir.iterdir() if x.match('*.genes.results')]
    samples = [str(x.name).replace('.genes.results', '') for x in samples]
    # assemble_expr_beds(samples, expr_dir, ref_anno, out_dir)
    # cluster_junctions(samples, splice_dir)
    # assemble_splicing_bed(splice_dir, ref_anno, out_dir)
    # TODO assemble exon/intron ratios
    # prepare_rembrandts_metadata(samples, stab_dir)
    # run_rembrandts(htseq_dir, stab_dir)


def message(text: str):
    """Print a message with some formatting to distinguish it from other software messages"""
    print(f'=== {text} {"=" * max(3, 75 - len(text))}', flush=True)


def get_fastq_paths(args) -> List[List[Path]]:
    """Get FASTQ paths from command line arguments and confirm that the files exist
    
    Returns a list containing one (single-end) or two (paired-end) lists of paths.
    """
    if not args.paired_end:
        if args.fastq_map is not None:
            fastq_map = pd.read_csv(args.fastq_map, sep='\t', dtype=str, names=['file', 'sample'])
            fastqs = [list(fastq_map.loc[fastq_map['sample'] == args.id, 'file'])]
        else:
            fastqs = [args.fastq]
    else:
        if args.fastq_map is not None:
            fastq_map = pd.read_csv(args.fastq_map, sep='\t', dtype=str, names=['file1', 'file2', 'sample'])
            fastqs1 = list(fastq_map.loc[fastq_map['sample'] == args.id, 'file1'])
            fastqs2 = list(fastq_map.loc[fastq_map['sample'] == args.id, 'file2'])
            fastqs = [fastqs1, fastqs2]
        else:
            fastqs = [args.fastq1, args.fastq2]
    paths = []
    for group in fastqs:
        paths.append([])
        for fastq in group:
            path = args.fastq_dir / fastq
            assert path.exists(), f'FASTQ file not found: {path}'
            paths[-1].append(path)
    return paths


def prepare_star_index(ref_seq: Path, ref_anno: Path, ref_dir: Path, read_length: int, threads: int = 1):
    """Generate STAR reference index"""
    index_dir = ref_dir / 'star_index'
    index_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        'STAR',
        '--runMode', 'genomeGenerate',
        '--outTmpDir', str(index_dir / 'tmp'),
        '--outFileNamePrefix', str(index_dir / 'STAR_'),  # prefix of Log.out
        '--genomeDir', str(index_dir),
        '--genomeFastaFiles', str(ref_seq),
        '--sjdbGTFfile', str(ref_anno),
        '--sjdbOverhang', str(read_length - 1),
        '--runThreadN', str(threads),
        # These two recommended by STAR author for low memory:
        # '--genomeSAsparseD', '3',
        # '--genomeSAindexNbases', '12',
    ]
    subprocess.run(cmd, check=True)


def prepare_rsem_reference(ref_seq: Path, ref_anno: Path, ref_dir: Path, threads: int = 1):
    """Generate RSEM reference index"""
    index_dir = ref_dir / 'rsem_index'
    index_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        'rsem-prepare-reference',
        str(ref_seq),
        str(index_dir / 'rsem_ref'),
        '--gtf', str(ref_anno),
        '--num-threads', str(threads),
    ]
    subprocess.run(cmd, check=True)


def run_star(fastqs: List[List[Path]], ref_dir: Path, prefix: str, threads: int = 1):
    """Run STAR alignment"""
    # Pass an unpacked tuple so paired-end isn't seen as one arg containing space
    files = tuple([','.join([str(x) for x in group]) for group in fastqs])
    cmd = [
        'STAR',
        '--runMode', 'alignReads',
        '--genomeDir', str(ref_dir / 'star_index'),
        '--readFilesIn', *files,
        '--readFilesCommand', 'zcat',
        '--twopassMode', 'Basic',
        '--quantMode', 'TranscriptomeSAM',
        '--outSAMstrandField', 'intronMotif',
        '--outSAMtype', 'BAM', 'SortedByCoordinate',
        '--outFileNamePrefix', prefix,
        '--runThreadN', str(threads),
    ]
    subprocess.run(cmd, check=True)


def run_rsem(bam: Path, paired: bool, ref_dir: Path, prefix: str, threads: int = 1):
    """Run RSEM to quantify total expression"""
    cmd = [
        'rsem-calculate-expression',
        '--paired-end' if paired else '',
        '--num-threads', str(threads),
        '--quiet',
        '--estimate-rspd',
        '--no-bam-output',
        '--alignments', str(bam),
        str(ref_dir / 'rsem_index' / 'rsem_ref'),
        prefix,
    ]
    subprocess.run(cmd)


def run_junctions(bam: Path, junc: Path):
    """Run regtools junctions for LeafCutter"""
    cmd = [
        'regtools', 'junctions', 'extract',
        '-a' '8',       # Minimum anchor length
        '-m' '50',      # Minimum intron length
        '-M' '500000',  # Maximum intron length
        '-s' '0',       # Strand specificity of RNA library prep
        str(bam),
        '-o', str(junc),
    ]
    subprocess.run(cmd)


def run_htseq(bam: Path, feature: str, gtf: Path, output: Path, threads: int = 1):
    """Run HTSeq-count to get exon or intron counts"""
    assert feature in {'exon', 'intron'}
    cmd = [
        'htseq-count',
        '--mode', 'intersection-strict' if feature == 'exon' else 'union',
        '--format', 'bam',
        '--type', feature,
        '--stranded', 'reverse',
        '--nprocesses', str(threads),  # Actually only works for multiple BAM files
        str(bam),
        str(gtf),
        '>', str(output),
    ]
    # print(' '.join(cmd))
    # subprocess.run(cmd, stdout=open(output, 'w'))
    subprocess.run(' '.join(cmd), shell=True)
    if os.stat(output).st_size == 0:
        os.remove(output)


def run_featureCounts(bam: Path, paired: bool, feature: str, gtf: Path, output: Path, threads: int = 1):
    """Run featureCounts from subread to get exon or intron counts"""
    assert feature in {'exon', 'intron'}
    cmd = [
        'featureCounts',
        str(bam),
        '-p' if paired else '',
        '-a', str(gtf),
        '-t', feature,
        '-T', str(threads),
        '-o', str(output),
    ]
    subprocess.run(' '.join(cmd), shell=True)
    if os.stat(output).st_size == 0:
        os.remove(output)


def load_tss(ref_anno: Path) -> pd.DataFrame:
    """Load TSS annotations
    
    Returns TSS as the first four columns of the BED format, meaning the coordinates are 0-based and chromEnd is just chromStart + 1.
    """
    anno = read_gtf(ref_anno)
    anno = anno.loc[anno['feature'] == 'gene', :]
    anno['chromEnd'] = np.where(anno['strand'] == '+', anno['start'], anno['end'])
    anno['chromStart'] = anno['chromEnd'] - 1  # BED coordinates are 0-based
    # anno = anno.rename(columns={'seqname': '#chrom', 'gene_id': 'name'})
    anno['#chrom'] = anno['seqname']
    anno = anno.sort_values(['#chrom', 'chromStart'])
    return anno[['#chrom', 'chromStart', 'chromEnd', 'gene_id']]


def assemble_expr_beds(samples: list, expr_dir: Path, ref_anno: Path, out_dir: Path):
    """Assemble RSEM log2 and tpm outputs into expression BED files"""
    for i, sample in enumerate(samples):
        fname = expr_dir / f'{sample}.genes.results'
        d = pd.read_csv(fname, sep='\t', index_col='gene_id')
        if i == 0:
            log2 = pd.DataFrame(index=d.index)
            tpm = pd.DataFrame(index=d.index)
        else:
            assert d.index.equals(log2.index)
        log2[sample] = np.log2(d['expected_count'] + 1)
        tpm[sample] = d['TPM']
    anno = load_tss(ref_anno)
    log2 = anno.merge(log2.reset_index(), on='gene_id', how='inner')
    tpm = anno.merge(tpm.reset_index(), on='gene_id', how='inner')
    log2 = log2.rename(columns={'gene_id': 'name'})
    tpm = tpm.rename(columns={'gene_id': 'name'})
    log2.to_csv(out_dir / 'expression.log2.bed', sep='\t', index=False, float_format='%g')
    tpm.to_csv(out_dir / 'expression.tpm.bed', sep='\t', index=False, float_format='%g')


def load_exons(ref_anno: Path) -> pd.DataFrame:
    """Load exon annotations
    
    Returns exons with start and end oriented on the gene's strand."""
    anno = read_gtf(ref_anno)
    anno = anno.loc[anno['feature'] == 'exon', :]
    anno['chrom'] = anno['seqname']
    anno['exonStart'] = np.where(anno['strand'] == '+', anno['start'], anno['end'])
    anno['exonEnd'] = np.where(anno['strand'] == '+', anno['end'], anno['start'])
    return anno[['gene_id', 'chrom', 'exonStart', 'exonEnd']]


def map_introns_to_genes(introns: list, exons: pd.DataFrame) -> pd.DataFrame:
    exons['exonStart'] = exons['exonStart'].astype(str)
    exons['exonEnd'] = exons['exonEnd'].astype(str)
    df = pd.DataFrame({'intron': introns})
    df[['chrom', 'chr_start', 'chr_end', 'clu', 'cluster', 'strand']] = df['intron'].str.split(r':|_', expand=True)
    df['start'] = np.where(df['strand'] == '+', df['chr_start'], df['chr_end'])
    df['end'] = np.where(df['strand'] == '+', df['chr_end'], df['chr_start'])
    start_matches = df.merge(
        exons[['chrom', 'exonEnd', 'gene_id']].rename(columns={'exonEnd': 'start'}),
        on=['chrom', 'start'],
        how='inner'
    )
    end_matches = df.merge(
        exons[['chrom', 'exonStart', 'gene_id']].rename(columns={'exonStart': 'end'}),
        on=['chrom', 'end'],
        how='inner'
    )
    df = pd.concat([start_matches, end_matches])
    clust_genes = df.groupby('cluster', group_keys=False).agg({'gene_id': pd.Series.unique})
    clust_genes = clust_genes.reset_index().explode('gene_id')
    return clust_genes


def cluster_junctions(samples: list, splice_dir: Path):
    class Options(object):
        pass
    options = Options()
    options.outprefix = 'clust'
    options.verbose = True
    options.rundir = str(splice_dir)
    options.maxintronlen = 100000
    options.minclureads = 30
    options.mincluratio = 0.001
    options.cluster = None
    options.checkchrom = False
    juncfiles = [str(splice_dir / f'{sample}.junc') for sample in samples]
    leafcutter_cluster_regtools_py3.main(options, juncfiles)


def assemble_splicing_bed(splice_dir: Path, ref_anno: Path, out_dir: Path):
    df = pd.read_csv(splice_dir / 'clust_perind_numers.counts.gz', sep=' ')
    samples = list(df.columns)
    df.index = df.index.rename('intron')
    df = df.reset_index()
    # df = df.rename(columns={'chrom': 'name'})
    df['cluster'] = df['intron'].str.extract(r'clu_(\d+)_', expand=False)
    for col in samples:
    #     df[col] = pd.eval(df[col])  # convert fraction string to float
        df[col] = df.groupby('cluster', group_keys=False).apply(lambda g: g[col] / g[col].sum())
    exons = load_exons(ref_anno)
    genes = map_introns_to_genes(df['intron'], exons)
    df = df.merge(genes, on='cluster', how='left')
    anno = load_tss(ref_anno)
    df = anno.merge(df, on='gene_id', how='inner')
    df['name'] = df['gene_id'] + ':' + df['intron']
    df = df[['#chrom', 'chromStart', 'chromEnd', 'name'] + samples]
    df.to_csv(out_dir / 'splicing.bed', sep='\t', index=False, float_format='%g')


def prepare_rembrandts_metadata(samples: list, stab_dir: Path):
    df = pd.DataFrame({
        'Label': [x for x in samples for _ in ('exonic', 'intronic')],
        'File': [f'{x}.{y}.counts.txt' for x in samples for y in ('consExons', 'Introns')],
        'ReadType': ['exonic', 'intronic'] * len(samples),
        'Batch': [1] * (len(samples) * 2),
    })
    output = stab_dir / 'metadata.txt'
    df.to_csv(output, sep='\t', index=False)


def run_rembrandts(htseq_dir: Path, stab_dir: Path):
    working_dir = Path(__file__).parent.parent / 'tools/REMBRANDTS'
    cmd = [
        'sh', 'REMBRANDTS.sh',
        'stability',
        str(stab_dir.resolve() / 'metadata.txt'),
        str(htseq_dir.resolve()),
        '0.99',     # stringency for genes to be included in analysis. Most examples use 0.99.
        'linear'  # biasMode (currently only linear is accepted)
    ]
    subprocess.run(cmd, cwd=working_dir)


if __name__ == '__main__':
    main()
