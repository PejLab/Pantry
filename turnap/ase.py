    # parser_extract.add_argument('--vcf', type=Path, metavar='PATH', help='Path of the VCF file with sample genotypes')

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
        (out_dir / 'phaser').mkdir(exist_ok=True)
        (out_dir / 'phaser_gene_ae').mkdir(exist_ok=True)


def extract(args):
    """Extract phenotypes for one sample"""
    message(f'TURNAP: Extracting phenotypes for sample {args.id}')
    out_dir = args.out_dir
    ref_dir = args.ref_dir if args.ref_dir is not None else out_dir / 'reference'
    bam_dir = out_dir / 'bam'
    bam = bam_dir / f'{args.id}.Aligned.sortedByCoord.out.bam'
    bam_trans = bam_dir / f'{args.id}.Aligned.toTranscriptome.out.bam'
    phaser_dir = out_dir / 'phaser'
    gene_ae_dir = out_dir / 'phaser_gene_ae'

    if not (phaser_dir / f'{args.id}.haplotypic_counts.txt').exists() or args.rerun_all:
        message('Running phASER')
        assert args.vcf is not None, 'VCF file must be provided to run phASER.'
        vcf = args.vcf
        assert vcf.exists()
        if not args.dry_run:
            run_phaser(bam, vcf, args.id, phaser_dir, args.threads)

    output = gene_ae_dir / f'{args.id}.gene_ae.txt'
    if not output.exists() or args.rerun_all:
        message('Running phASER Gene AE')
        counts = phaser_dir / f'{args.id}.haplotypic_counts.txt'
        if not args.dry_run:
            run_phaser_gene_ae(counts, ref_dir / 'genes.bed', output)


def combine(args):
    """Combine phenotypes for multiple samples"""
    message('TURNAP: Assembling phenotype tables')
    out_dir = args.out_dir
    ref_dir = args.ref_dir if args.ref_dir is not None else out_dir / 'reference'
    ref_anno = args.ref_anno
    expr_dir = out_dir / 'expression'
    gene_ae_dir = out_dir / 'phaser_gene_ae'
    splice_dir = out_dir / 'splicing'
    samples = [x for x in expr_dir.iterdir() if x.match('*.genes.results')]
    samples = [str(x.name).replace('.genes.results', '') for x in samples]
    assemble_ase(samples, gene_ae_dir, ref_dir / 'genes.bed', out_dir)


def run_phaser(bam: Path, vcf: Path, sample_id: str, phaser_dir: Path, threads: int = 1):
    """Run phASER to get ASE counts"""
    # subprocess.run(
    #     ['python', 'setup.py', 'build_ext', '--inplace'],
    #     cwd=str(Path(__file__).parent / 'phaser' / 'phaser')
    # )
    # phaser_py = str(Path(__file__).parent / 'phaser' / 'phaser' / 'phaser.py')
    phaser_py = '/home/dmunro/tools/phaser/phaser/phaser.py'
    cmd = [
        'python2', phaser_py,
        '--bam', str(bam),
        '--vcf', str(vcf),
        '--sample', sample_id,
        '--baseq', '10',
        '--mapq', '255',
        '--isize', '1e6',
        '--paired_end', '0',
        '--o', str(phaser_dir / sample_id),
        '--include_indels', '0',
        '--gw_phase_vcf', '1',
        # '--temp_dir', '$TMPDIR',
        '--threads', str(threads),
    ]
    subprocess.run(cmd)
    cmd = [
        'sed', '-i', 's/.Aligned.sortedByCoord.out//',
        str(phaser_dir / f'{sample_id}.haplotypic_counts.txt')
    ]
    subprocess.run(cmd)


def run_phaser_gene_ae(counts: Path, gene_bed: Path, output: Path):
    """Run phASER Gene AE to get gene-level haplotype-specific counts"""
    phaser_gene_ae_py = '/home/dmunro/tools/phaser/phaser_gene_ae/phaser_gene_ae.py'
    cmd = [
        'python', phaser_gene_ae_py,
        '--haplotypic_counts', str(counts),
        '--features', str(gene_bed),
        '--o', str(output),
    ]
    subprocess.run(cmd)


def assemble_ase(samples: list, gene_ae_dir: Path, gene_bed: Path, out_dir: Path):
    popdir = out_dir / 'phaser_pop'
    popdir.mkdir(parents=True, exist_ok=True)
    phaser_expr_matrix_py = '/home/dmunro/tools/phaser/phaser_pop/phaser_expr_matrix.py'
    cmd = [
        'python', phaser_expr_matrix_py,
        '--gene_ae_dir', str(gene_ae_dir),
        '--features', str(gene_bed),
        '--o', str(popdir / 'phaser_pop'),
        # rm {output.bed}.tbi {output.bedgw}.tbi
    ]
    subprocess.run(cmd)
    # cmd = [
    #     'zcat', str(popdir / 'phaser_pop.bed.gz'),
    #     '| sed "s/contig\tstart\tstop\tname'
    #     '| sed "s///g" >',
    #     str(out_dir / 'ASE.bed')
    #     ]
    # subprocess.run()


