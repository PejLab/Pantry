localrules:
    splicing_pheno_groups,

rule regtools_junctions:
    """Run regtools junctions for LeafCutter"""
    input:
        bam = project_dir / 'bam' / '{sample_id}.Aligned.sortedByCoord.out.bam',
        bai = project_dir / 'bam' / '{sample_id}.Aligned.sortedByCoord.out.bam.bai',
    output:
        project_dir / 'splicing' / '{sample_id}.junc'
    params:
        splice_dir = project_dir / 'splicing',
        min_anchor_len = 8,
        min_intron_len = 50,
        max_intron_len = 500000,
        strandedness = 0,
    shell:
        """
        mkdir -p {params.splice_dir}
        regtools junctions extract \
            {input.bam} \
            -a {params.min_anchor_len} \
            -m {params.min_intron_len} \
            -M {params.max_intron_len} \
            -s {params.strandedness} \
            -o {output}
        """

rule cluster_junctions:
    """Cluster splice junctions with common boundaries"""
    input:
        lambda w: expand(str(project_dir / 'splicing' / '{sample_id}.junc'), sample_id=samples),
    output:
        project_dir / 'splicing' / 'leafcutter_perind_numers.counts.gz',
    params:
        juncfile_list = project_dir / 'splicing' / 'juncfiles.txt',
        splice_dir = project_dir / 'splicing',
        max_intron_len = 100000,
        min_clust_reads = 30,
        min_clust_ratio = 0.001,
    shell:
        # TODO get proper path to script
        """
        printf '%s\\n' {input} > {params.juncfile_list}
        python3 TURNAP/src/leafcutter_cluster_regtools_py3.py \
            --juncfiles {params.juncfile_list} \
            --rundir {params.splice_dir} \
            --maxintronlen {params.max_intron_len} \
            --minclureads {params.min_clust_reads} \
            --mincluratio {params.min_clust_ratio} \
            --checkchrom False
        """

def load_exons(ref_anno: Path) -> pd.DataFrame:
    """Load exon annotations
    
    Returns exons with start and end oriented on the gene's strand.
    """
    anno = read_gtf(ref_anno)
    anno = anno.loc[anno['feature'] == 'exon', :]
    anno['chrom'] = anno['seqname']
    anno['exonStart'] = np.where(anno['strand'] == '+', anno['start'], anno['end'])
    anno['exonEnd'] = np.where(anno['strand'] == '+', anno['end'], anno['start'])
    return anno[['gene_id', 'chrom', 'exonStart', 'exonEnd']]

def map_introns_to_genes(introns: list, exons: pd.DataFrame) -> pd.DataFrame:
    """Map de novo splice junctions to genes based on known exon boundaries"""
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

rule assemble_splicing_bed:
    """Convert leafcutter output into splicing BED file"""
    input:
        counts = project_dir / 'splicing' / 'leafcutter_perind_numers.counts.gz',
        ref_anno = ref_anno,
    output:
        bed = project_dir / 'unnorm' / 'splicing.bed',
    run:
        df = pd.read_csv(input.counts, sep=' ')
        samples = list(df.columns)
        df.index = df.index.rename('intron')
        df = df.reset_index()
        # df = df.rename(columns={'chrom': 'name'})
        df['cluster'] = df['intron'].str.extract(r'clu_(\d+)_', expand=False)
        for col in samples:
        #     df[col] = pd.eval(df[col])  # convert fraction string to float
            df[col] = df.groupby('cluster', group_keys=False).apply(lambda g: g[col] / g[col].sum())
        exons = load_exons(input.ref_anno)
        genes = map_introns_to_genes(df['intron'], exons)
        df = df.merge(genes, on='cluster', how='left')
        anno = load_tss(input.ref_anno)
        df = anno.merge(df, on='gene_id', how='inner')
        df['name'] = df['gene_id'] + ':' + df['intron']
        df = df[['#chrom', 'chromStart', 'chromEnd', 'name'] + samples]
        df.to_csv(output.bed, sep='\t', index=False, float_format='%g')

rule normalize_splicing:
    """Quantile-normalize values for QTL mapping"""
    input:
        project_dir / 'unnorm' / 'splicing.bed',
    output:
        project_dir / 'splicing.bed.gz',
    params:
        bed = project_dir / 'splicing.bed',
    shell:
        """
        python3 TURNAP/src/normalize_phenotypes.py \
            --input {input} \
            --output {params.bed}
        bgzip {params.bed}
        """

rule splicing_pheno_groups:
    """Group phenotypes by gene for tensorQTL"""
    input:
        project_dir / 'splicing.bed.gz',
    output:
        project_dir / 'splicing.phenotype_groups.txt',
    run:
        df = pd.read_csv(input[0], sep='\t', usecols=['name'])
        df['group'] = df['name'].str.split(':', expand=False).str[0]
        df.to_csv(output[0], sep='\t', index=False, header=False)
