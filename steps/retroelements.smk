rule run_telescope:
    """Run telescope to get retroelement read counts"""
    input:
        bam = project_dir / 'bam' / '{sample_id}.Aligned.sortedByCoord.out.bam',
        retro_anno = retro_anno,
    output:
        project_dir / 'retroelements' / '{sample_id}-telescope_report.tsv',
    params:
        retro_dir = project_dir / 'retroelements',
    resources:
        walltime = 8,
    shell:
        """
        mkdir -p {params.retro_dir}
        telescope assign \
            {input.bam} \
            {input.retro_anno} \
            --outdir {params.retro_dir} \
            --exp_tag {wildcards.sample_id}
        """

def load_retros(retro_anno: Path) -> pd.DataFrame:
    """Load retroelement annotations
    
    For retroelements with multiple entries (corresponding to multiple repeats),
    returns the median of their starts as its 'TSS'.
    """
    anno = read_gtf(retro_anno)
    anno['chromEnd'] = np.where(anno['strand'] == '+', anno['start'], anno['end'])
    anno = anno[['gene_id', 'seqname', 'chromEnd']]
    anno = anno.groupby(['gene_id', 'seqname']).median().astype(int).reset_index()
    # assert len(anno['gene_id'].unique()) == len(anno['gene_id'])
    anno['chromStart'] = anno['chromEnd'] - 1  # BED coordinates are 0-based
    anno = anno.rename(columns={'seqname': '#chrom', 'gene_id': 'name'})
    anno = anno.sort_values(['#chrom', 'chromStart'])
    return anno[['#chrom', 'chromStart', 'chromEnd', 'name']]

rule assemble_retroelement_bed:
    """Assemble telescope output files into retroelement expression BED file"""
    input:
        telescope_files = lambda w: expand(str(project_dir / 'retroelements' / '{sample_id}-telescope_report.tsv'), sample_id=samples),
        retro_anno = retro_anno,
    output:
        bed = project_dir / 'retroelements.bed',
    params:
        samples = samples,
        retro_dir = project_dir / 'retroelements',
    run:
        for i, sample in enumerate(params.samples):
            fname = params.retro_dir / f'{sample}-telescope_report.tsv'
            d = pd.read_csv(fname, sep='\t', index_col='transcript', skiprows=1)
            d = d.rename(columns={'final_count': sample})[[sample]]
            if i == 0:
                df = d
            else:
                df = df.join(d, how='outer')
        df = df.fillna(0).astype(int)
        df = df.loc[np.sum(df > 0, axis=1) > 0, :]
        anno = load_retros(input.retro_anno)
        df.index = df.index.rename('name')
        df = anno.merge(df.reset_index(), on='name', how='inner')
        df.to_csv(output.bed, sep='\t', index=False, float_format='%g')
