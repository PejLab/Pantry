rule exons_introns_from_GTF:
    """Extract exons and introns from GTF"""
    input:
        ref_anno = ref_anno,
    output:
        exon_gtf = ref_dir / 'constit_exons.gtf',
        intron_gtf = ref_dir / 'introns.gtf',
    shell:
        """
        sh TURNAP/src/exons_introns_from_gtf.sh \
            {input.ref_anno} \
            {output.exon_gtf} \
            {output.intron_gtf}
        """

rule run_featureCounts:
    """Run featureCounts from subread to get exon or intron read counts"""
    input:
        bam = project_dir / 'bam' / '{sample_id}.Aligned.sortedByCoord.out.bam',
        gtf = ref_dir / '{feature_type}.gtf',
    output:
        project_dir / 'stability' / '{sample_id}.{feature_type}.counts.txt',
    params:
        stab_dir = project_dir / 'stability',
        paired_flag = '-p' if paired_end else '',
        feature_id = lambda w: {'constit_exons': 'exon', 'introns': 'intron'}[w.feature_type],
        # TODO add strandedness parameter
    resources:
        cpus = threads,
        walltime = 8,
    shell:
        """
        mkdir -p {params.stab_dir}
        featureCounts \
            {input.bam} \
            {params.paired_flag} \
            -a {input.gtf} \
            -t {params.feature_id} \
            -T {resources.cpus} \
            -o {output}
        """

rule assemble_stability_bed:
    """Assemble exon to intron read ratios into mRNA stability BED file"""
    input:
        exon = lambda w: expand(str(project_dir / 'stability' / '{sample_id}.constit_exons.counts.txt'), sample_id=samples),
        intron = lambda w: expand(str(project_dir / 'stability' / '{sample_id}.introns.counts.txt'), sample_id=samples),
        ref_anno = ref_anno,
    output:
        bed = project_dir / 'unnorm' / 'stability.bed',
    params:
        samples = samples,
        stab_dir = project_dir / 'stability',
    run:
        for i, sample in enumerate(params.samples):
            fname_ex = params.stab_dir / f'{sample}.constit_exons.counts.txt'
            d_ex = pd.read_csv(fname_ex, sep='\t', index_col='Geneid', skiprows=1)
            fname_in = params.stab_dir / f'{sample}.introns.counts.txt'
            d_in = pd.read_csv(fname_in, sep='\t', index_col='Geneid', skiprows=1)
            if i == 0:
                exon = pd.DataFrame(index=d_ex.index)
                intron = pd.DataFrame(index=d_in.index)
            else:
                assert d_ex.index.equals(exon.index)
                assert d_in.index.equals(intron.index)
            exon[sample] = d_ex.iloc[:, 5]
            intron[sample] = d_in.iloc[:, 5]
            exon.loc[exon[sample] < 10, sample] = np.nan
            intron.loc[intron[sample] < 10, sample] = np.nan
        genes = exon.index[np.isin(exon.index, intron.index)]
        # genes = set(exon.index).intersection(intron.index)
        assert exon.loc[genes, :].index.equals(intron.loc[genes, :].index)
        assert exon.columns.equals(intron.columns)
        df = exon.loc[genes, :] / intron.loc[genes, :]
        df = df[df.isnull().mean(axis=1) <= 0.5]
        anno = load_tss(input.ref_anno)
        df.index = df.index.rename('gene_id')
        df = anno.merge(df.reset_index(), on='gene_id', how='inner')
        df = df.rename(columns={'gene_id': 'name'})
        df.to_csv(output.bed, sep='\t', index=False, float_format='%g')

rule normalize_stability:
    """Quantile-normalize values for QTL mapping"""
    input:
        project_dir / 'unnorm' / 'stability.bed',
    output:
        project_dir / 'stability.bed.gz',
    params:
        bed = project_dir / 'stability.bed',
    shell:
        """
        python3 TURNAP/src/normalize_phenotypes.py \
            --input {input} \
            --output {params.bed}
        bgzip {params.bed}
        """
