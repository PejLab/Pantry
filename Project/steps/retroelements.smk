rule run_telescope:
    """Run telescope to get retroelement read counts"""
    input:
        bam = interm_dir / 'bam' / '{sample_id}.Aligned.sortedByCoord.out.bam',
        bai = interm_dir / 'bam' / '{sample_id}.Aligned.sortedByCoord.out.bam.bai',
        retro_anno = retro_anno,
    output:
        interm_dir / 'retroelements' / '{sample_id}-telescope_report.tsv',
    params:
        retro_dir = interm_dir / 'retroelements',
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

rule assemble_retroelement_bed:
    """Assemble telescope output files into retroelement expression BED file"""
    input:
        telescope_files = lambda w: expand(str(interm_dir / 'retroelements' / '{sample_id}-telescope_report.tsv'), sample_id=samples),
        samples = samples_file,
        retro_anno = retro_anno,
    output:
        bed = interm_dir / 'unnorm' / 'retroelements.bed',
    params:
        unnorm_dir = interm_dir / 'unnorm',
        project_dir = project_dir,
        retro_dir = interm_dir / 'retroelements',
    shell:
        """
        mkdir -p {params.unnorm_dir}
        python3 {params.project_dir}/src/assemble_bed.py \
            --type retroelements \
            --input-dir {params.retro_dir} \
            --samples {input.samples} \
            --ref_anno {input.retro_anno} \
            --output {output.bed}
        """

rule normalize_retroelements:
    """Quantile-normalize values for QTL mapping"""
    input:
        bed = interm_dir / 'unnorm' / 'retroelements.bed',
        samples = samples_file,
    output:
        output_dir / 'retroelements.bed.gz',
    params:
        project_dir = project_dir,
        bed = output_dir / 'retroelements.bed',
    shell:
        """
        python3 {params.project_dir}/src/normalize_phenotypes.py \
            --input {input.bed} \
            --samples {input.samples} \
            --output {params.bed}
        bgzip {params.bed}
        """