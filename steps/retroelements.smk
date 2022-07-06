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

rule assemble_retroelement_bed:
    """Assemble telescope output files into retroelement expression BED file"""
    input:
        telescope_files = lambda w: expand(str(project_dir / 'retroelements' / '{sample_id}-telescope_report.tsv'), sample_id=samples),
        samples = samples_file,
        retro_anno = retro_anno,
    output:
        bed = project_dir / 'unnorm' / 'retroelements.bed',
    params:
        # samples = samples,
        retro_dir = project_dir / 'retroelements',
    shell:
        """
        python3 TURNAP/src/assemble_bed.py \
            --type retroelements \
            --input-dir {params.retro_dir} \
            --samples {input.samples} \
            --ref_anno {input.retro_anno} \
            --output {output.bed}
        """

rule normalize_retroelements:
    """Quantile-normalize values for QTL mapping"""
    input:
        bed = project_dir / 'unnorm' / 'retroelements.bed',
        samples = samples_file,
    output:
        project_dir / 'retroelements.bed.gz',
    params:
        bed = project_dir / 'retroelements.bed',
    shell:
        """
        python3 TURNAP/src/normalize_phenotypes.py \
            --input {input.bed} \
            --samples {input.samples} \
            --output {params.bed}
        bgzip {params.bed}
        """
