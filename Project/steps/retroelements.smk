rule run_telescope:
    """Run telescope to get retroelement read counts"""
    input:
        bam = lambda w: bam_map[w.sample_id],
        bai = lambda w: f'{bam_map[w.sample_id]}.bai',
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
        telescope_files = expand(str(interm_dir / 'retroelements' / '{sample_id}-telescope_report.tsv'), sample_id=samples),
        samples = samples_file,
        retro_anno = retro_anno,
    output:
        bed = interm_dir / 'unnorm' / 'retroelements.bed',
    params:
        unnorm_dir = interm_dir / 'unnorm',
        retro_dir = interm_dir / 'retroelements',
    shell:
        """
        mkdir -p {params.unnorm_dir}
        python3 scripts/assemble_bed.py \
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
        bed = output_dir / 'retroelements.bed',
    shell:
        """
        python3 scripts/normalize_phenotypes.py \
            --input {input.bed} \
            --samples {input.samples} \
            --output {params.bed}
        bgzip {params.bed}
        """
