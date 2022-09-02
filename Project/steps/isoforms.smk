# Uses RSEM outputs produced in expression.smk

localrules:
    isoforms_pheno_groups,

rule assemble_isoforms_bed:
    """Assemble RSEM isoform percentages into a BED file"""
    input:
        rsem = lambda w: expand(str(interm_dir / 'expression' / '{sample_id}.isoforms.results.gz'), sample_id=samples),
        samples = samples_file,
        ref_anno = ref_anno,
    output:
        interm_dir / 'unnorm' / 'isoforms.bed',
    params:
        unnorm_dir = interm_dir / 'unnorm',
        project_dir = project_dir,
        expr_dir = interm_dir / 'expression',
    shell:
        """
        mkdir -p {params.unnorm_dir}
        python3 {params.project_dir}/src/assemble_bed.py \
            --type isoforms \
            --input-dir {params.expr_dir} \
            --samples {input.samples} \
            --ref_anno {input.ref_anno} \
            --output {output}
        """

rule normalize_isoforms:
    """Quantile-normalize values for QTL mapping"""
    input:
        bed = interm_dir / 'unnorm' / 'isoforms.bed',
        samples = samples_file,
    output:
        output_dir / 'isoforms.bed.gz',
    params:
        project_dir = project_dir,
        bed = output_dir / 'isoforms.bed',
    shell:
        """
        python3 {params.project_dir}/src/normalize_phenotypes.py \
            --input {input.bed} \
            --samples {input.samples} \
            --output {params.bed}
        bgzip {params.bed}
        """

rule isoforms_pheno_groups:
    """Group phenotypes by gene for tensorQTL"""
    input:
        output_dir / 'isoforms.bed.gz',
    output:
        output_dir / 'isoforms.phenotype_groups.txt',
    shell:
        """
        zcat < {input} \
            | tail -n +2 \
            | cut -f4 \
            | awk '{{ g=$1; sub(/:.*$/, "", g); print $1 "\t" g }}' \
            > {output}
        """
