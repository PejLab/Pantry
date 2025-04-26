localrules:
    isoforms_pheno_groups,

rule kallisto_index:
    """Generate the index for kallisto."""
    input:
        ref_cdna,
    output:
        ref_dir / 'kallisto.idx',
    params:
        ref_dir = ref_dir,
    resources:
        mem_mb = 32000,
    shell:
        """
        mkdir -p {params.ref_dir}
        kallisto index --index {output} {input}
        """

rule kallisto:
    """Quantify expression from FASTQ files using kallisto."""
    input:
        index = ref_dir / 'kallisto.idx',
        fastq = fastq_inputs,
    output:
        h5 = interm_dir / 'expression' / '{sample_id}' / 'abundance.h5',
        tsv = interm_dir / 'expression' / '{sample_id}' / 'abundance.tsv',
        json = interm_dir / 'expression' / '{sample_id}' / 'run_info.json',
    params:
        expr_dir = interm_dir / 'expression',
        out_dir = str(interm_dir / 'expression' / '{sample_id}'),
        single_end_flag = kallisto_single_end_flag,
        # TODO add strandedness parameter
    threads: 16
    resources:
        runtime = '16h',
    shell:
        """
        mkdir -p {params.expr_dir}
        kallisto quant \
            --index {input.index} \
            --output-dir {params.out_dir} \
            {params.single_end_flag} \
            --threads {threads} \
            {input.fastq}
        """

rule assemble_expression_bed:
    """Assemble kallisto outputs into expression BED files
    
    Currently only tpm is used for downstream analysis.
    """
    input:
        kallisto = expand(str(interm_dir / 'expression' / '{sample_id}' / 'abundance.tsv'), sample_id=samples),
        samples = samples_file,
        ref_anno = ref_anno,
    output:
        iso = output_dir / 'unnorm' / 'isoforms.bed',
        gene = output_dir / 'unnorm' / 'expression.bed',
    params:
        unnorm_dir = output_dir / 'unnorm',
        expr_dir = interm_dir / 'expression',
    resources:
        mem_mb = 32000,
    shell:
        """
        mkdir -p {params.unnorm_dir}
        python3 scripts/assemble_bed.py \
            --type expression \
            --input-dir {params.expr_dir} \
            --samples {input.samples} \
            --ref_anno {input.ref_anno} \
            --output {output.iso} \
            --output2 {output.gene}
        """

rule normalize_expression:
    """Quantile-normalize values for QTL mapping"""
    input:
        bed = output_dir / 'unnorm' / 'expression.bed',
        samples = samples_file,
    output:
        output_dir / 'expression.bed.gz',
    params:
        bed = str(output_dir / 'expression.bed'),
    shell:
        """
        python3 scripts/normalize_phenotypes.py \
            --input {input.bed} \
            --samples {input.samples} \
            --output {params.bed}
        bgzip {params.bed}
        """

rule normalize_isoforms:
    """Quantile-normalize values for QTL mapping"""
    input:
        bed = output_dir / 'unnorm' / 'isoforms.bed',
        samples = samples_file,
    output:
        output_dir / 'isoforms.bed.gz',
    params:
        bed = str(output_dir / 'isoforms.bed'),
    shell:
        """
        mkdir -p {output_dir}
        python3 scripts/normalize_phenotypes.py \
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
            | awk '{{ g=$1; sub(/__.*$/, "", g); print $1 "\t" g }}' \
            > {output}
        """
