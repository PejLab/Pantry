localrules:
    expression_isoforms_pheno_groups,

rule kallisto_index:
    """Generate the index for kallisto."""
    input:
        ref_cdna,
    output:
        ref_dir / 'kallisto.idx',
    params:
        ref_dir = ref_dir,
    shell:
        """
        mkdir -p {params.ref_dir}
        kallisto index --index {output} {input}
        """

rule kallisto:
    """Quantify expression from FASTQ files using kallisto."""
    input:
        index = ref_dir / 'kallisto.idx',
        # bam = interm_dir / 'bam' / '{sample_id}.Aligned.toTranscriptome.out.bam',
        fastq = fastqs_kallisto,
    output:
        h5 = interm_dir / 'kallisto' / '{sample_id}' / 'abundance.h5',
        tsv = interm_dir / 'kallisto' / '{sample_id}' / 'abundance.tsv',
        json = interm_dir / 'kallisto' / '{sample_id}' / 'run_info.json',
    params:
        kallisto_dir = interm_dir / 'kallisto',
        out_dir = str(interm_dir / 'kallisto' / '{sample_id}'),
        single_end_flag = '' if paired_end else '--single',
        # TODO add strandedness parameter
    resources:
        cpus = threads,
    shell:
        """
        mkdir -p {params.kallisto_dir}
        kallisto quant \
            --index {input.index} \
            --output-dir {params.out_dir} \
            {params.single_end_flag} \
            --threads {resources.cpus} \
            {input.fastq}
        """

rule assemble_expression_bed:
    """Assemble kallisto est_counts and tpm outputs into expression BED files
    
    Currently only tpm is used for downstream analysis.
    """
    input:
        kallisto = expand(str(interm_dir / 'expression' / '{sample_id}' / 'abundance.tsv'), sample_id=samples),
        samples = samples_file,
        ref_anno = ref_anno,
    output:
        iso = interm_dir / 'unnorm' / 'expression.isoforms.bed',
        gene = interm_dir / 'unnorm' / 'expression.genes.bed',
    params:
        unnorm_dir = interm_dir / 'unnorm',
        project_dir = project_dir,
        expr_dir = interm_dir / 'expression',
    shell:
        """
        mkdir -p {params.unnorm_dir}
        python3 {params.project_dir}/src/assemble_bed.py \
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
        bed = interm_dir / 'unnorm' / 'expression.{level}.bed',
        samples = samples_file,
    output:
        output_dir / 'expression.{level}.bed.gz',
    params:
        project_dir = project_dir,
        bed = str(output_dir / 'expression.{level}.bed'),
    shell:
        """
        python3 {params.project_dir}/src/normalize_phenotypes.py \
            --input {input.bed} \
            --samples {input.samples} \
            --output {params.bed}
        bgzip {params.bed}
        """

rule expression_isoforms_pheno_groups:
    """Group phenotypes by gene for tensorQTL"""
    input:
        output_dir / 'expression.isoforms.bed.gz',
    output:
        output_dir / 'expression.isoforms.phenotype_groups.txt',
    shell:
        """
        zcat < {input} \
            | tail -n +2 \
            | cut -f4 \
            | awk '{{ g=$1; sub(/:.*$/, "", g); print $1 "\t" g }}' \
            > {output}
        """
