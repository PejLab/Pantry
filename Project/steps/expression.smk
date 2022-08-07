rule rsem_index:
    """Generate the index for RSEM."""
    input:
        ref_genome = ref_genome,
        ref_anno = ref_anno,
    output:
        ref_dir / 'rsem_index' / 'rsem_ref.transcripts.fa'
    params:
        index_dir = ref_dir / 'rsem_index',
        ref_prefix = ref_dir / 'rsem_index' / 'rsem_ref',
    resources:
        cpus = threads,
    shell:
        """
        mkdir -p {params.index_dir}
        rsem-prepare-reference \
            {input.ref_genome} \
            {params.ref_prefix} \
            --gtf {input.ref_anno} \
            --num-threads {resources.cpus}
        """

rule rsem:
    """Quantify expression from a BAM file."""
    input:
        ref = ref_dir / 'rsem_index' / 'rsem_ref.transcripts.fa',
        bam = interm_dir / 'bam' / '{sample_id}.Aligned.toTranscriptome.out.bam',
    output:
        genes = interm_dir / 'expression' / '{sample_id}.genes.results.gz',
        isoforms = interm_dir / 'expression' / '{sample_id}.isoforms.results.gz',
    params:
        ref_prefix = ref_dir / 'rsem_index' / 'rsem_ref',
        expr_dir = interm_dir / 'expression',
        out_prefix = str(interm_dir / 'expression' / '{sample_id}'),
        paired_end_flag = '--paired-end' if paired_end else '',
        # TODO add strandedness parameter
    resources:
        cpus = threads,
        walltime = 4,
    shell:
        """
        mkdir -p {params.expr_dir}
        rsem-calculate-expression \
            {params.paired_end_flag} \
            --num-threads {resources.cpus} \
            --quiet \
            --estimate-rspd \
            --no-bam-output \
            --alignments {input.bam} \
            {params.ref_prefix} \
            {params.out_prefix}
        gzip {params.out_prefix}.genes.results
        gzip {params.out_prefix}.isoforms.results
        rm -r {params.out_prefix}.stat
        """

rule assemble_expression_bed:
    """Assemble RSEM log2 and tpm outputs into expression BED files"""
    input:
        rsem = lambda w: expand(str(interm_dir / 'expression' / '{sample_id}.genes.results.gz'), sample_id=samples),
        samples = samples_file,
        ref_anno = ref_anno,
    output:
        log2 = interm_dir / 'unnorm' / 'expression.log2.bed',
        tpm = interm_dir / 'unnorm' / 'expression.tpm.bed',
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
            --output {output.log2} \
            --output2 {output.tpm}
        """


rule normalize_expression:
    """Quantile-normalize values for QTL mapping"""
    input:
        bed = interm_dir / 'unnorm' / 'expression.tpm.bed',
        samples = samples_file,
    output:
        output_dir / 'expression.bed.gz',
    params:
        project_dir = project_dir,
        bed = output_dir / 'expression.bed',
    shell:
        """
        python3 {params.project_dir}/src/normalize_phenotypes.py \
            --input {input.bed} \
            --samples {input.samples} \
            --output {params.bed}
        bgzip {params.bed}
        """
