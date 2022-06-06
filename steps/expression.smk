import numpy as np
import pandas as pd

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
        bam = project_dir / 'bam' / '{sample_id}.Aligned.toTranscriptome.out.bam',
    output:
        genes = project_dir / 'expression' / '{sample_id}.genes.results.gz',
        isoforms = project_dir / 'expression' / '{sample_id}.isoforms.results.gz',
    params:
        ref_prefix = ref_dir / 'rsem_index' / 'rsem_ref',
        expr_dir = project_dir / 'expression',
        out_prefix = str(project_dir / 'expression' / '{sample_id}'),
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

rule assemble_expression:
    """Assemble RSEM log2 and tpm outputs into expression BED files"""
    input:
        rsem = lambda w: expand(str(project_dir / 'expression' / '{sample_id}.genes.results.gz'), sample_id=samples),
        ref_anno = ref_anno,
    output:
        log2 = project_dir / 'expression.log2.bed',
        tpm = project_dir / 'expression.tpm.bed',
    params:
        samples = samples,
        expr_dir = project_dir / 'expression',
    run:
        for i, sample in enumerate(params.samples):
            fname = params.expr_dir / f'{sample}.genes.results.gz'
            d = pd.read_csv(fname, sep='\t', index_col='gene_id')
            if i == 0:
                log2 = pd.DataFrame(index=d.index)
                tpm = pd.DataFrame(index=d.index)
            else:
                assert d.index.equals(log2.index)
            log2[sample] = np.log2(d['expected_count'] + 1)
            tpm[sample] = d['TPM']
        anno = load_tss(input.ref_anno)
        log2 = anno.merge(log2.reset_index(), on='gene_id', how='inner')
        tpm = anno.merge(tpm.reset_index(), on='gene_id', how='inner')
        log2 = log2.rename(columns={'gene_id': 'name'})
        tpm = tpm.rename(columns={'gene_id': 'name'})
        log2.to_csv(output.log2, sep='\t', index=False, float_format='%g')
        tpm.to_csv(output.tpm, sep='\t', index=False, float_format='%g')
