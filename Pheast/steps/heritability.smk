rule calculate_heritability_chr_batch:
    """Get heritability for phenotypes in one batch using plink and gcta64.
    
    Runs one batch of phenotypes from one chromosome. Batches are
    assigned using the last digit of the TSS position (i.e. the `end` column)
    so that phenotypes for the same gene are in the same batch.
    """
    input:
        bed = pheno_dir / '{modality}.bed.gz',
        geno = multiext(geno_prefix, '.bed', '.bim', '.fam'),
        covar = interm_dir / 'covar' / '{modality}.covar.plink.tsv',
    output:
        hsq = interm_dir / 'heritability' / '{modality}.hsq.{chrom}.batch{batch}.tsv',
    params:
        geno_prefix = geno_prefix,
        grm_dir = str(interm_dir / 'heritability' / 'grm_{modality}_{chrom}'),
        tmp_dir = str(interm_dir / 'heritability' / 'tmp_{modality}_{chrom}'),
    resources:
        mem_mb = 8000,
        walltime = 24,
    shell:
        """
        mkdir -p {params.grm_dir}
        mkdir -p {params.tmp_dir}
        python3 scripts/heritability.py \
            --bed {input.bed} \
            --geno {params.geno_prefix} \
            --covar {input.covar} \
            --chrom {wildcards.chrom} \
            --batch {wildcards.batch} \
            --grm-dir {params.grm_dir} \
            --tmp-dir {params.tmp_dir} \
            --output {output.hsq}
        """

rule combine_heritability_chr:
    """Combine per-chromosome heritability stats into one file"""
    input:
        lambda w: expand(interm_dir / 'heritability' / '{modality}.hsq.{chrom}.batch{batch}.tsv', modality=w.modality, chrom=config['modalities'][w.modality]['chroms'], batch=range(10)),
    output:
        hsq = output_dir / 'heritability' / '{modality}.hsq.tsv',
    params:
        hsq_dir = output_dir / 'heritability',
    shell:
        """
        mkdir -p {params.hsq_dir}
        head -n 1 {input[0]} > {output.hsq}
        for FILE in {input}; do
            tail -n +2 $FILE >> {output.hsq}
        done
        """
