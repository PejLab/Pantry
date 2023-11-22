# rule calculate_heritability_chr:
#     """Get heritability for phenotypes on one chromosome using plink and gcta64.
    
#     Splits into 10 batches run in parallel within the same job. Batches are
#     assigned using the last digit of the TSS position (i.e. the `end` column)
#     so that phenotypes for the same gene are in the same batch.
#     """
#     input:
#         bed = pheno_dir / '{pheno}.bed.gz',
#         geno = multiext(geno_prefix, '.bed', '.bim', '.fam'),
#         covar = interm_dir / 'covar' / '{pheno}.covar.plink.tsv',
#     output:
#         hsq = expand(interm_dir / 'heritability' / '{{pheno}}.hsq.{{chrom}}.batch{batch}.tsv', batch=range(10)),
#     params:
#         geno_prefix = geno_prefix,
#         grm_dir = str(interm_dir / 'heritability' / 'grm_{pheno}_{chrom}'),
#         tmp_dir = str(interm_dir / 'heritability' / 'tmp_{pheno}_{chrom}'),
#         out_template = str(interm_dir / 'heritability' / '{pheno}.hsq.{chrom}.batch{}.tsv'),
#     threads: 10,
#     resources:
#         mem_mb = 32000,
#         walltime = 24,
#     shell:
#         """
#         mkdir -p {params.grm_dir}
#         mkdir -p {params.tmp_dir}
#         parallel -j {threads} --halt now,fail=1 \
#             python3 scripts/heritability.py \
#             --bed {input.bed} \
#             --geno {params.geno_prefix} \
#             --covar {input.covar} \
#             --chrom {wildcards.chrom} \
#             --batch {{}} \
#             --grm-dir {params.grm_dir} \
#             --tmp-dir {params.tmp_dir} \
#             --output {params.out_template} \
#             ::: {{0..9}}
#         rm -r {params.grm_dir}
#         rm -r {params.tmp_dir}
#         """

rule calculate_heritability_chr_batch:
    """Get heritability for phenotypes in one batch using plink and gcta64.
    
    Runs one batch of phenotypes from one chromosome. Batches are
    assigned using the last digit of the TSS position (i.e. the `end` column)
    so that phenotypes for the same gene are in the same batch.
    """
    input:
        bed = pheno_dir / '{pheno}.bed.gz',
        geno = multiext(geno_prefix, '.bed', '.bim', '.fam'),
        covar = interm_dir / 'covar' / '{pheno}.covar.plink.tsv',
    output:
        hsq = interm_dir / 'heritability' / '{pheno}.hsq.{chrom}.batch{batch}.tsv',
    params:
        geno_prefix = geno_prefix,
        grm_dir = str(interm_dir / 'heritability' / 'grm_{pheno}_{chrom}'),
        tmp_dir = str(interm_dir / 'heritability' / 'tmp_{pheno}_{chrom}'),
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
        lambda w: expand(interm_dir / 'heritability' / '{pheno}.hsq.{chrom}.batch{batch}.tsv', pheno=w.pheno, chrom=config['phenotypes'][w.pheno]['chroms'], batch=range(10)),
    output:
        hsq = output_dir / 'heritability' / '{pheno}.hsq.tsv',
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
