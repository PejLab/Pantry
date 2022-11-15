rule calculate_heritability_chr:
    """Get heritability for phenotypes on one chromosome using plink and gcta64."""
    input:
        bed = pheno_dir / '{pheno}.bed.gz',
        geno = multiext(geno_prefix, '.bed', '.bim', '.fam'),
    output:
        interm_dir / 'heritability' / '{pheno}.hsq.{chrom}.tsv',
    params:
        geno_prefix = geno_prefix,
        grm_dir = lambda w: interm_dir / 'heritability' / f'grm_{w.pheno}_{w.chrom}',
        tmp_dir = lambda w: interm_dir / 'heritability' / f'tmp_{w.pheno}_{w.chrom}',
    resources:
        walltime = 12,
    shell:
        """
        mkdir -p {params.grm_dir}
        mkdir -p {params.tmp_dir}
        python3 scripts/heritability.py \
            --bed {input.bed} \
            --geno {params.geno_prefix} \
            --chrom {wildcards.chrom} \
            --grm-dir {params.grm_dir} \
            --tmp-dir {params.tmp_dir} \
            --output {output}
        rm -r {params.grm_dir}
        rm -r {params.tmp_dir}
        """

rule combine_heritability_chr:
    """Combine per-chromosome heritability stats into one file"""
    input:
        lambda w: expand(interm_dir / 'heritability' / '{pheno}.hsq.{chrom}.tsv', pheno=w.pheno, chrom=config['phenotypes'][w.pheno]['chroms']),
    output:
        output_dir / 'heritability' / '{pheno}.hsq.tsv',
    params:
        hsq_dir = output_dir / 'heritability',
    shell:
        """
        mkdir -p {params.hsq_dir}
        head -n 1 {input[0]} > {output}
        for FILE in {input}; do
            tail -n +2 $FILE >> {output}
        done
        """
