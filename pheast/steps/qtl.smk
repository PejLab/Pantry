def groups_arg(wildcards, input):
    """Pass the phenotype groups file as a Python arg if applicable"""
    if modalities[wildcards.modality]['grouped']:
        return f'--groups {input.groups}'
    else:
        return ''

def groups_input(wildcards):
    """Include the phenotype groups file as an input if applicable"""
    if modalities[wildcards.modality]['grouped']:
        return pheno_dir / f'{wildcards.modality}.phenotype_groups.txt'
    else:
        return []

def qtl_filter_args(wildcards):
    """Pass QTL genotype filter settings to tensorQTL wrapper scripts."""
    args = f'--min_allele_count {qtl_min_allele_count}'
    if qtl_maf_threshold is not None:
        args += f' --maf_threshold {qtl_maf_threshold}'
    return args

rule qtl_genotype_qc:
    """Report genotype variants below the QTL allele-count/MAF threshold."""
    input:
        geno = multiext(geno_prefix, '.bed', '.bim', '.fam'),
    output:
        output_dir / 'qtl' / 'genotype_qc.tsv',
    params:
        geno_prefix = geno_prefix,
        qtl_dir = output_dir / 'qtl',
        filter_args = qtl_filter_args,
    shell:
        """
        mkdir -p {params.qtl_dir}
        python3 scripts/genotype_qc.py \
            {params.geno_prefix} \
            {output} \
            {params.filter_args}
        """

rule tensorqtl_cis:
    """Map cis-QTLs, determining significance using permutations.
    Outputs the top association per phenotype.
    """
    input:
        geno = multiext(geno_prefix, '.bed', '.bim', '.fam'),
        bed = pheno_dir / '{modality}.bed.gz',
        bedi = pheno_dir / '{modality}.bed.gz.tbi',
        covar = interm_dir / 'covar' / '{modality}.covar.tsv',
        groups = groups_input,
    output:
        output_dir / 'qtl' / '{modality}.cis_qtl.txt.gz',
    params:
        geno_prefix = geno_prefix,
        qtl_dir = output_dir / 'qtl',
        groups_arg = groups_arg,
        filter_args = qtl_filter_args,
    resources:
        mem_mb = 32000,
        runtime = '12h',
    shell:
        ## Cluster environments may require cuda to be loaded, e.g.:
        # module load cuda
        """
        mkdir -p {params.qtl_dir}
        python3 scripts/run_tensorqtl.py \
            {params.geno_prefix} \
            {input.bed} \
            {output} \
            --covariates {input.covar} \
            {params.groups_arg} \
            {params.filter_args} \
            --mode cis
        """

rule tensorqtl_cis_independent:
    """Use stepwise regression to identify multiple conditionally independent cis-QTLs per phenotype."""
    input:
        geno = multiext(geno_prefix, '.bed', '.bim', '.fam'),
        bed = pheno_dir / '{modality}.bed.gz',
        bedi = pheno_dir / '{modality}.bed.gz.tbi',
        covar = interm_dir / 'covar' / '{modality}.covar.tsv',
        groups = groups_input,
        cis = output_dir / 'qtl' / '{modality}.cis_qtl.txt.gz',
    output:
        output_dir / 'qtl' / '{modality}.cis_independent_qtl.txt.gz',
    params:
        geno_prefix = geno_prefix,
        groups_arg = groups_arg,
        filter_args = qtl_filter_args,
    resources:
        mem_mb = 32000,
        runtime = '20h',
    shell:
        ## Cluster environments may require cuda to be loaded, e.g.:
        # module load cuda
        """
        python3 scripts/run_tensorqtl.py \
            {params.geno_prefix} \
            {input.bed} \
            {output} \
            --covariates {input.covar} \
            --cis_output {input.cis} \
            {params.groups_arg} \
            {params.filter_args} \
            --mode cis_independent
        """

rule tensorqtl_trans:
    """Map trans-QTL associations without permutation testing.
    Outputs all variant pairs with TSS distance > 5 Mb, configured MAF threshold, and p < 1e-5.
    """
    input:
        geno = multiext(geno_prefix, '.bed', '.bim', '.fam'),
        bed = pheno_dir / '{modality}.bed.gz',
        bedi = pheno_dir / '{modality}.bed.gz.tbi',
        covar = interm_dir / 'covar' / '{modality}.covar.tsv',
    output:
        output_dir / 'qtl' / '{modality}.trans_qtl.txt.gz',
    params:
        geno_prefix = geno_prefix,
        qtl_dir = output_dir / 'qtl',
        filter_args = qtl_filter_args,
    resources:
        mem_mb = 32000,
        runtime = '12h',
    shell:
        ## Cluster environments may require cuda to be loaded, e.g.:
        # module load cuda
        """
        mkdir -p {params.qtl_dir}
        python3 scripts/run_tensorqtl.py \
            {params.geno_prefix} \
            {input.bed} \
            {output} \
            --covariates {input.covar} \
            {params.filter_args} \
            --mode trans
        """
