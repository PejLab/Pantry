def groups_arg(wildcards, input):
    """Pass the phenotype groups file as a Python arg if applicable"""
    if phenotypes[wildcards.pheno]['grouped']:
        return f'--groups {input.groups}'
    else:
        return ''

def groups_input(wildcards):
    """Include the phenotype groups file as an input if applicable"""
    if phenotypes[wildcards.pheno]['grouped']:
        return project_dir / f'{wildcards.pheno}.phenotype_groups.txt'
    else:
        return []

rule tensorqtl_perm:
    """Map cis-QTLs, determining significance using permutations.
    Outputs the top association per phenotype.
    """
    input:
        geno = multiext('data/genotype/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup', '.bed', '.bim', '.fam'),
        bed = project_dir / '{pheno}.bed.gz',
        bedi = project_dir / '{pheno}.bed.gz.tbi',
        covar = 'data/genotype/GEUVADIS.445_samples.covariates.txt',
        groups = groups_input,
    output:
        project_dir / 'qtl' / '{pheno}.cis_qtl.txt.gz',
    params:
        geno_prefix = 'data/genotype/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup',
        qtl_dir = project_dir / 'qtl',
        groups_arg = groups_arg,
    resources:
        walltime = 12,
        partition = '--partition=gpu',
    shell:
        """
        module load cuda
        mkdir -p {params.qtl_dir}
        python3 TURNAP/src/run_tensorqtl.py \
            {params.geno_prefix} \
            {input.bed} \
            {output} \
            --covariates {input.covar} \
            {params.groups_arg} \
            --mode cis
        """

rule tensorqtl_independent:
    """Use stepwise regression to identify multiple conditionally independent cis-QTLs per phenotype."""
    input:
        geno = multiext('data/genotype/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup', '.bed', '.bim', '.fam'),
        bed = project_dir / '{pheno}.bed.gz',
        bedi = project_dir / '{pheno}.bed.gz.tbi',
        covar = 'data/genotype/GEUVADIS.445_samples.covariates.txt',
        groups = groups_input,
        cis = project_dir / 'qtl' / '{pheno}.cis_qtl.txt.gz',
    output:
        project_dir / 'qtl' / '{pheno}.cis_independent_qtl.txt.gz',
    params:
        geno_prefix = 'data/genotype/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup',
        groups_arg = groups_arg,
    resources:
        walltime = 20,
        partition = '--partition=gpu',
    shell:
        """
        module load cuda
        python3 TURNAP/src/run_tensorqtl.py \
            {params.geno_prefix} \
            {input.bed} \
            {output} \
            --covariates {input.covar} \
            --cis_output {input.cis} \
            {params.groups_arg} \
            --mode cis_independent
        """
