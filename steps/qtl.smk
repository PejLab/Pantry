rule tensorqtl_perm:
    """Map cis-QTLs, determining significance using permutations.
    Outputs the top association per phenotype.
    """
    input:
        geno = multiext('data/genotype/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup', '.bed', '.bim', '.fam')
        bed = project_dir / '{pheno}.bed.gz',
        bedi = project_dir / '{pheno}.bed.gz.tbi',
        covar = 'data/genotype/GEUVADIS.445_samples.covariates.txt',
    output:
        project_dir / 'qtl' / '{pheno}.cis_qtl.txt.gz',
    params:
        geno_prefix = 'data/genotype/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup',
    resources:
        walltime = 12,
        partition = '--partition=gpu',
    shell:
        """
        module load cuda
        python3 src/run_tensorqtl.py \
            {params.geno_prefix} \
            {input.bed} \
            {output} \
            --covariates {input.covar} \
            --mode cis
        """

rule tensorqtl_independent:
    """Use stepwise regression to identify multiple conditionally independent cis-QTLs per phenotype."""
    input:
        geno = multiext('data/genotype/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup', '.bed', '.bim', '.fam'),
        bed = project_dir / '{pheno}.bed.gz',
        bedi = project_dir / '{pheno}.bed.gz.tbi',
        covar = 'data/genotype/GEUVADIS.445_samples.covariates.txt',
        cis = project_dir / 'qtl' / '{pheno}.cis_qtl.txt.gz',
    output:
        project_dir / 'qtl' / '{pheno}.cis_independent_qtl.txt.gz',
    params:
        geno_prefix = 'data/genotype/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup',
    resources:
        walltime = 20,
        partition = '--partition=gpu',
    shell:
        """
        module load cuda
        python3 src/run_tensorqtl.py \
            {params.geno_prefix} \
            {input.bed} \
            {output} \
            --covariates {input.covar} \
            --cis_output {input.cis} \
            --mode cis_independent
        """
