rule tensorqtl_perm:
    """Map cis-eQTLs, determining significance using permutations.
    Outputs the top association per gene.
    """
    input:
        geno = multiext('data/genotype/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup', '.bed', '.bim', '.fam')
        bed = '{tissue}/{tissue}.expr.iqn.filtered.bed.gz',
        bedi = '{tissue}/{tissue}.expr.iqn.filtered.bed.gz.tbi',
        covar = 'data/genotype/GEUVADIS.445_samples.covariates.txt',
    output:
        '{tissue}/{tissue}.cis_qtl.txt.gz'
    params:
        geno_prefix = '{tissue}/geno',
    resources:
        walltime = 12,
        # partition = '--partition=gpu',
    shell:
        # module load cuda
        """
        python3 src/run_tensorqtl.py \
            {params.geno_prefix} \
            {input.bed} \
            {output} \
            --covariates {input.covar} \
            --mode cis
        """


rule tensorqtl_independent:
    """Use stepwise regression to identify multiple conditionally independent cis-eQTLs per gene."""
    input:
        geno = multiext('data/genotype/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup', '.bed', '.bim', '.fam')
        bed = '{tissue}/{tissue}.expr.iqn.filtered.bed.gz',
        bedi = '{tissue}/{tissue}.expr.iqn.filtered.bed.gz.tbi',
        covar = 'data/genotype/GEUVADIS.445_samples.covariates.txt',
        cis = '{tissue}/{tissue}.cis_qtl.txt.gz'
    output:
        '{tissue}/{tissue}.cis_independent_qtl.txt.gz'
    params:
        geno_prefix = '{tissue}/geno',
    resources:
        walltime = 20,
        # partition = '--partition=gpu',
    shell:
        # module load cuda
        """
        python3 src/run_tensorqtl.py \
            {params.geno_prefix} \
            {input.bed} \
            {output} \
            --covariates {input.covar} \
            --cis_output {input.cis} \
            --mode cis_independent
        """
