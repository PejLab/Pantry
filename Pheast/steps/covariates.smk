localrules:
    plink_covariates,

rule prune_for_covar:
    """Prune genotypes to compute covariate PCs.
    --indep-pairwise parameters are based on GTEx methods.
    """
    input:
        multiext(geno_prefix, '.bed', '.bim', '.fam'),
    output:
        interm_dir / 'covar' / 'geno_pruned.vcf.gz',
    params:
        geno_prefix = geno_prefix,
        pruned_dir = interm_dir / 'covar',
        pruned_prefix = interm_dir / 'covar' / 'geno_pruned',
    shell:
        # --geno 0.05 filters variants with >5% missing values (the rest will be imputed)
        """
        mkdir -p {params.pruned_dir}
        plink2 \
            --bfile {params.geno_prefix} \
            --geno 0.05 \
            --maf 0.05 \
            --indep-pairwise 200 100 0.1 \
            --out {params.pruned_prefix}
        plink2 \
            --bfile {params.geno_prefix} \
            --extract {params.pruned_prefix}.prune.in \
            --export vcf bgz id-paste=iid \
            --out {params.pruned_prefix}
        """

rule covariates:
    """Compute genotype and expression PCs and combine."""
    input:
        vcf = interm_dir / 'covar' / 'geno_pruned.vcf.gz',
        bed = pheno_dir / '{pheno}.bed.gz',
    output:
        interm_dir / 'covar' / '{pheno}.covar.tsv',
    params:
        n_geno_pcs = 5,
        n_pheno_pcs = 20,
    shell:
        """
        Rscript scripts/covariates.R \
            {input.vcf} \
            {input.bed} \
            {params.n_geno_pcs} \
            {params.n_pheno_pcs} \
            {output}
        """

rule plink_covariates:
    """Convert covariates to PLINK format."""
    input:
        covar = interm_dir / 'covar' / '{pheno}.covar.tsv',
    output:
        plink = interm_dir / 'covar' / '{pheno}.covar.plink.tsv',
    run:
        covar = pd.read_csv(input.covar, sep='\t', index_col=0)
        covar.index.name = None
        covar = covar.T
        covar.index.name = 'IID'
        covar = covar.reset_index()
        covar.insert(0, 'FID', 0)  # Family must match that in the genotype files.
        covar.to_csv(output.plink, sep='\t', index=False, header=False)
