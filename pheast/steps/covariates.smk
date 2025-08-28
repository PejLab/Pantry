localrules:
    plink_covariates,

rule prune_for_covar:
    """Prune genotypes to compute covariate PCs.
    --indep-pairwise parameters are based on GTEx methods.
    """
    input:
        multiext(geno_prefix, '.bed', '.bim', '.fam'),
    output:
        multiext(str(interm_dir / 'covar' / 'geno_pruned'), '.bed', '.bim', '.fam'),
    params:
        geno_prefix = geno_prefix,
        pruned_dir = interm_dir / 'covar',
        pruned_prefix = interm_dir / 'covar' / 'geno_pruned',
    shell:
        # --geno 0.05 filters variants with >5% missing values (the rest will be imputed).
        # Default is 0 so that samples with many more genotyped variants than others don't
        # result in other samples having mostly missing values after pruning.
        """
        mkdir -p {params.pruned_dir}
        plink2 \
            --bfile {params.geno_prefix} \
            --geno 0.00 \
            --maf 0.05 \
            --indep-pairwise 200 100 0.1 \
            --out {params.pruned_prefix}
        plink2 \
            --bfile {params.geno_prefix} \
            --extract {params.pruned_prefix}.prune.in \
            --make-bed \
            --out {params.pruned_prefix}
        """

rule covariates:
    """Compute genotype and expression PCs and combine."""
    input:
        geno = multiext(str(interm_dir / 'covar' / 'geno_pruned'), '.bed', '.bim', '.fam'),
        bed = pheno_dir / '{modality}.bed.gz',
    output:
        interm_dir / 'covar' / '{modality}.covar.tsv',
    params:
        pruned_prefix = interm_dir / 'covar' / 'geno_pruned',
        n_geno_pcs = 5,
        n_pheno_pcs = 20,
    resources:
        mem_mb = 32000,
    shell:
        """
        Rscript scripts/covariates.R \
            {params.pruned_prefix} \
            {input.bed} \
            {params.n_geno_pcs} \
            {params.n_pheno_pcs} \
            {output}
        """

rule plink_covariates:
    """Convert covariates to PLINK format."""
    input:
        covar = interm_dir / 'covar' / '{modality}.covar.tsv',
        fam = f'{geno_prefix}.fam',
    output:
        plink = interm_dir / 'covar' / '{modality}.covar.plink.tsv',
    run:
        covar = pd.read_csv(input.covar, sep='\t', index_col=0, dtype=str)
        covar.index.name = None
        covar = covar.T
        covar.index.name = 'IID'
        covar = covar.reset_index()
        ## Get FIDs from genotypes:
        fam = pd.read_csv(input.fam, sep=r'\s+', header=None, dtype=str)
        # In some versions, dtype doesn't apply to index, so set index later:
        fam = fam.set_index(1).to_dict()[0]
        ## Insert FIDs as first column of covar, joining by IID:
        covar.insert(0, 'FID', covar['IID'].map(fam))
        covar.to_csv(output.plink, sep='\t', index=False, header=False)
