TWAS_N_BATCHES = 64

rule twas_compute_weights_batch:
    """Use FUSION to compute TWAS weights from expression and genotypes."""
    input:
        geno = multiext(geno_prefix, '.bed', '.bim', '.fam'),
        bed = pheno_dir / '{pheno}.bed.gz',
        # covar = interm_dir / 'covar' / '{pheno}.covar.tsv',
    output:
        # expand(interm_dir / 'twas' / 'weights_{{pheno}}' / '{gene}.wgt.RDat', gene=gene_list),
        interm_dir / 'twas' / 'hsq_{pheno}' / '{batch_start}_{batch_end}.hsq',
    params:
        geno_prefix = geno_prefix,
        twas_interm_dir = interm_dir / 'twas',
    shell:
        """
        sh scripts/fusion_compute_weights.sh \
            {params.geno_prefix} \
            {input.bed} \
            {wildcards.pheno} \
            {wildcards.batch_start} \
            {wildcards.batch_end} \
            {params.twas_interm_dir}
        """

def twas_batch_hsq_input(wildcards):
    """Get start and end indices of BED file for TWAS batches."""
    bed = pheno_dir / f'{wildcards.pheno}.bed.gz'
    n = int(subprocess.check_output(f'zcat {bed} | wc -l', shell=True))
    batch_size = math.ceil(n / TWAS_N_BATCHES)
    for i in range(TWAS_N_BATCHES):
        start = i * batch_size + 1
        end = min((i + 1) * batch_size, n)
        yield interm_dir / 'twas' / f'hsq_{wildcards.pheno}' / f'{start}_{end}.hsq'

rule twas_assemble_weights:
    """Combine TWAS weights from all batches/genes."""
    input:
        # expand(interm_dir / 'twas' / 'weights_{{pheno}}' / '{gene}.wgt.RDat', gene=gene_list),
        # expand(interm_dir / 'twas' / 'hsq_{{pheno}}' / '{batch}.hsq', batch=[f'{i}_{i+9}' for i in range(1, 100, 10)]),
        twas_batch_hsq_input,
    output:
        output_dir / 'twas' / '{pheno}.weights.tsv',
    params:
        weights_dir = str(interm_dir / 'twas' / 'weights_{pheno}'),
    shell:
        """
        ls {params.weights_dir}/*.wgt.RDat > {params.weights_dir}.list
        Rscript scripts/fusion_twas/utils/FUSION.profile_wgt.R \
            {params.weights_dir}.list \
            > {output}
        """
