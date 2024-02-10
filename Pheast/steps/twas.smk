localrules:
    twas_pos_file,

TWAS_N_BATCHES = 64
MAX_BATCH_SIZE = 300
twas_geno_prefix = interm_dir / 'twas' / 'geno' if 'twas_snps' in config else geno_prefix
twas_snps = Path(config['twas_snps']) if 'twas_snps' in config else ''

rule twas_geno:
    """Subset genotypes to SNPs used in TWAS LD reference panel."""
    input:
        geno = multiext(geno_prefix, '.bed', '.bim', '.fam'),
        snps = twas_snps,
    output:
        multiext(str(interm_dir / 'twas' / 'geno'), '.bed', '.bim', '.fam'),
    params:
        in_geno_prefix = geno_prefix,
        out_geno_prefix = interm_dir / 'twas' / 'geno',
    shell:
        """
        plink \
            --bfile {params.in_geno_prefix} \
            --extract {input.snps} \
            --make-bed \
            --out {params.out_geno_prefix}
        """

rule twas_compute_weights_batch:
    """Use FUSION to compute TWAS weights from expression and genotypes.
    
    Outputs also include the actual weights
    ('{interm_dir}/twas/weights_{modality}/{gene}.wgt.RDat'), but the genes for
    which we end up with outputs are not known until the rule is run.
    """
    input:
        geno = multiext(str(twas_geno_prefix), '.bed', '.bim', '.fam'),
        bed = pheno_dir / '{modality}.bed.gz',
        covar = interm_dir / 'covar' / '{modality}.covar.plink.tsv',
        gcta = 'scripts/fusion_twas/gcta_nr_robust',
        gemma = 'scripts/fusion_twas/gemma',
    output:
        # expand(interm_dir / 'twas' / 'weights_{{modality}}' / '{gene}.wgt.RDat', gene=gene_list),
        interm_dir / 'twas' / 'hsq_{modality}' / '{batch_start}_{batch_end}.hsq',
    params:
        twas_geno_prefix = twas_geno_prefix,
        twas_interm_dir = interm_dir / 'twas',
    resources:
        walltime = 8,
    shell:
        """
        sh scripts/fusion_compute_weights.sh \
            {params.twas_geno_prefix} \
            {input.bed} \
            {input.covar} \
            {wildcards.modality} \
            {wildcards.batch_start} \
            {wildcards.batch_end} \
            {params.twas_interm_dir}
        """

def twas_batch_hsq_input(wildcards):
    """Get start and end indices of BED file for TWAS batches."""
    bed = pheno_dir / f'{wildcards.modality}.bed.gz'
    n = int(subprocess.check_output(f'zcat < {bed} | wc -l', shell=True)) - 1
    batch_size = min(math.ceil(n / TWAS_N_BATCHES), MAX_BATCH_SIZE)
    # Given the necessary batch size, items might fit into fewer batches:
    n_batches = math.ceil(n / batch_size)
    for i in range(n_batches):
        start = i * batch_size + 1
        end = min((i + 1) * batch_size, n)
        yield interm_dir / 'twas' / f'hsq_{wildcards.modality}' / f'{start}_{end}.hsq'

rule twas_assemble_summary:
    """Summarize weights from all batches/genes.
    
    The hsq 'input' files aren't actually used, but signify the batch was run.
    The true input files are the per-gene weights, but the genes for which we
    end up with outputs are not known until the rule is run.
    """
    input:
        twas_batch_hsq_input,
    output:
        file_list = interm_dir / 'twas' / '{modality}.list',
        profile = interm_dir / 'twas' / '{modality}.profile',
        summary = interm_dir / 'twas' / '{modality}.profile.err',
    params:
        twas_interm_dir = interm_dir / 'twas',
    shell:
        """
        # Avoid using relative path in case intermediate dir is a symlink:
        scripts_dir="$(realpath scripts)"
        cd {params.twas_interm_dir}
        ls {wildcards.modality}/*.wgt.RDat > {wildcards.modality}.list
        Rscript $scripts_dir/fusion_twas/utils/FUSION.profile_wgt.R \
            {wildcards.modality}.list \
            > {wildcards.modality}.profile \
            2> {wildcards.modality}.profile.err
        """

rule twas_pos_file:
    input:
        file_list = interm_dir / 'twas' / '{modality}.list',
        bed = pheno_dir / '{modality}.bed.gz',
    output:
        pos_file = interm_dir / 'twas' / '{modality}.pos',
    params:
        n_samples = len(samples),
    shell:
        """
        echo 'WGT\tID\tCHR\tP0\tP1\tN' > {output.pos_file}
        cut -d'/' -f2 {input.file_list} \
            | sed 's/.wgt.RDat//' \
            | paste {input.file_list} - \
            | join -1 2 -2 4 - <(zcat {input.bed} | cut -f1-4 | sort -k4) \
            | awk '{{OFS="\t"; print $2, $1, $3, $4, $5, {params.n_samples}}}' \
            >> {output.pos_file}
        """

rule twas_compress_output:
    """Combine TWAS weights and summary files into a single archive.

    Uses the same format as these:
    http://gusevlab.org/projects/fusion/#single-tissue-gene-expression
    """
    input:
        file_list = interm_dir / 'twas' / '{modality}.list',
        profile = interm_dir / 'twas' / '{modality}.profile',
        summary = interm_dir / 'twas' / '{modality}.profile.err',
        pos_file = interm_dir / 'twas' / '{modality}.pos',
    output:
        output_dir / 'twas' / '{modality}.tar.bz2',
    params:
        twas_interm_dir = interm_dir / 'twas',
        twas_output_dir = output_dir / 'twas',
    shell:
        """
        mkdir -p {params.twas_output_dir}
        tar -cjf {output} \
            -C {params.twas_interm_dir} \
            {wildcards.modality}.list \
            {wildcards.modality}.profile \
            {wildcards.modality}.profile.err \
            {wildcards.modality}.pos \
            {wildcards.modality}/
        """
