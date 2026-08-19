localrules:
    twas_pos_file,

twas_geno_prefix = interm_dir / 'twas' / 'geno' if 'twas_snps' in config else geno_prefix
twas_snps = Path(config['twas_snps']) if 'twas_snps' in config else ''

def twas_model_count(modality):
    with gzip.open(pheno_dir / f'{modality}.bed.gz', 'rt') as handle:
        return sum(1 for _ in handle) - 1

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
    ('{interm_dir}/twas/{modality}/{gene}.wgt.RDat'), but the genes for
    which we end up with outputs are not known until the rule is run.
    """
    input:
        geno = multiext(str(twas_geno_prefix), '.bed', '.bim', '.fam'),
        bed = pheno_dir / '{modality}.bed.gz',
        covar = interm_dir / 'covar' / '{modality}.covar.plink.tsv',
        worker = 'scripts/fit_twas_weights.py',
        fusion = 'scripts/fusion_twas/FUSION.compute_weights.R',
        gcta = 'scripts/fusion_twas/gcta_nr_robust',
    output:
        interm_dir / 'twas' / 'status_{modality}' / '{batch_start}_{batch_end}.tsv',
    params:
        twas_geno_prefix = twas_geno_prefix,
        twas_interm_dir = interm_dir / 'twas',
    resources:
        runtime = '8h',
    threads: config['twas_threads']
    shell:
        """
        python {input.worker} \
            --geno {params.twas_geno_prefix} \
            --bed {input.bed} \
            --covar {input.covar} \
            --modality {wildcards.modality} \
            --batch-start {wildcards.batch_start} \
            --batch-end {wildcards.batch_end} \
            --output-dir {params.twas_interm_dir} \
            --status {output} \
            --threads {threads} \
            --fusion-script {input.fusion} \
            --gcta {input.gcta}
        """

def twas_batch_status_input(wildcards):
    """Get the status files for all model-sized TWAS batches."""
    n_models = twas_model_count(wildcards.modality)
    if n_models < 1:
        raise ValueError(f'No TWAS models found for modality {wildcards.modality}')
    for start, end in twas_batch_ranges(n_models, config['twas_models_per_job']):
        yield interm_dir / 'twas' / f'status_{wildcards.modality}' / f'{start}_{end}.tsv'

rule twas_assemble_summary:
    """Summarize weights from all batches/genes.
    
    Batch status tables define which phenotypes produced weights and retain
    compact diagnostics for phenotypes that FUSION intentionally skipped.
    """
    input:
        status = twas_batch_status_input,
        assembler = 'scripts/assemble_twas_weight_list.py',
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
        python {input.assembler} \
            --status {input.status} \
            --weights-dir {params.twas_interm_dir} \
            --output {output.file_list}
        cd {params.twas_interm_dir}
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
            | sort -k 2b,2 \
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
            --files-from={input.file_list}
        """
