localrules:
    laddr_manifest,
    laddr_config,
    pheno_groups_latent,

LADDR_PROJECT_DIR = laddr_config['project_dir']
LADDR_DATASET = laddr_config['dataset_name']
LADDR_MANIFEST = LADDR_PROJECT_DIR / 'coverage_manifest.tsv'
LADDR_CONFIG = LADDR_PROJECT_DIR / 'config.yaml'
LADDR_BW_DIR = LADDR_PROJECT_DIR / 'covg_bigwig' / LADDR_DATASET
LADDR_BINS_READY = LADDR_PROJECT_DIR / '.bins_ready'
LADDR_COVERAGE_DONE = LADDR_PROJECT_DIR / 'covg_norm' / LADDR_DATASET / '.coverage_done'
LADDR_FIT_DONE = LADDR_PROJECT_DIR / 'models' / '.fit_done'

def laddr_pheno_inputs(wildcards=None):
    return [output_dir / f'{name}.bed.gz' for name in laddr_config['regress_modalities']]

rule bam_to_bigwig:
    """Convert BAM alignments to base-resolution bigWig coverage."""
    input:
        bam = interm_dir / 'bam' / '{sample_id}.bam',
        bai = interm_dir / 'bam' / '{sample_id}.bam.bai',
    output:
        bw = LADDR_BW_DIR / '{sample_id}.bw',
    params:
        outdir = LADDR_BW_DIR,
        bin_size = laddr_config['bam_coverage_bin_size'],
    threads: laddr_config['bam_coverage_threads']
    shell:
        """
        mkdir -p {params.outdir}
        bamCoverage \
            -b {input.bam} \
            -o {output.bw} \
            -of bigwig \
            --binSize {params.bin_size} \
            -p {threads}
        """

rule laddr_manifest:
    """Write coverage manifest for LaDDR."""
    input:
        expand(str(LADDR_BW_DIR / '{sample_id}.bw'), sample_id=samples),
    output:
        manifest = LADDR_MANIFEST,
    params:
        dataset = LADDR_DATASET,
    run:
        Path(str(output.manifest)).parent.mkdir(parents=True, exist_ok=True)
        with open(output.manifest, 'w') as f:
            for sample_id in samples:
                f.write(f'{params.dataset}\t{sample_id}\t{params.dataset}/{sample_id}.bw\n')

rule laddr_config:
    """Write config used for LaDDR commands."""
    input:
        manifest = LADDR_MANIFEST,
        phenos = laddr_pheno_inputs,
    output:
        config = LADDR_CONFIG,
    params:
        project_dir = LADDR_PROJECT_DIR,
        manifest_rel = LADDR_MANIFEST.name,
        min_samples_expressed = laddr_config['min_samples_expressed'],
        binning_method = laddr_config['binning_method'],
        batch_size = laddr_config['batch_size'],
        max_bin_width = laddr_config['max_bin_width'],
        adaptive_max_samples = laddr_config['adaptive_max_samples'],
        adaptive_bins_per_gene = laddr_config['adaptive_bins_per_gene'],
        adaptive_min_mean_total_covg = laddr_config['adaptive_min_mean_total_covg'],
        adaptive_max_corr = laddr_config['adaptive_max_corr'],
        model_var_explained_max = laddr_config['model_var_explained_max'],
        model_n_pcs_max = laddr_config['model_n_pcs_max'],
        use_existing_bins_flag = '--use-existing-bins' if laddr_config['use_existing_bins'] else '',
    shell:
        """
        mkdir -p {params.project_dir}
        python3 scripts/write_laddr_config.py \
            --output {output.config} \
            --gtf {ref_anno} \
            --manifest {params.manifest_rel} \
            --coverage-directory covg_bigwig \
            --pheno-files {input.phenos:q} \
            {params.use_existing_bins_flag} \
            --min-samples-expressed {params.min_samples_expressed} \
            --binning-method {params.binning_method} \
            --batch-size {params.batch_size} \
            --max-bin-width {params.max_bin_width} \
            --adaptive-max-samples {params.adaptive_max_samples} \
            --adaptive-bins-per-gene {params.adaptive_bins_per_gene} \
            --adaptive-min-mean-total-covg {params.adaptive_min_mean_total_covg} \
            --adaptive-max-corr {params.adaptive_max_corr} \
            --model-var-explained-max {params.model_var_explained_max} \
            --model-n-pcs-max {params.model_n_pcs_max}
        """

if laddr_config['use_existing_bins']:
    rule laddr_prepare_existing_bins:
        """Attach existing LaDDR info and bins directories."""
        input:
            config = LADDR_CONFIG,
        output:
            ready = LADDR_BINS_READY,
        params:
            project_dir = LADDR_PROJECT_DIR,
            info_dir = laddr_config['existing_info_dir'],
            bins_dir = laddr_config['existing_bins_dir'],
        shell:
            """
            mkdir -p {params.project_dir}
            ln -sfn {params.info_dir:q} {params.project_dir}/info
            ln -sfn {params.bins_dir:q} {params.project_dir}/gene_bins
            test -s {params.project_dir}/info/n_batches.txt
            test -d {params.project_dir}/gene_bins
            touch {output.ready}
            """
else:
    rule laddr_setup:
        """Prepare gene metadata and batch layout for LaDDR."""
        input:
            config = LADDR_CONFIG,
            manifest = LADDR_MANIFEST,
            bw = expand(str(LADDR_BW_DIR / '{sample_id}.bw'), sample_id=samples),
        output:
            genes = LADDR_PROJECT_DIR / 'info' / 'genes.tsv',
            exons = LADDR_PROJECT_DIR / 'info' / 'exons.tsv.gz',
            n_batches = LADDR_PROJECT_DIR / 'info' / 'n_batches.txt',
            median_covg = LADDR_PROJECT_DIR / 'info' / 'median_coverage.txt',
        params:
            project_dir = LADDR_PROJECT_DIR,
        shell:
            """
            laddr setup -p {params.project_dir} -c {input.config}
            """

    rule laddr_binning:
        """Generate bins for all LaDDR batches."""
        input:
            config = LADDR_CONFIG,
            genes = LADDR_PROJECT_DIR / 'info' / 'genes.tsv',
            exons = LADDR_PROJECT_DIR / 'info' / 'exons.tsv.gz',
            n_batches = LADDR_PROJECT_DIR / 'info' / 'n_batches.txt',
            bw = expand(str(LADDR_BW_DIR / '{sample_id}.bw'), sample_id=samples),
        output:
            ready = LADDR_BINS_READY,
        params:
            project_dir = LADDR_PROJECT_DIR,
        shell:
            """
            laddr binning -p {params.project_dir} -c {input.config}
            test -d {params.project_dir}/gene_bins
            touch {output.ready}
            """

rule laddr_coverage:
    """Generate normalized, binned coverage and regress out explicit phenotypes."""
    input:
        config = LADDR_CONFIG,
        manifest = LADDR_MANIFEST,
        bins_ready = LADDR_BINS_READY,
        bw = expand(str(LADDR_BW_DIR / '{sample_id}.bw'), sample_id=samples),
        phenos = laddr_pheno_inputs,
    output:
        done = LADDR_COVERAGE_DONE,
    params:
        project_dir = LADDR_PROJECT_DIR,
        dataset = LADDR_DATASET,
    shell:
        """
        laddr coverage -p {params.project_dir} -c {input.config} -d {params.dataset}
        touch {output.done}
        """

rule laddr_fit:
    """Fit latent phenotype models from normalized coverage."""
    input:
        config = LADDR_CONFIG,
        coverage_done = LADDR_COVERAGE_DONE,
    output:
        done = LADDR_FIT_DONE,
    params:
        project_dir = LADDR_PROJECT_DIR,
    shell:
        """
        laddr fit -p {params.project_dir} -c {input.config}
        touch {output.done}
        """

rule laddr_transform:
    """Generate latent phenotypes from fitted models."""
    input:
        config = LADDR_CONFIG,
        coverage_done = LADDR_COVERAGE_DONE,
        fit_done = LADDR_FIT_DONE,
    output:
        phenos = LADDR_PROJECT_DIR / 'phenotypes' / f'latent_phenos.{LADDR_DATASET}.tsv.gz',
    params:
        project_dir = LADDR_PROJECT_DIR,
        dataset = LADDR_DATASET,
    shell:
        """
        laddr transform -p {params.project_dir} -c {input.config} -d {params.dataset}
        """

rule assemble_latent_bed:
    """Convert latent phenotypes to BED format."""
    input:
        latent = LADDR_PROJECT_DIR / 'phenotypes' / f'latent_phenos.{LADDR_DATASET}.tsv.gz',
        ref_anno = ref_anno,
    output:
        bed = output_dir / 'unnorm' / 'latent_residual.bed',
    params:
        unnorm_dir = output_dir / 'unnorm',
    resources:
        mem_mb = 32000,
    shell:
        """
        mkdir -p {params.unnorm_dir}
        python3 scripts/assemble_bed.py latent \
            --input {input.latent} \
            --ref-anno {input.ref_anno} \
            --output {output.bed}
        """

rule normalize_latent:
    """Quantile-normalize latent phenotypes for QTL mapping."""
    input:
        bed = output_dir / 'unnorm' / 'latent_residual.bed',
        samples = samples_file,
    output:
        output_dir / 'latent_residual.bed.gz',
    params:
        bed = output_dir / 'latent_residual.bed',
    resources:
        mem_mb = 32000,
    shell:
        """
        mkdir -p {output_dir}
        python3 scripts/normalize_phenotypes.py \
            --input {input.bed} \
            --samples {input.samples} \
            --output {params.bed}
        bgzip {params.bed}
        """

rule pheno_groups_latent:
    """Group latent phenotypes by gene for tensorQTL."""
    input:
        output_dir / 'latent_residual.bed.gz',
    output:
        output_dir / 'latent_residual.phenotype_groups.txt',
    shell:
        """
        zcat < {input} \
            | tail -n +2 \
            | cut -f4 \
            | awk '{{ g=$1; sub(/__.*$/, "", g); print $1 "\t" g }}' \
            > {output}
        """
