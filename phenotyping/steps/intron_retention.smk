from math import ceil

localrules:
    make_majiq_experiments_tsv,
    pheno_groups_intron_retention,

majiq_dir = interm_dir / 'intron_retention'
majiq_build_dir = majiq_dir / 'build'
majiq_psicov_dir = majiq_dir / 'psicov'
majiq_params = modality_groups['intron_retention']
majiq_annotation_gff3 = Path(majiq_params['annotation_gff3']).expanduser() if 'annotation_gff3' in majiq_params else ref_dir / 'majiq' / 'annotation.gff3'
majiq_license = Path(majiq_params.get('license', 'scripts/intron_retention/majiq_license_academic_official.lic')).expanduser()
majiq_threads = majiq_params.get('threads', 8)
majiq_psicov_batch_size = majiq_params.get('psicov_batch_size', 25)
majiq_min_experiments = majiq_params.get('min_experiments', 0.01)
majiq_batches = list(range(ceil(len(samples) / majiq_psicov_batch_size)))
majiq_sj = [
    directory(path)
    for path in expand(str(majiq_build_dir / '{sample_id}.sj'), sample_id=samples)
]


def majiq_batch_samples(batch):
    start = int(batch) * majiq_psicov_batch_size
    return samples[start : start + majiq_psicov_batch_size]


def majiq_batch_sj(wildcards):
    return expand(str(majiq_build_dir / '{sample_id}.sj'), sample_id=majiq_batch_samples(wildcards.batch))


def majiq_batch_prefixes(wildcards):
    return ' '.join(majiq_batch_samples(wildcards.batch))


if 'annotation_gff3' not in majiq_params:
    rule majiq_gff3_from_gtf:
        """Convert the reference GTF to GFF3 for MAJIQ."""
        input:
            ref_anno,
        output:
            majiq_annotation_gff3,
        params:
            outdir = ref_dir / 'majiq',
        resources:
            mem_mb = 32000,
        shell:
            """
            mkdir -p {params.outdir}
            python3 scripts/gtf_to_majiq_gff3.py \
                --input {input} \
                --output {output}
            """


rule make_majiq_experiments_tsv:
    """Write MAJIQ experiment groups and BAM paths."""
    input:
        bams = expand(str(interm_dir / 'bam' / '{sample_id}.bam'), sample_id=samples),
        bai = expand(str(interm_dir / 'bam' / '{sample_id}.bam.bai'), sample_id=samples),
    output:
        majiq_dir / 'experiments.tsv',
    params:
        group = majiq_params.get('group', 'pantry'),
    run:
        Path(output[0]).parent.mkdir(parents=True, exist_ok=True)
        with open(output[0], 'w') as out:
            out.write('group\tpath\n')
            for bam in input.bams:
                out.write(f'{params.group}\t{bam}\n')


rule majiq_build:
    """Build MAJIQ splice graph and per-sample junction coverage."""
    input:
        gff3 = majiq_annotation_gff3,
        experiments = rules.make_majiq_experiments_tsv.output,
    output:
        splicegraph = directory(majiq_build_dir / 'splicegraph.zarr'),
        sj = majiq_sj,
    params:
        outdir = majiq_build_dir,
        min_experiments = majiq_min_experiments,
        license = majiq_license,
    conda: 'majiq3'
    threads: majiq_threads
    resources:
        mem_mb = 64000,
        runtime = '12h',
    shell:
        """
        mkdir -p {params.outdir}
        rm -rf {output.splicegraph} {output.sj}
        majiq build \
            --license {params.license} \
            --all-introns \
            --min-experiments {params.min_experiments} \
            --nthreads {threads} \
            --overwrite \
            {input.gff3} \
            {input.experiments} \
            {params.outdir}
        """


rule majiq_psi_coverage:
    """Calculate MAJIQ PSI coverage in sample batches."""
    input:
        splicegraph = rules.majiq_build.output.splicegraph,
        sj = majiq_batch_sj,
    output:
        psicov = directory(majiq_psicov_dir / 'batch_{batch}.psicov.zarr'),
    params:
        prefixes = majiq_batch_prefixes,
        license = majiq_license,
    conda: 'majiq3'
    threads: majiq_threads
    resources:
        mem_mb = 32000,
    shell:
        """
        mkdir -p {majiq_psicov_dir}
        rm -rf {output.psicov}
        majiq psi-coverage \
            --license {params.license} \
            --nthreads {threads} \
            --prefixes {params.prefixes} \
            --overwrite \
            {input.splicegraph} \
            {output.psicov} \
            {input.sj}
        """


rule majiq_psi:
    """Estimate MAJIQ PSI values across all batches."""
    input:
        splicegraph = rules.majiq_build.output.splicegraph,
        psicov = expand(str(majiq_psicov_dir / 'batch_{batch}.psicov.zarr'), batch=majiq_batches),
    output:
        majiq_dir / 'majiq.psi.tsv',
    params:
        license = majiq_license,
    conda: 'majiq3'
    threads: majiq_threads
    resources:
        mem_mb = 64000,
    shell:
        """
        majiq psi \
            --license {params.license} \
            --splicegraph {input.splicegraph} \
            --output-tsv {output} \
            --quantiles 0.025 0.5 0.975 \
            --nthreads {threads} \
            --overwrite \
            {input.psicov}
        """


rule extract_retained_introns:
    """Extract retained-intron rows from MAJIQ PSI output."""
    input:
        majiq_dir / 'majiq.psi.tsv',
    output:
        majiq_dir / 'retained_intron_psi.tsv.gz',
    resources:
        mem_mb = 32000,
    shell:
        """
        python3 scripts/intron_retention/extract_ir_psi.py \
            --input {input} \
            --output {output}
        """


rule assemble_intron_retention_bed:
    """Convert MAJIQ retained-intron PSI values into BED format."""
    input:
        ir = rules.extract_retained_introns.output,
        ref_anno = ref_anno,
    output:
        bed = output_dir / 'unnorm' / 'intron_retention.bed',
    params:
        unnorm_dir = output_dir / 'unnorm',
    resources:
        mem_mb = 32000,
    shell:
        """
        mkdir -p {params.unnorm_dir}
        python3 scripts/assemble_bed.py intron-retention \
            --input {input.ir} \
            --ref-anno {input.ref_anno} \
            --output {output.bed}
        """


rule normalize_intron_retention:
    """Quantile-normalize retained-intron PSI values for QTL mapping."""
    input:
        bed = output_dir / 'unnorm' / 'intron_retention.bed',
        samples = samples_file,
    output:
        output_dir / 'intron_retention.bed.gz',
    params:
        bed = output_dir / 'intron_retention.bed',
    resources:
        mem_mb = 32000,
    shell:
        """
        python3 scripts/normalize_phenotypes.py \
            --input {input.bed} \
            --samples {input.samples} \
            --output {params.bed}
        bgzip {params.bed}
        """


rule pheno_groups_intron_retention:
    """Group retained-intron phenotypes by gene for tensorQTL."""
    input:
        output_dir / 'intron_retention.bed.gz',
    output:
        output_dir / 'intron_retention.phenotype_groups.txt',
    shell:
        """
        zcat < {input} \
            | tail -n +2 \
            | cut -f4 \
            | awk '{{ g=$1; sub(/__.*$/, "", g); print $1 "\t" g }}' \
            > {output}
        """
