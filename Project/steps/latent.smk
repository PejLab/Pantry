MODALITIES_FOR_LATENT = ['expression', 'isoforms', 'splicing', 'alt_TSS', 'alt_polyA', 'stability']
N_BATCHES = 1000 # After getting gene bins, count unique genes, divide by batch size, and round up to set N_BATCHES

localrules:
    get_chrom_lengths,
    latent_pheno_groups,

rule get_chrom_lengths:
    """Get chromosome lengths from genome FASTA index for bedtools"""
    input:
        lambda w: f'{ref_genome}.fai',
    output:
        ref_dir / 'latent' / 'chr_lengths.genome',
    shell:
        'cut -f1,2 {input} > {output}'

rule collapse_annotation:
    """Merge the exons for all isoforms of a gene into one set of non-overlapping exon regions."""
    input:
        ref_anno,
    output:
        ref_dir / 'latent' / 'collapsed.gtf',
    shell:
        'python scripts/latent-rna/scripts/collapse_annotation.py {input} {output} --collapse_only'

rule get_gene_bins:
    """Divide gene features into bins"""
    input:
        gtf = ref_dir / 'latent' / 'collapsed.gtf',
        chrom = ref_dir / 'latent' / 'chr_lengths.genome',
    output:
        ref_dir / 'latent' / 'gene_bins.bed.gz',
    params:
        n_bins = 10,
    shell:
        """
        python3 scripts/latent-rna/scripts/get_gene_bins.py \
            --gtf {input.gtf} \
            --chromosomes {input.chrom} \
            --n-bins {params.n_bins} \
            --output {output}
        """

rule bedtools_coverage:
    """Get RNA-Seq read coverage for feature bins"""
    input:
        bam = interm_dir / 'bam' / '{sample_id}.bam',
        bed = ref_dir / 'latent' / 'gene_bins.bed.gz',
        chrom = ref_dir / 'latent' / 'chr_lengths.genome',
    output:
        interm_dir / 'latent' / 'covg_sample' / '{sample_id}.txt',
    params:
        covg_sample_dir = interm_dir / 'latent' / 'covg_sample',
    shell:
        # -split is necessary I think to avoid counting coverage between spliced exons
        """
        mkdir -p {params.covg_sample_dir}
        bedtools coverage -split -sorted -counts \
            -a {input.bed} \
            -b {input.bam} \
            -g {input.chrom} \
            | cut -f7 > {output}
        """

rule prepare_coverage_full:
    """Assemble per-sample coverage count files, normalize, and bin into batches for PCA
    
    The flag output can be listed in the config so coverage for multiple datasets can be
    generated in preparation for model fitting across all datasets.
    """
    input:
        covg_sample = lambda w: expand(str(interm_dir / 'latent' / 'covg_sample' / '{sample_id}.txt'), sample_id=samples),
        bins = ref_dir / 'latent' / 'gene_bins.bed.gz',
    output:
        covg = expand(str(interm_dir / 'latent' / 'covg_batch_full' / 'covg_{batch}.npy'), batch=range(N_BATCHES)),
        regions = expand(str(interm_dir / 'latent' / 'covg_batch_full' / 'covg_{batch}.regions.tsv.gz'), batch=range(N_BATCHES)),
        flag = touch(str(output_dir / 'latent_full.coverage_done')),
    params:
        covg_batch_dir = str(interm_dir / 'latent' / 'covg_batch_full'),
        covg_file_list = str(interm_dir / 'latent' / 'covg_sample_files.txt'),
        batch_size = 200,
    shell:
        """
        mkdir -p {params.covg_batch_dir}
        printf '%s\\n' {input.covg_sample} > {params.covg_file_list}
        python scripts/latent-rna/latent_RNA.py prepare \
            --inputs {params.covg_file_list} \
            --regions {input.bins} \
            --output-dir {params.covg_batch_dir} \
            --batch-size {params.batch_size}
        """

rule prepare_coverage_residual:
    """Assemble per-sample coverage count files, normalize, regress out explicit phenotypes, and bin into batches for PCA
    
    The flag output can be listed in the config so coverage for multiple datasets can be
    generated in preparation for model fitting across all datasets.
    """
    input:
        covg_sample = lambda w: expand(str(interm_dir / 'latent' / 'covg_sample' / '{sample_id}.txt'), sample_id=samples),
        bins = ref_dir / 'latent' / 'gene_bins.bed.gz',
        phenos = expand(output_dir / '{modality}.bed.gz', modality=MODALITIES_FOR_LATENT),
    output:
        covg = expand(str(interm_dir / 'latent' / 'covg_batch_residual' / 'covg_{batch}.npy'), batch=range(N_BATCHES)),
        regions = expand(str(interm_dir / 'latent' / 'covg_batch_residual' / 'covg_{batch}.regions.tsv.gz'), batch=range(N_BATCHES)),
        flag = touch(str(output_dir / 'latent_residual.coverage_done')),
    params:
        covg_batch_dir = str(interm_dir / 'latent' / 'covg_batch_residual'),
        covg_file_list = str(interm_dir / 'latent' / 'covg_sample_files.txt'),
        phenos_list = str(interm_dir / 'latent' / 'pheno_files.txt'),
        batch_size = 200,
    shell:
        """
        mkdir -p {params.covg_batch_dir}
        printf '%s\\n' {input.covg_sample} > {params.covg_file_list}
        printf '%s\\n' {input.phenos} > {params.phenos_list}
        python scripts/latent-rna/latent_RNA.py prepare \
            --inputs {params.covg_file_list} \
            --regions {input.bins} \
            --pheno-paths-file {params.phenos_list} \
            --output-dir {params.covg_batch_dir} \
            --batch-size {params.batch_size}
        """

rule get_latent_phenotypes:
    """Apply PCA models to generate latent phenotypes"""
    input:
        covg_batch = expand(str(interm_dir / 'latent' / 'covg_batch_{{version}}' / 'covg_{batch}.npy'), batch=range(N_BATCHES)),
        models = expand(str(ref_dir / 'latent' / 'models_{{version}}' / 'models_{batch}.pickle'), batch=range(N_BATCHES)),
    output:
        phenos = interm_dir / 'latent' / 'latent_{version}.tsv.gz',
    params:
        covg_batch_dir = str(interm_dir / 'latent' / 'covg_batch_{version}'),
        models_dir = str(ref_dir / 'latent' / 'models_{version}'),
    shell:
        """
        python scripts/latent-rna/latent_RNA.py transform \
            --batch-covg-dir {params.covg_batch_dir} \
            --models-dir {params.models_dir} \
            --output {output}
        """

rule assemble_latent_bed:
    """Add gene info to latent phenotypes to generate BED file"""
    input:
        phenos = interm_dir / 'latent' / 'latent_{version}.tsv.gz',
        ref_anno = ref_anno,
    output:
        bed = interm_dir / 'unnorm' / 'latent_{version}.bed',
    params:
        unnorm_dir = interm_dir / 'unnorm',
    shell:
        """
        mkdir -p {params.unnorm_dir}
        python3 scripts/assemble_bed.py \
            --type latent \
            --input {input.phenos} \
            --ref_anno {input.ref_anno} \
            --output {output.bed}
        """

rule normalize_latent:
    """Quantile-normalize values for QTL mapping"""
    input:
        bed = interm_dir / 'unnorm' / 'latent_{version}.bed',
        samples = samples_file,
    output:
        output_dir / 'latent_{version}.bed.gz',
    params:
        bed = str(output_dir / 'latent_{version}.bed'),
    shell:
        """
        mkdir -p {output_dir}
        python3 scripts/normalize_phenotypes.py \
            --input {input.bed} \
            --samples {input.samples} \
            --output {params.bed}
        bgzip {params.bed}
        """

rule latent_pheno_groups:
    """Group phenotypes by gene for tensorQTL"""
    input:
        output_dir / 'latent_{version}.bed.gz',
    output:
        output_dir / 'latent_{version}.phenotype_groups.txt',
    shell:
        """
        zcat < {input} \
            | tail -n +2 \
            | cut -f4 \
            | awk '{{ g=$1; sub(/:.*$/, "", g); print $1 "\t" g }}' \
            > {output}
        """
