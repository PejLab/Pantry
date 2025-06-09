localrules:
    map_edit_sites_to_genes,

rule query_editing_level:
    input:
        edit_sites = edit_sites_bed,
        ref_genome = ref_genome,
        bam = interm_dir / 'star_out' / '{sample_id}.Aligned.sortedByCoord.out.bam',
    output:
        levels = interm_dir / 'RNA_editing' / 'edit_levels' / '{sample_id}.rnaeditlevel.tsv.gz',
    params:
        edit_levels_dir = interm_dir / 'RNA_editing' / 'edit_levels',
    resources:
        runtime = '8h',
    shell:
        """
        mkdir -p {params.edit_levels_dir}
        python3 scripts/RNA_editing/query_editing_level.py \
            --edit_sites {input.edit_sites} \
            --ref {input.ref_genome} \
            --bam {input.bam} \
            --output {output.levels}
        """

rule shared_sample_sites:
    input:
        levels = expand(interm_dir / 'RNA_editing' / 'edit_levels' / '{sample_id}.rnaeditlevel.tsv.gz', sample_id=samples),
        samples = samples_file,
    output:
        matrix = interm_dir / 'RNA_editing' / f'edMat.{edit_sites_min_coverage}cov.{edit_sites_min_samples}samps.tsv',
    params:
        edit_levels_dir = interm_dir / 'RNA_editing' / 'edit_levels',
        min_cov = edit_sites_min_coverage,
        min_sam = edit_sites_min_samples,
    shell:
        """
        python3 scripts/RNA_editing/shared_samples_sites_matrix.py \
            --path_to_edit_files {params.edit_levels_dir} \
            --samples_file {input.samples} \
            --output_file {output.matrix} \
            --min_coverage {params.min_cov} \
            --min_samples {params.min_sam}
        """

rule map_edit_sites_to_genes:
    input:
        edit_sites = edit_sites_bed,
        ref_anno = ref_anno,
    output:
        tsv = ref_dir / 'edit_sites_to_genes.tsv',
    shell:
        """
        Rscript scripts/RNA_editing/map_edit_sites_to_genes.R {input.edit_sites} {input.ref_anno} {output.tsv}
        """

rule assemble_RNA_editing_bed:
    """Convert RNA editing matrix into BED file"""
    input:
        matrix = interm_dir / 'RNA_editing' / f'edMat.{edit_sites_min_coverage}cov.{edit_sites_min_samples}samps.tsv',
        edit_sites_to_genes = ref_dir / 'edit_sites_to_genes.tsv',
        ref_anno = ref_anno,
    output:
        bed = output_dir / 'unnorm' / 'RNA_editing.bed',
    params:
        unnorm_dir = output_dir / 'unnorm',
    resources:
        mem_mb = 32000,
    shell:
        """
        mkdir -p {params.unnorm_dir}
        python3 scripts/assemble_bed.py \
            --type RNA_editing \
            --input {input.matrix} \
            --ref_anno {input.ref_anno} \
            --edit_sites_to_genes {input.edit_sites_to_genes} \
            --output {output.bed}
        """

rule normalize_RNA_editing:
    """Quantile-normalize values for QTL mapping"""
    input:
        bed = output_dir / 'unnorm' / 'RNA_editing.bed',
        samples = samples_file,
    output:
        output_dir / 'RNA_editing.bed.gz',
    params:
        bed = output_dir / 'RNA_editing.bed',
    shell:
        """
        mkdir -p {output_dir}
        python3 scripts/normalize_phenotypes.py \
            --input {input.bed} \
            --samples {input.samples} \
            --output {params.bed}
        bgzip {params.bed}
        """

rule RNA_editing_pheno_groups:
    """Group phenotypes by gene for tensorQTL"""
    input:
        output_dir / 'RNA_editing.bed.gz',
    output:
        output_dir / 'RNA_editing.phenotype_groups.txt',
    shell:
        """
        zcat < {input} \
            | tail -n +2 \
            | cut -f4 \
            | awk '{{ g=$1; sub(/__.*$/, "", g); print $1 "\t" g }}' \
            > {output}
        """
