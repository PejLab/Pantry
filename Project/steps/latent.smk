localrules:
    latent_pheno_groups,

rule get_gene_bins:
    """Divide gene features into bins"""
    input:
        ref_anno = ref_anno,
        chrom = ref_dir / 'star_index' / 'chrNameLength.txt',
    output:
        ref_dir / 'gene_bins.bed.gz',
    params:
        project_dir = project_dir,
        n_bins = 10,
        bed = ref_dir / 'gene_bins.bed',
    shell:
        """
        python3 {params.project_dir}/src/get_gene_bins.py \
            -g {input.ref_anno} \
            -c {input.chrom} \
            -n {params.n_bins} \
            -o {params.bed}
        bgzip {params.bed}
        """

rule bedtools_coverage:
    """Get RNA-Seq read coverage for feature bins"""
    input:
        bam = interm_dir / 'bam' / '{sample_id}.Aligned.sortedByCoord.out.bam',
        bed = ref_dir / 'gene_bins.bed.gz',
        chrom = ref_dir / 'star_index' / 'chrNameLength.txt',
    output:
        interm_dir / 'latent' / '{sample_id}.bed.gz',
    params:
        latent_dir = interm_dir / 'latent',
    shell:
        # -split is necessary I think to avoid counting coverage between spliced exons
        """
        mkdir -p {params.latent_dir}
        bedtools coverage -split -sorted -counts \
            -a {input.bed} \
            -b {input.bam} \
            -g {input.chrom} \
            | bgzip -c > {output}
        """

rule assemble_latent_bed:
    """Run PCA on feature bin coverage and create BED file"""
    input:
        beds = expand(str(interm_dir / 'latent' / '{sample_id}.bed.gz'), sample_id=samples),
        ref_anno = ref_anno,
    output:
        interm_dir / 'unnorm' / 'latent.bed',
    params:
        unnorm_dir = interm_dir / 'unnorm',
        project_dir = project_dir,
        bedfile_list = interm_dir / 'latent' / 'bedfiles.txt',
        var_expl = 0.95,
    shell:
        """
        mkdir -p {params.unnorm_dir}
        printf '%s\\n' {input.beds} > {params.bedfile_list}
        python3 {params.project_dir}/src/get_PC_features.py \
            -i {params.bedfile_list} \
            -g {input.ref_anno} \
            -v {params.var_expl} \
            -o {output}
        """

rule normalize_latent:
    """Quantile-normalize values for QTL mapping"""
    input:
        bed = interm_dir / 'unnorm' / 'latent.bed',
        samples = samples_file,
    output:
        output_dir / 'latent.bed.gz',
    params:
        project_dir = project_dir,
        bed = output_dir / 'latent.bed',
    shell:
        """
        python3 {params.project_dir}/src/normalize_phenotypes.py \
            --input {input.bed} \
            --samples {input.samples} \
            --output {params.bed}
        bgzip {params.bed}
        """

rule latent_pheno_groups:
    """Group phenotypes by gene for tensorQTL"""
    input:
        output_dir / 'latent.bed.gz',
    output:
        output_dir / 'latent.phenotype_groups.txt',
    shell:
        """
        zcat < {input} \
            | tail -n +2 \
            | cut -f4 \
            | awk '{{ g=$1; sub(/:.*$/, "", g); print $1 "\t" g }}' \
            > {output}
        """

# rule get_chrom_lengths:
#     """Get chromosome lengths from genome FASTA header lines for pyBedGraph"""
#     input:
#         ref_genome,
#     output:
#         ref_dir / 'chr_lengths.genome',
#     shell:
#         """
#         grep '^>' {input} \
#             | sed 's/^>//' | tr ':' ' ' | tr ' ' '\t' | cut -f1,8 \
#             > {output}
#         """

# rule get_genome_coverage:
#     """Get RNA-Seq read coverage"""
#     input:
#         project_dir / 'bam' / '{sample_id}.Aligned.sortedByCoord.out.bam',
#     output:
#         project_dir / 'latent' / '{sample_id}.bedgraph',
#     params:
#         latent_dir = project_dir / 'latent',
#     shell:
#         # -split is necessary I think to avoid counting coverage between spliced exons
#         """
#         mkdir -p {params.latent_dir}
#         bedtools genomecov -trackline -bg -split \
#             -ibam {input} \
#             | bedtools sort \
#             > {output}
#         """

# rule coverage_for_gene_bins:
#     """Extract coverage for bins of gene feature regions"""
#     input:
#         bedgraph = project_dir / 'latent' / '{sample_id}.bedgraph',
#         regions = ref_dir / 'gene_bins.tsv.gz',
#         chr_lengths = ref_dir / 'chr_lengths.genome',
#     output:
#         project_dir / 'latent' / '{sample_id}.bin_covg.txt.gz',
#     shell:
#         """
#         python3 TURNAP/src/region_coverage.py \
#             -b {input.bedgraph} \
#             -r {input.regions} \
#             -c {input.chr_lengths} \
#             -o {output}
#         """
