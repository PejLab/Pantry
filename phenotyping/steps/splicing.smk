localrules:
    pheno_groups_splicing,

rule regtools_junctions:
    """Run regtools junctions for LeafCutter"""
    input:
        bam = interm_dir / 'bam' / '{sample_id}.bam',
        bai = interm_dir / 'bam' / '{sample_id}.bam.bai',
    output:
        interm_dir / 'splicing' / '{sample_id}.junc'
    params:
        splice_dir = interm_dir / 'splicing',
        min_anchor_len = 8,
        min_intron_len = 50,
        max_intron_len = 500000,
        strandedness = "XS",
    shell:
        """
        mkdir -p {params.splice_dir}
        regtools junctions extract \
            -a {params.min_anchor_len} \
            -m {params.min_intron_len} \
            -M {params.max_intron_len} \
            -s {params.strandedness} \
            -o {output} \
            {input.bam}
        """

rule cluster_junctions:
    """Cluster splice junctions with common boundaries"""
    input:
        expand(str(interm_dir / 'splicing' / '{sample_id}.junc'), sample_id=samples),
    output:
        interm_dir / 'splicing' / 'leafcutter_perind_numers.counts.gz',
    params:
        juncfile_list = interm_dir / 'splicing' / 'juncfiles.txt',
        splice_dir = interm_dir / 'splicing',
        max_intron_len = 100000,
        min_clust_reads = 30,
        min_clust_ratio = 0.001,
    shell:
        """
        printf '%s\\n' {input} > {params.juncfile_list}
        python3 scripts/leafcutter_cluster_regtools_py3.py \
            --juncfiles {params.juncfile_list} \
            --rundir {params.splice_dir} \
            --maxintronlen {params.max_intron_len} \
            --minclureads {params.min_clust_reads} \
            --mincluratio {params.min_clust_ratio}
        """

rule assemble_splicing_bed:
    """Convert leafcutter output into splicing BED file"""
    input:
        counts = interm_dir / 'splicing' / 'leafcutter_perind_numers.counts.gz',
        ref_anno = ref_anno,
    output:
        bed = output_dir / 'unnorm' / 'splicing.bed',
    params:
        unnorm_dir = output_dir / 'unnorm',
    resources:
        mem_mb = 32000,
    shell:
        """
        mkdir -p {params.unnorm_dir}
        python3 scripts/assemble_bed.py splicing \
            --input {input.counts} \
            --ref-anno {input.ref_anno} \
            --output {output.bed}
        """

rule normalize_splicing:
    """Quantile-normalize values for QTL mapping"""
    input:
        bed = output_dir / 'unnorm' / 'splicing.bed',
        samples = samples_file,
    output:
        output_dir / 'splicing.bed.gz',
    params:
        bed = output_dir / 'splicing.bed',
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

rule pheno_groups_splicing:
    """Group phenotypes by gene for tensorQTL"""
    input:
        output_dir / 'splicing.bed.gz',
    output:
        output_dir / 'splicing.phenotype_groups.txt',
    shell:
        """
        zcat < {input} \
            | tail -n +2 \
            | cut -f4 \
            | awk '{{ g=$1; sub(/__.*$/, "", g); print $1 "\t" g }}' \
            > {output}
        """
