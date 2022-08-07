localrules:
    splicing_pheno_groups,

rule regtools_junctions:
    """Run regtools junctions for LeafCutter"""
    input:
        bam = interm_dir / 'bam' / '{sample_id}.Aligned.sortedByCoord.out.bam',
        bai = interm_dir / 'bam' / '{sample_id}.Aligned.sortedByCoord.out.bam.bai',
    output:
        interm_dir / 'splicing' / '{sample_id}.junc'
    params:
        splice_dir = interm_dir / 'splicing',
        min_anchor_len = 8,
        min_intron_len = 50,
        max_intron_len = 500000,
        strandedness = 0,
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
        lambda w: expand(str(interm_dir / 'splicing' / '{sample_id}.junc'), sample_id=samples),
    output:
        interm_dir / 'splicing' / 'leafcutter_perind_numers.counts.gz',
    params:
        juncfile_list = interm_dir / 'splicing' / 'juncfiles.txt',
        project_dir = project_dir,
        splice_dir = interm_dir / 'splicing',
        max_intron_len = 100000,
        min_clust_reads = 30,
        min_clust_ratio = 0.001,
    shell:
        """
        printf '%s\\n' {input} > {params.juncfile_list}
        python3 {params.project_dir}/src/leafcutter_cluster_regtools_py3.py \
            --juncfiles {params.juncfile_list} \
            --rundir {params.splice_dir} \
            --maxintronlen {params.max_intron_len} \
            --minclureads {params.min_clust_reads} \
            --mincluratio {params.min_clust_ratio} \
            --checkchrom False
        """

rule assemble_splicing_bed:
    """Convert leafcutter output into splicing BED file"""
    input:
        counts = interm_dir / 'splicing' / 'leafcutter_perind_numers.counts.gz',
        ref_anno = ref_anno,
    output:
        bed = interm_dir / 'unnorm' / 'splicing.bed',
    params:
        unnorm_dir = interm_dir / 'unnorm',
        project_dir = project_dir,
    shell:
        """
        mkdir -p {params.unnorm_dir}
        python3 {params.project_dir}/src/assemble_bed.py \
            --type splicing \
            --input {input.counts} \
            --ref_anno {input.ref_anno} \
            --output {output.bed}
        """

rule normalize_splicing:
    """Quantile-normalize values for QTL mapping"""
    input:
        bed = interm_dir / 'unnorm' / 'splicing.bed',
        samples = samples_file,
    output:
        output_dir / 'splicing.bed.gz',
    params:
        project_dir = project_dir,
        bed = output_dir / 'splicing.bed',
    shell:
        """
        python3 {params.project_dir}/src/normalize_phenotypes.py \
            --input {input.bed} \
            --samples {input.samples} \
            --output {params.bed}
        bgzip {params.bed}
        """

rule splicing_pheno_groups:
    """Group phenotypes by gene for tensorQTL"""
    input:
        output_dir / 'splicing.bed.gz',
    output:
        output_dir / 'splicing.phenotype_groups.txt',
    shell:
        """
        gzcat {input} \
            | tail -n +2 \
            | cut -f4 \
            | awk '{{ g=$1; sub(/:.*$/, "", g); print $1 "\t" g }}' \
            > {output}
        """
