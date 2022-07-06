localrules:
    splicing_pheno_groups,

rule regtools_junctions:
    """Run regtools junctions for LeafCutter"""
    input:
        bam = project_dir / 'bam' / '{sample_id}.Aligned.sortedByCoord.out.bam',
        bai = project_dir / 'bam' / '{sample_id}.Aligned.sortedByCoord.out.bam.bai',
    output:
        project_dir / 'splicing' / '{sample_id}.junc'
    params:
        splice_dir = project_dir / 'splicing',
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
        lambda w: expand(str(project_dir / 'splicing' / '{sample_id}.junc'), sample_id=samples),
    output:
        project_dir / 'splicing' / 'leafcutter_perind_numers.counts.gz',
    params:
        juncfile_list = project_dir / 'splicing' / 'juncfiles.txt',
        splice_dir = project_dir / 'splicing',
        max_intron_len = 100000,
        min_clust_reads = 30,
        min_clust_ratio = 0.001,
    shell:
        # TODO get proper path to script
        """
        printf '%s\\n' {input} > {params.juncfile_list}
        python3 TURNAP/src/leafcutter_cluster_regtools_py3.py \
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
        counts = project_dir / 'splicing' / 'leafcutter_perind_numers.counts.gz',
        ref_anno = ref_anno,
    output:
        bed = project_dir / 'unnorm' / 'splicing.bed',
    shell:
        """
        python3 TURNAP/src/assemble_bed.py \
            --type splicing \
            --input {input.counts} \
            --ref_anno {input.ref_anno} \
            --output {output.bed}
        """

rule normalize_splicing:
    """Quantile-normalize values for QTL mapping"""
    input:
        bed = project_dir / 'unnorm' / 'splicing.bed',
        samples = samples_file,
    output:
        project_dir / 'splicing.bed.gz',
    params:
        bed = project_dir / 'splicing.bed',
    shell:
        """
        python3 TURNAP/src/normalize_phenotypes.py \
            --input {input.bed} \
            --samples {input.samples} \
            --output {params.bed}
        bgzip {params.bed}
        """

rule splicing_pheno_groups:
    """Group phenotypes by gene for tensorQTL"""
    input:
        project_dir / 'splicing.bed.gz',
    output:
        project_dir / 'splicing.phenotype_groups.txt',
    shell:
        """
        gzcat {input} \
            | tail -n +2 \
            | cut -f4 \
            | awk '{{ g=$1; sub(/:.*$/, "", g); print $1 "\t" g }}' \
            > {output}
        """
