rule exons_introns_from_GTF:
    """Extract exons and introns from GTF"""
    input:
        ref_anno = ref_anno,
    output:
        exon_gtf = ref_dir / 'constit_exons.gtf',
        intron_gtf = ref_dir / 'introns.gtf',
    shell:
        """
        sh TURNAP/src/exons_introns_from_gtf.sh \
            {input.ref_anno} \
            {output.exon_gtf} \
            {output.intron_gtf}
        """

rule run_featureCounts:
    """Run featureCounts from subread to get exon or intron read counts"""
    input:
        bam = project_dir / 'bam' / '{sample_id}.Aligned.sortedByCoord.out.bam',
        gtf = ref_dir / '{feature_type}.gtf',
    output:
        project_dir / 'stability' / '{sample_id}.{feature_type}.counts.txt',
    params:
        stab_dir = project_dir / 'stability',
        paired_flag = '-p' if paired_end else '',
        feature_id = lambda w: {'constit_exons': 'exon', 'introns': 'intron'}[w.feature_type],
        # TODO add strandedness parameter
    resources:
        cpus = threads,
        walltime = 8,
    shell:
        """
        mkdir -p {params.stab_dir}
        featureCounts \
            {input.bam} \
            {params.paired_flag} \
            -a {input.gtf} \
            -t {params.feature_id} \
            -T {resources.cpus} \
            -o {output}
        """

rule assemble_stability_bed:
    """Assemble exon to intron read ratios into mRNA stability BED file"""
    input:
        exon = lambda w: expand(str(project_dir / 'stability' / '{sample_id}.constit_exons.counts.txt'), sample_id=samples),
        intron = lambda w: expand(str(project_dir / 'stability' / '{sample_id}.introns.counts.txt'), sample_id=samples),
        samples = samples_file,
        ref_anno = ref_anno,
    output:
        bed = project_dir / 'unnorm' / 'stability.bed',
    params:
        # samples = samples,
        stab_dir = project_dir / 'stability',
    shell:
        """
        python3 TURNAP/src/assemble_bed.py \
            --type stability \
            --input-dir {params.stab_dir} \
            --samples {input.samples} \
            --ref_anno {input.ref_anno} \
            --output {output.bed}
        """

rule normalize_stability:
    """Quantile-normalize values for QTL mapping"""
    input:
        bed = project_dir / 'unnorm' / 'stability.bed',
        samples = samples_file,
    output:
        project_dir / 'stability.bed.gz',
    params:
        bed = project_dir / 'stability.bed',
    shell:
        """
        python3 TURNAP/src/normalize_phenotypes.py \
            --input {input.bed} \
            --samples {input.samples} \
            --output {params.bed}
        bgzip {params.bed}
        """
