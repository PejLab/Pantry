rule exonic_intronic_from_GTF:
    """Extract exonic and intronic regions from GTF"""
    input:
        ref_anno = ref_anno,
    output:
        exon_gtf = ref_dir / 'exonic.gtf',
        intron_gtf = ref_dir / 'intronic.gtf',
    shell:
        """
        python3 scripts/exonic_intronic_from_gtf.py \
            {input.ref_anno} \
            {output.exon_gtf} \
            {output.intron_gtf}
        """

rule run_featureCounts:
    """Run featureCounts from subread to get exon or intron read counts
    
    Intron counts represent reads that could only come from pre-mRNA, so they
    need only overlap an intron by one base. Exon counts represent reads that
    could come from pre-mRNA or mature mRNA, so they must be fully within an
    exon.
    """
    input:
        bam = interm_dir / 'bam' / '{sample_id}.bam',
        gtf = ref_dir / '{feature_type}.gtf',
    output:
        interm_dir / 'stability' / '{sample_id}.{feature_type}.counts.txt',
    params:
        stab_dir = interm_dir / 'stability',
        paired_flag = lambda w: '-p' if fastq_map[w.sample_id][1] else '',
        feature_id = lambda w: {'exonic': 'exon', 'intronic': 'intron'}[w.feature_type],
        frac_overlap = lambda w: {'exonic': 1, 'intronic': 0}[w.feature_type],
        # TODO add strandedness parameter
    threads: 8
    resources:
        runtime = '16h',
    shell:
        """
        mkdir -p {params.stab_dir}
        featureCounts \
            {input.bam} \
            {params.paired_flag} \
            -a {input.gtf} \
            -t {params.feature_id} \
            --fracOverlap {params.frac_overlap} \
            -T {threads} \
            -o {output}
        """

rule assemble_stability_bed:
    """Assemble exon to intron read ratios into mRNA stability BED file"""
    input:
        exon = expand(str(interm_dir / 'stability' / '{sample_id}.exonic.counts.txt'), sample_id=samples),
        intron = expand(str(interm_dir / 'stability' / '{sample_id}.intronic.counts.txt'), sample_id=samples),
        samples = samples_file,
        ref_anno = ref_anno,
    output:
        bed = output_dir / 'unnorm' / 'stability.bed',
    params:
        unnorm_dir = output_dir / 'unnorm',
        stab_dir = interm_dir / 'stability',
    resources:
        mem_mb = 32000,
    shell:
        """
        mkdir -p {params.unnorm_dir}
        python3 scripts/assemble_bed.py stability \
            --input-dir {params.stab_dir} \
            --samples {input.samples} \
            --ref-anno {input.ref_anno} \
            --output {output.bed}
        """

rule normalize_stability:
    """Quantile-normalize values for QTL mapping"""
    input:
        bed = output_dir / 'unnorm' / 'stability.bed',
        samples = samples_file,
    output:
        output_dir / 'stability.bed.gz',
    params:
        bed = output_dir / 'stability.bed',
    shell:
        """
        mkdir -p {output_dir}
        python3 scripts/normalize_phenotypes.py \
            --input {input.bed} \
            --samples {input.samples} \
            --output {params.bed}
        bgzip {params.bed}
        """
