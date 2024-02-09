rule exons_introns_from_GTF:
    """Extract exons and introns from GTF"""
    input:
        ref_anno = ref_anno,
    output:
        exon_gtf = ref_dir / 'constit_exons.gtf',
        intron_gtf = ref_dir / 'introns.gtf',
    shell:
        """
        sh scripts/exons_introns_from_gtf.sh \
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
        paired_flag = '-p' if paired_end else '',
        feature_id = lambda w: {'constit_exons': 'exon', 'introns': 'intron'}[w.feature_type],
        frac_overlap = lambda w: {'constit_exons': 1, 'introns': 0}[w.feature_type],
        # TODO add strandedness parameter
    threads: 8
    resources:
        walltime = 16,
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
        exon = expand(str(interm_dir / 'stability' / '{sample_id}.constit_exons.counts.txt'), sample_id=samples),
        intron = expand(str(interm_dir / 'stability' / '{sample_id}.introns.counts.txt'), sample_id=samples),
        samples = samples_file,
        ref_anno = ref_anno,
    output:
        bed = output_dir / 'unnorm' / 'stability.bed',
    params:
        unnorm_dir = output_dir / 'unnorm',
        stab_dir = interm_dir / 'stability',
    shell:
        """
        mkdir -p {params.unnorm_dir}
        python3 scripts/assemble_bed.py \
            --type stability \
            --input-dir {params.stab_dir} \
            --samples {input.samples} \
            --ref_anno {input.ref_anno} \
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
