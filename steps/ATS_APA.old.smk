mc_dir = Path('TURNAP/src/mountainClimber/src')
mc_env = '../envs/mountainClimber.yml'

localrules:
    get_chrom_lengths,

rule get_chrom_lengths:
    """Get chromosome lengths from genome FASTA header lines"""
    input:
        ref_genome,
    output:
        ref_dir / 'chr_lengths.genome',
    shell:
        """
        grep '^>' {input} \
            | sed 's/^>//' | tr ':' ' ' | tr ' ' '\t' | cut -f1,8 \
            > {output}
        """

rule get_transcript_annotations:
    """The GTF entries must contain gene_name and transcript_id fields
    
    And only transcript entries are used, so filter for smaller file size.
    """
    input:
        ref_anno,
    output:
        ref_dir / 'transcripts.gtf',
    shell:
        """
        awk '$3=="transcript"' {input} \
            | grep gene_name | grep transcript_id \
            > {output}
        """

rule get_genome_coverage:
    """Get RNA-Seq read coverage for mountainClimber"""
    input:
        project_dir / 'bam' / '{sample_id}.Aligned.sortedByCoord.out.bam',
    output:
        project_dir / 'ATS_APA' / '{sample_id}.bedgraph',
    params:
        ATS_APA_dir = project_dir / 'ATS_APA',
    shell:
        """
        mkdir -p {params.ATS_APA_dir}
        bedtools genomecov -trackline -bg -split \
            -ibam {input} \
            | bedtools sort \
            > {output}
        """

rule mountainClimber_junction_counts:
    """Get splice junction counts for mountainClimber"""
    input:
        project_dir / 'bam' / '{sample_id}.Aligned.sortedByCoord.out.bam',
    output:
        project_dir / 'ATS_APA' / '{sample_id}.jxn.bed',
    params:
        script = mc_dir / 'get_junction_counts.py',
    conda:
        mc_env,
    shell:
        # conda activate py2
        """
        python2 {params.script} \
            -i {input} \
            -s fr-unstrand \
            -o {output}
        """

rule mountainClimber_TU:
    """Call transcription units de novo"""
    input:
        coverage = project_dir / 'ATS_APA' / '{sample_id}.bedgraph',
        junctions = project_dir / 'ATS_APA' / '{sample_id}.jxn.bed',
        chr_lengths = ref_dir / 'chr_lengths.genome',
    output:
        project_dir / 'ATS_APA' / '{sample_id}.tu.bed',
    params:
        script = mc_dir / 'mountainClimberTU.py',
        strandedness = 0,  # 1, -1, or 0
    conda:
        mc_env,
    shell:
        # conda activate py2
        """
        python2 {params.script} \
            -b {input.coverage} \
            -j {input.junctions} \
            -s {params.strandedness} \
            -g {input.chr_lengths} \
            -o {output}
        """

rule mountainClimber_merge_TUs:
    input:
        tu = lambda w: expand(str(project_dir / 'ATS_APA' / '{sample_id}.tu.bed'), sample_id=samples),
        gtf = ref_dir / 'transcripts.gtf',
    output:
        project_dir / 'ATS_APA' / 'tus_merged.annot.transcripts_singleGenes.bed'
    params:
        script = mc_dir / 'merge_tus.py',
        strandedness = 'n',  # 'y' or 'n'
        out_prefix = project_dir / 'ATS_APA' / 'tus_merged',
    conda:
        mc_env,
    shell:
        # conda activate py2
        """
        python2 {params.script} \
            -i {input.tu} \
            -s {params.strandedness} \
            -g {input.gtf} \
            -o {params.out_prefix}
        """

rule mountainClimber_CP:
    """Identify change points in read coverage in each TU"""
    input:
        coverage = project_dir / 'ATS_APA' / '{sample_id}.bedgraph',
        tus_merged = project_dir / 'ATS_APA' / 'tus_merged.annot.transcripts_singleGenes.bed',
        junctions = project_dir / 'ATS_APA' / '{sample_id}.jxn.bed',
        ref_genome = ref_genome,
    output:
        project_dir / 'ATS_APA' / '{sample_id}.cp.bed',
    params:
        script = mc_dir / 'mountainClimberCP.py',
    conda:
        mc_env,
    shell:
        # conda activate py2
        """
        python2 {params.script} \
            -i {input.coverage} \
            -g {input.tus_merged} \
            -j {input.junctions} \
            -x {input.ref_genome} \
            -o {output}
        """

rule mountainClimber_RU:
    """Calculate relative usage of changepoints"""
    input:
        project_dir / 'ATS_APA' / '{sample_id}.cp.bed',
    output:
        project_dir / 'ATS_APA' / '{sample_id}.ru.bed',
    params:
        script = mc_dir / 'mountainClimberCP.py',
    conda:
        mc_env,
    shell:
        # conda activate py2
        """
        python2 {params.script} -i {input} -o {output}
        """
