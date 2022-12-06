rule star_index:
    """Generate the index for STAR."""
    input:
        ref_genome = ref_genome,
        ref_anno = ref_anno,
    output:
        # Among other generated files:
        ref_dir / 'star_index' / 'SAindex',
        ref_dir / 'star_index' / 'chrNameLength.txt',
    params:
        index_dir = ref_dir / 'star_index',
        overhang = read_length - 1,
        genomeSAindexNbases = int(np.log2(genome_size) / 2 - 1),
    threads: 16
    resources:
        mem_mb = 60000,
        walltime = 4,
    shell:
        """
        mkdir -p {params.index_dir}
        STAR --runMode genomeGenerate \
            --outTmpDir {params.index_dir}/tmp \
            --outFileNamePrefix {params.index_dir}/STAR_ \
            --genomeDir {params.index_dir} \
            --genomeFastaFiles {input.ref_genome} \
            --sjdbGTFfile {input.ref_anno} \
            --sjdbOverhang {params.overhang} \
            --genomeSAindexNbases {params.genomeSAindexNbases} \
            --runThreadN {threads}
        """

def fastq_star_param(wildcards):
    """Get a string listing the fastq input files, with two lists if paired_end.
    
    This is the string supplied directly to the STAR command.
    """
    files = fastq_map[wildcards.sample_id]
    if paired_end:
        return ' '.join([','.join(files[0]), ','.join(files[1])])
    else:
        return ','.join(files)

rule star_align:
    """Align RNA-Seq reads for a sample using STAR."""
    input:
        fastq = fastq_inputs,
        index = ref_dir / 'star_index' / 'SAindex',
    output:
        interm_dir / 'bam' / '{sample_id}.Aligned.sortedByCoord.out.bam',
    params:
        fastq_list = fastq_star_param,
        index_dir = ref_dir / 'star_index',
        bam_dir = interm_dir / 'bam',
        prefix = str(interm_dir / 'bam' / '{sample_id}.'),
    threads: 16
    resources:
        mem_mb = 60000,
        walltime = 8,
    shell:
        """
        mkdir -p {params.bam_dir}
        STAR --runMode alignReads \
            --genomeDir {params.index_dir} \
            --readFilesIn {params.fastq_list} \
            --readFilesCommand 'zcat < ' \
            --twopassMode Basic \
            --outSAMstrandField intronMotif \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix {params.prefix} \
            --runThreadN {threads}
        """

rule shrink_bam:
    """Remove SEQ and QUAL fields from BAM file to reduce size."""
    input:
        interm_dir / 'bam' / '{sample_id}.Aligned.sortedByCoord.out.bam',
    output:
        interm_dir / 'bam' / '{sample_id}.bam',
    shell:
        """
        samtools view -h {input} \
            | awk -v OFS="\t" '{{if (substr($0, 1, 1) != "@") {{$10="*"; $11="*"}}; print ;}}' \
            | samtools view -h -b \
            > {output} \
            && rm {input}
        """

rule bam_to_fastq_paired_end:
    """Extract reads from BAM file for tool requiring unmapped reads as input"""
    input:
        bam = lambda w: bam_map[w.sample_id],
    output:
        fastq1 = interm_dir / 'fastq' / '{sample_id}.PE1.fastq.gz',
        fastq2 = interm_dir / 'fastq' / '{sample_id}.PE2.fastq.gz',
    params:
        fastq_dir = interm_dir / 'fastq',
        add_threads = lambda w, threads: threads - 1,
    threads: 8
    shell:
        """
        mkdir -p {params.fastq_dir}
        samtools fastq \
            -1 {output.fastq1} \
            -2 {output.fastq2} \
            --threads {params.add_threads} \
            {input.bam}
        """

rule bam_to_fastq_single_end:
    """Extract reads from BAM file for tool requiring unmapped reads as input"""
    input:
        bam = lambda w: bam_map[w.sample_id],
    output:
        fastq = interm_dir / 'fastq' / '{sample_id}.SE.fastq.gz',
    params:
        fastq_dir = interm_dir / 'fastq',
        add_threads = lambda w, threads: threads - 1,
    threads: 8
    shell:
        """
        mkdir -p {params.fastq_dir}
        samtools fastq \
            -o {output.fastq} \
            --threads {params.add_threads} \
            {input.bam}
        """
