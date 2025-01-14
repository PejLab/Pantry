rule star_index:
    """Generate the index for STAR."""
    input:
        ref_genome = ref_genome,
        ref_anno = ref_anno,
    output:
        # Among other generated files:
        ref_dir / f'star_index_{read_length}' / 'SAindex',
        ref_dir / f'star_index_{read_length}' / 'chrNameLength.txt',
    params:
        index_dir = ref_dir / f'star_index_{read_length}',
        overhang = read_length - 1,
        genomeSAindexNbases = int(np.log2(genome_size) / 2 - 1),
    threads: 16
    resources:
        mem_mb = 60000,
        runtime = '4h',
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
        index = ref_dir / f'star_index_{read_length}' / 'SAindex',
    output:
        interm_dir / 'star_out' / '{sample_id}.Aligned.sortedByCoord.out.bam',
    params:
        fastq_list = fastq_star_param,
        index_dir = ref_dir / f'star_index_{read_length}',
        star_out_dir = interm_dir / 'star_out',
        prefix = str(interm_dir / 'star_out' / '{sample_id}.'),
    threads: 16
    resources:
        mem_mb = 60000,
        runtime = '8h',
    shell:
        """
        mkdir -p {params.star_out_dir}
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
        interm_dir / 'star_out' / '{sample_id}.Aligned.sortedByCoord.out.bam',
    output:
        interm_dir / 'bam' / '{sample_id}.bam',
    params:
        bam_dir = interm_dir / 'bam',
    shell:
        """
        mkdir -p {params.bam_dir}
        samtools view -h {input} \
            | awk -v OFS="\t" '{{if (substr($0, 1, 1) != "@") {{$10="*"; $11="*"}}; print ;}}' \
            | samtools view -h -b \
            > {output}
        """
