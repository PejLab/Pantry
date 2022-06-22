rule star_index:
    """Generate the index for STAR."""
    input:
        ref_genome = ref_genome,
        ref_anno = ref_anno,
    output:
        # Among other generated files:
        ref_dir / 'star_index' / 'SAindex',
    params:
        index_dir = ref_dir / 'star_index',
        unnorm_dir = project_dir / 'unnorm', # TODO: find better rule to make this
        overhang = read_length - 1,
    resources:
        mem_mb = 60000,
        cpus = threads,
        walltime = 4,
    shell:
        """
        mkdir -p {params.index_dir}
        mkdir -p {params.unnorm_dir}
        STAR --runMode genomeGenerate \
            --outTmpDir {params.index_dir}/tmp \
            --outFileNamePrefix {params.index_dir}/STAR_ \
            --genomeDir {params.index_dir} \
            --genomeFastaFiles {input.ref_genome} \
            --sjdbGTFfile {input.ref_anno} \
            --sjdbOverhang {params.overhang} \
            --runThreadN {resources.cpus}
        """

def fastqs(sample_id: str) -> list:
    """Get the list of FASTQ file paths for a sample.

    If paired_end is True, return a list of two lists of corresponding paths.
    """
    paths = [[], []] if paired_end else []
    with open(fastq_map, 'r') as f:
        for line in f.read().splitlines():
            if paired_end:
                fastq1, fastq2, s_id = line.split('\t')
                if s_id == sample_id:
                    paths[0].append(str(fastq_dir / fastq1))
                    paths[1].append(str(fastq_dir / fastq2))
            else:
                fastq, s_id = line.split('\t')
                if s_id == sample_id:
                    paths.append(str(fastq_dir / fastq))
    return paths

def fastq_input(wildcards):
    """Get the FASTQ file/s for a sample.

    This is the input list for the star_align rule, so if paired_end is True,
    concatenate into one list of files.
    """
    files = fastqs(wildcards.sample_id)
    return files[0] + files[1] if paired_end else files

def fastq_param(wildcards, input):
    """Get a string listing the fastq input files, with two lists if paired_end.
    
    This is the string supplied directly to the STAR command.
    """
    if paired_end:
        files = fastqs(wildcards.sample_id)
        return ' '.join([','.join(files[0]), ','.join(files[1])])
    else:
        return ','.join(input.fastq)

rule star_align:
    """Align RNA-Seq reads for a sample using STAR."""
    input:
        fastq = fastq_input,
        index = ref_dir / 'star_index' / 'SAindex',
    output:
        bam1 = project_dir / 'bam' / '{sample_id}.Aligned.sortedByCoord.out.bam',
        bam2 = project_dir / 'bam' / '{sample_id}.Aligned.toTranscriptome.out.bam',
    params:
        fastq_list = fastq_param,
        index_dir = ref_dir / 'star_index',
        bam_dir = project_dir / 'bam',
        prefix = str(project_dir / 'bam' / '{sample_id}.'),
    resources:
        mem_mb = 60000,
        cpus = threads,
        walltime = 8,
    shell:
        """
        mkdir -p {params.bam_dir}
        STAR --runMode alignReads \
            --genomeDir {params.index_dir} \
            --readFilesIn {params.fastq_list} \
            --readFilesCommand zcat \
            --twopassMode Basic \
            --quantMode TranscriptomeSAM \
            --outSAMstrandField intronMotif \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix {params.prefix} \
            --runThreadN {resources.cpus}
        """

rule index_bam:
    """Index a BAM file."""
    input:
        project_dir / 'bam' / '{basename}.bam',
    output:
        project_dir / 'bam' / '{basename}.bam.bai',
    params:
        add_threads = threads - 1,
    resources:
        cpus = threads,
    shell:
        # It expects the number of *additional* threads to use beyond the first.
        """
        samtools index -@ {params.add_threads} {input}
        """
