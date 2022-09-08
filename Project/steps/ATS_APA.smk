localrules:
    ATS_APA_pheno_groups,

TXREVISE_N_BATCHES = 10

rule txrevise_extract_tags:
    """Extract transcript tags from the GTF file
    
    txrevise rules are adapted from https://github.com/kauralasoo/txrevise/blob/master/scripts/Snakefile
    """
    input:
        ref_anno,
    output:
        tags = ref_dir / 'txrevise' / 'transcript_tags.txt',
    shell:
        """
        python3 src/txrevise/extractTranscriptTags.py \
            --gtf <(gzip -c {input}) \
            > {output.tags}
        """

rule txrevise_prepare_annotations:
    """Prepare annotations"""
    input:
        gtf = ref_anno,
        tags = ref_dir / 'txrevise' / 'transcript_tags.txt',
    output:
        ref_dir / 'txrevise' / 'txrevise_annotations.rds',
    shell:
        """
        Rscript src/txrevise/prepareAnnotations.R \
            --gtf {input.gtf} \
            --tags {input.tags} \
            --out {output}
        """

rule txrevise_construct_events:
    """Construct events for annotations"""
    input:
        ref_dir / 'txrevise' / 'txrevise_annotations.rds',
    output:
        ref_dir / 'txrevise' / 'batch' / 'txrevise.grp_1.upstream.{batch}_{n_batches}.gff3',
        ref_dir / 'txrevise' / 'batch' / 'txrevise.grp_2.upstream.{batch}_{n_batches}.gff3',
        ref_dir / 'txrevise' / 'batch' / 'txrevise.grp_1.contained.{batch}_{n_batches}.gff3',
        ref_dir / 'txrevise' / 'batch' / 'txrevise.grp_2.contained.{batch}_{n_batches}.gff3',
        ref_dir / 'txrevise' / 'batch' / 'txrevise.grp_1.downstream.{batch}_{n_batches}.gff3',
        ref_dir / 'txrevise' / 'batch' / 'txrevise.grp_2.downstream.{batch}_{n_batches}.gff3',
    params:
        batch_str = "'{batch} {n_batches}'",
        outdir = ref_dir / 'txrevise' / 'batch',
    shell:
        """
        Rscript src/txrevise/constructEvents.R \
            --annot {input} \
            --batch {params.batch_str} \
            --out {params.outdir}
        """

rule txrevise_merge_gff_files:
    """Merge txrevise output files"""
    input:
        gff = expand(str(ref_dir / 'txrevise' / 'batch' / 'txrevise.{{group}}.{{position}}.{batch}_{n_batches}.gff3'),
            batch = [i for i in range(1, TXREVISE_N_BATCHES + 1)],
            n_batches = TXREVISE_N_BATCHES),
    output:
        gff = ref_dir / 'txrevise' / 'txrevise.{group}.{position}.gff3',
    shell:
        'cat {input.gff} | grep -v "^#" > {output.gff}'

rule txrevise_make_one:
    """Iterate over groups and positions"""
    input:
        gff = expand(str(ref_dir / 'txrevise' / 'txrevise.{group}.{position}.gff3'),
            group = ['grp_1', 'grp_2'],
            position = ['upstream', 'contained', 'downstream']),
        # metadata = 'processed/{annotation}_{kind}/txrevise_{kind}_phenotype_metadata.tsv.gz'
    output:
        ref_dir / 'txrevise' / 'completed.txt',
    shell:
        "echo 'Done!' > {output}"

rule ATS_APA_transcript_seqs:
    """Generate transcript sequences from txrevise annotations for kallisto."""
    input:
        ref_genome = ref_genome,
        gff3 = ref_dir / 'txrevise' / 'txrevise.{type}.gff3',
    output:
        ref_dir / 'txrevise' / 'txrevise.{type}.fa',
    shell:
        """
        gffread \
            -w {output} \
            -g {input.ref_genome} \
            {input.gff3}
        """

rule ATS_APA_kallisto_index:
    """Generate the index for kallisto using annotations from txrevise."""
    input:
        ref_dir / 'txrevise' / 'txrevise.{type}.fa',
    output:
        ref_dir / 'txrevise' / 'txrevise.{type}.idx',
    shell:
        'kallisto index --index {output} {input}'

rule ATS_APA_kallisto:
    """Quantify expression using annotations from txrevise."""
    input:
        index = ref_dir / 'txrevise' / 'txrevise.{type}.idx',
        fastq = fastqs_kallisto,
    output:
        h5 = interm_dir / 'ATS_APA' / '{type}' / '{sample_id}' / 'abundance.h5',
        tsv = interm_dir / 'ATS_APA' / '{type}' / '{sample_id}' / 'abundance.tsv',
        json = interm_dir / 'ATS_APA' / '{type}' / '{sample_id}' / 'run_info.json',
    params:
        atsapa_type_dir = str(interm_dir / 'ATS_APA' / '{type}'),
        out_dir = str(interm_dir / 'ATS_APA' / '{type}' / '{sample_id}'),
        single_end_flag = '' if paired_end else '--single',
        # TODO add strandedness parameter
    resources:
        cpus = threads,
    shell:
        """
        mkdir -p {params.atsapa_type_dir}
        kallisto quant \
            --index {input.index} \
            --output-dir {params.out_dir} \
            {params.single_end_flag} \
            --threads {resources.cpus} \
            {input.fastq}
        """

rule assemble_ATS_APA_bed:
    """Assemble kallisto log2 and tpm outputs into ATS/APA BED files"""
    input:
        kallisto = expand(str(interm_dir / 'ATS_APA' / '{{type}}' / '{sample_id}' / 'abundance.tsv'), sample_id=samples),
        samples = samples_file,
        ref_anno = ref_anno,
    output:
        interm_dir / 'unnorm' / 'ATS_APA.{type}.bed',
    params:
        unnorm_dir = str(interm_dir / 'unnorm'),
        project_dir = project_dir,
        atsapa_type_dir = str(interm_dir / 'ATS_APA' / '{type}'),
    shell:
        """
        mkdir -p {params.unnorm_dir}
        python3 {params.project_dir}/src/assemble_bed.py \
            --type ATS_APA \
            --input-dir {params.atsapa_type_dir} \
            --samples {input.samples} \
            --ref_anno {input.ref_anno} \
            --output {output}
        """

rule normalize_ATS_APA:
    """Quantile-normalize values for QTL mapping"""
    input:
        bed = interm_dir / 'unnorm' / 'ATS_APA.{type}.bed',
        samples = samples_file,
    output:
        output_dir / 'ATS_APA.{type}.bed.gz',
    params:
        project_dir = project_dir,
        bed = str(output_dir / 'ATS_APA.{type}.bed'),
    shell:
        """
        python3 {params.project_dir}/src/normalize_phenotypes.py \
            --input {input.bed} \
            --samples {input.samples} \
            --output {params.bed}
        bgzip {params.bed}
        """

rule ATS_APA_pheno_groups:
    """Group phenotypes by gene for tensorQTL"""
    input:
        output_dir / 'ATS_APA.{type}.bed.gz',
    output:
        output_dir / 'ATS_APA.{type}.phenotype_groups.txt',
    shell:
        """
        zcat < {input} \
            | tail -n +2 \
            | cut -f4 \
            | awk '{{ g=$1; sub(/\..*$/, "", g); print $1 "\t" g }}' \
            > {output}
        """
