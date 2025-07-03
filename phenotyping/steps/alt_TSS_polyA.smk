localrules:
    alt_TSS_polyA_pheno_groups,

TXREVISE_N_BATCHES = 100

rule preprocess_gtf_for_txrevise:
    """Preprocess the GTF file for txrevise"""
    input:
        ref_anno,
    output:
        ref_dir / 'txrevise' / 'gene_annotations.gtf',
    shell:
        """
        python3 scripts/txrevise/preprocess_gtf.py \
            --input {input} \
            --output {output}
        """

rule txrevise_extract_tags:
    """Extract transcript tags from the GTF file
    
    txrevise rules are adapted from https://github.com/kauralasoo/txrevise/blob/master/scripts/Snakefile
    """
    input:
        ref_dir / 'txrevise' / 'gene_annotations.gtf',
    output:
        tags = ref_dir / 'txrevise' / 'transcript_tags.txt',
    shell:
        """
        python3 scripts/txrevise/extractTranscriptTags.py \
            --gtf <(gzip -c {input}) \
            > {output.tags}
        """

rule txrevise_prepare_annotations:
    """Prepare annotations"""
    input:
        gtf = ref_dir / 'txrevise' / 'gene_annotations.gtf',
        tags = ref_dir / 'txrevise' / 'transcript_tags.txt',
    output:
        ref_dir / 'txrevise' / 'txrevise_annotations.rds',
    shell:
        """
        Rscript scripts/txrevise/prepareAnnotations.R \
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
    resources:
        runtime = '16h',
    shell:
        """
        Rscript scripts/txrevise/constructEvents.R \
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

rule alt_TSS_polyA_transcript_seqs:
    """Generate transcript sequences from txrevise annotations for kallisto."""
    input:
        ref_genome = ref_genome,
        gff3 = ref_dir / 'txrevise' / 'txrevise.{group}.{position}.gff3',
    output:
        ref_dir / 'txrevise' / 'txrevise.{group}.{position}.fa.gz',
    shell:
        'gffread {input.gff3} -g {input.ref_genome} -w - | bgzip -c > {output}'

rule alt_TSS_polyA_kallisto_index:
    """Generate the index for kallisto using annotations from txrevise."""
    input:
        ref_dir / 'txrevise' / 'txrevise.{group}.{position}.fa.gz',
    output:
        ref_dir / 'txrevise' / 'txrevise.{group}.{position}.idx',
    resources:
        mem_mb = 32000,
    shell:
        'kallisto index --index {output} {input}'

rule alt_TSS_polyA_kallisto:
    """Quantify expression using annotations from txrevise."""
    input:
        index = ref_dir / 'txrevise' / 'txrevise.{group}.{position}.idx',
        fastq = fastq_inputs,
    output:
        h5 = interm_dir / 'alt_TSS_polyA' / '{group}.{position}' / '{sample_id}' / 'abundance.h5',
        tsv = interm_dir / 'alt_TSS_polyA' / '{group}.{position}' / '{sample_id}' / 'abundance.tsv',
        json = interm_dir / 'alt_TSS_polyA' / '{group}.{position}' / '{sample_id}' / 'run_info.json',
    params:
        alt_group_pos_dir = str(interm_dir / 'alt_TSS_polyA' / '{group}.{position}'),
        out_dir = str(interm_dir / 'alt_TSS_polyA' / '{group}.{position}' / '{sample_id}'),
        single_end_flag = kallisto_single_end_flag,
        # TODO add strandedness parameter
    threads: 16
    resources:
        runtime = '16h',
    shell:
        """
        mkdir -p {params.alt_group_pos_dir}
        kallisto quant \
            --index {input.index} \
            --output-dir {params.out_dir} \
            {params.single_end_flag} \
            --threads {threads} \
            {input.fastq}
        """

rule assemble_alt_TSS_polyA_bed:
    """Assemble kallisto log2 and tpm outputs into alt TSS/polyA BED files"""
    input:
        kallisto = lambda w: expand(
            str(interm_dir / 'alt_TSS_polyA' / 'grp_{grp}.{position}' / '{sample_id}' / 'abundance.tsv'),
            grp=[1, 2],
            position=dict(TSS='upstream', polyA='downstream')[w.type],
            sample_id=samples
        ),
        samples = samples_file,
        ref_anno = ref_anno,
    output:
        output_dir / 'unnorm' / 'alt_{type}.bed',
    params:
        unnorm_dir = str(output_dir / 'unnorm'),
        alt_group1_pos_dir = lambda w: str(interm_dir / 'alt_TSS_polyA' / f"grp_1.{dict(TSS='upstream', polyA='downstream')[w.type]}"),
        alt_group2_pos_dir = lambda w: str(interm_dir / 'alt_TSS_polyA' / f"grp_2.{dict(TSS='upstream', polyA='downstream')[w.type]}"),
    resources:
        mem_mb = 32000,
    shell:
        """
        mkdir -p {params.unnorm_dir}
        python3 scripts/assemble_bed.py \
            --type alt_TSS_polyA \
            --input-dir {params.alt_group1_pos_dir} \
            --input-dir2 {params.alt_group2_pos_dir} \
            --samples {input.samples} \
            --ref_anno {input.ref_anno} \
            --output {output}
        """

rule normalize_alt_TSS_polyA:
    """Quantile-normalize values for QTL mapping"""
    input:
        bed = output_dir / 'unnorm' / 'alt_{type}.bed',
        samples = samples_file,
    output:
        output_dir / 'alt_{type}.bed.gz',
    params:
        bed = str(output_dir / 'alt_{type}.bed'),
    shell:
        """
        mkdir -p {output_dir}
        python3 scripts/normalize_phenotypes.py \
            --input {input.bed} \
            --samples {input.samples} \
            --output {params.bed}
        bgzip {params.bed}
        """

rule alt_TSS_polyA_pheno_groups:
    """Group phenotypes by gene for tensorQTL"""
    input:
        output_dir / 'alt_{type}.bed.gz',
    output:
        output_dir / 'alt_{type}.phenotype_groups.txt',
    shell:
        """
        zcat < {input} \
            | tail -n +2 \
            | cut -f4 \
            | awk '{{ g=$1; sub(/__.*$/, "", g); print $1 "\\t" g }}' \
            > {output}
        """
