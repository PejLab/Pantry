rule calculate_heritability_chr:
    """Get heritability for phenotypes on one chromosome using plink and gcta64."""
    input:
        bed = project_dir / '{pheno}.bed.gz',
        geno = multiext('data/genotype/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup', '.bed', '.bim', '.fam'),
    output:
        project_dir / 'heritability' / '{pheno}_hsq.{chrom}.tsv',
    params:
        geno_prefix = 'data/genotype/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup',
        grm_dir = lambda w: project_dir / 'heritability' / f'grm_{w.chrom}',
        tmp_dir = lambda w: project_dir / 'heritability' / f'tmp_{w.chrom}',
    shell:
        """
        mkdir -p {params.grm_dir}
        mkdir -p {params.tmp_dir}
        python3 TURNAP/src/heritability.py \
            --bed {input.bed} \
            --geno {params.geno_prefix} \
            --chrom {wildcards.chrom} \
            --grm-dir {params.grm_dir} \
            --tmp-dir {params.tmp_dir} \
            --output {output}
        rm -r {params.grm_dir}
        rm -r {params.tmp_dir}
        """

CHROMS=[f'chr{i}' for i in range(1, 23)] + ['chrX']

rule combine_heritability_chr:
    """Combine per-chromosome heritability stats into one file"""
    input:
        expand(project_dir / 'heritability' / '{{pheno}}_hsq.{chrom}.tsv', chrom=CHROMS),
    output:
        project_dir / 'heritability' / '{pheno}_hsq.tsv',
    shell:
        """
        head -n 1 {input[0]} > {output}
        for FILE in {input}; do
            tail -n +2 $FILE >> {output}
        done
        """
