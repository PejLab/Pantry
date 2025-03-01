# Template for running FUSION TWAS using the weights produced by Pantry Pheast
# This example is based on GWAS data from https://zenodo.org/records/3629742

import pandas as pd

modalities = ['alt_polyA', 'alt_TSS', 'expression', 'isoforms', 'splicing', 'stability']
trait_df = pd.read_csv('../data/gwas_metadata.txt', sep='\t')
traits = trait_df['Tag'].tolist()
gwas_n = {traits[i]: n for i, n in enumerate(trait_df['Sample_Size'].tolist())}
chroms = [f'chr{i + 1}' for i in range(22)]

wildcard_constraints:
    modality="[^.]+"

rule all:
    input:
        # expand('output/{modality}/{modality}.{trait}.tsv', modality=modalities, trait=traits),
        'output/twas_hits.tsv',

rule assoc_test_trait:
    """Run all chromosomes for a trait together to submit fewer jobs"""
    input:
        sumstats = '../data/sumstats/{trait}.sumstats',
        weights = 'weights/{modality}.pos',
        ref_ld_chr = expand('../data/LDREF_b38ids_chr/1000G.EUR.{chrom}.{ext}', chrom=chroms, ext=['bed', 'bim', 'fam']),
    output:
        expand('intermediate/{{modality}}/{{trait}}/{{modality}}.{{trait}}.{chrom}.tsv', chrom=chroms)
    params:
        output_dir = 'intermediate/{modality}/{trait}',
        weights_dir = 'weights',
        ref_ld_prefix = '../data/LDREF_b38ids_chr/1000G.EUR.',
        chroms = ' '.join(chroms),
        coloc_p = '5e-8', # Hits with p below this are tested for coloc, recommended to be TWAS threshold. A threshold of 5e-8 / 6 is recommended for 6 modalities, but this will run coloc for everything below 5e-8 in case the values are useful.
        gwas_n = lambda w: gwas_n[w.trait],
    group: 'twas'
    threads: 4
    resources:
        mem_mb = 32000,
    shell:
        """
        mkdir -p {params.output_dir}
        parallel -j {threads} \
            Rscript ../scripts/FUSION.assoc_test.R \
                --sumstats {input.sumstats} \
                --weights {input.weights} \
                --weights_dir {params.weights_dir} \
                --ref_ld_chr {params.ref_ld_prefix} \
                --chr {{}} \
                --coloc_P {params.coloc_p} \
                --GWASN {params.gwas_n} \
                --separate_human_MHC \
                --out {params.output_dir}/{wildcards.modality}.{wildcards.trait}.{{}}.tsv \
            ::: {params.chroms}
        """

rule assemble_assoc_tests:
    input:
        expand('intermediate/{{modality}}/{{trait}}/{{modality}}.{{trait}}.{chrom}.tsv', chrom=chroms)
    output:
        out = 'output/{modality}/{modality}.{trait}.tsv'
    params:
        output_dir = 'output/{modality}',
    group: 'twas'
    shell:
        """
        mkdir -p {params.output_dir}
        head -n1 {input[0]} | cut -f3- > {output.out}
        for fname in {input}; do
            tail -n+2 $fname | cut -f3- >> {output.out}
        done
        """

rule assemble_hits_all_traits_modalities:
    input:
        expand('output/{modality}/{modality}.{trait}.tsv', modality=modalities, trait=traits)
    output:
        out = 'output/twas_hits.tsv'
    params:
        traits = ' '.join(traits),
        modalities = ' '.join(modalities),
        n_modalities = len(modalities),
    shell:
        """
        head -n1 {input[0]} \
            | sed 's/^/TRAIT\tMODALITY\t/' \
            > {output.out}
        for trait in {params.traits}; do
            for modality in {params.modalities}; do
                tail -n+2 output/$modality/$modality.$trait.tsv \
                    | awk '$18 < 5e-8 / {params.n_modalities}' \
                    | sed "s/^/$trait\t$modality\t/" \
                    >> {output.out}
            done
        done
        """
