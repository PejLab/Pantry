# These are the commands used to generate a small test dataset, which consists
# of 445 Geuvadis samples, subsetted to chr1 genotypes and phenotypes for the
# first 200 genes per modality. Keeping all samples allows for LD pruning for
# covariates and sufficient power to detect some cis-QTLs.

set -euo pipefail

# phenoprefix=../../../gdml/laddr/pheast/geuvadis/input/phenotypes/gtex-residual-Geuvadis-cross_latent
# phenodir=../../../gdml/pantry/phenos/geuvadis
phenodir=../../../gdml/etc/rna-editing/geuvadis
geno="../../../gdml/data/geuvadis/genotype/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup.nochr"
genonew="input/GEUVADIS.445_samples.GRCh38.chr1_120Mb"

for modality in expression stability; do
mkdir -p input/phenotypes
# cp ../../../gdml/data/geuvadis/samples.txt input/
cp "${phenodir}/input/samples.txt" input/

## Modalities with non-grouped phenotypes (one per gene)
for modality in expression stability; do
    zcat "${phenodir}/output/${modality}.bed.gz" | head -201 | bgzip -c > "input/phenotypes/${modality}.bed.gz"
    tabix -p bed "input/phenotypes/${modality}.bed.gz"
done

## Modalies with grouped phenotypes
for modality in alt_polyA alt_TSS isoforms RNA_editing splicing; do
    python3 scripts/test_subset_grouped_pheno.py \
        --in-bed "${phenodir}/output/${modality}.bed.gz" \
        --in-groups "${phenodir}/output/${modality}.phenotype_groups.txt" \
        --out-bed "input/phenotypes/${modality}.bed" \
        --out-groups "input/phenotypes/${modality}.phenotype_groups.txt" \
        --chrom chr1 \
        --n-groups 200

    bgzip -f "input/phenotypes/${modality}.bed"
    tabix -p bed "input/phenotypes/${modality}.bed.gz"
done

plink2 --bfile "$geno" --chr 1 --make-bed --from-bp 1 --to-bp 120000000 --output-chr chrM --out "$genonew"

tar -czvf test_input_pheast.tar.gz input
