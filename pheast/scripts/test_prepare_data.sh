# These are the commands used to generate a small test dataset, which consists
# of 445 Geuvadis samples, subsetted to chr1 genotypes and phenotypes for the
# first 100 genes per modality.

phenodir=../../Geuvadis/output

#######################
## Subset phenotypes ##
#######################

mkdir -p input/phenotypes

cp ../../Geuvadis/input/samples.txt input/

## Test a modality with non-grouped phenotypes (one per gene)

zcat $phenodir/expression.bed.gz | grep -P '^(#chr)|(1\t)' | head -201 | bgzip -c > input/phenotypes/expression.bed.gz
tabix -p bed input/phenotypes/expression.bed.gz

## Test a modality with grouped phenotypes

python3 scripts/test_subset_grouped_pheno.py \
    --in-bed $phenodir/alt_polyA.bed.gz \
    --in-groups $phenodir/alt_polyA.phenotype_groups.txt \
    --out-bed input/phenotypes/alt_polyA.bed \
    --out-groups input/phenotypes/alt_polyA.phenotype_groups.txt \
    --chrom 1 \
    --n-groups 200
bgzip input/phenotypes/alt_polyA.bed
tabix -p bed input/phenotypes/alt_polyA.bed.gz

##############################
## Subset genotypes to chr1 ##
##############################

## Keep all samples to allow for LD pruning for covariates and to detect QTLs
geno="../../Geuvadis/input/genotype/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup.nochr"
genonew="input/GEUVADIS.445_samples.GRCh38.chr1"
plink2 --bfile $geno \
    --chr 1 \
    --make-bed \
    --out $genonew

####################
## Create archive ##
####################

tar -czvf test_input_pheast.tar.gz input
