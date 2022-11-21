# These are the commands used to generate a small test dataset, which consists
# of 445 Geuvadis samples, subsetted to chr1 genotypes and phenotypes for the
# first 100 genes per phenotype category.

phenodir=$1

#######################
## Subset phenotypes ##
#######################

mkdir -p input/phenotypes

for pheno in expression retroelements stability
do
    # echo $pheno
    zcat $phenodir/$pheno.bed.gz | grep -P '^(#chr)|(1\t)' | head -101 | bgzip -c > input/phenotypes/$pheno.bed.gz
    tabix -p bed input/phenotypes/$pheno.bed.gz
done

for pheno in alt_polyA alt_TSS isoforms latent splicing
do
    # echo $pheno
    python3 scripts/test_subset_grouped_pheno.py \
        --in-bed $phenodir/$pheno.bed.gz \
        --in-groups $phenodir/$pheno.phenotype_groups.txt \
        --out-bed input/phenotypes/$pheno.bed \
        --out-groups input/phenotypes/$pheno.phenotype_groups.txt \
        --chrom 1 \
        --n-groups 100
    bgzip input/phenotypes/$pheno.bed
    tabix -p bed input/phenotypes/$pheno.bed.gz
done

##############################
## Subset genotypes to chr1 ##
##############################

## Keep all samples to allow for LD pruning for covariates and to detect QTLs
geno="../../Geuvadis/input/genotype/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup.nochr"
genonew="input/GEUVADIS.445_samples.GRCh38.chr1"
mkdir -p ../Pheast/input
plink2 --bfile $geno \
    --chr 1 \
    --make-bed \
    --out $genonew
