# These are the commands used to generate the small test dataset, which
# consists of 16 randomly chosen Geuvadis samples, subsetted to the first
# 2 Mb of chr1, with reference files subsetted accordingly.

#################
## Subset data ##
#################

## Get sample subset
geno="../../Geuvadis/input/genotype/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup.nochr"
mkdir -p input
shuf -n 16 $geno.fam | cut -f2 | sort > input/samples.txt

## Subset bam files
mkdir -p input/bam
while read f; do
    echo $f
    samtools view ../../Geuvadis/intermediate/bam/$f.Aligned.sortedByCoord.out.bam "1:1-2000000" -b > input/bam/$f.bam
done < input/samples.txt

## Make BAM map
awk '{print $1".bam\t"$1}' input/samples.txt > input/bam_map.txt

## Subset genotypes (Add 1 Mb to include the full cis-window for the test region genes)
## Keep all samples to allow for LD pruning for covariates
genonew="../Pheast/input/GEUVADIS.445_samples.GRCh38.chr1_0-2Mb"
mkdir -p ../Pheast/input
plink2 --bfile $geno \
    --chr 1 \
    --from-mb 0 \
    --to-mb 3 \
    --make-bed \
    --out $genonew

########################
## Subset annotations ##
########################

## Subset gene annotations
mkdir -p input/ref
gtf="../../Geuvadis/input/human_ref/Homo_sapiens.GRCh38.106.gtf"
gtfnew="input/ref/Homo_sapiens.GRCh38.106.chr1_0-2Mb.gtf"
head -5 "$gtf" > "$gtfnew"
awk '($1 == "1") && ($4 < 2000000) && ($5 < 2000000)' "$gtf" >> "$gtfnew"

## Subset retroelement annotations
retro="../../Geuvadis/input/human_ref/retro.hg38.v1.nochr.gtf"
retronew="input/ref/retro.hg38.v1.nochr.chr1_0-2Mb.gtf"
awk '($1 == "1") && ($4 < 2000000) && ($5 < 2000000)' "$retro" > "$retronew"

## Subset genome
fa="../../Geuvadis/input/human_ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
fanew="input/ref/Homo_sapiens.GRCh38.dna.primary_assembly.chr1_0-2Mb.fa"
# 60 bases per line, so first 2 Mb is 33334 lines
head -33335 "$fa" > "$fanew"
