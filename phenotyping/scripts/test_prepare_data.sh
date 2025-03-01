# These are the commands used to generate the small test dataset, which
# consists of 16 randomly chosen Geuvadis samples, subsetted to the first
# 2 Mb of chr1, with reference files subsetted accordingly.

#################
## Subset data ##
#################

## Get sample subset
geno="../../Geuvadis/input/genotype/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup.nochr"
mkdir -p input
# shuf -n 16 $geno.fam | cut -f2 | sort > input/samples.txt
## For reproducibility, these were used:
## 'HG00108', 'HG00120', 'HG00185', 'HG00231', 'HG00240', 'HG00275', 'HG00328', 'HG00339',
## 'NA07357', 'NA11843', 'NA12872', 'NA19129', 'NA19175', 'NA20516', 'NA20518', 'NA20521'


python3 scripts/test_subset_fastq.py

######################
## Subset reference ##
######################

## Subset gene annotations
mkdir -p input/ref
gtf="../../Geuvadis/input/human_ref/Homo_sapiens.GRCh38.106.gtf"
gtfnew="input/ref/Homo_sapiens.GRCh38.106.chr1_0-2Mb.gtf"
head -5 "$gtf" > "$gtfnew"
awk '($1 == "1") && ($4 < 2000000) && ($5 < 2000000)' "$gtf" >> "$gtfnew"

## Subset genome
fa="../../Geuvadis/input/human_ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
fanew="input/ref/Homo_sapiens.GRCh38.dna.primary_assembly.chr1_0-2Mb.fa"
# 60 bases per line, so first 2 Mb is 33334 lines
head -33335 "$fa" > "$fanew"

####################
## Create archive ##
####################

tar -czvf test_input_phenotyping.tar.gz input
