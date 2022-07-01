# These are the commands used to generate the small test dataset, which
# consists of 10 randomly chosen Geuvadis samples, subsetted to the first
# Mb of chr1, with reference files subsetted accordingly.

## Subset fastq and make new fastq mapping file for them
python3 src/test_subset_fastq.py

## Subset gene annotations
mkdir -p test/ref
gtf="../data/human_ref/Homo_sapiens.GRCh38.106.gtf"
gtfnew="test/ref/Homo_sapiens.GRCh38.106.chr1_0-1Mb.gtf"
head -5 "$gtf" > "$gtfnew"
awk '($1 == "1") && ($4 < 1000000) && ($5 < 1000000)' "$gtf" >> "$gtfnew"

## Subset retroelement annotations
retro="../data/human_ref/retro.hg38.v1.nochr.gtf"
retronew="test/ref/retro.hg38.v1.nochr.chr1_0-1Mb.gtf"
awk '($1 == "1") && ($4 < 1000000) && ($5 < 1000000)' "$retro" > "$retronew"

## Subset genome
fa="../data/human_ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
fanew="test/ref/Homo_sapiens.GRCh38.dna.primary_assembly.chr1_0-1Mb.fa"
# 60 bases per line, so first Mb is in 16667 lines
head -16668 "$fa" > "$fanew"
