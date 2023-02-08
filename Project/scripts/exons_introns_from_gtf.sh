# Adapted from https://github.com/csglab/CRIES

# This is the source annotation file
# Note that the annotation file should be in Ensembl format (downloaded from Ensembl FTP)
# Specifically, the chromosome names should lack the "chr" prefix
# Also, the first tag in column 9 must be gene_id, and the third tag must be transcript_id

## Changes by Daniel Munro:
## - Removed grep filters for 'transcript_source "ensembl_havana"' to allow for
##   use of non-ensembl annotations, and because for some species like rat,
##   there are very few havana calls, and for other species there are none.

input=$1
output_exons=$2
output_introns=$3

# 1. Identify and write constitutive exons, i.e. those that appear in all Ensembl/Havana isoforms of a gene
cat $input \
    | awk -v FS='\t' '$3=="exon" { exonName=$1":"$4":"$5":"$7; split($9, fields, ";"); geneName=fields[1]; transcriptName=fields[3]; printf("%s\t%s\t%s\n",exonName,geneName,transcriptName); }' \
    | sort \
    | uniq \
    | awk -v FS='\t' '{ eCount[$1]++; tCount[$3]++; exonHost[$1]=$2; if(tCount[$3]==1) gCount[$2]++; } END { for(i in eCount) if(eCount[i]==gCount[exonHost[i]]) { split(i,fields,":"); printf("%s\t.\texon\t%s\t%s\t.\t%s\t.\t%s;\n",fields[1],fields[2],fields[3],fields[4],exonHost[i]); } }' \
    > $output_exons

# 1. Find all exons
# 2. Write the coordinate of intronic regions, i.e. the regions that separate two adjacent exons of the same gene
cat $input \
    | awk -v FS='\t' '$3=="exon" { exonName=$1":"$4":"$5":"$7; split($9, fields, ";"); geneName=fields[1]; transcriptName=fields[3]; printf("%s\t%s\t%s\n",exonName,geneName,transcriptName); }' \
    | sort \
    | uniq \
    | awk -v FS='\t' '{ eCount[$1]++; tCount[$3]++; exonHost[$1]=$2; if(tCount[$3]==1) gCount[$2]++; } END { for(i in eCount) { split(i,fields,":"); printf("%s\t.\texon\t%s\t%s\t.\t%s\t.\t%s;\n",fields[1],fields[2],fields[3],fields[4],exonHost[i]); } }' \
    | bedtools sort -i stdin \
    | awk -v FS='\t' '{ if( last_exon[$9]==1 && (last_exon_end[$9]+1)<($4-1) ) printf("%s\t%s\tintron\t%i\t%i\t%s\t%s\t%s\t%s\n",$1,$2,last_exon_end[$9]+1,$4-1,$6,$7,$8,$9); last_exon[$9]=1; last_exon_end[$9]=$5; }' \
    > $output_introns
