REF_SEQ=$1
REF_ANNO=$2
OUT_CHR=$3
OUT_GTF=$4

## Get chromosome lengths from genome FASTA header lines
grep '^>' $REF_SEQ | sed 's/^>//' | tr ':' ' ' | tr ' ' '\t' | cut -f1,8 > $OUT_CHR

## The GTF entries must contain gene_name and transcript_id fields
## (and only transcript entries are used)
awk '$3=="transcript"' $REF_ANNO | grep gene_name | grep transcript_id > $OUT_GTF
