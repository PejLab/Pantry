BAM=$1
SAMPLE=$2
REF_SEQ=$3
CODE_DIR=$4
OUT_DIR=$5
ATS_APA_DIR="$OUT_DIR/ATS_APA"
MC_DIR="$CODE_DIR/src/mountainClimber/src"

# Get genome coverage
# Add sort step so chr order is 1, 10, etc. instead of 1, 2, etc. so that merge step called by mountainClimberTU.py works.
bedtools genomecov -trackline -bg -split \
    -ibam $BAM \
    | bedtools sort \
    > $ATS_APA_DIR/$SAMPLE.bedgraph

conda activate py2

# 
python2 $MC_DIR/get_junction_counts.py \
    -i $BAM \
    -s fr-unstrand \
    -o $ATS_APA_DIR/$SAMPLE.jxn.bed

# Call transcription units de novo
python2 $MC_DIR/mountainClimberTU.py \
    -b $ATS_APA_DIR/$SAMPLE.bedgraph \
    -j $ATS_APA_DIR/$SAMPLE.jxn.bed \
    -s 0 \
    -g $OUT_DIR/reference/chr_lengths.genome \
    -o $ATS_APA_DIR/$SAMPLE.tu.bed


python2 $MC_DIR/merge_tus.py \
    -i $ATS_APA_DIR/$SAMPLE.tu.bed \
    -s n \
    -g $OUT_DIR/reference/transcripts.gtf \
    -o $ATS_APA_DIR/tus_merged

# Identify change points in read coverage in each TU
python2 $MC_DIR/mountainClimberCP.py \
    -i $ATS_APA_DIR/$SAMPLE.bedgraph \
    -g $ATS_APA_DIR/tus_merged.annot.transcripts_singleGenes.bed \
    -j $ATS_APA_DIR/$SAMPLE.jxn.bed \
    -o $ATS_APA_DIR/$SAMPLE.cp.bed \
    -x $REF_SEQ
python2 $MC_DIR/mountainClimberRU.py \
    -i $ATS_APA_DIR/$SAMPLE.cp.bed \
    -o $ATS_APA_DIR/$SAMPLE.ru.bed
