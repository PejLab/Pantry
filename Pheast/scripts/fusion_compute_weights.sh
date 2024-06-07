#!/bin/sh
# Adapted from https://github.com/gusevlab/fusion_twas/blob/master/examples/GEUV.compute_weights.sh

GENO=$1
BED=$2
COVAR=$3
PHENO=$4
B_START=$5
B_END=$6
OUTDIR=$7

# I don't think set -e can be used because plink has error when no cis-window variants present, in which case this script should skip it and continue.
# set -e

echo $GENO $BED $COVAR $PHENO $B_START $B_END $OUTDIR

# PATH TO DIRECTORY CONTAINING LDREF DATA (FROM FUSION WEBSITE or https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2)
# LDREF="scripts/fusion_twas/LDREF"
# THIS IS USED TO RESTRICT INPUT SNPS TO REFERENCE IDS ONLY
mkdir -p $OUTDIR

# Setting family ID to 0, which is current plink default when coverting from VCF
# zcat < $BED | head -n1 | tr '\t' '\n' | tail -n+5 | awk '{ print 0,$1 }' > $OUTDIR/$PHENO.ID
# Instead of assuming family ID is 0, get it from $GENO.fam plink file:
zcat < $BED | head -n1 | tr '\t' '\n' | tail -n+5 \
	| awk 'NR==FNR{a[$2]=$0;next}{print a[$1]}' $GENO.fam - \
	| awk '{ print $1,$2 }' \
	> $OUTDIR/$PHENO.ID

NR="${B_START}_${B_END}"
mkdir -p $OUTDIR/tmp_$PHENO/$NR
mkdir -p $OUTDIR/hsq_$PHENO
# THIS IS DIRECTORY WHERE THE OUTPUT WILL GO:
mkdir -p $OUTDIR/$PHENO

if [ -f "$OUTDIR/tmp_$PHENO/$NR.hsq" ]; then
	rm $OUTDIR/tmp_$PHENO/$NR.hsq
fi
# Create batch hsq file even if empty to indicate the batch has been run
touch $OUTDIR/tmp_$PHENO/$NR.hsq

# Loop through each phenotype in the batch (and NR starts at 1)
zcat < $BED | awk -vs=$B_START -ve=$B_END 'NR >= s + 1 && NR <= e + 1' | while read PARAM; do

	# Get the gene positions +/- 500kb
	CHR=`echo $PARAM | awk '{ print $1 }'`
	P0=`echo $PARAM | awk '{ p = $3 - 0.5e6; if (p < 1) p = 1; print p }'`
	P1=`echo $PARAM | awk '{ print $3 + 0.5e6 }'`
	GNAME=`echo $PARAM | awk '{ print $4 }'`

	G_TMP="$OUTDIR/tmp_$PHENO/$NR/$GNAME"

	echo $GNAME $CHR $P0 $P1

	# Pull out the current phenotype
	echo $PARAM | tr ' ' '\n' | tail -n+5 | paste $OUTDIR/$PHENO.ID - > $G_TMP.pheno

	# Get the locus genotypes for all samples and set current gene expression as the phenotype
	# Omit filtering to LD reference SNPs. Instead, Pantry does optional pre-filtering.
		# --extract $LDREF/1000G.EUR.$CHR.bim \
		# --force-intersect
	plink --bfile $GENO \
		--pheno $G_TMP.pheno \
		--keep $G_TMP.pheno \
		--make-bed \
		--out $G_TMP \
		--chr $CHR \
		--from-bp $P0 \
		--to-bp $P1

	G_OUT="$OUTDIR/$PHENO/$GNAME"

	Rscript scripts/fusion_twas/FUSION.compute_weights.R \
		--bfile $G_TMP \
		--covar $COVAR \
		--tmp $G_TMP.tmp \
		--out $G_OUT \
		--verbose 0 \
		--save_hsq \
		--PATH_plink plink \
		--PATH_gcta ./scripts/fusion_twas/gcta_nr_robust \
		--PATH_gemma ./scripts/fusion_twas/gemma \
		--models blup,lasso,top1,enet

	# Append heritability output to hsq file
	if [ -f "$G_OUT.hsq" ]; then
		cat $G_OUT.hsq >> $OUTDIR/tmp_$PHENO/$NR.hsq
	fi

	# Clean-up just in case
	rm -f $G_OUT.hsq $G_TMP.tmp.*

	# Remove all intermediate files
	rm $G_TMP.*

done

# Wait until now to move hsq to final location so  that terminated jobs don't leave a partial output file
mv $OUTDIR/tmp_$PHENO/$NR.hsq $OUTDIR/hsq_$PHENO/$NR.hsq
