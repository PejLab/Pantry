#!/bin/bash

####################### define executables
deseq="./src/_R/DESeq.R"
rambrandts="./src/_R/REMBRANDTS.R"

####################### identify the input arguments
jobid=$1
metadata=$2
refdir=$3
stringency=$4
fitmode=$5

if [ "$jobid" = "" ]; then
	echo -e "\nUsage: bash REMBRANDTS.sh <jobID> <metadata.txt> <inputDir> <stringency> <biasMode>\n"
	exit
fi

echo "Job ID: "$jobid
echo "Input metadata file: "$metadata
echo "Reference directory for HTSeq-Count files: "$refdir
echo "Stringency for filtering measurements: "$stringency
echo "The mode of identifying bias parameters: "$fitmode

if [ -e "$metadata" ]; then
	echo "Metadata file found."
else
	echo "ERROR: Metadata file was not found."
	exit
fi


####################### define temporary path
tmp_folder="./tmp/"$jobid
mkdir -p $tmp_folder

####################### run DESeq
Rscript $deseq $jobid $metadata $refdir


####################### define output path
out_folder="./out/"$jobid
mkdir -p $out_folder
mkdir -p $out_folder"/sampleScatterplots"

####################### run REMBRANDTS

Rscript $rambrandts $jobid $stringency $fitmode
