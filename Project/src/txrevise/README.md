# _txrevise_
_txrevise_ R package provides utilites to pre-process Ensembl transcript annotations to quantify differences in transcript strucuture (alternative promoters, alternative splicing, alternative poly-adenylation) either between experimental conditions or genotypes (e.g. for transcript usage quntitative trait loci (tuQTL) mapping).

If you use txrevise in your research, please cite the following paper: [Alasoo, Kaur, et al. "Genetic effects on promoter usage are highly context-specific and contribute to complex traits." Elife 8 (2019): e41673](https://doi.org/10.7554/eLife.41673)


## Pre-computed event annotations
Running _txrevise_ on the latest version of Ensembl can be quite timeconsuming. Thus, to make it easier to get started, we have pre-computed alternatve transcription events in the GFF3 format for both GRCh37 and GRCh38 reference genomes.

We have constructed two types of events. In the main event files we have masked the alternative internal exons in the promoter and 3' end events (labeled upstream and downstream). This ensures that we can reliably distinguish promoter usage and 3' end usage from alternative splicing, which are likely to be driven by distinct molecular mechanisms. The main annotation files are:

- [Homo_sapiens.GRCh37.87.version_1.tar.gz](https://zenodo.org/record/1302499/files/Homo_sapiens.GRCh37.87.version_1.tar.gz)
- [Homo_sapiens.GRCh38.92.version_1.tar.gz](https://zenodo.org/record/1302499/files/Homo_sapiens.GRCh38.92.version_1.tar.gz)
- [Homo_sapiens.GRCh38.96.version_1.tar.gz](https://zenodo.org/record/3232932/files/Homo_sapiens.GRCh38.96.version_1.tar.gz)

We have also constructed an alternative set of "raw" annotation files where the alternative internal exons in promoter and 3' end events have not masked. This maximises QTL discovery, because we can now also detect splicing event near promoters and 3 ends, but it comes at the expense that we are no longer able to distinguish between different molecular mechanisms. The raw annotation files are:

-   [Homo_sapiens.GRCh37.87.raw_events.version_1.tar.gz](https://zenodo.org/record/1302499/files/Homo_sapiens.GRCh37.87.raw_events.version_1.tar.gz)
-   [Homo_sapiens.GRCh38.92.raw_events.version_1.tar.gz](https://zenodo.org/record/1302499/files/Homo_sapiens.GRCh38.92.raw_events.version_1.tar.gz)


## Constructing transcription events

### Automatically
We have written a [Snakefile](https://github.com/kauralasoo/txrevise/blob/master/scripts/Snakefile) for [snakemake](http://snakemake.readthedocs.io/en/stable/) that
 - automates all of the steps
 - can be easily [executed in parallel](http://snakemake.readthedocs.io/en/latest/executable.html) on multiple cores or on a compute cluster if desired
 - can use a [Singularity](https://sylabs.io/guides/3.6/user-guide/) container, preventing R package version conflicts and ensuring reproducible results
Thus, using the Snakefile is the recommended method. Building the Singlarity container is described [here](https://github.com/kauralasoo/txrevise/blob/master/build_container.sh) and running the snakemake pipeline [here](https://github.com/kauralasoo/txrevise/blob/master/scripts/run_snakemake.sh).

### Manually
This section contains step-by-step instructions for manually constructing transcriptional events based on Ensembl transcript annotations.

#### Dependencies
Make sure that you have R 3.5 installed together with the following packages:

 - [optparse](https://cran.r-project.org/package=optparse)
 - [dplyr](https://cran.r-project.org/package=dplyr)
 - [purrr](https://cran.r-project.org/package=purrr)
 - [data.table](https://cran.r-project.org/package=data.table)
 - [ggplot2](https://cran.r-project.org/package=ggplot2)
 - [devtools](https://cran.r-project.org/package=devtools)
 - [rtracklayer](https://bioconductor.org/packages/rtracklayer/)
 - [GenomicFeatures](https://bioconductor.org/packages/GenomicFeatures/)
 - txrevise

You can install _txrevise_ directly from GitHub using the following command:

	library("devtools")
	install_github("kauralasoo/txrevise")

#### Step 1: Download the GTF file
First, you need to download the GTF file from the Ensembl website. For example, if you want to use Ensembl version 96, you should download the following file:

	wget ftp://ftp.ensembl.org/pub/release-96/gtf/homo_sapiens/Homo_sapiens.GRCh38.96.gtf.gz

#### Step 2: Extract tanscript tags from the GTF file
Ensembl GTF file contains a tags field marking protein coding transcript that are truncated at 5' or 3' ends. Txrevise uses this information to extend truncated transcripts. Unfortunately, the `import` function from rtracklayer does not import the tags field correctly. Thus, we need to first extract the transcript tags manually using a custom Python 3 script that comes with txrevise.

	python scripts/extractTranscriptTags.py --gtf Homo_sapiens.GRCh38.96.gtf.gz > Homo_sapiens.GRCh38.96.transcript_tags.txt

#### Step 3: Prepare transcript annotations for event construction
Next, we need to convert the transcript annotations (in GTF format) to a binary representation that can be efficiently used by txrevise.

	Rscript scripts/prepareAnnotations.R --gtf Homo_sapiens.GRCh38.96.gtf.gz --tags Homo_sapiens.GRCh38.96.transcript_tags.txt --out Homo_sapiens.GRCh38.96.txrevise_annotations.rds

#### Step 4: Construct transcription events
Finally, we can use the `constructEvents.R` script to construct alternative transcription events. Since the implementation of the various algorithms in txrevise have not been optimised for efficiency, processing the full Ensembl GTF file can take several days. However, different genes can trivially processed in parallel. To simplify parallel execution, `constructEvents.R` has the `--batch` option that takes as an input two integers separated by a space. For example, '1 200' would mean that all genes are split into 200 batches and and only genes in the first are be processed. To process all genes, you simply need to iterate from '1 200' to '200 200'.

	Rscript scripts/constructEvents.R --annot Homo_sapiens.GRCh38.96.txrevise_annotations.rds --batch '1 2000' --out txrevise_events  --fill TRUE

#### Step 5: Merge output files
Each run of `constructEvents.R` produces up to six output files: alternative promoter, internal exon and 3' end events (labeled as upstream, contained and downstream) for two possible sets of shared exons (grp_1 and grp_2). See the [vignette](http://htmlpreview.github.io/?https://github.com/kauralasoo/txrevise/blob/master/inst/doc/construct_events.html) for more details on how the events are constructed.

To be able to use these events with transcript quantification software such as [Salmon](http://salmon.readthedocs.io/en/latest/) or [kallisto](https://pachterlab.github.io/kallisto/), we first need to merge all files from different batches. The `grep -v "^#"` command removes the comment lines from the beginning of each individual GFF3 file.

	cat txrevise_events/txrevise.grp_1.upstream.* | grep -v "^#" > txrevise.grp_1.upstream.gff3
	cat txrevise_events/txrevise.grp_1.contained.* | grep -v "^#" > txrevise.grp_1.contained.gff3
	cat txrevise_events/txrevise.grp_1.downstream.* | grep -v "^#" > txrevise.grp_1.downstream.gff3

	cat txrevise_events/txrevise.grp_2.upstream.* | grep -v "^#" > txrevise.grp_2.upstream.gff3
	cat txrevise_events/txrevise.grp_2.contained.* | grep -v "^#" > txrevise.grp_2.contained.gff3
	cat txrevise_events/txrevise.grp_2.downstream.* | grep -v "^#" > txrevise.grp_2.downstream.gff3

#### Step 6: Quantify event expression
The GFF3 files constructed by _txrevise_ can be used with any transcript quantification algorithm (e.g. [Salmon](http://salmon.readthedocs.io/en/latest/) or [kallisto](https://pachterlab.github.io/kallisto/)). **Importantly, since the promoter, internal exon and 3' end events constructed by txrevise share some of the same exons, they should always quantified independently**. For example, if you have six GFF3 files from txrevise (txrevise.grp_1.upstream.gff3, txrevise.grp_1.contained.gff3, txrevise.grp_1.downstream.gff3, txrevise.grp_2.upstream.gff3, txrevise.grp_2.contained.gff3, txrevise.grp_2.downstream.gff3), you should also run the quantification algorithm six times with each of the GFF3 file independently.

## Extracting event sequences for quantification
Many transcript expression quantification tools (e.g. Salmon or kallisto) do not directly work with transcript annotations in GFF3 format and require the transcript sequences in FASTA format instead. The simplest way to convert GFF3 annotations into transcript (or event) sequences is to use the `gffread` tool from the [cufflinks](http://cole-trapnell-lab.github.io/cufflinks/) package:

	gffread -w <output.fa> -g <reference_genome.fa> <input.gff3>
