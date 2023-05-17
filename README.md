# Pantry

Pantry (PAN-TRanscriptome phenotYping) is a framework for generating molecular phenotypes from RNA-Seq. Existing tools have been developed to generate different transcriptomic phenotypes, and Pantry allows you to run them with full flexibility in a convenient and organized process. The goal is to make molecular QTL mapping and other analyses as easy to perform for a multitude of transcriptomic phenotypes as for basic gene expression only.

A project is specified with a config file containing general parameters, e.g. paths to FASTQ and reference annotation file, single-end or paired-end, and strandedness. It also lists the phenotypes to generate.

## Installation

It is recommended to install Pantry by cloning or downloading the repository, since you will then need to copy the project template directory for your own analysis. If using a conda environment, hold off on installing the Pantry Python package until the conda environment has been created and activated.

```sh
git clone https://github.com/daniel-munro/Pantry.git
```

## Code

The `Project/` directory is a template for a project. To run Pantry on a dataset, copy the contents of `Project/` to a new directory that you can write to. Your project directory now contains all the data processing code, so that if you modify it or add custom phenotypes, you have a record of it that stays with the results and is not automatically changed by package updates or other projects.

The commands and code to run the steps are given in [Snakemake](https://snakemake.github.io/) files, located by default in `Project/steps/`. This package includes sensible defaults for a set of tools, but you can edit the commands and add new phenotypes.

The goal of a project is to generate phenotype tables for a set of samples. Each tool could involve any combination of:

- Setup like indexing or annotation filtering
- Processing each sample's BAM file
- Intermediate steps requiring data from all samples, e.g. clustering
- Final assembly of all results into a phenotype table

Each file in `steps/` contains all of these steps for one phenotype category, or multiple categories produced in a similar way. Snakemake looks at the inputs and outputs for each rule and figures out the order and number of times to run each one.

## Config file

The config file is a file in YAML format specifying general parameters, input files, and phenotypes to generate. The project template directory includes a config, which is for a small example dataset and a small subset of the human reference:

```yaml
## Raw RNA-Seq data
paired_end: True
read_length: 75
fastq_dir: input/fastq
fastq_map: input/fastq_map.txt
...
```

## Input files

Some phenotypes are extracted from raw sequences (e.g. using kallisto), while others are extracted from aligned reads. Pantry input sequences must be in FASTQ format. If you only have BAM files, you can convert them to FASTQ using `samtools`, e.g.:

```shell
samtools fastq -1 pair1.fq -2 pair2.fq -0 /dev/null -s /dev/null -n in.bam
```

Pantry will generate its own BAM files from the FASTQ files using the reference files. These new BAM files will not contain the sequences themselves so they will be smaller than typical BAM files. If you are sure your existing BAM files are compatible with the reference files, you can symlink to them in the intermediate directory, named as expected by the pipeline, and update the timestamps to avoid re-running the alignment step. Note that the BAM files output by STAR are deleted by default after being converted to shrunken BAM files, so it's recommended to use symlinks instead of the actual files, and name the symlinks either as the STAR output BAM files or the final, shrunken BAM files.

### FASTQ sample map

A file containing FASTQ file paths and the name of the sample they map to, separated by tabs, with no header. The paths should be relative to the `fastq_dir` specified in the config file.

Example for single-end reads:

```
batch1/sampleA_run1.fastq.gz	sampleA
batch1/sampleA_run2.fastq.gz	sampleA
batch2/sampleA_run3.fastq.gz	sampleA
batch1/sampleB_run1.fastq.gz	sampleB
...
```

Example for paired-end reads:

```
batch1/sampleA_run1_1.fastq.gz	batch1/sampleA_run1_2.fastq.gz	sampleA
batch1/sampleA_run2_1.fastq.gz	batch1/sampleA_run2_2.fastq.gz	sampleA
batch2/sampleA_run3_1.fastq.gz	batch2/sampleA_run3_2.fastq.gz	sampleA
batch1/sampleB_run1_1.fastq.gz	batch1/sampleB_run1_2.fastq.gz	sampleB
...
```

### FASTQ files

FASTQ files must be compressed (decompressible with `zcat`) unless you modify the `STAR` command to handle uncompressed files. Multiple FASTQ files (or paired-end file pairs) can be provided for each sample and the reads will be combined into one BAM file per sample.

### Reference genome and annotations

- The FASTA file for a reference genome, e.g. `Homo_sapiens.GRCh38.dna.primary_assembly.fa`.
- A GTF file containing the gene, exon, and other annotations, compatible with the supplied reference genome. For example, `Homo_sapiens.GRCh38.106.gtf`. Any desired quality filtering should be done beforehand, e.g. including only transcripts validated by both Ensembl and HAVANA.

## Running from the command line

Phenotypes are usually computed by running one subprocess per phenotype category per sample. Within the specified output directory, one directory will be created for each phenotype category, containing intermediate (often per-sample) files from each program. These data will then be combined into one BED file per phenotype group so that QTL mapping can be run separately for each.

All this is done using Snakemake, so general guides to using Snakemake can be found online to learn its features. For example, you can specify a profile that determines how steps get run, and is different from the project config file described above. Here is an example profile config for use on a computing cluster with slurm scheduling:

`~/.config/snakemake/slurm/config.yaml`:

```yaml
use-conda: true
cluster: "sbatch -t {resources.walltime}:00:00 --mem={resources.mem_mb} -c {threads} {resources.partition} --mail-type=FAIL --mail-user=myemail@address.edu"
default-resources: [walltime=1, mem_mb=4000, partition=""]
# partition should be e.g. "--partition=gpu"
```

Resources are specified within some of the snakemake rules, which are plugged into this command and automatically submitted as cluster jobs.

## Output files

### BED output files

The [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) files include four annotation columns, `chr`, `start`, `end`, and `phenotype_id`, followed by one column per sample containing the phenotype values. There may be one phenotype per gene (e.g. expression) or multiple (e.g. splicing), but in either case the coordinates indicate the transcription start site of the phenotype's gene. This ensures that the same cis-window variants are tested for all phenotypes of the same gene.

```
#chr	start	end	phenotype_id	HG00315	HG00106	NA18489
1	29569	29570	ENSG00000227232:1:17055:17233:clu_1_-	0.671533	0.654321	0.673716
1	29569	29570	ENSG00000227232:1:17055:17606:clu_1_-	0.0340633	0.037037	0.0241692
1	29569	29570	ENSG00000227232:1:17368:17606:clu_1_-	0.294404	0.308642	0.302115
1	778668	778669	ENSG00000228327:1:729804:729898:clu_2_-	0.275362	0.203125	0.247748
1	778668	778669	ENSG00000228327:1:729804:733307:clu_2_-	0.362319	0.265625	0.166667
1	778668	778669	ENSG00000228327:1:729804:736713:clu_2_-	0	0.234375	0.0765766
1	778668	778669	ENSG00000228327:1:733213:733307:clu_2_-	0	0	0.238739
1	778668	778669	ENSG00000228327:1:736619:736713:clu_2_-	0.224638	0.296875	0.175676
1	778668	778669	ENSG00000228327:1:736619:740129:clu_2_-	0.137681	0	0.0945946
...
```

### Phenotype groups

For phenotypes categories in which multiple phenotypes are produced per gene (e.g., splice junctions), a file is produced that specifies which gene each phenotype belongs to. This is used by tensorQTL for grouped testing:

```
ENSG00000227232:1:17055:17233:clu_4_-	ENSG00000227232
ENSG00000227232:1:17055:17606:clu_4_-	ENSG00000227232
...
ENSG00000227232:1:18379:24738:clu_5_-	ENSG00000227232
ENSG00000279457:1:187577:187755:clu_6_-	ENSG00000279457
ENSG00000279457:1:187577:188130:clu_6_-	ENSG00000279457
...
```

## Running

We suggest these steps to run Pantry.

### 1. Run the included test data

This will require installation of all programs used in the snakefiles. [Miniconda](https://docs.conda.io/en/latest/miniconda.html) is recommended for easy installation and management of all these programs. A conda environment specification is provided in `conda_env.yaml`:

```sh
conda env create -n pantry --file conda_env.yaml
conda activate pantry
pip3 install -e Pantry
```

Once you think you have everything installed, try running on the included test data, which is small:

```sh
cd Project
snakemake -j1
```

### 2. Run with your data

Snakemake has features to handle many execution needs such as threads, computational resources, and automatic cluster job submission.

## Pheast

Pantry includes another template directory, `Pheast`, which can run downstream analyses, such as mapping QTLs and generating TWAS models, on all of the phenotypes generated by Pantry. At this stage, details and files for raw sequencing data, reference data, and phenotyping tools are no longer needed, while genotypes and multiple downstream analyses are now employed. The Pheast stage is similar in structure to the phenotyping stage, with a config file, snakefiles, and scripts. The phenotyping outputs (BED files and phenotype groups) can be used in place or transfered to a different computer to run Pheast.
