=================
TURNAP User Guide
=================

TURNAP (Tool for Uniform RNA Phenotyping) is a framework for generating molecular phenotypes from RNA-Seq. Existing tools have been developed to generate different transcriptomic phenotypes, and TURNAP allows you to run them with full flexibility in a convenient and organized process. The goal is to make molecular QTL mapping and other analyses as easy to perform for a multitude of transcriptomic phenotypes as for basic gene expression only.

A project is specified with a config file containing general parameters, e.g. paths to FASTQ and reference annotation file, single-end or paired-end, and strandedness. It also lists the phenotypes to generate.

Installation
============

Code
====

The commands and code to run the steps are given in `Snakemake <https://snakemake.github.io/>`_ files, located by default in ``steps/``. This package includes sensible defaults for a set of tools, but you can edit the commands and add new phenotypes.

The goal of a project is to generate phenotype tables for a set of samples. Each tool could involve any combination of:

- Setup like indexing or annotation filtering
- Processing each sample's BAM file
- Intermediate steps requiring data from all samples, e.g. clustering
- Final assembly of all results into a phenotype table

Each file in ``steps/`` contains all of these steps for one phenotype. Snakemake looks at the inputs and outputs for each rule and figures out the order and number of times to run each one.

Config file
===========

The config file is a file in YAML format specifying general parameters, input files, and phenotypes to generate. This is the default config:

.. code-block:: yaml

   # Raw RNA-Seq data
   fastq_map: data/geuvadis/fastq_map.txt
   fastq_dir: data/geuvadis
   paired_end: True
   read_length: 75

   # Reference files
   ref_genome: data/human_ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa
   ref_anno: data/human_ref/Homo_sapiens.GRCh38.106.gtf
   retro_anno: data/human_ref/retro.hg38.v1.nochr.gtf

   # Project
   project_dir: test
   ref_dir: test/reference
   threads: 4

   # Phenotypes
   phenotypes:
   expression:
      files:
      - expression.log2.bed
      - expression.tpm.bed
   splicing:
      files:
      - splicing.bed
   stability:
      files:
      - stability.bed
   retroelements:
      files:
      - retroelements.bed
   latent:
      files:
      - latent.bed

Input files
===========

FASTQ sample map
----------------

A file containing FASTQ file paths and the name of the sample they map to, separated by tabs, with no header.

Example for single-end reads::

   batch1/sampleA_run1.fastq.gz	sampleA
   batch1/sampleA_run2.fastq.gz	sampleA
   batch2/sampleA_run3.fastq.gz	sampleA
   batch1/sampleB_run1.fastq.gz	sampleB
   ...

Example for paired-end reads::

   batch1/sampleA_run1_1.fastq.gz	batch1/sampleA_run1_2.fastq.gz	sampleA
   batch1/sampleA_run2_1.fastq.gz	batch1/sampleA_run2_2.fastq.gz	sampleA
   batch2/sampleA_run3_1.fastq.gz	batch2/sampleA_run3_2.fastq.gz	sampleA
   batch1/sampleB_run1_1.fastq.gz	batch1/sampleB_run1_2.fastq.gz	sampleB
   ...

FASTQ files
-----------

FASTQ files must be compressed (decompressible with ``zcat``) unless you modify the ``STAR`` command to handle uncompressed files. Multiple FASTQ files (or paired-end file pairs) can be provided for each sample and all reads will be used.

Reference genome sequence
-------------------------

The FASTA file for a reference genome, e.g. ``Homo_sapiens.GRCh38.dna.primary_assembly.fa``.

Reference annotation file
-------------------------

A GTF file containing the gene, exon, and other annotations, compatible with the supplied reference genome. For example, ``Homo_sapiens.GRCh38.106.gtf``.

Running from the command line
=============================

Phenotypes are computed by running one subprocess per phenotype category per sample. Within the specified output directory, one directory will be created for each phenotype category, containing intermediate (often per-sample) files from each program.

These data will then be combined into one BED file per phenotype group so that QTL mapping can be run separately for each.

Output files
============

BED output files
----------------

The `BED <https://genome.ucsc.edu/FAQ/FAQformat.html#format1>`_ files include four annotation columns, ``chrom``, ``chromStart``, ``chromEnd``, and ``name``, followed by one column per sample containing the phenotype values. There may be one phenotype per gene (e.g. expression) or multiple (e.g. splicing), but in either case the coordinates indicate the transcription start site of the phenotype's gene. This ensures that the same cis-window variants are tested for all phenotypes of the same gene.

::

   #chrom	chromStart	chromEnd	name	HG00315	HG00106	NA18489
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

Running
=======

We suggest these steps to run TURNAP.

1. Run the included test data
-----------------------------

This will require installation of all programs used in the snakefiles. `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ is recommended for easy installation and management of all these programs. Once you think you have everything installed, try running on the included test data, which is small:

.. code-block:: shell

   snakemake -s TURNAP/Snakefile -c1

1. Run with your data
---------------------

Snakemake has features to handle many execution needs such as threads, computational resources, and automatic cluster job submission.


.. toctree::
   :maxdepth: 2
   :caption: Contents:


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
