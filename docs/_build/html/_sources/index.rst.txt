TURNAP User Guide
=================

Tool for Unified RNA Phenotypes

Installation
------------



Input files
-----------

FASTQ sample map
^^^^^^^^^^^^^^^^

A file containing FASTQ file paths and the name of the sample they map to, separated by tabs, with no header. For example:

   batch1/sample1_run1.fastq.gz   sample1
   batch1/sample1_run2.fastq.gz   sample1
   batch2/sample1_run3.fastq.gz   sample1
   batch1/sample2_run1.fastq.gz   sample2
   ...

FASTQ files
^^^^^^^^^^^

[need to be zipped?]
Multiple FASTQ files can be provided for each sample and all reads will be used.
[how to specify paired-end files]

Running from the command line
-----------------------------

Phenotypes are computed by running one subprocess per phenotype category per sample. Within the specified output directory, one directory will be created for each phenotype category, containing per-sample output files from each program. For example, ``expression/`` will contain ...

These data will then be combined into one BED file per phenotype group so that QTL mapping can be run separately for each.


Output files
------------

BED output files
^^^^^^^^^^^^^^^^

The `BED <https://genome.ucsc.edu/FAQ/FAQformat.html#format1>`_ files include four annotation columns, ``chrom``, ``chromStart``, ``chromEnd``, and ``name``, followed by one column per sample containing the phenotype values. There may be one phenotype per gene (e.g. expression) or multiple (e.g. splicing), but in either case the coordinates indicate the transcription start site of the phenotype's gene. This ensures that the same cis-window variants are tested for all phenotypes of the same gene.

.. toctree::
   :maxdepth: 2
   :caption: Contents:


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
