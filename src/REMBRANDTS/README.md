# REMBRANDTS

### REMoving Bias from Rna-seq ANalysis of Differential Transcript Stability

REMBRANDTS is a package for analysis of RNA-seq data across multiple samples in order to obtain unbiased estimates of differential mRNA stability. It uses DESeq to obtain estimates of differential pre-mRNA and mature mRNA abundance across samples, and then estimates a gene-specific bias function that is then subtracted from &Delta;exon–&Delta;intron to provide unbiased differential mRNA stability measures.

### Requirements

1. Unix-compatible OS
2. [R](https://www.r-project.org/) version 3.2.3 or later
3. R [gplots](https://cran.r-project.org/web/packages/gplots/index.html) library
4. R [DESeq](http://bioconductor.org/packages/release/bioc/html/DESeq.html) library

### Installation

REMBRANDTS is ready to use once the [Requirements](#requirements) are in place.

### Running REMBRANDTS


#### Input files

For running REMBRANDTS, you need the following files:

* **Metadata file**. The metadata file is a tab-delimited table with each row corresponding to either the exonic or intronic reads of one sample. The first line contains column headers, which must be exactly the same as follows:
	1. **Label**. This column contains the sample labels. For example, both the exonic and intronic read sets that come from sample #1 should be labeled as `Sample1`.
	2. **File**. This column contains the path to the HTSeq-count output files. The path can be either absolute, or relative to the "inputDir" provided as argument to REMBRANDTS.sh.
	3. **ReadType**. This column contains the read type for each of the HTSeq-Count files. Should be either `intronic` or `exonic`.
	4. **Batch**. In case batch-specific normalization is required, this column can be used to specify different batches as integer numbers 1 to *n*, where *n* is the total number of batches. If batch-specific normalization is not required, this number should always be 1.

  An example metadata table is shown below:
  
  | Label   | File                                    | ReadType | Batch |
  | ------- | --------------------------------------- | -------- | ----- |
  | Sample1 | ./htseq/Sample1.htseqCount.exonic.tab   | exonic   | 1     |
  | Sample1 | ./htseq/Sample1.htseqCount.intronic.tab | intronic | 1     |
  | Sample2 | ./htseq/Sample2.htseqCount.exonic.tab   | exonic   | 1     |
  | Sample2 | ./htseq/Sample2.htseqCount.intronic.tab | intronic | 1     |

* **Read count files**. In the read count files, each row corresponds to one gene, with the first column representing the gene ID and the second column representing the total number of reads mapped to that gene (either intronic or exonic reads). These files can be generated using HTSeq-count. A complete workflow for generating read count files that are compatible with REMBRENDTS is described [here](https://github.com/csglab/CRIES).

#### Usage

To run REMBRANDTS, use the following command:

```bash
bash ./REMBRANDTS.sh <jobID> <metadata.txt> <inputDir> <stringency> <biasMode>
```
* `jobID`: A job name that is used to create the output directory.
* `inputDir`: The directory relative to which the read count file paths are determined.
* `stringency`: The stringency for determining the cutoff for genes to be included in the analysis. REMBRANDTS determines a read count cutoff that results in an overall Pearson correlation (&rho;) between &Delta;exon and &Delta;intron equal to &rho;<sub>min</sub>+(&rho;<sub>max</sub>–&rho;<sub>min</sub>)×stringency.
* `biasMode`: Currently, only `linear` is accepted.

#### Output

REMBRANDTS creates the following output files in `./out/<jobID>/`

* `exonic.filtered.mx.txt`: The estimated log<sub>2</sub> of abundance of exonic fragments (&Delta;exon), relative to the average of all samples.
* `exonic.filtered.correl.heatmap.jpg`: The heatmap of Pearson similarities of samples with respect to the above estimates.
* `intronic.filtered.mx.txt`: The estimated log<sub>2</sub> of abundance of intronic fragments (&Delta;intron), relative to the average of all samples.
* `intronic.filtered.correl.heatmap.jpg`: The heatmap of Pearson similarities of samples with respect to the above estimates.
* `stability.filtered.mx.txt`: The unbiased estimates of differential mRNA stability (&Delta;exon–&Delta;intron–*bias*), relative to the average of all samples.
* `stability.filtered.correl.heatmap.jpg`: The heatmap of Pearson similarities of samples with respect to the above estimates.
* `scatterplot.jpg`: The scatterplot of &Delta;exon vs. &Delta;intron for all filtered genes in all samples.
* `sampleScatterplots/scatterplot.<Label>.jpg`: The scatterplot of &Delta;exon–&Delta;intron vs &Delta;intron for each sample, before and after removing the bias term.


#### Example

Three example datasets are provided in the `./examples/` folder. You can run REMBRANDTS on these examples using these commands:
```bash
bash ./REMBRANDTS.sh Human_tissue_stability ./examples/Tissues.SRP056969/table.txt ./examples/Tissues.SRP056969 0.99 linear
```
```bash
bash ./REMBRANDTS.sh AD_stability ./examples/AD.GSE53697/table.txt ./examples/AD.GSE53697 0.7 linear
```
```bash
bash ./REMBRANDTS.sh Mouse_mixed_stability ./examples/Mouse.PMID26098447/table.txt ./examples/Mouse.PMID26098447 0.99 linear
```
```bash
bash ./REMBRANDTS.sh Shen_2012_GSE29278_stability ./examples/Shen_2012_GSE29278_counts/table.txt ./examples/Shen_2012_GSE29278_counts 0.99 linear
```
```bash
bash ./REMBRANDTS.sh Furlow_2015_GSE45162_stability ./examples/Furlow_2015_GSE45162_counts/table.txt ./examples/Furlow_2015_GSE45162_counts 0.99 linear
```
These commands will replicate the stability estimates presented in [Alkallas et al. (Nat Commun, 2017)](https://www.nature.com/articles/s41467-017-00867-z).

## Citation
Alkallas R, Fish L, Goodarzi H, Najafabadi HS (2017). Inference of RNA decay rate from transcriptional profiling highlights the regulatory programs of Alzheimer's disease. Nat Commun [8:909](https://www.nature.com/articles/s41467-017-00867-z)
