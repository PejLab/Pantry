## Map edit sites to genes for grouped cis-QTL mapping

library(tidyverse)
library(GenomicRanges)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript map_edit_sites_to_genes.R <edit_site_file> <gtf_file> <output_file>")
}

edit_site_file <- args[1]
gtf_file <- args[2]
output_file <- args[3]

edit_ref <- read_table(edit_site_file, col_names = c("chrom", "start", "end", "strand"), col_types = "cii--c") |>
  mutate(site_id = str_glue("{chrom}_{end}"))

genes <- rtracklayer::readGFF(gtf_file) |>
  filter(type == "gene")

# (Add 1 to start position to convert BED 0-based to 1-based coordinates)
edit_ref_grange <- GRanges(
  seqnames = edit_ref$chrom,
  ranges = IRanges(start = edit_ref$start + 1, end = edit_ref$end),
  strand = edit_ref$strand
)

genes_grange <- GRanges(
  seqnames = genes$seqid,
  ranges = IRanges(start = genes$start, end = genes$end),
  strand = genes$strand
)

# Find strand-specific overlaps
hits <- findOverlaps(edit_ref_grange, genes_grange, ignore.strand = FALSE)

# Check if any overlaps were found
if (length(hits) == 0) {
  stop("No overlapping regions found between edit sites and genes. This likely indicates an input data problem, such as a mismatch between chromosome names in the input files.")
}

pairs <- tibble(
  edit_site = edit_ref$site_id[queryHits(hits)],
  gene_id = genes$gene_id[subjectHits(hits)]
)

write_tsv(pairs, output_file, col_names = FALSE)
