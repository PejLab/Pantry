suppressPackageStartupMessages(library(snpStats))
library(impute)

load_geno <- function(filename) {
    geno <- read.plink(filename)
    # Convert to 0,1,2 coding
    geno_mat <- as(geno$genotypes, "numeric")
    t(geno_mat)
}

get_PCs <- function(df, n_pcs) {
    if (sum(is.na(df)) > 0) {
        df <- impute.knn(df)$data # Expects samples as columns
    }
    df <- df[apply(df, 1, var) != 0, ]
    pca <- prcomp(t(df), center = TRUE, scale = TRUE)
    n_pcs <- min(n_pcs, ncol(pca$x) - 1)
    pcs <- round(pca$x[, 1:n_pcs, drop = FALSE], 6)
    pcs_df <- data.frame(ID = colnames(pcs), t(pcs))
    colnames(pcs_df) <- c("ID", rownames(pcs))  # Column names starting with digits get 'fixed' and must be changed back
    pcs_df
}

args <- commandArgs(trailingOnly = TRUE)
GENO_PREFIX <- args[1]
BED_FILE <- args[2]
N_GENO_PCS <- as.integer(args[3])
N_PHENO_PCS <- as.integer(args[4])
OUT_FILE <- args[5]

pheno <- read.delim(BED_FILE, check.names = FALSE, row.names = 4)[, -(1:3)]
if (ncol(pheno) < 2) stop("Computing covariate PCs requires more than 1 sample.")
pheno_pcs <- get_PCs(pheno, N_PHENO_PCS)
pheno_pcs$ID <- paste("pheno", pheno_pcs$ID, sep = "_")

geno <- load_geno(GENO_PREFIX)
geno <- geno[, colnames(pheno)]
geno_pcs <- get_PCs(geno, N_GENO_PCS)
geno_pcs$ID <- paste("geno", geno_pcs$ID, sep = "_")

stopifnot(identical(colnames(geno_pcs), colnames(pheno_pcs)))
covars <- rbind(geno_pcs, pheno_pcs)
write.table(covars, OUT_FILE, sep = "\t", quote = FALSE, row.names = FALSE)
