import argparse

import pandas as pd
from tensorqtl import genotypeio


def effective_maf_threshold(n_samples, min_allele_count, maf_threshold):
    if maf_threshold is not None:
        return maf_threshold
    if min_allele_count <= 0:
        return 0
    return min_allele_count / (2 * n_samples)


parser = argparse.ArgumentParser(description="Summarize genotype variants below the QTL allele-count/MAF threshold.")
parser.add_argument("geno_prefix", help="Prefix of plink (bed, bim, fam) genotype files.")
parser.add_argument("output", help="Output TSV file.")
parser.add_argument("--min_allele_count", type=int, default=10)
parser.add_argument("--maf_threshold", type=float)
args = parser.parse_args()

pr = genotypeio.PlinkReader(args.geno_prefix)
genotype_df = pr.load_genotypes()

called = genotype_df.where(genotype_df >= 0)
n_called = called.notna().sum(axis=1)
allele_number = 2 * n_called
alt_allele_count = called.sum(axis=1)
minor_allele_count = pd.concat([alt_allele_count, allele_number - alt_allele_count], axis=1).min(axis=1)
maf = minor_allele_count / allele_number
missingness = 1 - n_called / genotype_df.shape[1]

threshold = effective_maf_threshold(genotype_df.shape[1], args.min_allele_count, args.maf_threshold)
if args.maf_threshold is not None or args.min_allele_count <= 0:
    below = maf < threshold
else:
    below = (minor_allele_count < args.min_allele_count) | (maf < threshold)

summary = pd.Series({
    "total_variants": genotype_df.shape[0],
    "n_samples": genotype_df.shape[1],
    "min_allele_count_threshold": args.min_allele_count,
    "maf_threshold": threshold,
    "variants_below_threshold": below.sum(),
    "fraction_below_threshold": below.mean(),
    "variants_below_min_allele_count": (minor_allele_count < args.min_allele_count).sum(),
    "variants_below_maf_threshold": (maf < threshold).sum(),
    "min_minor_allele_count": minor_allele_count.min(),
    "median_minor_allele_count": minor_allele_count.median(),
    "min_maf": maf.min(),
    "median_maf": maf.median(),
    "max_missingness": missingness.max(),
    "median_missingness": missingness.median(),
})

summary.rename_axis("stat").reset_index(name="value").to_csv(
    args.output, sep="\t", index=False, float_format="%.6g"
)
