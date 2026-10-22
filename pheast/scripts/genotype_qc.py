import argparse

import numpy as np
import pandas as pd
from pandas_plink import read_plink


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
parser.add_argument("--chunk_size", type=int, default=10000, help="Variants decoded per chunk (default: 10000).")
args = parser.parse_args()
if args.chunk_size <= 0:
    parser.error("--chunk_size must be positive")

# Use the lazy BED array directly: PlinkReader also builds large position
# dictionaries that this summary does not need. No genotype imputation is done.
bim, fam, bed = read_plink(args.geno_prefix, verbose=False)
n_variants, n_samples = bed.shape
del bim, fam
minor_allele_count = np.empty(n_variants, dtype=np.float64)
maf = np.empty(n_variants, dtype=np.float64)
missingness = np.empty(n_variants, dtype=np.float64)
for start in range(0, n_variants, args.chunk_size):
    stop = min(start + args.chunk_size, n_variants)
    # Single-threaded computation bounds concurrent decoder allocations.
    genotypes = bed[start:stop].compute(scheduler="single-threaded")
    n_called = np.count_nonzero(~np.isnan(genotypes), axis=1)
    allele_number = 2 * n_called
    allele_count = np.nansum(genotypes, axis=1, dtype=np.float64)
    # Minor counts are invariant to the allele flip performed by PlinkReader.
    minor_allele_count[start:stop] = np.minimum(allele_count, allele_number - allele_count)
    with np.errstate(divide="ignore", invalid="ignore"):
        maf[start:stop] = minor_allele_count[start:stop] / allele_number
    missingness[start:stop] = 1 - n_called / n_samples
    del genotypes
    print(f"Processed {stop:,}/{n_variants:,} variants", flush=True)

# Pandas reductions preserve the original skip-NaN behavior and exact medians.
minor_allele_count = pd.Series(minor_allele_count, copy=False)
maf = pd.Series(maf, copy=False)
missingness = pd.Series(missingness, copy=False)

threshold = effective_maf_threshold(n_samples, args.min_allele_count, args.maf_threshold)
if args.maf_threshold is not None or args.min_allele_count <= 0:
    below = maf < threshold
else:
    below = (minor_allele_count < args.min_allele_count) | (maf < threshold)

summary = pd.Series({
    "total_variants": n_variants,
    "n_samples": n_samples,
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
