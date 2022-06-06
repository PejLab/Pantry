"""Normalize phenotypes in a BED file for direct use in QTL mapping"""

import argparse
import numpy as np
import pandas as pd
import scipy.stats as stats


def normalize_quantiles(df):
    """
    Quantile normalization to the average empirical distribution
    Note: replicates behavior of R function normalize.quantiles from library("preprocessCore")
    Author: Francois Aguet

    Reference:
     [1] Bolstad et al., Bioinformatics 19(2), pp. 185-193, 2003

    Adapted from https://github.com/andrewdyates/quantile_normalize
    """
    M = df.values.copy()

    Q = M.argsort(axis=0)
    m,n = M.shape

    # compute quantile vector
    quantiles = np.zeros(m)
    for i in range(n):
        quantiles += M[Q[:,i],i]
    quantiles = quantiles / n

    for i in range(n):
        # Get equivalence classes; unique values == 0
        dupes = np.zeros(m, dtype=int)  # replaced deprecated np.int with int -D. Munro
        for j in range(m-1):
            if M[Q[j,i],i]==M[Q[j+1,i],i]:
                dupes[j+1] = dupes[j]+1

        # Replace column with quantile ranks
        M[Q[:,i],i] = quantiles

        # Average together equivalence classes
        j = m-1
        while j >= 0:
            if dupes[j] == 0:
                j -= 1
            else:
                idxs = Q[j-dupes[j]:j+1,i]
                M[idxs,i] = np.median(M[idxs,i])
                j -= 1 + dupes[j]
        assert j == -1

    return pd.DataFrame(M, index=df.index, columns=df.columns)


def inverse_normal_transform(M):
    """
    Transform rows to a standard normal distribution
    Author: Francois Aguet
    """
    R = stats.mstats.rankdata(M, axis=1)  # ties are averaged
    if isinstance(M, pd.DataFrame):
        Q = pd.DataFrame(stats.norm.ppf(R/(M.shape[1]+1)), index=M.index, columns=M.columns)
    else:
        Q = stats.norm.ppf(R/(M.shape[1]+1))
    return Q


parser = argparse.ArgumentParser(description="Normalize phenotypes in a BED file for direct use in QTL mapping")
parser.add_argument("-i", "--input", help="Input BED file", required=True)
parser.add_argument("-o", "--output", help="Output BED file", required=True)
# parser.add_argument("--upper-quartile", help="Prior to quantile transformation, scale values per sample such that the upper quartile value is the same across samples", action="store_true", required=True)
args = parser.parse_args()

df = pd.read_csv(args.input, sep="\t", dtype={"#chrom": str, "chromStart": int, "chromEnd": int, "name": str})
data = df.drop(["#chrom", "chromStart", "chromEnd", "name"], axis=1)

# Remove rows with no variation:
print(args.input, df.shape, data.shape)
no_var = data.apply(lambda x: x.nunique() == 1, axis=1)
df = df[~no_var]
data = data[~no_var]
print(args.input, df.shape, data.shape)
if no_var.sum() > 0:
    print(f"{args.output}: Removed {no_var.sum()} rows with no variation")

iqn = normalize_quantiles(data)
iqn = inverse_normal_transform(iqn)
# Replace data values in df with transformed values:
df.loc[:, data.columns] = iqn
df.to_csv(args.output, sep="\t", index=False, float_format="%g")
