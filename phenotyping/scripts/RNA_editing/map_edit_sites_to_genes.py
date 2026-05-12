#!/usr/bin/env python3
"""Map RNA editing sites to one overlapping gene on the same strand."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from gtfparse import read_gtf


def parse_args():
    parser = argparse.ArgumentParser(
        description="Map BED edit sites to same-strand overlapping genes by nearest TSS"
    )
    parser.add_argument("edit_site_file", type=Path, help="BED file of edit sites")
    parser.add_argument("gtf_file", type=Path, help="Reference annotation GTF")
    parser.add_argument("output_file", type=Path, help="Output two-column TSV: site, gene_id")
    return parser.parse_args()


def load_edit_sites(edit_site_file: Path) -> pd.DataFrame:
    sites = pd.read_csv(edit_site_file, sep="\t", header=None, comment="#")
    if sites.shape[1] >= 6:
        sites = sites.iloc[:, [0, 1, 2, 5]]
    else:
        sites = sites.iloc[:, [0, 1, 2, 3]]
    sites.columns = ["chrom", "start0", "end", "strand"]
    sites["chrom"] = sites["chrom"].astype(str)
    sites["start"] = sites["start0"] + 1
    sites["site_id"] = sites["chrom"].astype(str) + "_" + sites["end"].astype(str)
    return sites[["chrom", "start", "end", "strand", "site_id"]]


def load_genes(gtf_file: Path) -> pd.DataFrame:
    genes = read_gtf(gtf_file)
    # Newer versions return a polars DF by default, but not all versions allow
    # return type to be specified, so this handles older and newer versions:
    if type(genes).__module__ == "polars.dataframe.frame":
        genes = genes.to_pandas()
    genes = genes.loc[genes["feature"] == "gene", :].copy()
    genes["seqname"] = genes["seqname"].astype(str)
    genes["tss"] = np.where(genes["strand"] == "+", genes["start"], genes["end"])
    genes = genes.rename(columns={"seqname": "chrom"})
    return genes[["chrom", "start", "end", "strand", "gene_id", "tss"]]


def map_group(sites: pd.DataFrame, genes: pd.DataFrame) -> list[tuple[str, int, str]]:
    sites = sites.sort_values(["start", "end", "site_id"])
    genes = genes.sort_values(["start", "end", "gene_id"])

    mappings = []
    active = []
    gene_iter = genes.itertuples(index=False)
    next_gene = next(gene_iter, None)

    for site in sites.itertuples(index=False):
        while next_gene is not None and next_gene.start <= site.end:
            active.append(next_gene)
            next_gene = next(gene_iter, None)

        active = [gene for gene in active if gene.end >= site.start]
        overlapping = [gene for gene in active if gene.start <= site.end]
        if not overlapping:
            continue

        best_gene = min(overlapping, key=lambda gene: (abs(site.start - gene.tss), gene.gene_id))
        mappings.append((site.site_id, site.start, best_gene.gene_id))

    return mappings


def main():
    args = parse_args()

    sites = load_edit_sites(args.edit_site_file)
    genes = load_genes(args.gtf_file)

    mappings = []
    for (chrom, strand), site_group in sites.groupby(["chrom", "strand"], sort=False):
        gene_group = genes.loc[(genes["chrom"] == chrom) & (genes["strand"] == strand), :]
        if gene_group.empty:
            continue
        mappings.extend(map_group(site_group, gene_group))

    if not mappings:
        raise ValueError(
            "No overlapping regions found between edit sites and genes. This likely "
            "indicates an input data problem, such as a mismatch between chromosome "
            "names in the input files."
        )

    output = pd.DataFrame(mappings, columns=["edit_site", "site_pos", "gene_id"])
    output = output.sort_values(["gene_id", "site_pos"])
    output[["edit_site", "gene_id"]].to_csv(
        args.output_file,
        sep="\t",
        header=False,
        index=False,
    )


if __name__ == "__main__":
    main()
