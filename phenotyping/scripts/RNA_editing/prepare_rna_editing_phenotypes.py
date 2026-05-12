#!/usr/bin/env python3
"""Prepare RNA editing phenotypes from a covered site matrix."""

import argparse
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform


def parse_args():
    parser = argparse.ArgumentParser(
        description="Map RNA editing sites to genes and optionally cluster correlated sites"
    )
    parser.add_argument("--edit-matrix", type=Path, required=True,
                        help="Input site x sample matrix with num/den editing ratios")
    parser.add_argument("--site-to-gene", type=Path, required=True,
                        help="Two-column TSV mapping edit site IDs to gene IDs")
    parser.add_argument("--output-matrix", type=Path, required=True,
                        help="Output RNA editing phenotype x sample matrix")
    parser.add_argument("--output-site-map", type=Path, required=True,
                        help="Output TSV mapping source sites to output phenotype IDs")
    parser.add_argument("--cluster", action="store_true",
                        help="Cluster correlated sites within each gene")
    parser.add_argument("--correlation-threshold", type=float, default=0.9,
                        help="Minimum Pearson r for sites to cluster together")
    parser.add_argument("--max-missing-fraction", type=float, default=0.4,
                        help="Drop thinned phenotypes with more than this fraction missing")
    return parser.parse_args()


def ratios_to_fractions(matrix: pd.DataFrame) -> pd.DataFrame:
    fractions = pd.DataFrame(index=matrix.index, columns=matrix.columns, dtype=float)
    for sample_id in matrix.columns:
        num_den = matrix[sample_id].str.split("/", expand=True)
        num = pd.to_numeric(num_den[0], errors="coerce")
        den = pd.to_numeric(num_den[1], errors="coerce")
        fractions[sample_id] = (num + 0.5) / (den + 0.5)
        fractions.loc[den == 0, sample_id] = np.nan
    return fractions


def filter_and_impute_sites(edit_levels: pd.DataFrame, max_missing_fraction: float) -> pd.DataFrame:
    max_missing = edit_levels.shape[1] * max_missing_fraction
    edit_levels = edit_levels.loc[edit_levels.isna().sum(axis=1) <= max_missing].copy()
    return edit_levels.apply(lambda row: row.fillna(row.mean()), axis=1)


def cluster_sites(edit_levels: pd.DataFrame, threshold: float) -> dict[str, list[str]]:
    if edit_levels.shape[0] == 1:
        return {"Cluster_1": edit_levels.index.tolist()}

    corr = edit_levels.T.corr().fillna(0).clip(lower=-1.0, upper=1.0)
    dist = 1 - corr
    dist[dist < 0] = 0
    dist_matrix = squareform(dist.values, checks=False)

    clusters = defaultdict(list)
    linkage_matrix = linkage(dist_matrix, method="average")
    cluster_labels = fcluster(linkage_matrix, t=1 - threshold, criterion="distance")
    for site_id, label in zip(edit_levels.index, cluster_labels):
        clusters[f"Cluster_{label}"].append(site_id)
    return dict(clusters)


def main():
    args = parse_args()

    matrix = pd.read_csv(args.edit_matrix, sep="\t", index_col=0, dtype=str)
    site_to_gene = pd.read_csv(
        args.site_to_gene,
        sep="\t",
        header=None,
        names=["site", "gene_id"],
    )

    edit_levels = ratios_to_fractions(matrix)
    edit_levels = filter_and_impute_sites(edit_levels, args.max_missing_fraction)
    site_to_gene = site_to_gene[site_to_gene["site"].isin(edit_levels.index)]

    output_rows = []
    site_map_rows = []

    for gene_id, gene_sites in site_to_gene.groupby("gene_id", sort=True):
        site_ids = gene_sites["site"].tolist()
        gene_matrix = edit_levels.loc[edit_levels.index.intersection(site_ids)]
        if gene_matrix.empty:
            continue

        if args.cluster:
            clusters = cluster_sites(gene_matrix, args.correlation_threshold)
        else:
            clusters = {site_id: [site_id] for site_id in gene_matrix.index}

        for cluster_id, cluster_sites_ids in sorted(clusters.items()):
            if args.cluster:
                phenotype_id = f"{gene_id}__RNAedit_{cluster_id}"
            else:
                phenotype_id = f"{gene_id}__{cluster_id}"

            cluster_mean = gene_matrix.loc[cluster_sites_ids].mean(axis=0)
            output_rows.append(pd.DataFrame([cluster_mean], index=[phenotype_id]))
            for site_id in cluster_sites_ids:
                site_map_rows.append((site_id, gene_id, phenotype_id))

    if not output_rows:
        raise ValueError("No RNA editing phenotypes remained after gene mapping and filtering")

    phenotypes = pd.concat(output_rows)
    phenotypes.index.name = "phenotype_id"
    phenotypes.to_csv(args.output_matrix, sep="\t", float_format="%g")

    site_map = pd.DataFrame(
        site_map_rows,
        columns=["site", "gene_id", "phenotype_id"],
    )
    site_map.to_csv(args.output_site_map, sep="\t", index=False)


if __name__ == "__main__":
    main()
