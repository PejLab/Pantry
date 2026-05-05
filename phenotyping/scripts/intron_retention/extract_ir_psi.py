#!/usr/bin/env python3
"""Extract retained-intron PSI values from a MAJIQ PSI TSV."""

import argparse
from pathlib import Path

import pandas as pd


def normalize_gene_id(value):
    if pd.isna(value):
        return value
    return str(value).split(".", 1)[0]


def parse_bool(value):
    if isinstance(value, bool):
        return value
    if pd.isna(value):
        return False
    text = str(value).strip().lower()
    if text in {"true", "t", "1", "yes"}:
        return True
    if text in {"false", "f", "0", "no", ""}:
        return False
    raise ValueError(f"Cannot parse boolean value: {value!r}")


def read_majiq_tsv(path):
    return pd.read_csv(path, sep="\t", comment="#")


def explode_semicolon_columns(df):
    semicolon_cols = [
        col
        for col in df.columns
        if df[col].map(lambda x: isinstance(x, str) and ";" in x).any()
    ]
    if not semicolon_cols:
        return df

    rows = []
    for _, row in df.iterrows():
        split_values = {}
        max_len = 1
        for col in semicolon_cols:
            value = row[col]
            parts = str(value).split(";") if isinstance(value, str) else [value]
            split_values[col] = parts
            max_len = max(max_len, len(parts))
        for col, parts in split_values.items():
            if len(parts) not in {1, max_len}:
                raise ValueError(
                    f"Cannot explode row with inconsistent semicolon field lengths in {semicolon_cols}"
                )
        for idx in range(max_len):
            new_row = row.copy()
            for col, parts in split_values.items():
                new_row[col] = parts[idx] if len(parts) == max_len else parts[0]
            rows.append(new_row)
    return pd.DataFrame(rows, columns=df.columns)


def psi_columns(df):
    suffixes = [" raw_psi_mean", " psi", " PSI"]
    for suffix in suffixes:
        cols = [col for col in df.columns if col.endswith(suffix)]
        if cols:
            return cols, suffix
    available = ", ".join(df.columns)
    raise ValueError(
        "Could not identify sample PSI columns. Expected columns ending in "
        "' raw_psi_mean', ' psi', or ' PSI'. Available columns: "
        f"{available}"
    )


def build_ir_id(row):
    ec_idx = row.get("ec_idx", row.name)
    return (
        f"IR:{row['gene_id_base']}:{row['seqid']}:{int(row['start'])}:"
        f"{int(row['end'])}:{row['strand']}:{ec_idx}"
    )


def extract_ir(input_path, output_path):
    df = explode_semicolon_columns(read_majiq_tsv(input_path))
    if "is_intron" not in df.columns:
        available = ", ".join(df.columns)
        raise ValueError(
            "MAJIQ PSI TSV does not contain an 'is_intron' column. "
            f"Available columns: {available}"
        )
    required = {"seqid", "start", "end", "strand", "gene_id"}
    missing = sorted(required - set(df.columns))
    if missing:
        available = ", ".join(df.columns)
        raise ValueError(
            f"MAJIQ PSI TSV is missing required columns {missing}. Available columns: {available}"
        )

    sample_cols, suffix = psi_columns(df)
    ir = df.loc[df["is_intron"].map(parse_bool)].copy()
    ir["gene_id_base"] = ir["gene_id"].map(normalize_gene_id)
    if "ec_idx" not in ir.columns:
        ir.insert(0, "ec_idx", ir.index)
    ir["majiq_ir_id"] = ir.apply(build_ir_id, axis=1)

    rename = {col: col[: -len(suffix)] for col in sample_cols}
    ir = ir.rename(columns=rename)
    sample_cols = [rename[col] for col in sample_cols]

    metadata = [
        col
        for col in [
            "majiq_ir_id",
            "ec_idx",
            "seqid",
            "start",
            "end",
            "strand",
            "gene_id",
            "gene_id_base",
            "gene_name",
            "event_type",
            "is_denovo",
            "event_denovo",
            "ref_exon_start",
            "ref_exon_end",
            "other_exon_start",
            "other_exon_end",
        ]
        if col in ir.columns
    ]
    out = ir[metadata + sample_cols]
    out.to_csv(output_path, sep="\t", index=False, compression="infer")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Extract retained-intron rows from a MAJIQ PSI TSV and write a "
            "metadata-plus-sample PSI table."
        )
    )
    parser.add_argument(
        "--input",
        required=True,
        help="MAJIQ PSI TSV produced by `majiq psi --output-tsv`.",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output retained-intron PSI table; compression is inferred from the suffix.",
    )
    args = parser.parse_args()
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    extract_ir(args.input, args.output)


if __name__ == "__main__":
    main()
