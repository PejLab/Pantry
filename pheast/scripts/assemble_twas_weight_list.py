#!/usr/bin/env python3
"""Assemble the ordered FUSION weight list from batch status tables."""

import argparse
import csv
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--status", required=True, nargs="+", type=Path)
    parser.add_argument("--weights-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def successful_weights(status_paths):
    rows = []
    for status_path in status_paths:
        with status_path.open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row["outcome"] == "success":
                    rows.append((int(row["bed_index"]), Path(row["weight_path"])))
                elif row["outcome"] != "skipped":
                    raise ValueError(f"Unexpected outcome in {status_path}: {row['outcome']}")
    rows.sort(key=lambda item: item[0])
    indices = [item[0] for item in rows]
    if len(indices) != len(set(indices)):
        raise ValueError("Duplicate BED indices in TWAS batch status tables")
    return [path for _, path in rows]


def main():
    args = parse_args()
    paths = successful_weights(args.status)
    weights_dir = args.weights_dir.resolve()
    with args.output.open("w") as handle:
        for path in paths:
            handle.write(f"{path.resolve().relative_to(weights_dir)}\n")


if __name__ == "__main__":
    main()
