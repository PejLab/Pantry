#!/usr/bin/env python3
"""Fit a deterministic batch of FUSION models concurrently."""

import argparse
import csv
import gzip
import os
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path


STATUS_FIELDS = [
    "bed_index",
    "phenotype_id",
    "chromosome",
    "position",
    "outcome",
    "weight_path",
    "heritability",
    "heritability_se",
    "heritability_p",
    "diagnostic",
]

SKIP_PATTERNS = [
    (re.compile(r"heritability .*skipping gene", re.IGNORECASE), "insufficient heritability"),
    (re.compile(r"likely GCTA could not converge, skipping gene", re.IGNORECASE), "GCTA did not converge"),
    (re.compile(r"Only one SNP available, skipping this gene", re.IGNORECASE), "only one cis-SNP"),
]


@dataclass(frozen=True)
class Model:
    bed_index: int
    chromosome: str
    position: int
    phenotype_id: str
    values: tuple[str, ...]


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--geno", required=True, type=Path)
    parser.add_argument("--bed", required=True, type=Path)
    parser.add_argument("--covar", required=True, type=Path)
    parser.add_argument("--modality", required=True)
    parser.add_argument("--batch-start", required=True, type=int)
    parser.add_argument("--batch-end", required=True, type=int)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--status", required=True, type=Path)
    parser.add_argument("--threads", required=True, type=int)
    parser.add_argument("--fusion-script", required=True, type=Path)
    parser.add_argument("--gcta", required=True, type=Path)
    parser.add_argument("--gemma", required=True, type=Path)
    parser.add_argument("--plink", default="plink")
    parser.add_argument("--rscript", default="Rscript")
    return parser.parse_args()


def classify_missing_weight(log_text):
    for pattern, diagnostic in SKIP_PATTERNS:
        if pattern.search(log_text):
            return "skipped", diagnostic
    return "failed", "FUSION exited without producing a weight file"


def classify_command_failure(error, log_text):
    """Recognize narrowly defined scientific skips from failed commands."""
    command = [str(item) for item in error.cmd]
    empty_cis_window = (
        error.returncode == 12
        and "Error: All variants excluded." in log_text
        and "--make-bed" in command
        and "--from-bp" in command
        and "--to-bp" in command
    )
    if empty_cis_window:
        return "skipped", "no variants in cis-window"
    return "failed", f"Command exited with status {error.returncode}"


def read_batch(bed_path, start, end):
    models = []
    with gzip.open(bed_path, "rt", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        samples = header[4:]
        for bed_index, row in enumerate(reader, start=1):
            if bed_index < start:
                continue
            if bed_index > end:
                break
            if len(row) != len(header):
                raise ValueError(
                    f"BED row {bed_index} has {len(row)} columns; expected {len(header)}"
                )
            models.append(
                Model(
                    bed_index=bed_index,
                    chromosome=row[0],
                    position=int(row[2]),
                    phenotype_id=row[3],
                    values=tuple(row[4:]),
                )
            )
    if not models:
        raise ValueError(f"No phenotype rows found for batch {start}-{end}")
    return samples, models


def read_sample_ids(fam_path, samples):
    sample_ids = {}
    with fam_path.open() as handle:
        for line in handle:
            fields = line.split()
            sample_ids[fields[1]] = (fields[0], fields[1])
    missing = [sample for sample in samples if sample not in sample_ids]
    if missing:
        raise ValueError(f"Phenotype samples missing from genotypes: {missing}")
    return [sample_ids[sample] for sample in samples]


def run_logged(command, log_handle):
    result = subprocess.run(
        [str(item) for item in command],
        stdout=log_handle,
        stderr=subprocess.STDOUT,
        text=True,
    )
    if result.returncode != 0:
        raise subprocess.CalledProcessError(result.returncode, command)


def validate_tools(args):
    """Fail once, before model dispatch, when a shared executable is unusable."""
    commands = [
        [args.plink, "--help"],
        [args.rscript, "--version"],
        [args.gcta, "--version"],
        [args.gemma, "-h"],
    ]
    for command in commands:
        subprocess.run(
            [str(item) for item in command],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=True,
        )


def parse_hsq(path):
    if not path.exists():
        return "", "", ""
    fields = path.read_text().strip().split()
    if len(fields) < 4:
        raise ValueError(f"Malformed FUSION heritability output: {path}")
    return fields[-3], fields[-2], fields[-1]


def diagnostic_tail(log_path, lines=3):
    text = log_path.read_text(errors="replace")
    nonempty = [line.strip() for line in text.splitlines() if line.strip()]
    return " | ".join(nonempty[-lines:])


def result_row(model, outcome, weight_path="", hsq=("", "", ""), diagnostic=""):
    return {
        "bed_index": model.bed_index,
        "phenotype_id": model.phenotype_id,
        "chromosome": model.chromosome,
        "position": model.position,
        "outcome": outcome,
        "weight_path": str(weight_path),
        "heritability": hsq[0],
        "heritability_se": hsq[1],
        "heritability_p": hsq[2],
        "diagnostic": diagnostic,
    }


def fit_model(model, sample_ids, args, batch_temp, weights_dir):
    safe_id = re.sub(r"[^A-Za-z0-9_.-]", "_", model.phenotype_id)
    work_dir = batch_temp / f"{model.bed_index:09d}_{safe_id}"
    if work_dir.exists():
        shutil.rmtree(work_dir)
    work_dir.mkdir(parents=True)
    log_path = work_dir / "model.log"
    pheno_path = work_dir / "phenotype.tsv"
    locus_prefix = work_dir / "locus"
    staging_prefix = work_dir / "weight"
    fusion_tmp = work_dir / "fusion_tmp"
    final_weight = weights_dir / f"{model.phenotype_id}.wgt.RDat"
    final_weight.unlink(missing_ok=True)
    start = max(1, model.position - 500_000)
    end = model.position + 500_000

    with pheno_path.open("w") as handle:
        for (fid, iid), value in zip(sample_ids, model.values):
            handle.write(f"{fid}\t{iid}\t{value}\n")

    try:
        with log_path.open("w") as log_handle:
            run_logged(
                [
                    args.plink,
                    "--bfile", args.geno,
                    "--pheno", pheno_path,
                    "--keep", pheno_path,
                    "--make-bed",
                    "--out", locus_prefix,
                    "--chr", model.chromosome,
                    "--from-bp", start,
                    "--to-bp", end,
                ],
                log_handle,
            )
            run_logged(
                [
                    args.rscript, args.fusion_script,
                    "--bfile", locus_prefix,
                    "--covar", args.covar,
                    "--tmp", fusion_tmp,
                    "--out", staging_prefix,
                    "--verbose", "0",
                    "--save_hsq",
                    "--noclean",
                    "--PATH_plink", args.plink,
                    "--PATH_gcta", args.gcta,
                    "--PATH_gemma", args.gemma,
                    "--models", "blup,lasso,top1,enet",
                ],
                log_handle,
            )
    except subprocess.CalledProcessError as error:
        log_text = log_path.read_text(errors="replace")
        outcome, diagnostic = classify_command_failure(error, log_text)
        if outcome == "skipped":
            row = result_row(model, outcome, diagnostic=diagnostic)
            shutil.rmtree(work_dir)
            return row
        return result_row(
            model,
            outcome,
            diagnostic=(
                f"{diagnostic}: {diagnostic_tail(log_path)}; "
                f"temporary files retained at {work_dir}"
            ),
        )
    except Exception as error:
        return result_row(
            model,
            "failed",
            diagnostic=f"{error}; temporary files retained at {work_dir}",
        )

    staged_weight = Path(f"{staging_prefix}.wgt.RDat")
    hsq_path = Path(f"{staging_prefix}.hsq")
    try:
        hsq = parse_hsq(hsq_path)
    except ValueError as error:
        return result_row(
            model,
            "failed",
            diagnostic=f"{error}; temporary files retained at {work_dir}",
        )

    if staged_weight.exists() and staged_weight.stat().st_size > 0:
        weights_dir.mkdir(parents=True, exist_ok=True)
        os.replace(staged_weight, final_weight)
        row = result_row(model, "success", final_weight, hsq, "weight fitted")
        shutil.rmtree(work_dir)
        return row

    log_text = log_path.read_text(errors="replace")
    outcome, diagnostic = classify_missing_weight(log_text)
    if outcome == "failed":
        diagnostic = f"{diagnostic}; temporary files retained at {work_dir}"
    row = result_row(model, outcome, hsq=hsq, diagnostic=diagnostic)
    if outcome == "skipped":
        shutil.rmtree(work_dir)
    return row


def write_status(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(f"{path}.tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=STATUS_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def main():
    args = parse_args()
    if args.threads < 1:
        raise ValueError("--threads must be at least 1")
    validate_tools(args)
    args.output_dir = args.output_dir.resolve()
    args.status = args.status.resolve()
    samples, models = read_batch(args.bed, args.batch_start, args.batch_end)
    sample_ids = read_sample_ids(Path(f"{args.geno}.fam"), samples)
    batch_name = f"{args.batch_start}_{args.batch_end}"
    batch_temp = args.output_dir / f"tmp_{args.modality}" / batch_name
    weights_dir = args.output_dir / args.modality
    batch_temp.mkdir(parents=True, exist_ok=True)

    worker = lambda model: fit_model(model, sample_ids, args, batch_temp, weights_dir)
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        rows = list(executor.map(worker, models))

    failures = [row for row in rows if row["outcome"] == "failed"]
    if failures:
        for row in failures:
            print(
                f"FAILED {row['phenotype_id']}: {row['diagnostic']}",
                file=sys.stderr,
            )
        raise RuntimeError(f"{len(failures)} of {len(rows)} TWAS models failed")

    write_status(args.status, rows)
    shutil.rmtree(batch_temp)


if __name__ == "__main__":
    main()
