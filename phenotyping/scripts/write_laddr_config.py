"""Write a LaDDR config file for Pantry outputs."""

import argparse
from pathlib import Path
import yaml


def parse_args():
    parser = argparse.ArgumentParser(description="Write a LaDDR config YAML")
    parser.add_argument("--output", type=Path, required=True, help="Output config YAML path")
    parser.add_argument("--gtf", type=Path, required=True, help="Reference GTF path")
    parser.add_argument("--manifest", type=Path, required=True, help="Coverage manifest path")
    parser.add_argument("--coverage-directory", type=Path, default=Path("covg_bigwig"), help="Coverage directory path in LaDDR project")
    parser.add_argument("--pheno-files", nargs="*", default=[], help="Pantry normalized BED files to regress out")
    parser.add_argument("--use-existing-bins", action="store_true", help="Set binning.use_existing=true")
    parser.add_argument("--min-samples-expressed", type=float, default=0.5)
    parser.add_argument("--binning-method", default="adaptive_diffvar")
    parser.add_argument("--batch-size", type=int, default=40)
    parser.add_argument("--max-bin-width", type=int, default=1024)
    parser.add_argument("--adaptive-max-samples", type=int, default=256)
    parser.add_argument("--adaptive-bins-per-gene", type=int, default=256)
    parser.add_argument("--adaptive-min-mean-total-covg", type=float, default=128.0)
    parser.add_argument("--adaptive-max-corr", type=float, default=0.8)
    parser.add_argument("--model-var-explained-max", type=float, default=0.8)
    parser.add_argument("--model-n-pcs-max", type=int, default=16)
    return parser.parse_args()


def main():
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    cfg = {
        "input": {
            "gtf": str(args.gtf.resolve()),
            "coverage": {
                "method": "manifest",
                "directory": str(args.coverage_directory),
                "manifest": str(args.manifest),
            },
            "min_samples_expressed": args.min_samples_expressed,
            "pheno_paths": [str(Path(p).resolve()) for p in args.pheno_files],
        },
        "binning": {
            "use_existing": args.use_existing_bins,
            "method": args.binning_method,
            "batch_size": args.batch_size,
            "max_bin_width": args.max_bin_width,
            "adaptive": {
                "max_samples": args.adaptive_max_samples,
                "bins_per_gene": args.adaptive_bins_per_gene,
                "min_mean_total_covg": args.adaptive_min_mean_total_covg,
                "max_corr": args.adaptive_max_corr,
            },
        },
        "model": {
            "use_existing": False,
            "var_explained_max": args.model_var_explained_max,
            "n_pcs_max": args.model_n_pcs_max,
            "use_fpca": False,
        },
    }

    with open(args.output, "w") as f:
        yaml.safe_dump(cfg, f, sort_keys=False)


if __name__ == "__main__":
    main()
