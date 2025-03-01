# Run in python script to use random_tiebreak option
# Based on examples at https://github.com/broadinstitute/tensorqtl

import argparse
import pandas as pd
from rpy2.robjects.packages import importr
import tensorqtl
from tensorqtl import genotypeio, cis, trans
import torch

parser = argparse.ArgumentParser(description="Run tensorQTL from command line but with extra options")
parser.add_argument("geno_prefix")
parser.add_argument("expression")
parser.add_argument("output", help="Output file path, or, for cis_nominal mode, output files prefix (with no directories)")
parser.add_argument("--covariates", help="Covariates file")
parser.add_argument("--cis_output", help="For cis_independent mode, output of tensorQTL in cis mode")
parser.add_argument("--groups", help="File with phenotype groups if applicable")
parser.add_argument("--window", type=int, default=1000000, help="cis-window size, default 1000000")
parser.add_argument("--output_dir", help="For cis_nominal mode, output directory")
parser.add_argument("--mode", required=True, help="Run mode: currently cis, cis_independent, cis_nominal, genome_wide, or trans. genome_wide runs `map_trans` without removing cis-variants, returning all associations with p<1e-5. trans does the same but removes all pairs with TSS distance < 5 Mb and all variants with MAF < 0.05.")
args = parser.parse_args()

# Check for rpy2 and qvalue packages first, since otherwise tensorQTL will do
# all the eQTL mapping and then fail.
_ = importr("qvalue")

if torch.cuda.is_available():
    print(f'  * using GPU ({torch.cuda.get_device_name(torch.cuda.current_device())})')
else:
    print('  * WARNING: using CPU!')
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

if args.output_dir is not None:
    assert args.mode == 'cis_nominal'

pheno, pheno_pos = tensorqtl.read_phenotype_bed(args.expression)
if args.covariates is not None:
    covar = pd.read_csv(args.covariates, sep="\t", index_col=0).T
else:
    covar = None
pr = genotypeio.PlinkReader(args.geno_prefix)
genotype_df = pr.load_genotypes()
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]
if args.groups is not None:
    groups = pd.read_csv(args.groups, sep="\t", index_col=0, header=None).squeeze('columns')
else:
    groups = None

# Remove phenotypes on chromosomes not present in genotypes to avoid error
chroms = variant_df['chrom'].unique()
keep = pheno_pos['chr'].isin(chroms)
pheno = pheno.loc[keep]
pheno_pos = pheno_pos.loc[keep]

if args.mode == "cis":
    d = cis.map_cis(genotype_df, variant_df, pheno, pheno_pos, covar, group_s=groups, window=args.window, random_tiebreak=True)
    tensorqtl.calculate_qvalues(d, qvalue_lambda=0.85)
    d.to_csv(args.output, sep="\t", float_format="%.6g")
elif args.mode == "cis_independent":
    cis_df = pd.read_csv(args.cis_output, sep="\t", index_col=0)
    d = cis.map_independent(genotype_df, variant_df, cis_df, pheno, pheno_pos, covar, group_s=groups, window=args.window, random_tiebreak=True)
    d.to_csv(args.output, sep="\t", index=False, float_format="%.6g")
elif args.mode == "genome_wide":
    d = trans.map_trans(genotype_df, pheno, covariates_df=covar, maf_threshold=0) # cis has default min MAF 0, trans has default 0.05
    d.to_csv(args.output, sep="\t", index=False, float_format="%.6g")
elif args.mode == "trans":
    d = trans.map_trans(genotype_df, pheno, covariates_df=covar)
    d = trans.filter_cis(d, pheno_pos, variant_df)
    d.to_csv(args.output, sep="\t", index=False, float_format="%.6g")
## trans-QTL permutation testing does not use phenotype groupings
# elif args.mode == "trans_perm":
#     trans_df = trans.map_trans(genotype_df, pheno, covariates_df=covar, return_r2=True)
#     trans_df = trans.filter_cis(trans_df, pheno_pos, variant_df)
#     perm_df = trans.map_permutations(genotype_df, covariates_df)
#     perm_output = trans.apply_permutations(perm_df, trans_df)
#     perm_output.to_csv(args.output, sep="\t", index=False, float_format="%.6g")
elif args.mode == "cis_nominal":
    assert args.output_dir is not None
    assert args.groups is None
    cis.map_nominal(genotype_df, variant_df, pheno, pheno_pos, args.output, covariates_df=covar,
                    group_s=None, window=args.window, output_dir=args.output_dir)
else:
    print("Mode not recognized.")
