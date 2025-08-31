"""Alter non-Ensembl GTF files to be compatible with txrevise

1. Rename gene_type/transcript_type to gene_biotype/transcript_biotype
2. Add gene_biotype/transcript_biotype fields to all non-gene entries that don't have them
"""

import argparse
import os
from gtfparse import read_gtf

def update_line(line, gene_type_dict):
    fields = line.strip().split("\t")
    fields[8] = fields[8].replace("gene_type", "gene_biotype")
    if fields[2] != "gene":
        fields[8] = fields[8].replace("transcript_type", "transcript_biotype")
        gene_type_present = "gene_biotype" in fields[8]
        transcript_type_present = "transcript_biotype" in fields[8]
        if (not gene_type_present) or (not transcript_type_present):
            gene_id = fields[8].split(";")[0].split(" ")[1].strip('"')
            assert gene_id in gene_type_dict, f"gene entry for gene_id {gene_id} not found"
        if not gene_type_present:
            fields[8] += f' gene_biotype "{gene_type_dict[gene_id]}";'
        if not transcript_type_present:
            fields[8] += f' transcript_biotype "{gene_type_dict[gene_id]}";'
    return "\t".join(fields) + "\n"

parser = argparse.ArgumentParser(description="Alter non-Ensembl GTF files to be compatible with txrevise")
parser.add_argument("--input", "-i", required=True, help="Input GTF file")
parser.add_argument("--output", "-o", required=True, help="Output GTF file")
args = parser.parse_args()

output_dir = os.path.dirname(args.output)
if output_dir:
    assert os.path.exists(output_dir), f"Output directory {output_dir} does not exist"

anno = read_gtf(args.input)
# Newer versions return a polars DF by default, but not all versions allow
# return type to be specified, so this handles older and newer versions:
if type(anno).__module__ == 'polars.dataframe.frame':
    anno = anno.to_pandas()

# In case there is no transcript_type/transcript_biotype field (e.g. RefSeq),
# store the gene_type/gene_biotype field and copy it to the
# transcript_type/transcript_biotype field for all of the gene's other entries.
gene_type_col = 'gene_type' if 'gene_type' in anno.columns else 'gene_biotype'
gene_type_dict = dict(
    anno.loc[anno['feature'] == 'gene', ['gene_id', gene_type_col]].drop_duplicates().values
)

with open(args.input, "r") as f:
    with open(args.output, "w") as f_out:
        for line in f:
            if not line.startswith("#"):
                line = update_line(line, gene_type_dict)
            f_out.write(line)

