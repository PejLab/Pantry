#!/usr/bin/env python3
"""Convert a GTF annotation to a gene-parented GFF3 for MAJIQ."""

import argparse
from pathlib import Path
from urllib.parse import quote


def parse_attributes(text):
    attrs = {}
    for field in text.strip().rstrip(";").split(";"):
        field = field.strip()
        if not field:
            continue
        if " " not in field:
            attrs[field] = ""
            continue
        key, value = field.split(" ", 1)
        attrs[key] = value.strip().strip('"')
    return attrs


def gff3_attributes(attrs):
    return ";".join(
        f"{quote(str(key), safe='._:-')}={quote(str(value), safe='._:-,')}"
        for key, value in attrs.items()
        if value is not None and value != ""
    )


def versioned_id(attrs, id_key, version_key):
    value = attrs[id_key]
    version = attrs.get(version_key)
    if version and not value.endswith(f".{version}"):
        return f"{value}.{version}"
    return value


def convert_gtf_to_gff3(input_gtf, output_gff3):
    records = []
    genes = set()
    transcripts = set()
    gene_info = {}
    gene_bounds = {}
    transcript_info = {}
    transcript_bounds = {}

    with open(input_gtf) as in_fh:
        for line in in_fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"Expected 9 GTF columns, got {len(fields)}: {line.rstrip()}")

            attrs = parse_attributes(fields[8])
            if "gene_id" not in attrs:
                continue
            gene_id = versioned_id(attrs, "gene_id", "gene_version")
            records.append((fields, attrs, gene_id))
            if fields[2] == "gene":
                genes.add(gene_id)
            if fields[2] == "transcript" and "transcript_id" in attrs:
                transcripts.add(versioned_id(attrs, "transcript_id", "transcript_version"))
            start = int(fields[3])
            end = int(fields[4])
            if gene_id not in gene_bounds:
                gene_bounds[gene_id] = [fields[0], fields[1], start, end, fields[6]]
            else:
                bounds = gene_bounds[gene_id]
                bounds[2] = min(bounds[2], start)
                bounds[3] = max(bounds[3], end)
            gene_info.setdefault(
                gene_id,
                {
                    "gene_name": attrs.get("gene_name"),
                    "gene_type": attrs.get("gene_type", attrs.get("gene_biotype")),
                },
            )
            if fields[2] != "gene" and "transcript_id" in attrs:
                transcript_id = versioned_id(attrs, "transcript_id", "transcript_version")
                if transcript_id not in transcript_bounds:
                    transcript_bounds[transcript_id] = [fields[0], fields[1], start, end, fields[6]]
                else:
                    bounds = transcript_bounds[transcript_id]
                    bounds[2] = min(bounds[2], start)
                    bounds[3] = max(bounds[3], end)
                transcript_info.setdefault(
                    transcript_id,
                    {
                        "gene_id": gene_id,
                        "gene_name": attrs.get("gene_name"),
                        "gene_type": attrs.get("gene_type", attrs.get("gene_biotype")),
                        "transcript_type": attrs.get("transcript_type", attrs.get("transcript_biotype")),
                        "transcript_name": attrs.get("transcript_name"),
                        "source": fields[1],
                    },
                )

    Path(output_gff3).parent.mkdir(parents=True, exist_ok=True)
    with open(output_gff3, "w") as out_fh:
        out_fh.write("##gff-version 3\n")
        emitted_genes = set()
        emitted_transcripts = set()
        for fields, attrs, gene_id in records:
            seqid, source, feature, start, end, score, strand, phase, attr_text = fields
            if feature != "gene" and gene_id not in genes and gene_id not in emitted_genes:
                bounds = gene_bounds[gene_id]
                info = gene_info[gene_id]
                synth_attrs = {
                    "ID": gene_id,
                    "gene_id": gene_id,
                    "gene_name": info.get("gene_name"),
                    "gene_type": info.get("gene_type"),
                }
                out_fh.write(
                    "\t".join(
                        [
                            bounds[0],
                            bounds[1],
                            "gene",
                            str(bounds[2]),
                            str(bounds[3]),
                            ".",
                            bounds[4],
                            ".",
                            gff3_attributes(synth_attrs),
                        ]
                    )
                    + "\n"
                )
                emitted_genes.add(gene_id)

            out_attrs = {
                "gene_id": gene_id,
                "gene_name": attrs.get("gene_name"),
                "gene_type": attrs.get("gene_type", attrs.get("gene_biotype")),
            }

            if feature == "gene":
                if gene_id in emitted_genes:
                    continue
                out_attrs = {"ID": gene_id, **out_attrs}
                emitted_genes.add(gene_id)
            else:
                if "transcript_id" not in attrs:
                    continue
                transcript_id = versioned_id(attrs, "transcript_id", "transcript_version")
                out_attrs.update(
                    {
                        "transcript_id": transcript_id,
                        "transcript_type": attrs.get("transcript_type", attrs.get("transcript_biotype")),
                        "transcript_name": attrs.get("transcript_name"),
                    }
                )
                if feature == "transcript":
                    if transcript_id in emitted_transcripts:
                        continue
                    out_attrs = {"ID": transcript_id, "Parent": gene_id, **out_attrs}
                    emitted_transcripts.add(transcript_id)
                else:
                    if transcript_id not in transcripts and transcript_id not in emitted_transcripts:
                        bounds = transcript_bounds[transcript_id]
                        info = transcript_info[transcript_id]
                        synth_attrs = {
                            "ID": transcript_id,
                            "Parent": info["gene_id"],
                            "gene_id": info["gene_id"],
                            "gene_name": info.get("gene_name"),
                            "gene_type": info.get("gene_type"),
                            "transcript_id": transcript_id,
                            "transcript_type": info.get("transcript_type"),
                            "transcript_name": info.get("transcript_name"),
                        }
                        out_fh.write(
                            "\t".join(
                                [
                                    bounds[0],
                                    info.get("source", bounds[1]),
                                    "transcript",
                                    str(bounds[2]),
                                    str(bounds[3]),
                                    ".",
                                    bounds[4],
                                    ".",
                                    gff3_attributes(synth_attrs),
                                ]
                            )
                            + "\n"
                        )
                        emitted_transcripts.add(transcript_id)
                    exon_number = attrs.get("exon_number")
                    suffix = exon_number if exon_number else f"{start}-{end}"
                    feature_id = f"{feature}:{transcript_id}:{suffix}"
                    out_attrs = {"ID": feature_id, "Parent": transcript_id, **out_attrs}
                    if f"{feature}_id" in attrs:
                        out_attrs[f"{feature}_id"] = attrs[f"{feature}_id"]
                    elif "exon_id" in attrs:
                        out_attrs["exon_id"] = attrs["exon_id"]
                    if exon_number:
                        out_attrs["exon_number"] = exon_number

            out_fh.write(
                "\t".join(
                    [
                        seqid,
                        source,
                        feature,
                        start,
                        end,
                        score,
                        strand,
                        phase,
                        gff3_attributes(out_attrs),
                    ]
                )
                + "\n"
            )


def main():
    parser = argparse.ArgumentParser(description="Convert GTF to MAJIQ-compatible GFF3")
    parser.add_argument("--input", required=True, help="Input GTF")
    parser.add_argument("--output", required=True, help="Output GFF3")
    args = parser.parse_args()
    convert_gtf_to_gff3(args.input, args.output)


if __name__ == "__main__":
    main()
