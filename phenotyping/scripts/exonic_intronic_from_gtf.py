"""
Extract consensus exonic regions and intronic gaps from a GTF file.

- For exonic regions, find all regions per gene that are covered by at least a
  certain fraction of isoforms.
- For intronic regions, find all regions between exons of the gene that overlap no
  exons in any of the gene's isoforms.
- Any regions covered by at least one exon but less than the minimum fraction
  specified will not be reported.

Usage:
        python exons_introns_from_gtf.py input.gtf consensus_exons.gtf intronic_gaps.gtf \
                [--min-exon-fraction 0.8] [--min-transcripts 1] [--verbose]
"""

import argparse
import sys
from collections import defaultdict
from pathlib import Path
import gzip
from typing import Dict, List, Tuple, Set, Optional


class GTFRecord:
    """Represents a GTF record with parsed attributes."""
    
    def __init__(self, line: str):
        parts = line.strip().split('\t')
        if len(parts) != 9:
            raise ValueError(f"Invalid GTF line: {line}")
        
        self.seqname = parts[0]
        self.source = parts[1]
        self.feature = parts[2]
        self.start = int(parts[3])
        self.end = int(parts[4])
        self.score = parts[5]
        self.strand = parts[6]
        self.frame = parts[7]
        self.attributes = self._parse_attributes(parts[8])
    
    def _parse_attributes(self, attr_str: str) -> Dict[str, str]:
        """Parse GTF attributes string into dictionary."""
        attrs = {}
        for item in attr_str.split(';'):
            item = item.strip()
            if ' ' in item:
                key, value = item.split(' ', 1)
                # Remove quotes from value
                value = value.strip('"')
                attrs[key] = value
        return attrs
    
    def get(self, key: str, default: str = None) -> Optional[str]:
        """Get attribute value with fallback to default."""
        return self.attributes.get(key, default)


class AnnotationParser:
    """Parse GTF into gene/transcript structures and derive regions."""

    def __init__(self, min_exon_fraction: float = 0.8, min_transcripts_per_gene: int = 1):
        self.min_exon_fraction = min_exon_fraction
        self.min_transcripts_per_gene = min_transcripts_per_gene

    def parse_gtf(self, gtf_file: Path) -> Tuple[Dict[str, Set[str]], Dict[str, List[Dict]]]:
        genes: Dict[str, Set[str]] = defaultdict(set)
        transcripts: Dict[str, List[Dict]] = defaultdict(list)
        open_func = gzip.open if gtf_file.suffix == '.gz' else open
        mode = 'rt' if gtf_file.suffix == '.gz' else 'r'
        with open_func(gtf_file, mode) as f:
            for line_num, line in enumerate(f, 1):
                if not line or line.startswith('#'):
                    continue
                try:
                    record = GTFRecord(line)
                except ValueError as e:
                    print(f"Warning: Skipping invalid line {line_num}: {e}", file=sys.stderr)
                    continue
                if record.feature != 'exon':
                    continue
                gene_id = record.get('gene_id')
                transcript_id = record.get('transcript_id')
                if not gene_id or not transcript_id:
                    continue
                exon = {
                    'seqname': record.seqname,
                    'start': record.start,
                    'end': record.end,
                    'strand': record.strand,
                    'gene_id': gene_id,
                    'transcript_id': transcript_id,
                }
                genes[gene_id].add(transcript_id)
                transcripts[transcript_id].append(exon)
        return genes, transcripts

    def find_consensus_exons(self, genes: Dict[str, Set[str]], transcripts: Dict[str, List[Dict]]) -> List[Dict]:
        """Compute per-segment coverage and merge qualifying segments."""
        consensus: List[Dict] = []
        for gene_id, tx_ids in genes.items():
            n_tx = len(tx_ids)
            if n_tx < self.min_transcripts_per_gene:
                continue
            seqname = None
            strand = None
            per_tx: Dict[str, List[Tuple[int, int]]] = {}
            for tx in tx_ids:
                exons = transcripts.get(tx)
                if not exons:
                    continue
                ints: List[Tuple[int, int]] = []
                for ex in exons:
                    seqname = ex['seqname'] if seqname is None else seqname
                    strand = ex['strand'] if strand is None else strand
                    ints.append((ex['start'], ex['end']))
                per_tx[tx] = ints
            if not per_tx:
                continue
            breakpoints: Set[int] = set()
            for ints in per_tx.values():
                for s, e in ints:
                    breakpoints.add(s)
                    breakpoints.add(e + 1)
            bps = sorted(breakpoints)
            if len(bps) < 2:
                continue
            idx = {p: i for i, p in enumerate(bps)}
            diff = [0] * (len(bps) + 1)
            for ints in per_tx.values():
                for s, e in ints:
                    diff[idx[s]] += 1
                    if (e + 1) in idx:
                        diff[idx[e + 1]] -= 1
            cov = [0] * len(bps)
            running = 0
            for i in range(len(bps)):
                running += diff[i]
                cov[i] = running
            current = None  # [start, end, min_fraction]
            for i in range(len(bps) - 1):
                seg_s = bps[i]
                seg_e = bps[i + 1] - 1
                if seg_e < seg_s:
                    continue
                fraction = cov[i] / n_tx
                if fraction >= self.min_exon_fraction:
                    if current is None:
                        current = [seg_s, seg_e, fraction]
                    else:
                        if seg_s == current[1] + 1:
                            current[1] = seg_e
                            if fraction < current[2]:
                                current[2] = fraction
                        else:
                            consensus.append({
                                'seqname': seqname,
                                'start': current[0],
                                'end': current[1],
                                'strand': strand,
                                'gene_id': gene_id,
                                'fraction': current[2],
                            })
                            current = [seg_s, seg_e, fraction]
                else:
                    if current is not None:
                        consensus.append({
                            'seqname': seqname,
                            'start': current[0],
                            'end': current[1],
                            'strand': strand,
                            'gene_id': gene_id,
                            'fraction': current[2],
                        })
                        current = None
            if current is not None:
                consensus.append({
                    'seqname': seqname,
                    'start': current[0],
                    'end': current[1],
                    'strand': strand,
                    'gene_id': gene_id,
                    'fraction': current[2],
                })
        return consensus

    def find_intronic_gaps(self, genes: Dict[str, Set[str]], transcripts: Dict[str, List[Dict]]) -> List[Dict]:
        introns: List[Dict] = []
        for gene_id, tx_ids in genes.items():
            if len(tx_ids) < self.min_transcripts_per_gene:
                continue
            seqname = None
            strand = None
            all_exons: List[Tuple[int, int]] = []
            for tx in tx_ids:
                for ex in transcripts.get(tx, []):
                    seqname = ex['seqname'] if seqname is None else seqname
                    strand = ex['strand'] if strand is None else strand
                    all_exons.append((ex['start'], ex['end']))
            if not all_exons:
                continue
            all_exons.sort()
            merged: List[Tuple[int, int]] = []
            cur_s, cur_e = all_exons[0]
            for s, e in all_exons[1:]:
                if s <= cur_e + 1:
                    if e > cur_e:
                        cur_e = e
                else:
                    merged.append((cur_s, cur_e))
                    cur_s, cur_e = s, e
            merged.append((cur_s, cur_e))
            for i in range(len(merged) - 1):
                gap_s = merged[i][1] + 1
                gap_e = merged[i + 1][0] - 1
                if gap_e >= gap_s:
                    introns.append({
                        'seqname': seqname,
                        'start': gap_s,
                        'end': gap_e,
                        'strand': strand,
                        'gene_id': gene_id,
                    })
        return introns
    
    def write_gtf(self, records: List[Dict], output_file: Path, feature_type: str):
        with open(output_file, 'w') as f:
            for record in records:
                attrs = [f'gene_id "{record["gene_id"]}"']
                if 'fraction' in record:
                    attrs.append(f'fraction "{record["fraction"]:.3f}"')
                attr_str = '; '.join(attrs) + ';'
                f.write(f'{record["seqname"]}\t.\t{feature_type}\t{record["start"]}\t{record["end"]}\t.\t{record["strand"]}\t.\t{attr_str}\n')


def main():
    parser = argparse.ArgumentParser(
        description='Derive consensus exons and intronic gaps from a GTF',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python exons_introns_from_gtf.py input.gtf consensus_exons.gtf introns.gtf
  python exons_introns_from_gtf.py input.gtf consensus_exons.gtf introns.gtf --min-exon-fraction 0.7
  python exons_introns_from_gtf.py input.gtf consensus_exons.gtf introns.gtf --min-transcripts 3
        """
    )
    parser.add_argument('input', type=Path, help='Input GTF (optionally gz)')
    parser.add_argument('output_exons', type=Path, help='Output consensus exons GTF')
    parser.add_argument('output_introns', type=Path, help='Output intronic gaps GTF')
    parser.add_argument('--min-exon-fraction', type=float, default=0.8,
                        help='Minimum transcript coverage fraction for a segment (default 0.8)')
    parser.add_argument('--min-transcripts', type=int, default=1,
                        help='Minimum transcripts per gene to consider (default 1)')
    parser.add_argument('--verbose', action='store_true', help='Verbose logging')
    
    args = parser.parse_args()
    
    if not args.input.exists():
        print(f"Error: Input file {args.input} does not exist", file=sys.stderr)
        sys.exit(1)
    
    if args.verbose:
        print(f"Input: {args.input}")
        print(f"Min exon fraction: {args.min_exon_fraction}")
        print(f"Min transcripts per gene: {args.min_transcripts}")
    
    # Parse GTF file
    parser_obj = AnnotationParser(
        min_exon_fraction=args.min_exon_fraction,
        min_transcripts_per_gene=args.min_transcripts
    )
    
    try:
        genes, transcripts = parser_obj.parse_gtf(args.input)
        
        if args.verbose:
            print(f"Found {len(genes)} genes with {sum(len(t) for t in genes.values())} transcripts")
        
        consensus_exons = parser_obj.find_consensus_exons(genes, transcripts)
        if args.verbose:
            print(f"Consensus exon regions: {len(consensus_exons)}")
        introns = parser_obj.find_intronic_gaps(genes, transcripts)
        if args.verbose:
            print(f"Intronic gap regions: {len(introns)}")
        parser_obj.write_gtf(consensus_exons, args.output_exons, 'exon')
        parser_obj.write_gtf(introns, args.output_introns, 'intron')
        if args.verbose:
            print(f"Wrote exons: {args.output_exons}")
            print(f"Wrote introns: {args.output_introns}")
    
    except Exception as e:
        print(f"Error processing file: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main() 
