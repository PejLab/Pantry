"""
Extract constitutive exons and introns from GTF files for RNA stability analysis.

Usage:
    python exons_introns_from_gtf.py input.gtf output_exons.gtf output_introns.gtf [options]
"""

import argparse
import sys
from collections import defaultdict, Counter
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
    """Parse GTF files from different annotation sources."""
    
    def __init__(self, min_exon_fraction: float = 0.8, min_transcripts_per_gene: int = 2):
        self.min_exon_fraction = min_exon_fraction
        self.min_transcripts_per_gene = min_transcripts_per_gene
    
    def parse_gtf(self, gtf_file: Path) -> Tuple[Dict, Dict]:
        """Parse GTF file and return gene and transcript structures."""
        
        genes = defaultdict(list)  # gene_id -> list of transcripts
        transcripts = defaultdict(list)  # transcript_id -> list of exons
        
        # Determine if file is gzipped
        open_func = gzip.open if gtf_file.suffix == '.gz' else open
        mode = 'rt' if gtf_file.suffix == '.gz' else 'r'
        
        with open_func(gtf_file, mode) as f:
            for line_num, line in enumerate(f, 1):
                if line.startswith('#'):
                    continue
                
                try:
                    record = GTFRecord(line)
                except ValueError as e:
                    print(f"Warning: Skipping invalid line {line_num}: {e}", file=sys.stderr)
                    continue
                
                if record.feature != 'exon':
                    continue
                
                # Extract gene_id and transcript_id with flexible parsing
                gene_id = self._extract_gene_id(record)
                transcript_id = self._extract_transcript_id(record)
                
                if not gene_id or not transcript_id:
                    continue
                
                # Store exon information
                exon_info = {
                    'seqname': record.seqname,
                    'start': record.start,
                    'end': record.end,
                    'strand': record.strand,
                    'gene_id': gene_id,
                    'transcript_id': transcript_id
                }
                
                genes[gene_id].append(transcript_id)
                transcripts[transcript_id].append(exon_info)
        
        return genes, transcripts
    
    def _extract_gene_id(self, record: GTFRecord) -> Optional[str]:
        """Extract gene_id from GTF record with format flexibility."""
        # Try common gene_id attribute names
        for attr in ['gene_id', 'geneId', 'gene_name', 'geneName']:
            if attr in record.attributes:
                return record.attributes[attr]
        return None
    
    def _extract_transcript_id(self, record: GTFRecord) -> Optional[str]:
        """Extract transcript_id from GTF record with format flexibility."""
        # Try common transcript_id attribute names
        for attr in ['transcript_id', 'transcriptId', 'transcript_name', 'transcriptName']:
            if attr in record.attributes:
                return record.attributes[attr]
        return None
    
    def find_constitutive_exons(self, genes: Dict, transcripts: Dict) -> List[Dict]:
        """Find constitutive exons based on fraction threshold."""
        constitutive_exons = []
        
        for gene_id, transcript_ids in genes.items():
            if len(transcript_ids) < self.min_transcripts_per_gene:
                continue
            
            # Get all exons for this gene
            gene_exons = []
            for transcript_id in transcript_ids:
                if transcript_id in transcripts:
                    gene_exons.extend(transcripts[transcript_id])
            
            # Group exons by coordinates
            exon_groups = defaultdict(list)
            for exon in gene_exons:
                key = (exon['seqname'], exon['start'], exon['end'], exon['strand'])
                exon_groups[key].append(exon['transcript_id'])
            
            # Find constitutive exons
            for (seqname, start, end, strand), transcript_list in exon_groups.items():
                unique_transcripts = set(transcript_list)
                fraction = len(unique_transcripts) / len(transcript_ids)
                
                if fraction >= self.min_exon_fraction:
                    constitutive_exons.append({
                        'seqname': seqname,
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'gene_id': gene_id,
                        'fraction': fraction
                    })
        
        return constitutive_exons
    
    def find_introns(self, genes: Dict, transcripts: Dict) -> List[Dict]:
        """Find intronic regions between exons."""
        introns = []
        
        for gene_id, transcript_ids in genes.items():
            if len(transcript_ids) < self.min_transcripts_per_gene:
                continue
            
            # Process each transcript separately
            for transcript_id in transcript_ids:
                if transcript_id not in transcripts:
                    continue
                
                transcript_exons = sorted(transcripts[transcript_id], 
                                        key=lambda x: (x['start'], x['end']))
                
                # Find introns between adjacent exons
                for i in range(len(transcript_exons) - 1):
                    exon1 = transcript_exons[i]
                    exon2 = transcript_exons[i + 1]
                    
                    # Ensure exons are on same chromosome and strand
                    if (exon1['seqname'] != exon2['seqname'] or 
                        exon1['strand'] != exon2['strand']):
                        continue
                    
                    # Intron is between exon1.end and exon2.start
                    intron_start = exon1['end'] + 1
                    intron_end = exon2['start'] - 1
                    
                    # Skip if intron is too small or too large
                    if intron_end < intron_start or (intron_end - intron_start + 1) > 1000000:
                        continue
                    
                    introns.append({
                        'seqname': exon1['seqname'],
                        'start': intron_start,
                        'end': intron_end,
                        'strand': exon1['strand'],
                        'gene_id': gene_id,
                        'transcript_id': transcript_id
                    })
        
        return introns
    
    def write_gtf(self, records: List[Dict], output_file: Path, feature_type: str):
        """Write records to GTF format."""
        with open(output_file, 'w') as f:
            for record in records:
                if feature_type == 'exon':
                    attr_str = f'gene_id "{record["gene_id"]}";'
                else:  # intron
                    attr_str = f'gene_id "{record["gene_id"]}"; transcript_id "{record["transcript_id"]}";'
                
                f.write(f'{record["seqname"]}\t.\t{feature_type}\t{record["start"]}\t{record["end"]}\t.\t{record["strand"]}\t.\t{attr_str}\n')


def main():
    parser = argparse.ArgumentParser(
        description='Extract constitutive exons and introns from GTF files for RNA stability analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with default parameters
  python exons_introns_from_gtf.py input.gtf exons.gtf introns.gtf
  
  # More permissive constitutive exon definition (70% of transcripts)
  python exons_introns_from_gtf.py input.gtf exons.gtf introns.gtf --min-exon-fraction 0.7
  
  # Require at least 3 transcripts per gene
  python exons_introns_from_gtf.py input.gtf exons.gtf introns.gtf --min-transcripts 3
        """
    )
    
    parser.add_argument('input', type=Path, help='Input GTF file')
    parser.add_argument('output_exons', type=Path, help='Output constitutive exons GTF file')
    parser.add_argument('output_introns', type=Path, help='Output introns GTF file')
    parser.add_argument('--min-exon-fraction', type=float, default=0.8,
                       help='Minimum fraction of transcripts an exon must appear in to be constitutive (default: 0.8)')
    parser.add_argument('--min-transcripts', type=int, default=2,
                       help='Minimum number of transcripts per gene to consider (default: 2)')
    parser.add_argument('--verbose', action='store_true', help='Print verbose output')
    
    args = parser.parse_args()
    
    if not args.input.exists():
        print(f"Error: Input file {args.input} does not exist", file=sys.stderr)
        sys.exit(1)
    
    if args.verbose:
        print(f"Processing {args.input}")
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
        
        # Find constitutive exons
        constitutive_exons = parser_obj.find_constitutive_exons(genes, transcripts)
        
        if args.verbose:
            print(f"Found {len(constitutive_exons)} constitutive exons")
        
        # Find introns
        introns = parser_obj.find_introns(genes, transcripts)
        
        if args.verbose:
            print(f"Found {len(introns)} introns")
        
        # Write output files
        parser_obj.write_gtf(constitutive_exons, args.output_exons, 'exon')
        parser_obj.write_gtf(introns, args.output_introns, 'intron')
        
        if args.verbose:
            print(f"Wrote constitutive exons to {args.output_exons}")
            print(f"Wrote introns to {args.output_introns}")
    
    except Exception as e:
        print(f"Error processing file: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main() 
