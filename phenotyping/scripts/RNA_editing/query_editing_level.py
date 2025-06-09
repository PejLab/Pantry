#!/usr/bin/env python3
# Quantify editing levels of known sites in a BAM file
# Based on Query_Editing_Level.GRCh37.20161110.pl from https://github.com/vsoch/mpileup

import argparse
import os
import subprocess
from parse_pileup_query import parse_pileup

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Query editing level of known sites in a BAM file')
    parser.add_argument('--edit_sites', required=True, help='Input file with editing sites')
    parser.add_argument('--ref', required=True, help='Reference genome path')
    parser.add_argument('--bam', required=True, help='BAM file path')
    parser.add_argument('--output', required=True, help='Output file path (can be *.gz)')
    args = parser.parse_args()

    # If .gz output is desired, get pre-zipped file name and gzip at the end
    outputfile = args.output
    if outputfile.endswith('.gz'):
        outputfile = outputfile[:-3]

    # Global variables
    minbasequal = 20  # MINIMUM BASE QUALITY SCORE
    minmapqual = 255  # MINIMUM READ MAPPING QUALITY SCORE. 255 FOR UNIQUE MAPPING WITH STAR
    sampath = "samtools"  # PATH TO THE SAMTOOLS EXECUTABLE
    offset = 33  # BASE QUALITY SCORE OFFSET - 33 FOR SANGER SCALE, 64 FOR ILLUMINA SCALE

    # Create temporary files
    bedtemp = f"{outputfile}.bed"
    piletemp = f"{outputfile}.pileup"

    # Create BED file from input
    with open(args.edit_sites, 'r') as infile, open(bedtemp, 'w') as outfile:
        for line in infile:
            fields = line.strip().split()
            if fields[0] != 'chromosome':
                outfile.write(f"{fields[0]}\t{fields[1]}\t{fields[2]}\n")

    # Run samtools mpileup
    mpileup_cmd = [
        sampath, "mpileup",
        "-A", "-B",
        "-d", "1000000",
        "-q", str(minmapqual),
        "-Q", str(minbasequal),
        "-f", args.ref,
        "-l", bedtemp,
        args.bam
    ]
    
    with open(piletemp, 'w') as outfile:
        subprocess.run(mpileup_cmd, stdout=outfile, check=True)

    # Process pileup file
    sitehash = {}
    with open(piletemp, 'r') as pileup_file:
        for line in pileup_file:
            fields = line.strip().split()
            chr_name, position, refnuc, coverage, pile, qual = fields
            location = f"{chr_name}_{position}"
            refnuccount, acount, tcount, ccount, gcount = parse_pileup(line, minbasequal, offset)
            counts = f"{refnuccount},{ccount},{gcount}"
            sitehash[location] = counts

    # Clean up temporary files
    os.remove(bedtemp)
    os.remove(piletemp)

    # Process input file and write output
    with open(args.edit_sites, 'r') as infile, open(outputfile, 'w') as outfile:
        outfile.write("#chrom\tposition\tcoverage\teditedreads\teditlevel\n")
        
        for line in infile:
            fields = line.strip().split()
            if fields[0] == 'chromosome':
                continue
                
            chr_name, position = fields[0], fields[2]  # 3rd column is 1-based coordinates
            location = f"{chr_name}_{position}"
            strand = fields[5]

            if location in sitehash:
                refcount, ccount, gcount = map(int, sitehash[location].split(','))
                newcov = newmismatch = 0
                
                if strand == '+':
                    newmismatch = gcount
                else:
                    newmismatch = ccount
                    
                newcov = refcount + newmismatch
                
                if newcov:
                    varfreq = f"{newmismatch/newcov:.3f}"
                    outfile.write(f"{fields[0]}\t{fields[2]}\t{newcov}\t{newmismatch}\t{varfreq}\n")
                else:
                    outfile.write(f"{fields[0]}\t{fields[2]}\t0\t0\tN/A\n")
            else:
                outfile.write(f"{fields[0]}\t{fields[2]}\t0\t0\tN/A\n")

    if args.output.endswith('.gz'):
        # Compress output file
        subprocess.run(['gzip', outputfile], check=True)

if __name__ == "__main__":
    main()
