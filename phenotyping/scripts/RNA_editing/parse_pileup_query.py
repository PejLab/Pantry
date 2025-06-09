#!/usr/bin/env python3
# library to parse a pileup file line and return counts of refnuc, A, T, G, C
# parse_pileup function requires a pileup line and minimum base quality as input
# Based on parse_pileup_query.pl from https://github.com/vsoch/mpileup

def parse_pileup(pileup_line, minbasequal, offset):
    """
    Parse a pileup line and return counts of reference nucleotide and A, T, G, C.
    
    Args:
        pileup_line (str): A tab-separated pileup line
        minbasequal (int): Minimum base quality threshold
        offset (int): Base quality offset
        
    Returns:
        tuple: Counts of (refnuc, A, T, C, G)
    """
    if len(pileup_line.split('\t')) < 6:
        raise ValueError("Invalid pileup line format")
    
    pilefields = pileup_line.strip().split('\t')
    pileup = list(pilefields[4])
    qualities = list(pilefields[5])
    
    acount = tcount = ccount = gcount = refnuccount = 0
    indel = 0
    readcount = 0
    indelstop = 0
    
    j = 0
    while j < len(pileup):
        if pileup[j] == '^':  # ignore read map qualities
            indel = 1
            indellength = 1
            indelstop = j + indellength
            
        if indel == 0:
            if pileup[j] in ['+', '-']:  # ignore indels
                indel = 1
                indellength = int(pileup[j+1]) + 1
                
                if j+2 < len(pileup) and pileup[j+2].isdigit():
                    indellength = (10 * int(pileup[j+1])) + int(pileup[j+2]) + 2
                    
                indelstop = j + indellength
                
        if indel == 1:
            if j == indelstop:
                indel = 0
                j += 1
                continue
                
        if indel == 0:  # count up the occurrence of each nucleotide
            if pileup[j] != '$':
                readcount += 1
                
                if ord(qualities[readcount-1]) >= (minbasequal + offset):
                    if pileup[j] in ['.', ',']:
                        refnuccount += 1
                    elif pileup[j].lower() == 'a':
                        acount += 1
                    elif pileup[j].lower() == 't':
                        tcount += 1
                    elif pileup[j].lower() == 'c':
                        ccount += 1
                    elif pileup[j].lower() == 'g':
                        gcount += 1
                        
        j += 1
        
    return (refnuccount, acount, tcount, ccount, gcount)

