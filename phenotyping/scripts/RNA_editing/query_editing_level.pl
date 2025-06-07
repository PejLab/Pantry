#!/usr/bin/env perl
# Perl script queries editing level of known sites in a BAM file
# Originally from https://github.com/vsoch/mpileup

use warnings;
use strict;
use Getopt::Long;
use Cwd;
require "./scripts/RNA_editing/parse_pileup_query.pl"; # NEED PARSE PILEUP LIBRARY


# Defining variables to hold the arguments
my ($inputfile, $genomepath, $bamfile, $outputfile);

# Parsing command-line options
GetOptions(
    'edit_sites=s'    => \$inputfile,
	'ref=s'	  		  => \$genomepath,
    'bam=s'			  => \$bamfile,
	'output=s'    	  => \$outputfile,
) or die "Error in command line arguments: $!\n";

# Added to handel gzip if running by snakemake (also won't make any problem for normal runs). It removes the ".gz" extension.
if ($outputfile =~ /\.gz$/) {$outputfile =~ s/\.gz$//;}; # Added by Me


# GLOBAL VARIABLES - PLEASE MODIFY THESE
my $minbasequal = 20; # MINIMUM BASE QUALITY SCORE
my $minmapqual = 255; # MINIMUM READ MAPPING QUALITY SCORE. 255 FOR UNIQUE MAPPING WITH STAR. ≥1 for reads mapped to less than 10 loci
my $sampath = "samtools"; # PATH TO THE SAMTOOLS EXECUTABLE
my $offset = 33; # BASE QUALITY SCORE OFFSET - 33 FOR SANGER SCALE, 64 FOR ILLUMINA SCALE

## END GLOBAL VARIABLES

my $bedtemp = join '', $outputfile, '.bed';
system("awk \'\$1\!\=\"chromosome\"\{print \$1\"\t\"\$2\"\t\"\$3\}\' $inputfile \> $bedtemp");
my $piletemp = join '', $outputfile, '.pileup';
system("$sampath mpileup -A -B -d 1000000 -q $minmapqual -Q $minbasequal -f $genomepath -l $bedtemp $bamfile \> $piletemp");

my %sitehash;
open (my $PILEUP, "<", $piletemp);
while(<$PILEUP>) {
	chomp;
	my ($chr, $position, $refnuc, $coverage, $pile, $qual) = split;
	my $location = join '_', $chr, $position;
	my ($refnuccount, $acount, $tcount, $ccount, $gcount) = &parse_pileup($_, $minbasequal, $offset);# parse each line of pileup
	my $counts = join ',', $refnuccount, $ccount, $gcount;
	$sitehash{$location} = $counts;
}
system("rm $bedtemp");
system("rm $piletemp");

open (my $INPUT , "<", $inputfile) or die "error opening inputfile: $!\n";
open (my $OUTPUT, ">", $outputfile);
print $OUTPUT "#chrom\tposition\tcoverage\teditedreads\teditlevel\n";

while (<$INPUT>) { # READ IN LIST OF KNOWN EDITED SITES AND QUERY EDITING STATUS
	chomp;
	my @fields = split;
	next if ($fields[0] eq 'chromosome');
	my ($chr, $position) = ($fields[0], $fields[2]); # 3rd column is 1-based coordinates
	my $location = join '_', $chr, $position;
	my ($strand) = ($fields[5]);

	if ($sitehash{$location}) { # PRINT OUT RESULT
		my ($refcount, $ccount, $gcount) = split(/\,/,$sitehash{$location});
		my ($newcov, $newmismatch) = (0,0);
		if ($strand eq '+') {
			$newmismatch = $gcount;
		} else {
			$newmismatch = $ccount;
		}
		$newcov = $refcount + $newmismatch;
		if ($newcov) {		
			my $varfreq = 0;
			$varfreq = sprintf("%.3f", $newmismatch/$newcov);
			print $OUTPUT "$fields[0]\t$fields[2]\t$newcov\t$newmismatch\t$varfreq\n";
		} else {
			print $OUTPUT "$fields[0]\t$fields[2]\t0\t0\tN/A\n";
		}
	} else {
		print $OUTPUT "$fields[0]\t$fields[2]\t0\t0\tN/A\n";
	}
}
system("gzip $outputfile") == 0 or die "gzip failed: $?\n";
close $INPUT;	
close $OUTPUT;

### END