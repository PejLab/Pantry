#!/usr/bin/env perl
# Perl script for making edit level matrix
# Originally from https://github.com/vargasliqin/GTEx_edQTL

use strict;
use warnings;
use Getopt::Long;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use Cwd;

my ($path_edit_files, $samples_file, $output_file, $mincov, $minsamps);

# Parsing command-line options
GetOptions(
    'path_to_edit_files=s'    => \$path_edit_files,
    'samples_file=s'          => \$samples_file,
    'output_file=s'           => \$output_file,
    'min_coverage=i'          => \$mincov,
    'min_samples=i'           => \$minsamps,
) or die "Error in command line arguments: $!\n";

# Adding "/" to thr end of path if it's not there already
$path_edit_files .= "/" unless $path_edit_files =~ /\/$/;

print "Path to Edit Files: $path_edit_files\n";
print "Samples File: $samples_file\n";
print "Output File: $output_file\n";
print "Min Coverage: $mincov\n";
print "Min Samples: $minsamps\n";

# Read samples from file
my @samples;
open(my $SAMPLES, "<", $samples_file) or die "Could not open samples file '$samples_file': $!\n";
while(<$SAMPLES>) {
    chomp;
    push @samples, $_;
}
close $SAMPLES;

my %sitehash;
my %totalhash;
my %lvlhash;

foreach my $sample (@samples) {
    my $file = $path_edit_files . $sample . ".rnaeditlevel.tsv.gz";
    unless (-e $file) {
        print "Warning: File not found for sample $sample: $file\n";
        next;
    }
    print "Analyzing: $sample\n";
    
    my $INPUT = IO::Uncompress::Gunzip->new($file) or die "Failed to open $file: $GunzipError\n";
    while(<$INPUT>) {
        chomp;
        my @fields = split;
        my ($chr, $pos, $cov, $edit, $lvl) = ($fields[0], $fields[1], $fields[2], $fields[3], $fields[4]);

        next if ($chr eq '#chrom');

        my $site = join '_', $chr, $pos;
        my $ratio = join '/', $edit, $cov;
        if ($cov >= $mincov) {
            $sitehash{$sample}{$site} = $ratio;
            if ($totalhash{$site}) {
                $totalhash{$site}++;
                $lvlhash{$site} = join ',', $lvlhash{$site}, $ratio;
            } else {
                $totalhash{$site} = 1;
                $lvlhash{$site} = $ratio;
            }
        }
    }
    close $INPUT;
}

# Prepare output file
unlink $output_file if -e $output_file;
print "Output file: $output_file\n";
open(my $OUTPUT, ">", $output_file) or die "Could not open file '$output_file': $!\n";

print $OUTPUT "site";
foreach my $sample (@samples) {
	print $OUTPUT "\t$sample";
}
print $OUTPUT "\n";
foreach my $site (keys %totalhash) {
	if ($totalhash{$site} >= $minsamps) {
		my @lvls = split(/\,/, $lvlhash{$site});
		
			print $OUTPUT "$site";
			foreach my $sample (@samples) {
				if ($sitehash{$sample}{$site}) {
					print $OUTPUT "\t$sitehash{$sample}{$site}";
				} else {
					print $OUTPUT "\t0/0";
				}
			}
			print $OUTPUT "\n";
		
	}
}
close $OUTPUT;
