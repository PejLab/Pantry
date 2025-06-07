#!/usr/bin/env perl
# library to parse a pileup file line and return counts of refnuc, A, T, G, C
# parse_pileup subroutine requires a pileup line and minimum base quality as input
# Originally from https://github.com/vsoch/mpileup


use strict;

sub parse_pileup {

if (@_ < 3) {

	die "need to provide 2 inputs to parse_pileup: pileup line, minbasequal, base quality offset\n"

}

my @pilefields = split(/\t/, $_[0]);

my $minbasequal = int($_[1]);

my $offset = int($_[2]);

my @pileup = split(//,$pilefields[4]);

my @qualities = split(//, $pilefields[5]);

my ($acount, $tcount, $ccount, $gcount, $refnuccount) = (0,0,0,0,0);

my $indel = 0;

my $readcount = 0;

my $indelstop;	

LOOP: for (my $j = 0; $j < @pileup; $j++) {  #loop over each individual read and assay the nucleotide called

	if ($pileup[$j] eq '^') { #ignore the read map qualities

		$indel = 1;

		my $indellength = 1;

		$indelstop = $j + $indellength;

	}

	if ($indel == 0) {		

		if ($pileup[$j] eq '+' || $pileup[$j] eq '-') { #ignore indels

			$indel = 1;						

			my $indellength = int($pileup[$j+1]) + 1;

			if ($pileup[$j+2] =~ m/\d/) {

				$indellength = (10*int($pileup[$j+1])) + int($pileup[$j+2]) + 2; 

			}

			$indelstop = $j + $indellength;

		}

	}

	if ($indel == 1) {

		if ($j == $indelstop) {

			$indel = 0;

			next LOOP;

		}

	}			

	if ($indel == 0) {	#count up the occurance of each nucleotide

		if ($pileup[$j] ne '$') {

			$readcount = $readcount + 1;			

			if (ord($qualities[$readcount-1]) >= ($minbasequal + $offset)) {

				if ($pileup[$j] eq '.' || $pileup[$j] eq ',') {

					$refnuccount = $refnuccount + 1;

				}

				elsif ($pileup[$j] =~ /a/i) {

					$acount = $acount + 1;

				}

				elsif ($pileup[$j] =~ /t/i) {

					$tcount = $tcount + 1;

				}

				elsif ($pileup[$j] =~ /c/i) {

					$ccount = $ccount + 1;

				}

				elsif ($pileup[$j] =~ /g/i) {

					$gcount = $gcount + 1;

				}

			}

		}

	}

}

my @returnarray = ($refnuccount, $acount, $tcount, $ccount, $gcount); #return counts

return @returnarray;

}

1;
