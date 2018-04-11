#!/usr/bin/perl -w
# count_fasta.pl
# AUTHOR: Joseph Fass (modified from script by Brad Sickler)
# LAST REVISED: November 2010
# 
# The Bioinformatics Core at UC Davis Genome Center
# http://bioinformatics.ucdavis.edu
# Copyright (c) 2009 The Regents of University of California, Davis Campus.
# All rights reserved.
####
# Grabbed 20160105 from http://wiki.bioinformatics.ucdavis.edu/index.php/Count_fasta.pl
# Edited by Ralf Mueller, 20160105 (including N5, N10, N90, N95, min, max)
# Edit 20160107: commented the histogram part, prepared for multiple fasta files, optimised
# output for use in spread sheet or R (produce header and csv-file)
####
use strict;
use POSIX;
use Getopt::Std;

my $usage = "\n\nusage: $0 [-i <interval size in # of residues>] <fasta file(s)>\n".
	"Produces a histogram of sequence lengths and other measures to standard out.\n\n".
	"-i #          specify bin size for histogram (default 100)\n\n";
our($opt_i); # histogram interval
getopts('i:') or die $usage;
if (!defined($opt_i) or !($opt_i =~ m/^[0-9]+$/)) {$opt_i = 100;}

if( ( $#ARGV + 1 ) < 1 ) {
	die $usage;
}

# Read in sequences from one or more fasta files
my @data_files;
for(my $i = 0; $i < ($#ARGV + 1); $i++){
	$data_files[$i] = $ARGV[$i];
}
print "file,tot.length,tot.number,min.length,max.length,N5.count,N5.length,N10.count,N10.length,N25.count,N25.length,".
		"N50.count,N50.length,N75.count,N75.length,N90.count,N90.length,N95.count,N95.length,gc.length,gc.percent";
foreach my $file (@data_files){
	my $Id;
	my %seq;
	open(FASTA, $file) or die"Can't open file $file\n";
	while (<FASTA>) {
		if (/^>(.*)$/)  { $Id = $1; }						# $1 ... first match of previous regex pattern
		elsif (/^(\S+)$/)	{ $seq{$Id} .= $1 if $Id; }		# \S ... non-whitespace character, one or more; .= ... concatinating assignment
	}
	close (FASTA);

	# Count the number of sequences in the file and create a histogram of the distribution
	my $n = 0;
	my $int = 0;
	my $totalLength = 0;
	my $gcCount = 0;
	my %len = ();
	my @seqLengths = ();
	foreach my $id (keys %seq) {
		push @seqLengths, length($seq{$id}); # record length for N50 calc's
		$n++;
		$int = floor( length($seq{$id})/$opt_i );
		$totalLength += length($seq{$id});
		$gcCount += ($seq{$id} =~ tr/gGcC/gGcC/);
		if( !defined($len{$int}) ) {
			$len{$int} = 1;
		} else {
			$len{$int}++;
		}
	}

	# Calculate N5, N10, N25, N50, N75, N90 and N95 and counts
	my $N5; my $N10; my $N25; my $N50; my $N75; my $N90; my $N95;
	my $N5count=0; my $N10count=0; my $N25count=0; my $N50count=0; my $N75count=0; my $N90count=0; my $N95count=0;
	my $frac_covered = $totalLength;
	my $min_length; my $max_length;
	@seqLengths = reverse sort { $a <=> $b } @seqLengths;
	$min_length = $seqLengths[-1];
	$max_length = $seqLengths[0];

	$N5 = $seqLengths[0];
	while ($frac_covered > $totalLength*95/100) {
		$N5 = shift(@seqLengths);
		$N5count++; $N10count++; $N25count++; $N50count++; $N75count++; $N90count++; $N95count++;
		$frac_covered -= $N5;
	}
	$N10 = $N5;
	while ($frac_covered > $totalLength*90/100) {
		$N10 = shift(@seqLengths);
		$N10count++; $N25count++; $N50count++; $N75count++; $N90count++; $N95count++;
		$frac_covered -= $N10;
	}
	$N25 = $N10;
	while ($frac_covered > $totalLength*75/100) {
		$N25 = shift(@seqLengths);
		$N25count++; $N50count++; $N75count++; $N90count++; $N95count++;
		$frac_covered -= $N25;
	}
	$N50 = $N25;
	while ($frac_covered > $totalLength*50/100) {
		$N50 = shift(@seqLengths);
		$N50count++; $N75count++; $N90count++; $N95count++;
		$frac_covered -= $N50;
	}
	$N75 = $N50;
	while ($frac_covered > $totalLength*25/100) {
		$N75 = shift(@seqLengths);
		$N75count++; $N90count++; $N95count++;
		$frac_covered -= $N75;
	}
	$N90 = $N75;
	while ($frac_covered > $totalLength*10/100) {
		$N90 = shift(@seqLengths);
		$N90count++; $N95count++;
		$frac_covered -= $N90;
	}
	$N95 = $N90;
	while ($frac_covered > $totalLength*5/100) {
		$N95 = shift(@seqLengths);
		$N95count++;
		$frac_covered -= $N95;
	}

	# Print out the results
	print "\n";
	my @ints = sort { $a <=> $b } keys(%len);
#	for(my $i=$ints[0]; $i <= $ints[-1]; $i++) {
#		$len{$i} = 0 if(!defined($len{$i}));
#		printf "%d:%d \t$len{$i}\n", ( ($i*$opt_i), ($i*$opt_i+$opt_i-1) );
#	}
	printf $file.",".$totalLength.",".$n.",".$min_length.",".$max_length.",".
		$N5count.",".$N5.",".$N10count.",".$N10.",".$N25count.",".$N25.",".$N50count.",".$N50.",".
		$N75count.",".$N75.",".$N90count.",".$N90.",".$N95count.",".$N95.",".$gcCount.",%.2f",($gcCount/$totalLength * 100);
}
print "\n"
