#!/usr/bin/perl
use strict;
use warnings;

####
# Parse fasta file (first and only arumgent) with leading '>' sign for contig numbers
# Store contig number and calculate total number of nucleotides for each contig
# Write output file with extension '_stats.tsv' as tab-separated file for use in R or Calc
####
# in file: fasta file with leading '>' for header line
# out file: stats for contigs file:
# 1. column: contig number; 2. column: bp in contig
####

$| = 1;									# turn on instant output to STDOUT (for progress counter?)

# functional variables
my $infile = shift @ARGV;
my $outfile = $infile."_stats.tsv";
my $contig = "";
my $line_length = 0;
my $contig_length = 0;

my $total_bp = 0;

open MYINFILE, $infile or die $!;
open (MYOUTFILE, ">$outfile");
print MYOUTFILE "contig\tlength\n";		# write header

my $i = 0;								# start with 0 - we have not read anything yet

while (<MYINFILE>) {
	my($line) = $_;
	chomp($_);
	if($line =~ m/^>/){					# contig number
		if($contig ne ""){				# write 1. contig and 1. contig_length to outfile on 2. contig and so on ...
			print MYOUTFILE $contig."\t".$contig_length."\n";
			$total_bp = ($total_bp + $contig_length);
			$contig_length = 0;			# flush $contig_length (and calculate sum of bp, before doing so)
		}
		$contig = (split />/,$_)[1];	# split input line at regex '>', assign second entry ([1]) - the contig number - to $contig
		++$i;
	}else{
		$line_length = length($_);		# read line length into variable
		$contig_length = ($contig_length + $line_length);
	}
}

$total_bp = ($total_bp + $contig_length);			# correct total number of nucleotides (last contig)
print MYOUTFILE $contig."\t".$contig_length."\n";	# write last contig

print "\nFile contains ".$i." contigs and a total of ".$total_bp." nucleotides\n\n";

close(MYOUTFILE);
close(MYINFILE);
