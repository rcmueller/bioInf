#!/usr/bin/perl
use strict;
use warnings;

####
# input file: txt files with alternating paired headers and reads, usually derived from CLC-exported sam file via
# grep '^read_pair' <input_file>.sam | awk '{print $1,$10}' | sed 's|read_pair_|>|' | sed 's| |\n|' > <output_file>_pe.txt
# output files: fasta files with same headers; first file with forward, second one with reverse reads
####

$| = 1;									# turn on instant output to STDOUT

# functional variables
my $infile = shift @ARGV;
my $outfile_f = shift @ARGV;
my $outfile_r = shift @ARGV;
my $fw_rv = "";						# switch for forward and reverse reads
my $header = "";

open MYINFILE, $infile or die $!;
open (MYOUTFILE_F, ">$outfile_f");
open (MYOUTFILE_R, ">$outfile_r");

my $i = 0;								# start with 0 - we have not read anything yet

while (<MYINFILE>) {
	my($line) = $_;
	if($line =~ m/^>/){				# read header is a trigger for fw-rv switch
		if ($line eq $header){
			print MYOUTFILE_R $line;
			$fw_rv = "R";
		}else{
			++$i;					# ... first PE read
			$header = $line;
			print MYOUTFILE_F $line;
			$fw_rv = "F";
		}
	}else{
		if ($fw_rv =~ m/F/){
			print MYOUTFILE_F $line;
		}elsif($fw_rv =~ m/R/){
			print MYOUTFILE_R $line;
		}
	}
}

print "\nSplit $i PE reads\n";		# flush line before printing (see escape sequences/characters, \e)

close(MYOUTFILE_F);
close(MYOUTFILE_R);
close(MYINFILE);
