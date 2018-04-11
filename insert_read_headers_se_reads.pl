#!/usr/bin/perl
use strict;
use warnings;

####
# input file: fasta file without headers ("No name")
# output files: fasta files with headers (single end)
####

# functional variables
my $infile = shift @ARGV;
my $outfile = shift @ARGV;
# cosmetic variables
my $file_size = -s $infile;				# get file size once
my $size_fraction = sprintf("%d", ($file_size / 1000));	# every 0.1% step, do progress calculations
my $current_pos = 0;
my $progress = 0;

open MYINFILE, $infile or die $!;
open (MYOUTFILE, ">$outfile");

my $i = 0;								# start with 0 - we have not read any line yet

while (<MYINFILE>) {
	my($line) = $_;
	if($line =~ m/^>No_name/) {
		++$i;
		$line =~ s/^>No_name/>$i/;
	} 
	print MYOUTFILE $line;
	# Progress counter, only calculate every 100th of file_size
	$current_pos = tell MYINFILE;		# actual byte position
	if($current_pos % $size_fraction == 0) {	# does not catch every 1000th of file size, as current position at end of current line not necessarily a multiple of size fraction
		# "print progress in percent and integer (current position in MB and integer of file size in MB and integer)"
		print STDERR ("Progress: ".(sprintf("%d", ($current_pos / $file_size * 100)))."% (".(sprintf("%d", ($current_pos / 1000000)))." of ".(sprintf("%d", ($file_size / 1000000)))." MB)\r");
	}

}

#print "\e[KReplaced $i headers\n";		# flush line before printing (see escape sequences/characters, \e)
print "\nReplaced $i headers\n";		# flush line before printing (see escape sequences/characters, \e)

close(MYOUTFILE);
close(MYINFILE);
