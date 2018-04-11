#!/usr/bin/perl
use strict;
use warnings;

####
# input file: txt files with alternating paired reads, no headers
# output files: fasta files with headers; first one with forward, second one with reverse reads
####

$| = 1;									# turn on instant output to STDOUT

# functional variables
my $infile = shift @ARGV;
my $outfile_f = shift @ARGV;
my $outfile_r = shift @ARGV;
# cosmetic variables
my $file_size = -s $infile;				# get file size once
my $size_fraction = sprintf("%d", ($file_size / 1000));	# every 0.1% step, do progress calculations
my $current_pos = 0;
my $progress = 0;

open MYINFILE, $infile or die $!;
open (MYOUTFILE_F, ">$outfile_f");
open (MYOUTFILE_R, ">$outfile_r");

my $i = 0;								# start with 0 - we have not read any line yet
my $j = 1;								# decoupling of line counter with read counter

while (<MYINFILE>) {
	my($line) = $_;
		++$i;
		if($i % 2 == 1){				# odd ("forward read")
			print MYOUTFILE_F ">".$j."_1\n";
			print MYOUTFILE_F $line;
		}elsif($i % 2 == 0){			# even ("reverse read")
			print MYOUTFILE_R ">".$j."_2\n";
			print MYOUTFILE_R $line;
			++$j;
		}

	# Progress counter, only calculate every 100th of file_size
	$current_pos = tell MYINFILE;		# actual byte position
	if($current_pos % $size_fraction == 0) {	# does not catch every 1000th of file size, as current position at end of current line not necessarily a multiple of size fraction
		# "print progress in percent and integer (current position in MB and integer of file size in MB and integer)"
		print ("Progress: ".(sprintf("%d", ($current_pos / $file_size * 100)))."% (".(sprintf("%d", ($current_pos / 1000000)))." of ".(sprintf("%d", ($file_size / 1000000)))." MB)\r");
	}

}

#print "\e[KReplaced $i headers\n";		# flush line before printing (see escape sequences/characters, \e)
print "\nReplaced $i headers\n";		# flush line before printing (see escape sequences/characters, \e)

close(MYOUTFILE_F);
close(MYOUTFILE_R);
close(MYINFILE);
