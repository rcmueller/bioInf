#!/usr/bin/perl
use strict;
use warnings;

#
# Iteratively split fasta seuqences after x nucleotides
# Append scaffold(s) with unique suffix(es) and write remaining sequence(s)
# 

#two arguments, length of region and input file
if (( $#ARGV +1 ) != 2 ){
  die "\n  Usage: ".$0." <length> <infile.fa>\n\n";
}

#read arguments, initialise variables
my $length = shift @ARGV;
my $infile = shift @ARGV;
my $scaffold;
my %seqHash;

#parse input file into hash
open FASTA, $infile or die "\n  Unable to open ".$infile.": ".$!."\n\n";
while ( <FASTA> ){
  if ( /^>(.*)$/ ) { $scaffold = $1; }
  elsif ( /^(\S+)$/ ) { $seqHash{$scaffold} .= $1 if $scaffold; }
}
close ( FASTA );

#while sequence for each scaffold is greater or equal than length of region
#print scaffold with suffix (incrementing integer) and first length characters of region
#assign remaining region to sequence string and increase suffix
foreach my $scaffold ( sort keys %seqHash ){
  my $suffix = 0;
  my $sequence = $seqHash{$scaffold};
  while ( length($sequence) > $length - 1 ){
    printf "%s.%s\n%.${length}s\n", $scaffold, $suffix, $sequence;
    $sequence = substr $sequence, $length ; 
    $suffix++;
  }
}
