#!/usr/bin/perl
use strict;
use warnings;

# 20171127, v0.3, write output file, e.g., 20171127-1202_ss_10_1000.txt
#                 fix: range includes start and stop position -> $l_region-1
#                 corr: var name $start -> $first_possible, $end -> $last_possible
# 20171123, v0.2, random loci and range sampling (Fisher-Yates)
# 20171122, v0.1, file parsing and format checking
#################################################################################
# Stochastic filtering of common mapping coordinate regions file
# Input file format (user input; without []): [scaffold]\t[begin]\t[end]
# Filter for loci containing >= 1000 (or user input) positions
# Random sample 10 (or user input) loci
# Of these, random sample [r_begin] with
#   [r_begin] >= [begin]
#   [r_begin] <= ([end] - 1000)
# Write stochastic sampling output: ss_10_1000.txt (numbers depend on user input)
##################################################################################


# Usage
sub usage(){
	die ("\n\tUsage: ".$0." <number of loci> <length of region> <input file>\n\n");
}


# Argument check
if (($#ARGV+1) != 3){
	usage();
}


# Globals
my $n_loci = shift @ARGV;
my $l_region = shift @ARGV;
my $infile = shift @ARGV;
my $pattern = qr/^.*?\t\d+?\t\d+?$/;
my %regions = ();
#---------------------------------------
my ($YY,$MM,$DD,$hh,$mm)=(localtime)[5,4,3,2,1];
    $YY+=1900;                    # YY starts at 1900 -> add 1900 years
    $MM++;                        # MM starts at 0 -> add one month
    $MM=sprintf("%02d",$MM);
    $DD=sprintf("%02d",$DD);
    $hh=sprintf("%02d",$hh);
    $mm=sprintf("%02d",$mm);
my $timestamp=$YY.$MM.$DD."-".$hh.$mm;
#---------------------------------------
my $outfilename=$timestamp."_ss_".$n_loci."_".$l_region.".txt";


# Debug
#print "\ninput:\t".$infile."\nloci:\t".$n_loci."\nlength:\t".$l_region."\n\n";


# Parse input file
open (INFILE,"<$infile") or die "\n\tUnable to open ".$infile.": ".$!."\n\n";
while(my $inline=<INFILE>){
  chomp($inline);
  if($inline =~ $pattern){  #$inline corresponds to [scaffold]\t[begin]\t[end], or throw warning
    my @range = split(/\t/, $inline);
    my $length = $range[2] - $range[1];
    if ($length >= $l_region-1){  #keep lines with [begin]-[end]>=$l_region-1
      $regions{$.}=[@range];
    }
  }else{
    print "\n\tWARNING: [line ".$.."] cannot parse [".$inline."]\n";
  }
}
close(INFILE);


# Random sample $n_loci keys
my @keys=(keys %regions);
if ($#keys < $n_loci){
  die "\n\tERROR: number of required loci [".$n_loci."] exceeds number of available loci [".$#keys."]\n\n";
}else{
  fisher_yates_shuffle(\@keys);
  @keys = @keys[$#keys-($n_loci-1) .. $#keys];
}


# Schwartzian transform (https://stackoverflow.com/questions/21871337/sorting-numeric-hash-keys-separated-by-hyphen)
my @sorted_keys=map{$_->[0]}
                sort{$a->[1] <=> $b->[1] || $a->[2] <=> $b->[2]}
                #map{[$_, /(\d+)/g]} keys %regions;
                map{[$_, /(\d+)/g]} @keys;
 

# Random sample $l_region from random $n_loci
open(OUTFILE, ">$outfilename") or die "\n\tUnable to write ".$outfilename.": ".$!."\n\n";
foreach my $entry (@sorted_keys){
  print "\n\trandom locus ($entry):\t\t@{$regions{$entry}}\n";
  my $first_possible=$regions{$entry}[1];
  my $last_possible=$regions{$entry}[2]-($l_region-1);
  #Debug
  print "\tpotential starting points:\t".$first_possible." ".$last_possible."\n";
  my @limited_range=($first_possible..$last_possible);
  #Debug
  #print join(", ", @limited_range)."\n";
  fisher_yates_shuffle(\@limited_range);
  my $random_start=pop(@limited_range);
  print "\trandomised range with l=".$l_region.":\t".$random_start." ".($random_start+$l_region-1)."\n";
  print OUTFILE ($regions{$entry}[0]."\t".$random_start."\t".($random_start+$l_region-1)."\n");
}
close (OUTFILE);
print "\n\n\tGenerated output file: ".$outfilename."\n\n";


# Fisher-Yates-Shuffle: https://stackoverflow.com/questions/8963228/how-can-i-take-n-elements-at-random-from-a-perl-array#8964210
sub fisher_yates_shuffle{
  my $to_shuffle = shift;
  my $i = @$to_shuffle;
  while (--$i){
    my $j = int rand($i+1);
    @$to_shuffle[$i,$j] = @$to_shuffle[$j,$i];
  }
}
