#!/usr/bin/perl
use strict;
use warnings;

###########################################################
# Parse MEGAN output file (csv; read-name, taxon-path)
# Extract OTU (scaffold_ORF) and Phylum tab separated
# Extend Proteobacteria to Class level and cleanup
# messy MEGAN output (Superphylum, subphylum, etc.)
###########################################################

# Arguments check
if (($#ARGV + 1) != 2){
	die "\n\tUsage: ".$0." <input_file> <output_file>\n\n";
}

# Globals
my $infile = shift @ARGV;
my $outfile = shift @ARGV;
my %taxtable;

# File handling
open IN, $infile or die "\n\tUnable to open ".$infile.": ".$!."\n\n";
open (OUT, ">$outfile") or die "\n\tUnable to write ".$outfile.": ".$!."\n\n";

# Core
while (my $line = <IN>){
	chomp $line;
	$line =~ s/\"//g;							#Remove quotation marks
	$line =~ s/ <phylum>//g;					#Remove nasty spelling
	$line =~ s/unclassified/Unclassified/g;		#Correct nasty spelling
	my @splitline = split(/\t/,$line);			#Separate line by tab (read-name \t taxon-path)
	my @tax = split(/;/,$splitline[1]);			#Separate taxon-path by semicolons
	if (exists($tax[3])){						#Only consider classified
		my $superphyla = "(.*?\/.*? group|Alveolata|Euglenozoa|Viridiplantae|Proteobacteria)";
		if ($tax[3] =~ m/$superphyla/){			#Discard superphyla columns and extend Proteobacteria
			if (exists($tax[4])){				#Only consider classified
				my $subphylum = "delta\/epsilon subdivisions";
				if ($tax[4] =~ m/$subphylum/){	#Discard subphylum columns
					if (exists($tax[5])){		#Only consider classified
						$taxtable{$splitline[0]} = $tax[5];
					}else{
						$taxtable{$splitline[0]} = "NA";
					}
				}else{
				$taxtable{$splitline[0]} = $tax[4];
				}
			}else{
				$taxtable{$splitline[0]} = "NA";
			}
		}else{
			$taxtable{$splitline[0]} = $tax[3];
		}
	}else{
		$taxtable{$splitline[0]} = "NA";
	}
}

# Generate output file
print OUT "\tphylum_class\n";

## Returns hash in random order
#foreach my $key (keys %taxtable) {
#   print OUT $key."\t".$taxtable{$key}."\n";
#}

### Schwartzian transform: (http://stackoverflow.com/questions/21871337/sorting-numeric-hash-keys-separated-by-hyphen)
# From the end, we 1) first store the original string, plus the first and second number inside an anonymous array ref.
# The result is a list of array refs -- a cache -- which we 2) pass on to sort, where they are sorted based first on the
# first number, and if they are the same, on the second. This is achieved by using || inside the sort code block.
# Lastly we 3) restore the original string and discard the array refs.

my @sorted_keys =	map  { $_->[0] }									#3)
					sort { $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2] }	#2)
					map  { [$_, /(\d+)/g] } keys %taxtable;				#1) (\d refers to digits)

# Returns hash in sorted order
foreach my $key(@sorted_keys) {
    print OUT $key."\t".$taxtable{$key}."\n";
}

close IN;
close OUT;
