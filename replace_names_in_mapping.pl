#!/usr/bin/perl
#use Text::CSV;    ## wanted to try the csv-parser for that problem (to change field #7) but change #1 instead

open MYINFILE, "Soil_Microbs_C4rs_mapping.csv" or die $!;
#open MYINFILE, "Soil_Microbs_S4rs_mapping.csv" or die $!;
open (MYOUTFILE, '>data.txt');

my $i = 0;       # leave out header line and start with "1" in data lines

while (<MYINFILE>) {
	my($line) = $_;
	if($line !=~ m/^\"Name.*/) {      # ignore header line
		$line =~ s/^\"[0-9]*\"/\"$i\"/;
		$line =~ s/^\"[0-9]* mapping\"/\"$i\"/;
		++$i;
	}
	print MYOUTFILE $line;
}

print "Found and replaced $i times\n";

close(MYOUTFILE);
close(MYINFILE);
