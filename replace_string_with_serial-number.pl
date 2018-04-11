#!/usr/bin/perl
open MYINFILE, "3_finished.fasta" or die $!;
open (MYOUTFILE, '>3_finished_fixed.fa');

my $i = 1;

while (<MYINFILE>) {
	my($line) = $_;
	if($line =~ m/^>gi/) {
		$line =~ s/^>gi.*/>$i/;
		++$i;
	}
	print MYOUTFILE $line;
}

print "Found and replaced $i times\n";

close(MYOUTFILE);
close(MYINFILE);
