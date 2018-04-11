#!/usr/bin/perl
use strict;
use warnings;

#######################################
# Parse fasta file and extract sequence
# frames with sliding window approach
#---------------20150515----------------
# Put each contig in one line
# Slide through sequence, extract frames
# Default frame size: 5000kb
# Default increment (step): 1nt
#---------------Structure---------------
# I.		Arguments handling
# II.		Global variables
# III.	Sub routines (functions)
# IV.		Main
########################################


########################################
# I.		Arguments handling
########################################

sub usage(){
	die("\n\tusage:\t".$0." [OPTION] -i <assembly_file_name.fa>\n\n",
	"\t[OPTION] may contain one or more of the following:\n",
	"\t-w <win size>\tsliding window size [nt] (default: 5000)\n",
	"\t-s <step>\tpositions to slide [nt] (default: 1)\n\n");
}
	
if(scalar @ARGV<1){
	usage();
}

my $win_size=5000;							# set default window size
my $step=1;											# set default number of nucleotides to move forward
my $assembly_file;							# initialize fasta file variable (input file)
#---------------------------------------
my ($YY,$MM,$DD,$hh,$mm)=(localtime)[5,4,3,2,1];
	$YY+=1900;										# YY starts at 1900 -> add 1900 years
	$MM++;												# MM starts at 0 -> add one month
	$MM=sprintf("%02d",$MM);
	$DD=sprintf("%02d",$DD);
	$hh=sprintf("%02d",$hh);
	$mm=sprintf("%02d",$mm);
my $timestamp=$YY.$MM.$DD.$hh.$mm;
#---------------------------------------
my $out_dir="winextract_".$timestamp;	# initialze default directory name for sequences

while(scalar @ARGV>0){					# read arguments (order does not matter)
	my $curr_val=shift @ARGV;
	if($curr_val=~/^-[wsi]$/){
		if($curr_val=~/-w/){
			$win_size=shift @ARGV;
		}elsif($curr_val=~/-s/){
			$step=shift @ARGV;
		}elsif($curr_val=~/-i/){
			$assembly_file=shift @ARGV;
		}
	}else{
		usage();
	}
}



########################################
# II.		Global variables
########################################

my $start_time=time();						# actual time (seconds since 1.1.1970)
my $exec_time=time()-$start_time;			# execution time of the script
my %seq_hash=();							# store sequence number (key) to sequence (value)
my $seq_value="";							# sequence collector (concatenate single lines of sequences in one)
my @seq_no_array=();						# store sequence numbers in array
my $seq_counter=0;							# counter for sequence numbers
my $key;									# keys in sequence hash
$|=1;										# turn on instant output to STDOUT
my $parse_line=0;							# user info (parsed lines counter)
#---------------------------------------
my $prog_name=substr $0, 2;					# get the program name without leading './'



########################################
# III.		Sub routines (functions)
########################################

sub parse_fasta(){							# parsing of the fasta file (join sequence lines and store in hash)
	$exec_time=time()-$start_time;			# execution time
	print("\n".$exec_time."s:\tparsing ".$assembly_file." ");
	open(INFILE,"<$assembly_file") or die "\tUnable to open ".$assembly_file.": ".$!."\n\n";
	while(my $inline=<INFILE>){				# read line by line
		chomp($inline);						# remove new line character
		$parse_line++;						# user info: parsed lines counter
		if($parse_line%1000==0){			# modulus (remainder of division) = 0
			print(".");
			if($parse_line%100000==0){
				print("\n\tread ".$parse_line." lines ");
			}
		}
		if($inline =~/^>.*?([0-9]+).*$/){	# ('?': don't be greedy); contains sequence number, memorize it ($1) ...
			$seq_no_array[$seq_counter]=$1;	# ... and write it into sequence number array
			$seq_counter++;					# increase array address counter
			if($seq_counter>1){				# only on second found sequence number:
				$seq_hash{$seq_no_array[$seq_counter-2]}=$seq_value;	# flush sequence(s) as value to hash (key=$seq_counter-2)
																		# (-1 for previous run, -1 for increased counter one line before)
				$seq_value="";				# reset sequence collector
			}
		}elsif($inline =~/^[ACGTUMRWSYKVHDBN]+$/i){	# see http://droog.gs.washington.edu/parc/images/iupac.html
			$seq_value=$seq_value.$inline;	# if it is not a sequence number, it must be part of the sequence itself
		}else{
			print("\nWarning: invalid fasta format \@ line ".$parse_line.", sequence ".$seq_no_array[$seq_counter-1].":\n".$inline."\n");
		}
	}
	$seq_hash{$seq_no_array[$seq_counter-1]}=$seq_value;	# Flush last collected sequence lines into last hash key-value pair
	close (INFILE);							# from here on we do not need to read the input file anymore
	$exec_time=time()-$start_time;			# execution time
	print("\n".$exec_time."s:\tdone parsing\n");
}


sub write_seq_files(){						# write parsed sequence files into sequence dir
	$exec_time=time()-$start_time;			# execution time
	print($exec_time."s:\twriting sequence files into ./".$out_dir."\n");
	if(!-d $out_dir){						# create directory for parsed sequences
		mkdir($out_dir,0777) || die "Unable to create ".$out_dir.": ".$!."\n\n";
	}
	for($seq_counter=0;$seq_counter<scalar(@seq_no_array);$seq_counter++){	# sorted keys already known, print out key-value pairs
		my @vars;																				# helper array, hold values from seq_hash
		## map { $vars[$seq_counter] = $seq_hash{$seq_counter+1} } keys %seq_hash;	# 20170614: don't know how this could even work before (on Nitrobin4_1.fa)
		map { $vars[$seq_counter] = $seq_hash{$seq_no_array[$seq_counter]} } keys %seq_hash;	# write values of seq_hash key to array
		my @vars2=split(//,$vars[$seq_counter]);				# split at every character and write to new array
		for(my $i=0;$i<scalar(@vars2)-$win_size+1;$i=$i+$step){	# through whole array-win_size
			print(sprintf("%.2f", $i*100/(scalar(@vars2)-$win_size))." %\r");
			my $window=sprintf("%07d", $i/$step);
			my $seq_file=$out_dir."/seq".$seq_no_array[$seq_counter]."_win".$win_size."_inc".$step."_fr".$window.".fa";
			open(OUTFILE, ">$seq_file") or die "Unable to write ".$seq_file.": ".$!."\n\n";
			print OUTFILE (">".$window."\n");
			for(my $j=$i;$j<$i+$win_size;$j++){		# from actual position in array ($i) to win_size
				print OUTFILE ($vars2[$j]);									# write each hash value into file, key = file name
			}
			print OUTFILE ("\n");
			close(OUTFILE);
		}
#		print substr ($seq_hash{$seq_no_array[$seq_counter]}, $step, $win_size."\n");	# write each hash value into file, key = file name
	}
	$exec_time=time()-$start_time;			# execution time
	print("\n".$exec_time."s:\tdone writing\n\n");
}



########################################
# III.		Main
########################################

parse_fasta();
write_seq_files();
