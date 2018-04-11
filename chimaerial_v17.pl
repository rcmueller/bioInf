#!/usr/bin/perl
use strict;
use warnings;
use threads;

#######################################
# Parse fasta file (one sequence/line)
# Calculate GC-stats
#---------------20141211----------------
# Improvement of progress calc in '-d'
# mode (did not work in v16 on hpc2)
# Error fixing in regex (lower case)
# Error fixing $outfile creation if '-d'
#---------------Structure---------------
# I.		Arguments handling
# II.		Global variables
# III.		Sub routines (functions)
# IV.		Main
########################################

########################################
# I.		Arguments handling
########################################
sub usage(){
	die("\n\tusage:\t".$0." [OPTION] -i <assembly_file_name.fa>\n",
	"\tor:\t".$0." [OPTION] -d <sequence_file_directory>\n\n",
	"\t[OPTIONS] may contain one or more of the following:\n",
	"\t-w <window size>\tsize of the sliding window\n",
	"\t-n <number of windows>\tnumber of windows for standard deviation calculation\n",
	"\t-t <stdev threshold>\tthreshold [%] for maximum allowed standard deviation in sequence\n",
	"\t-l \t\t\twrite logfile(s)\n",
	"\t-v \t\t\tbe verbose (produces more output)\n\n");
}
	
if(scalar @ARGV<1){
	usage();
}

my $win_size=32;							# set default window size
my $no_of_win=5;							# set default number of windows for stdevp calculation
my $gc_thresh=0.05;							# set default threshold of stdevp
my $assembly_file;							# initialize fasta file variable (input file)
#---------------------------------------
my ($YY,$MM,$DD,$hh,$mm)=(localtime)[5,4,3,2,1];
	$YY+=1900;								# YY starts at 1900 -> add 1900 years
	$MM++;									# MM starts at 0 -> add one month
	$MM=sprintf("%02d",$MM);
	$DD=sprintf("%02d",$DD);
	$hh=sprintf("%02d",$hh);
	$mm=sprintf("%02d",$mm);
my $timestamp=$YY.$MM.$DD.$hh.$mm;
#---------------------------------------
my $seq_dir="gc_chimaerial_".$timestamp;	# initialze default directory name for sequences
my $file_parse="off";						# set default for file parsing mode
my $dir_parse="off";						# set default for directory parsing mode
my $logging="off";							# set default for logging
my $verbose="off";							# set default for verbose output

while(scalar @ARGV>0){						# read arguments (order does not matter)
	my $curr_val=shift @ARGV;
	if($curr_val=~/^-[wntidlv]$/){
		if($curr_val=~/-w/){
			$win_size=shift @ARGV;
		}elsif($curr_val=~/-n/){
			$no_of_win=shift @ARGV;
		}elsif($curr_val=~/-t/){
			$gc_thresh=(shift @ARGV)/100;	# percent -> fraction
		}elsif($curr_val=~/-l/){
			$logging="on";
		}elsif($curr_val=~/-v/){
			$verbose="on";
		}elsif($curr_val=~/-i/){
			$assembly_file=shift @ARGV;
			$file_parse="on";				# turn on file parsing
			$dir_parse="on";				# turn on diretory parsing
		}elsif($curr_val=~/-d/){
			$seq_dir=shift @ARGV;
			$dir_parse="on";				# turn on directory parsing
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
my $num_of_threads=`nproc`;					# get the number of processing units available
my $prog_name=substr $0, 2;					# get the program name without leading './'
my $chimalarm="off";						# default for chimera detection (needed in gc_stats for file handling)
my $temp_chim_dir="temp_".$timestamp;		# dir for temp. chimera output files

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
	print($exec_time."s:\twriting sequence files into ./".$seq_dir);
	if(!-d $seq_dir){						# create directory for parsed sequences
		mkdir($seq_dir,0777) || die "Unable to create ".$seq_dir.": ".$!."\n\n";
	}
	for($seq_counter=0;$seq_counter<scalar(@seq_no_array);$seq_counter++){	# sorted keys already known, print out key-value pairs
		my $seq_file=$seq_dir."/".$seq_no_array[$seq_counter];
		open(OUTFILE, ">$seq_file") or die "Unable to write ".$seq_file.": ".$!."\n\n";
		print OUTFILE ($seq_hash{$seq_no_array[$seq_counter]}."\n");	# write each hash value into file, key = file name
		close(OUTFILE);
	}
	$exec_time=time()-$start_time;			# execution time
	print("\n".$exec_time."s:\tdone writing\n");
}

sub initThreads(){							# initialize threading
	my @initThreads;
	for(my $thr_counter = 1;$thr_counter<=$num_of_threads;$thr_counter++){
		push(@initThreads,$thr_counter);
	}
	return @initThreads;
}

sub gc_stats(){								# calculate GC-statistics for each sequence file
	my $seq_stats_counter=1;				# user info; increment 1 (every sequence file)
	my $seq_stats_progress=0;				# user info; increment x% (up to approximately 100%)
	my $id=threads->tid();
	my $log_dir="gc_stats_log";
	my $gc_log=$log_dir."/".$timestamp."_pid_".$$."_thread_".$id."-".$prog_name.".log";
	if($file_parse eq "off"){				# work around if in '-d' mode as we do not know @seq_no_array
		opendir(INDIR,$seq_dir) or die "\nUnable to open ".$seq_dir.": ".$!."\n\n";	# get files
		@seq_no_array=readdir(INDIR);	# also includes '.', '..' and sub dirs, but close enough for progress calc ;)
		closedir(INDIR);
	}
	opendir(INDIR,$seq_dir) or die "\nUnable to open ".$seq_dir.": ".$!."\n\n";# open seq_dir holding the sequences
	while(defined(my $seq_file=readdir(INDIR))){				# loop through every file
		my $seq_file_path=$seq_dir."/".$seq_file;
		my $gc_stats_file=$seq_dir."/".$timestamp."_gc_stats"."/".$seq_file."_stats";
		my $gc_stats_dir=$seq_dir."/".$timestamp."_gc_stats";
		if(!-d $gc_stats_dir){
			mkdir($gc_stats_dir,0777) || warn "\tWarning: unable to create ".$gc_stats_dir.": ".$!."\n\n";
		}
		if($seq_file=~/^[0-9]+$/ && !-e $gc_stats_file){		# only seq_file and no '_stats' extension, obsolete with v12?
			open(INFILE,"<$seq_file_path") or die "\nUnable to open ".$seq_file_path.": ".$!."\n\n";
			my $inline=<INFILE>;
			chomp($inline);
			open(SEQSTATS,">>$gc_stats_file") or die "\nUnable to open ".$gc_stats_file.": ".$!."\n\n"; # gc_stats file
			if($logging eq "on"){								# do we want logging?
				if(!-d $log_dir){
					mkdir($log_dir,0777) || die "Unable to create ".$log_dir.": ".$!."\n\n";
				}
				open(LOGFILE,">>$gc_log") or die "\nUnable to open ".$gc_log.": ".$!."\n\n";		# debug/log file
				print LOGFILE ("\n\nProcessing ".$seq_file_path);
			}
			if($verbose eq "on"){
				print ("\t... processing ".$seq_file_path."\n");
			}
			my @nucleotide=split(//,$inline);					# read single characters into nucleotide array;
			my $win_counter=0;
			my $nuc_counter=0;
			my $divisor;
			my $gc=0;
			my $gc_ratio;
			my @gc_stats_array=();
			while($win_counter<scalar(@nucleotide)+1-$win_size){# step (increment=1) through sequence until EOF-win_size
				while($nuc_counter<$win_size+$win_counter){		# step (increment=1) through window until win_size
					$divisor=$win_size;
					if($nucleotide[$nuc_counter]=~/[ACGTUMRWSYKVHDB]/i){# valid nucleotides w/o N (IUPAC ambiguity code)
						if($nucleotide[$nuc_counter]=~/[CGS]/i){		# C or G or (C|G)
							$gc++;
						}elsif($nucleotide[$nuc_counter]=~/[MRYK]/i){	# (A|C) or (A|G) or (C|T) or (G|T)
							$gc=$gc+(0.5);
						}elsif($nucleotide[$nuc_counter]=~/[VB]/i){		# (A|C|G) or (C|G|T)
							$gc=$gc+(2/3);
						}elsif($nucleotide[$nuc_counter]=~/[HD]/i){		# (A|C|T) or (A|G|T)
							$gc=$gc+(1/3);
						}												# no increase of $g for A or T or U or W=(A|T)
					}else{
						$divisor--;										# N=unknown nucleotide
					}
					$nuc_counter++;
				}
				$gc_ratio=$gc/$divisor;
				$gc_ratio=sprintf("%.2f", $gc_ratio);			# round to two digits
				push(@gc_stats_array,$gc_ratio);
				my @n_array=();									# statistics start here
				if(($win_counter+1)%$no_of_win==0){				# modulus (remainder of division) = 0 (every $no_of_win)
					my $total_1 = 0;
					for(my $n=$no_of_win;$n>0;$n--){			# push last $no_of_win numbers into an array
						push(@n_array,$gc_stats_array[$win_counter+1-$n]);
						$total_1=$gc_stats_array[$win_counter+1-$n]+$total_1;	# sum up
					}
					my $mean_1=$total_1/$no_of_win;
					$mean_1=sprintf("%.3f", $mean_1);			# round to three digits
					my $total_2 = 0;							# total 2 and mean 2 (of squares) as prep for stdevp
					for(my $n=$no_of_win;$n>0;$n--){
						$total_2=(($mean_1-$gc_stats_array[$win_counter+1-$n])**2)+$total_2;
					}
					my $mean_2=$total_2/$no_of_win;
					my $stdevp=sqrt($mean_2);					# standard deviation based on the entire population
					$mean_2=sprintf("%.3f", $mean_2);			# round to three digits
					$stdevp=sprintf("%.3f", $stdevp);			# round to three digits
					if($stdevp>=$gc_thresh){					# threshold violation (program core)
						$chimalarm="on";						# set alarm flag
						my $temp_chim_file=$temp_chim_dir."/".$id;
						if(!-d $temp_chim_dir){
							mkdir($temp_chim_dir,0777) || warn "\tWarning: unable to create ".$temp_chim_dir.": ".$!."\n\n";
						}
						open(TEMPCHIM,">>$temp_chim_file") or die "\nUnable to open ".$temp_chim_file.": ".$!."\n\n";	
						print TEMPCHIM ($seq_file." ".($stdevp*100)." ".$win_counter." ".$inline."\n");
					}
					if($logging eq "on"){
						print LOGFILE ("\n\npos. ".($win_counter+1-$no_of_win)." - ".($win_counter).":");
						print LOGFILE ("\nmean\t= ".($mean_1*100)."%\nSTDEVP\t= ".($stdevp*100)."%\n\(values\t=");
						for(my $m=0;$m<scalar(@n_array);$m++){	# debug: write array contents into log file
							print LOGFILE (" ".$n_array[$m]);
						}
						print LOGFILE ("\)");
					}
				}
				$gc=0;											# reset gc counter
				$win_counter++;									# ... next window
				$nuc_counter=$win_counter;						# also set new start for nuc_counter
			}
			for(my $j=0;$j<$win_counter;$j++){
				print SEQSTATS ($j+1 ." ".$gc_stats_array[$j]."\n");# write into gc stats file
			}
			if($logging eq "on"){
				close(LOGFILE);
			}
			close(SEQSTATS);
			close(INFILE);
		}
		if($id==$num_of_threads && $verbose eq "off"){			# only show percentage for last thread and if not verbose
			if($seq_stats_counter/scalar(@seq_no_array)>=$seq_stats_progress/100){	# progress calc (estimate)
				print("done:\t".$seq_stats_progress."%\r");
				$seq_stats_progress=$seq_stats_progress+5;
			}
			$seq_stats_counter++;
		}
	}
	closedir(INDIR);
	if($chimalarm eq "on"){
		close(TEMPCHIM);
	}
    threads->exit();
}

sub outfile_generation(){
	my $outfile;
	if($file_parse eq "on"){
		$outfile=$assembly_file."_".$timestamp."-w".$win_size."-n".$no_of_win."-t".$gc_thresh*100 .".out";
	}elsif($file_parse eq "off"){
		$outfile=$timestamp."-w".$win_size."-n".$no_of_win."-t".$gc_thresh*100 .".out";
	}
	my $no_of_chimeras=0;
	open(OUTFILE,">$outfile") or die "Unable to write ".$outfile.": ".$!."\n\n";
	print OUTFILE ("sequence_number stdev position sequence\n");		# header preparation of outfile
	opendir(TEMPCHIM,$temp_chim_dir) or die "\nUnable to open ".$temp_chim_dir.": ".$!."\n\n";
	while(defined(my $infile=readdir(TEMPCHIM))){					# loop through every file
		my $infile_path=$temp_chim_dir."/".$infile;
		open(INFILE,"<$infile_path") or die "\tUnable to open ".$infile_path.": ".$!."\n\n";
		while(<INFILE>){
			chomp;
			print OUTFILE "$_\n";
			$no_of_chimeras++;
		}
		close INFILE;
	}
	close OUTFILE;
	closedir TEMPCHIM;
	system('/bin/rm','-rf',$temp_chim_dir);
	$exec_time=time()-$start_time;			# execution time
	print("\n".$exec_time."s:\t".$prog_name." finished\n\t".$no_of_chimeras." chimera(s) found\n\tSee ".$outfile." for details\n\n");
}

########################################
# III.		Main
########################################
if($file_parse eq "on"){
	parse_fasta();
	write_seq_files();
}
if($dir_parse eq "on"){
	my @threads=initThreads();					# create an array of threads and loop through the threads array:
	$exec_time=time()-$start_time;			# execution time
	print($exec_time."s:\tcalculating GC stats, number of threads: ".$num_of_threads);
	foreach(@threads){							# Tell each thread to perform our 'doOperation()' subroutine.
		$_ = threads->create(\&gc_stats);
	}
	foreach(@threads){							# Tell main program to keep running until all threads have finished.
		$_->join();
	}
	$exec_time=time()-$start_time;			# execution time
	print("\n".$exec_time."s:\tdone calculating");
}else{
	usage();									# if neither -i nor -d is given as input
}
if(-d $temp_chim_dir){
	outfile_generation();
}else{
	$exec_time=time()-$start_time;			# execution time
	print("\n".$exec_time."s:\t".$prog_name." finished\n\tNo Chimera found\n\n");
}
