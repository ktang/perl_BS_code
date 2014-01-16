#!/usr/bin/perl -w
use strict;
opendir(CDIR, ".");
my @dirs = readdir CDIR;
my $log_file = "prepare_brat_input.log";
foreach my $dir(@dirs){
	
	#if($dir =~ /Lane1_JKZ1_Col0/){
	#	next;
	#}
    if((-d $dir) && ($dir =~ /lane/)){
		#next if ($dir eq '.' || $dir eq '..');
        print STDERR "Make brat input under $dir ...\n";

			`perl get_paired_reads.pl $dir $dir >> $log_file`;
	    
	 
	}
}
