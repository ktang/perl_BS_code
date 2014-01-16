#!/usr/bin/perl -w
use strict;
opendir(CDIR, ".");
my @dirs = readdir CDIR;
my $log_file = "remove_adaptor.log";
foreach my $dir(@dirs){
	
	#if($dir =~ /Lane1_JKZ1_Col0/){
	#	next;
	#}
    if((-d $dir) && ($dir =~ /lane/)){
		#next if ($dir eq '.' || $dir eq '..');
        print STDERR "Remove adaptor under $dir ...\n";

			`perl remove_adaptor.pl $dir $dir >> $log_file`;
	    
	 
	}
}
