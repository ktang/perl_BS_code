#!/usr/bin/perl -w

use strict;

opendir(CDIR, ".");
my @dirs = readdir CDIR;
foreach my $dir(@dirs){
	
	#if($dir =~ /Lane1_JKZ1_Col0/){
	#	next;
	#}
    if((-d $dir) && ($dir =~ /lane/)){
		#next if ($dir eq '.' || $dir eq '..');
        print STDERR "Calculate base composition under $dir ...\n";

	    chdir($dir);
	    #print STDERR "checking quality under $dir ...\n";
	
	    `perl ../base_composition.pl *_noadapt.fa >>base_composition_stats.txt`;
	    chdir('..');
	}
}
