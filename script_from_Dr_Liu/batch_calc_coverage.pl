#!/usr/bin/perl -w
use strict;
opendir(CDIR, ".");
my @dirs = readdir CDIR;
foreach my $dir(@dirs){
    if((-d $dir) && ($dir =~ /lane/)){
		print STDERR "calculating coverage under $dir\n";
		chdir($dir);
		opendir(SDIR, ".");
		my @files = grep {/brat\.txt/} readdir SDIR;
		my $file_str = join(" ", @files);
		`perl ../calculate_coverage.pl $file_str >>coverage_stats.txt`;
	    chdir('..');
	    
	 
	}
}
