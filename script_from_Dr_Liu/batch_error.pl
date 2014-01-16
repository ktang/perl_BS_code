#!/usr/bin/perl -w
use strict;

opendir(CDIR, ".");
my @dirs = readdir CDIR;
foreach my $dir(@dirs){
    if((-d $dir) && ($dir =~ /lane/)){
		print STDERR "calculating error under $dir\n";
		chdir($dir);
		opendir(SDIR, ".");
		my @files = grep {/_forw\.txt/} readdir SDIR;
		
		foreach my $file(@files){
		    my $pre = "NONE";
		    if($file =~ /(s_\d+)/){
			    $pre = $1;
		    }
			my $rev_file = $pre . "_rev.txt";
			
		    my $out = $pre . "_conversion_error.txt";	
		    `perl ../conversion_error.pl $file >> $out`;
			`perl ../conversion_error.pl $rev_file >> $out`;
		}
	    chdir('..');
	    
	 
	}
}
