#!/usr/bin/perl -w

use strict;

my $qual_folder = "quality_h20";
if(@ARGV > 0){
	$qual_folder = $ARGV[0];
}
opendir(CDIR, ".");
my @dirs = readdir CDIR;
foreach my $dir(@dirs){
	
	#if($dir =~ /Lane1_JKZ1_Col0/){
	#	next;
	#}
    if((-d $dir) && ($dir =~ /lane/)){
		next if ($dir eq '.' || $dir eq '..');
        opendir(LDIR, $dir);
		my @files = readdir LDIR;
        foreach my $file(@files){
			if($file =~ /matrix$/ || $file =~ /png$/ || $file =~ /quality$/ ||
				$file =~ /pdf$/ || $file =~ /segments$/){
	
	         `mv $dir/$file $qual_folder`;
		    }
		}
	    
	}
}
