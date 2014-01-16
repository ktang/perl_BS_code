#!/usr/bin/perl -w

use strict;

my $debug = 1;

my $usage ="$0 <input> <output>";
die $usage unless (@ARGV == 2);
my $input = shift or die;
my $output = shift or die;

die unless(-e $input);
die if (-e $output);

open (IN, $input)  or die "cannot open $input:$!";
open (OUT, ">$output") or die "cannot open $output:$!";
	
my $i = 0;
while(<IN>){
	$i++;
	if ($i % 4 == 2){
		chomp;
		print OUT join("\t", ($_, 0, 0)), "\n";
	}
}

close (IN);
close (OUT);
exit;