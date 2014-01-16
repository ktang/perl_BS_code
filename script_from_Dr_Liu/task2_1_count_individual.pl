#!/usr/bin/perl -w

use strict;

my $usage = "$0 <input> <output>";

die $usage unless (@ARGV == 2);

my ($input, $output) = @ARGV[0..1];

die "wrong input or output" unless (-e $input and !(-e $output));
		 
open(IN, $input) or die "cannot open $input";
open(OUT, ">$output") or die "cannot open output $output";
		
print OUT join ("\t", ("depth", "number", "methylated_number")), "\n";
		
my (%covs, %meths);
		
<IN>;#skip head;
while(<IN>){
	chomp;
	my @a = split "\t";
	$covs{$a[2]} ++;
	if ($a[6] == 1){
		$meths{$a[2]}++;
	}
}
		
foreach my $key (sort {$a <=> $b} keys %covs){
	my ($cov, $met) = (0, 0);
	$cov = $covs{$key};
	if (defined $meths{$key}){
		$met = $meths{$key};
	}
	print OUT join ("\t", ($key, $cov, $met)), "\n";
}
exit;