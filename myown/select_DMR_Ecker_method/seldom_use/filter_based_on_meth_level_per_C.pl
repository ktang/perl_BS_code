#!/usr/bin/perl -w

use strict;

#@F[25..26]
#wt_MethLevel_perC       mut_MethLevel_perC

#my $usage = "$0 \n <input> <index4/7> STDOUT\n\n";
#die $usage unless (@ARGV == 2);

my $index = 25;


my $usage = "$0 \n <input> <fold> STDOUT\n\n";
die $usage unless (@ARGV == 2);

my $input = shift or die "input";
my $fold  = shift or die "fold";
die unless (-e $input);

open(IN, $input) or die;

my $head = <IN>;
chomp $head;
my @pts = split "\t", $head;
die  unless ($pts[$index] =~ /_perC$/ and $pts[$index + 1] =~ /_perC$/ );
#print join("\t", (@pts, "p_CG", "p_CHG", "p_CHH")), "\n";

print join("\t", @pts), "\n";

while(<IN>){
	chomp;
	my @a = split "\t";
	if( $a[$index] / ($a[$index + 1] + 0.00000001) >= $fold){
		print join("\t", @a), "\n";
	}
	
}

close(IN);

exit;