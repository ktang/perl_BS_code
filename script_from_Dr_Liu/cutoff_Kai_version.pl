#!/usr/bin/perl -w
# calcualte methylation status using Lister et al. 2009 method
# do not separate different position types (CG, CHG, CHH)
use strict;
use Math::CDF qw(:all);
my $FDR = 0.01;
my $print_values = 0;
my $usage = "$0 <max_depth> <error rate>";
die $usage unless(@ARGV == 2);

my ($meth_max, $rate) = @ARGV[0..1];

	for my $i(2..$meth_max){
		
		
			for (my $j=min($i,1500) ; $j >= 1; $j--){
				my $prob = pbinom($j, $i, $rate);
				$prob = $prob - pbinom($j-1, $i, $rate);
				my ($num_unC, $num_mC) = ($i - $j,  $j);
				my $left = $prob * $num_unC;
				my $right = $FDR * $num_mC;
				if($left >= $right){
					my $last_left = (pbinom($j + 1, $i, $rate) - pbinom($j, $i, $rate) ) * ($i -$j -1);
					my $last_rigth =  $FDR * ($j + 1);
					print  join("\t", ($i, $j + 1,$last_left, $last_rigth, $left , $right)), "\n";
					last;
				}
			}
	}

			
	sub min{
		my ($a, $b) = @_;
		if($a > $b) {return $b}
		else {return $a}
	}