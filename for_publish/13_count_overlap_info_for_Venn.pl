#!/usr/bin/perl -w


use strict;
use File::Spec;

my $debug = 0;


my $usage = "$0 \n <input>  <number of overlaping samples>\n\n";
die $usage unless(@ARGV == 2);

my $input = shift or die;
my $sample_index = shift or die; 

die unless (-e $input);

open(IN, $input) or die;

my $h = <IN>;
chomp $h;
my @h_a = split "\t", $h;

my @labels = @h_a[($sample_index * -1) .. -1];
#perl -e '@a = (2,5,6,9,23,56); @b = @a[-3..-1] ;print join("\n", @b), "\n" '
#9
#23
#56

#chr1_122579_127669      NOT

my $total = 0;

my $power = 2**$sample_index;

my @numbers = (0) x $power ;

while (<IN>){
	$total++;
	chomp;
	my @a = split "\t";
	
	my $index = 0;
	
	for my $i (1..$sample_index){
		my $base = 2 ** ($i - 1);
		if($a[-$i] =~/chr\d_\d/){
			$index += $base;
		}else{
			die $_ unless ( $a[-$i]  eq "NOT" );
		}
		
	}
	
	$numbers[$index]++;	
}

close IN;

print STDERR $input, "\t", $total, "\n";
print STDERR join("\t", @labels), "\t", "number" , "\n";

for my $i (0.. ($power - 1)){
	my $temp = sprintf "%03b", $i;
	
#	print join("\t", (   $temp, $numbers[$i] )), "\n";
	
	
	my @ts = split //, $temp;
	print STDERR join("\t",( @ts, $numbers[$i]) ), "\n";
}


exit;