#!/usr/bin/perl -w


use strict;
use File::Spec;

my $debug = 0;

my $usage = "$0 \n <in1> <label1> <in2> <label2> STDOUT \n\n";
die $usage unless(@ARGV == 4);

my $input1 = shift or die "input";
my $pre1 = shift or die "shift";

my $input2 = shift or die "input";
my $pre2 = shift or die "shift";



die unless (-e $input1);
die unless (-e $input2);

my %pos;
my $num1 = read_first($input1, \%pos);



open(IN, $input2) or die "cannot open $input2";

my $num2 = 0;

my $both = 0;

while (<IN>){
	$num2 ++;
	chomp;
	my ($chr, $pos) = split "\t";
	if(defined $pos{$chr}->[$pos]){
		$both ++;
	}
}
close(IN);

print  "$pre1\t$num1\n";
print  "$pre2\t$num2\n\n";


print  $pre1, "_only\t" , $num1-$both, "\n";

print "both\t$both\n";


print  $pre2, "_only\t" , $num2-$both, "\n";

print "\n\n";

exit;

#my $num1 = read_first($input1, \%pos);

sub read_first{
	my ($file, $ref) = @_;
	die unless (-e $file);
	open(IN, $file) or die;
	
	my $i = 0;
	while (<IN>){
		chomp;
		$i++;
		my ($chr, $pos) = split "\t";
		$ref->{$chr}->[$pos] = 1;
	}
	close(IN);
	return $i;
}