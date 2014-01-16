#!/usr/bin/perl -w
# purpose:
# Sep12, 2013: now is just want get number for Venn(3set) the all_overlap number
# as the input is Jacobsen_cell method so the overlap is easy


use strict;
use File::Spec;

my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}
#my $usage = "$0 \n <outdir> <cutoff> <number_of_files> <isMeth> <label> [<isMeth> <label> ...] \n\n";
#die $usage unless(@ARGV >= 7);

my $usage = "$0 \n <number_of_files> <input1>  <input2> [ <input3> ...] \n\n";
die $usage unless(@ARGV >= 3);


#my $postfix = shift or die;
#my $outdir = shift or die;
#my $dep_cutoff = shift or die;
my $num_of_files = shift or die;
#my $last_index = $num_of_files - 1;

my @inputs;
for my $i(0..($num_of_files - 1)){
	my $tmp = shift or die;
	push @inputs, $tmp;
	
}
print STDERR  join("\n", @inputs), "\n\n";

my %h;# record the coor

foreach my $file(@inputs){
	open(IN, $file) or die "$file";
	my $head = <IN>;
	while (<IN>) {
		chomp;
		my @a = split "\t";
		$h{$a[0]}++;
	}
}

my $all = 0;
foreach my $k (keys %h){
	die $k if($h{$k} > $num_of_files);
	$all++ if($h{$k} == $num_of_files);
}

print STDERR "all_overlap: $all \n\n";

exit;
