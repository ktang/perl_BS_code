#!/usr/bin/perl -w

#this overlap script is  for Jacobsen method 1:1-100
# according to Zhaobo's request
use strict;
use File::Spec;

my $debug = 0;

my $p = 0;
if($debug){
	print STDERR "debug = $debug\n\n";
}
my $win_len = 100;

my $usage = "$0 \n <input> <bench_file> <output>\n\n";

die $usage unless (@ARGV == 3);
#my $dir = shift or die $usage;
my $input = shift or die;
my $bench_file = shift or die;
my $output = shift or die "output";
die if(-e $output);

open( OUT, ">>$output" ) or die;

my @nums;
my %bench_records;
open(BENCH, $bench_file) or die;
my $h = <BENCH>;
while (<BENCH>){
	chomp;
	my @a = split "\t";
	$bench_records{$a[0]} = 1;
}
close BENCH;

open(IN, $input) or die;
$h = <IN>;
print OUT $h;
my $l = 0;
my $overlap = 0;
while(<IN>){
	$l++;
	chomp;
	my @a = split "\t";
	my ($chr, $s, $e) ;
	if ($a[0] =~ /(\S+):(\d+)-(\d+)/){
		($chr, $s, $e)  = ($1, $2, $3);
	}else{
		die $_;
	}
	my $former = $chr . ":" . ($s-$win_len) . "-" . ( $s-1);
	my $next   = $chr . ":" . ($s+$win_len) . "-" . ($s+$win_len*2 - 1);
	if($debug){
		exit if ($l >= 10);
		print STDERR join("\t", ($a[0], $former, $next)), "\n";
	}
	if (defined $bench_records{$a[0]} or defined $bench_records{$former} or defined $bench_records{$next} ){
		$overlap ++;
		print OUT $_, "\n";
	}
	
}
close IN;
close OUT;
exit;

sub overlap{
	my ($first, $second) = @_;
	open (IN1, $first)  or die "cannot open $first: $!";
	open (IN2, $second) or die "cannot open $second: $!";

	my %intervals;
	my $num_f2 = 0;
	my $h = <IN2>;
	while (<IN2>){
		$num_f2++;
		chomp;
		my @a = split /\t/ ;
		$intervals{$a[0]} = 1;
	}

	my $num_f1 = 0;
	my $num_overlap = 0;
	$h = <IN1>;
	while (<IN1>){
		$num_f1 ++;
		chomp;
		my @a = split /\t/ ;
		if ( defined $intervals{$a[0]} ){$num_overlap++}
	}
	my $per = sprintf ("%.1f",100 * $num_overlap / $num_f1);
	return "$num_overlap/$num_f1 = $per" . "%";
}