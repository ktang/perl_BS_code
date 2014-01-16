#!/usr/bin/perl -w
# calcualte methylation status using Lister et al. 2009 method
# do not separate different position types (CG, CHG, CHH)
use strict;
use Math::CDF qw(:all);
my $FDR = 0.01;
#my $prin  _values = 1;
my $usage = "$0 <meth file> <error rate> <max_pos_id>";
die $usage unless(@ARGV >= 2);
my ($meth_file, $rate) = @ARGV[0..1];
my $meth_max = 0;
my @accu_meth;
my @total_meth;
my @meth_level;
my $max_pos_id = 42859753; #don't want to include chrC, chrM, pUPC19
if(@ARGV > 2){
	$max_pos_id = $ARGV[2];
}

my $cutoff_file = "cutoff_".$rate.".txt";
if (-e $cutoff_file) {die "$cutoff_file exists\n"}
open(MF, $meth_file) or die "Can't open $meth_file: $!";
while(<MF>){
	next if (/SampleID/);
	chomp;
	my ($sample_id, $pos_id, $depth, $num_C, $percent, $type) = split /\t/;
	next if($pos_id > $max_pos_id);
				if($meth_max < $depth){
			$meth_max = $depth;
		}
}
close MF;



print STDERR "File: $meth_file, Conversion Error: $rate\n";	
print STDERR "Max depth: $meth_max\n";

my $perl_cutoff = "/Users/tang58/Kai_BS/cutoff_Kai_version.pl";

my $cmd = "perl $perl_cutoff $meth_max $rate > $cutoff_file";

print STDERR "$cmd\n";

system("$cmd");
		
# find cutoff at each depth level
my @cutoff;
 
open (IN, $cutoff_file) or die "cannot open $cutoff_file:$!";
while (<IN>){
	chomp;
	my @a = split "\t";
	$cutoff[$a[0]] = $a[1];
}

close IN;

my %num_pos;
 
open(MF, $meth_file) or die "Can't open $meth_file: $!";
while(<MF>){
	chomp;
	if (/SampleID/){
		print;
		print "\tisMeth\n";
		next;
	}
	
	my ($sample_id, $pos_id, $depth, $num_C, $percent, $type) = split /\t/;
	next if($pos_id > $max_pos_id);
	my $isMeth = 0;
		if($depth >= 2 && $num_C >= 1 && $num_C >= $cutoff[$depth]){
			$isMeth = 1;
			$num_pos{$type}++;
		}
		print join("\t", ($sample_id, $pos_id, $depth, $num_C, $percent, $type, 
		    $isMeth)), "\n";
	
}
close MF;

print STDERR "Total mC position in ", $meth_file, "\n";
my $total = 0;
foreach my $t(sort keys %num_pos){
    print STDERR "$t: ", $num_pos{$t}, "\t";
	$total += $num_pos{$t};
}
print STDERR "Total: $total\n";
