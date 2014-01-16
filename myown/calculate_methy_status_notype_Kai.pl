#!/usr/bin/perl -w
# calcualte methylation status using Lister et al. 2009 method
# do not separate different position types (CG, CHG, CHH)
use strict;
use Math::CDF qw(:all);
my $FDR = 0.01;
my $print_values = 0;
my $usage = "$0 <meth file> <error rate> <cutoff_file> <max_pos_id>";
die $usage unless(@ARGV >= 3);
my ($meth_file, $rate) = @ARGV[0..1];
my $cutoff_file = $ARGV[2];
my $meth_max = 0;
my @accu_meth;
my @total_meth;
my @meth_level;
my $max_pos_id = 42859753; #don't want to include chrC, chrM, pUPC19
if(@ARGV > 3){
	$max_pos_id = $ARGV[3];
}

print STDERR "max_pos_id is: $max_pos_id\n";

print STDERR "reading $meth_file..\t";

open(MF, $meth_file) or die "Can't open $meth_file: $!";
while(<MF>){
	next if (/SampleID/);
	chomp;
	my ($sample_id, $pos_id, $depth, $num_C, $percent, $type) = split /\t/;
	next if($pos_id > $max_pos_id);
	
	$meth_level[$depth]->[$num_C]++;
	$total_meth[$depth]++;
		
	if($meth_max < $depth){
		$meth_max = $depth;
	}
}
close MF;
print STDERR "Done\n";

print STDERR "File: $meth_file, Conversion Error: $rate\n";			
print STDERR "Max depth: $meth_max\n";

print STDERR "calculate accumulated number..\t";
# calculate accumulated number of positions with <= k unconverted Cs
	foreach my $i(0..$meth_max){
		next if($i == 0 or !defined $total_meth[$i]);

		if(!defined $meth_level[$i]->[0]){
			$meth_level[$i]->[0] = 0;
		}
		$accu_meth[$i]->[0] = $meth_level[$i]->[0];
		foreach my $j(1..$i){
			if(!defined $meth_level[$i]->[$j]){
					$accu_meth[$i]->[$j] = $accu_meth[$i]->[$j-1]			
			}
			else{
					$accu_meth[$i]->[$j] = $accu_meth[$i]->[$j-1] +
											$meth_level[$i]->[$j];
			}
		}
	}
print STDERR "Done\n";

print STDERR "depth = 2 : $total_meth[2]\n";
print STDERR "0: ",$meth_level[2]->[0] ,"\t1: ",$meth_level[2]->[1], "\t2: ", $meth_level[2]->[2], "\n";
print STDERR "accum:0: ",$accu_meth[2]->[0] ,"\t1: ",$accu_meth[2]->[1], "\t2: ", $accu_meth[2]->[2], "\n";

undef @meth_level;# empty memory

print STDERR "Empth meth_leverl\nfind cutoff at each depth level..\t";

# find cutoff at each depth level
my @cutoff;
	foreach my $i(0..$meth_max){
		next if(!defined $total_meth[$i]);
		if($i < 2){
			$cutoff[$i] = $i+1;
			#next;
		}else{
			$cutoff[$i] = 4;
			
			my $start = int ($i * $rate + 0.999999);
			foreach my $j($start..$i){
			
				my $prob = pbinom($j, $i, $rate);
				$prob = $prob - pbinom($j-1, $i, $rate);
				my ($num_unC, $num_mC) = (0,0);
				$num_unC = $accu_meth[$i]->[$j-1];
				if(!defined $num_unC){
					$num_unC = 0;
				}
				$num_mC = $total_meth[$i] - $num_unC;
				my $left = $prob * $num_unC;
				my $right = $FDR * $num_mC;
				if($left < $right){
					$cutoff[$i] = $j;
					last;
				}
			}
		}
	}
print STDERR "Done\n";
undef @total_meth;
print STDERR "Empty total_meth\n";

open (CUT, ">>$cutoff_file") or die "cannot open $cutoff_file: $!";

print CUT "cutoff for $meth_file:\n";
print CUT join("\t", ("depth", "cutoff")),"\n";

for my $index (0..$meth_max){
	if (!defined $cutoff[$index]){
		$cutoff[$index] = 0;
	}
	print CUT join("\t", ($index , $cutoff[$index])), "\n";
}

close(CUT);


my %num_pos;
open(MF, $meth_file) or die "Can't open $meth_file: $!";
#seek(MF, 0, 0) or die "Can't seek to beginning of file: $!";

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
