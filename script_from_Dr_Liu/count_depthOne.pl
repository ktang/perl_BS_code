#!/usr/bin/perl -w

use strict;

my $max = 42859753;

my (%one, %mC);

die unless (@ARGV == 1);

my $input = $ARGV[0];

open (IN, $input) or die "cannot open $input:$!";

while(<IN>){
	chomp;
	next if (/SampleID/);
	my @a = split "\t";
	if ($a[1] <= $max){
		if($a[2] == 1){
			$one{$a[5]}++;
			if($a[3] == 1){
				$mC{$a[5]}++;
			}
		}
	}else{
		last;
	}
}
close IN;

print "in $input\n";
#print join("\t", ("num_depthOne", "num_mC_depthOne") ), "\n";
#print join("\t", ($one, $mC) ), "\n";

my ($t1, $tm ) = (0, 0);
foreach my $key (sort keys %one){
	print "$key: ",$one{$key}, "\t";
	$t1 += $one{$key};
}

print "num_depthOne_total: ", $t1, "\n";

foreach my $key (sort keys %mC){
	print "$key: ", $mC{$key}, "\t";
	$tm += $mC{$key}
}

print "num_mC_depthOne_total: ", $tm, "\n";

#foreach my $t(sort keys %num_pos){
#    print STDERR "$t: ", $num_pos{$t}, "\t";#
#	$total += $num_pos{$t};
#}