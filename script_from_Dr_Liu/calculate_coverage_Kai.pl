#!/usr/bin/perl -w
# given brat output, calculate the depth of coverage at each genome position
# and percentage of genome that was covered by at least one or two reads.
use strict;
my %cover;
my %chr_len = ("chr1"=>30427671, "chr2"=>19698289, "chr3"=>23459830, "chr4"=>18585056, "chr5"=>26975502, "chrM"=>366924, "chrC"=>154478, "pUC19"=>2686);

my $debug = 0;
my $count = 0;
my $constant = 3000000;
if ($debug) {$constant = 300 ;}
print STDERR "debug\t$debug\n";

print STDERR join("\t", @ARGV), "\n";

while(<>){
	$count++;
	if ($count % $constant == 0) {print STDERR "$count\tDONE\n"}
	chomp;
	my @temp = split /\t/;
	
	if(@temp == 6){ ## single map
	    my ($id, $seq, $chr, $strand, $start, $mismatches) = @temp;
	    foreach my $i(($start+1)..($start+length($seq))){
		
		    if($strand eq "+"){
		        $cover{$chr . "W"}->[$i]++;
		    }elsif($strand eq "-"){
			    $cover{$chr . "C"}->[$i]++;
		    }else{
				die "strand direction is not + or - for read $id";
			}
		}
	}elsif(@temp == 9){ # pair map
	    my ($id, $seq1, $seq2, $chr, $strand, $start1, $start2, $mis1, $mis2) = @temp;
		my ($str1, $str2) = ("", "");
		if($strand eq "+"){
			($str1, $str2) = ("W", "C");
		}elsif($strand eq "-"){
			($str1, $str2) = ("C", "W");
		}else{
			die "strand is not + or - for read $id";
		}
		foreach my $i(($start1+1)..($start1+length($seq1))){
		        $cover{$chr . $str1}->[$i]++;
		}
		foreach my $i(($start2+1)..($start2+length($seq2))){
			    $cover{$chr . $str2}->[$i]++;
		}

	}else{
		die "number of items is not 6 or 9 for read ", $temp[0];
	}
}
my (%cover_1W, %cover_2W, %cover_1C, %cover_2C, %cover_depth);
my $total_len = 0;
foreach my $chr(keys %chr_len){
	$total_len += $chr_len{$chr};
	foreach my $i(1..$chr_len{$chr}){
		if(defined $cover{$chr . "W"}->[$i]){ 
			$cover_depth{$chr . "W"} += $cover{$chr . "W"}->[$i];
			if($cover{$chr . "W"}->[$i] >= 1){
			    $cover_1W{$chr}++;
			    $cover_1W{"total"}++;
			}
			if($cover{$chr . "W"}->[$i] >= 2){
				$cover_2W{$chr}++;
				$cover_2W{"total"}++;
			}
		}
		if(defined $cover{$chr . "C"}->[$i]){
			$cover_depth{$chr . "C"} += $cover{$chr . "C"}->[$i]; 
			if($cover{$chr . "C"}->[$i] >= 1){
			    $cover_1C{$chr}++;
			    $cover_1C{"total"}++;
			}
			if($cover{$chr . "C"}->[$i] >= 2){
				$cover_2C{$chr}++;
				$cover_2C{"total"}++;
			}
		}
	}
}
print "Coverage >= 1:\n";
print "On Watson strand:\n";
foreach my $chr(sort keys %chr_len){
	print $chr, "\t", $cover_1W{$chr}/$chr_len{$chr} * 100, "%\n";
}
print "Whole genome", "\t", $cover_1W{"total"}/$total_len * 100, "%\n";
print "On Crick strand:\n";
foreach my $chr(sort keys %chr_len){
	print $chr, "\t", $cover_1C{$chr}/$chr_len{$chr} * 100, "%\n";
}
print "Whole genome", "\t", $cover_1C{"total"}/$total_len * 100, "%\n";
print "Coverage >= 2:\n";
print "On Watson strand:\n";
foreach my $chr(sort keys %chr_len){
	print $chr, "\t", $cover_2W{$chr}/$chr_len{$chr} * 100, "%\n";
}
print "Whole genome", "\t", $cover_2W{"total"}/$total_len * 100, "%\n";
print "On Crick strand:\n";
foreach my $chr(sort keys %chr_len){
	print $chr, "\t", $cover_2C{$chr}/$chr_len{$chr} * 100, "%\n";
}
print "Whole genome", "\t", $cover_2C{"total"}/$total_len * 100, "%\n";		
print "Depth coverage:\n";
print "On Watson strand:\n";
foreach my $chr(sort keys %chr_len){
	print $chr, "\t", $cover_depth{$chr . "W"} / $chr_len{$chr}, "\n";
}
print "On Crick strand:\n";
foreach my $chr(sort keys %chr_len){
	print $chr, "\t", $cover_depth{$chr . "C"} /$chr_len{$chr}, "\n";
}

