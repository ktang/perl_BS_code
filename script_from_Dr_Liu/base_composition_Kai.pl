#!/usr/bin/perl -w

# Given a sequence file, calculate base composition

use strict;
use Bio::SeqIO;
use Bio::Tools::SeqStats;
use Number::Format;

my $usage = "$0 <seq file> [seq file 2...]";
die $usage unless (@ARGV > 0);

foreach my $file(@ARGV){
my %all; # stores counts for all sequences

my $format = Number::Format->new();
my $seqin = Bio::SeqIO->new(-file=>$file);
my $numSeq = 0;
my $total = 0;
while(my $seq = $seqin->next_seq){
	my $seqStats = Bio::Tools::SeqStats->new(-seq=>$seq);
	my $ref = $seqStats->count_monomers();
	foreach my $base (keys %$ref){
		$all{$base} += $ref->{$base};
	}
	$total += $seq->length;
	$numSeq++;
}

print "Base composition for $numSeq sequences in $file:\ttotal: $total\t";
foreach my $base(sort keys %all){
	print "\t", $base, ": ", $all{$base}, ", ",$format->round((($all{$base})/$total) * 100), "%";
}
print "\n";
}