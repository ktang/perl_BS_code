#!/usr/bin/perl -w
# take illumina data in fastq format and output sequence data in fasta format
# without considering quality
# seqs of length < 20 are discarded
# seqs contain "N" are discarded
use strict;
#use Bio::SeqIO;
#use Bio::Seq;

my $min_len = 20;

my $usage = "$0 <fastq file> <output seq file>";
die $usage unless (@ARGV >=2);
my ($input, $output) = @ARGV[0..1];
my $qual_file;

#my $out = Bio::SeqIO->new(-file=>">$output", -format=>"fasta");
my $out;
open($out, ">$output") or die "Can't open $output: $!";

open(IN, $input) or die "Cannot open $input:$!";

my $line = 0;
my ($seq_name, $seq_str);
my ($num_contain_N, $num_short, $num_kept) = (0,0,0);

while(<IN>){
	$line++;
	chomp;
	my $n = $line % 4;
	if($n == 1){
		if(/\:(\d+)\:(\d+)\:(\d+)\:(\d+)\#/){
			$seq_name = join("_", ($1, $2, $3, $4));
		}else{
			die "On line $line, seq name $_ does not match pattern\n";
		}
	}
	if($n == 2){
		$seq_str = $_;
		#my $seq = Bio::Seq->new(-seq => $seq_str, -id=>$seq_name);
		#if($yes_or_no eq "Y"){
			if(length($seq_str) < $min_len){
				$num_short++;
			}elsif($seq_str =~ /N/){
				$num_contain_N++;
			}else{
		        #$out->write_seq($seq);
				print $out ">", $seq_name, "\n";
				print $out $seq_str, "\n";
				$num_kept++;
			}
		#}
	}
}
close IN;
close $out;
print "In file $input:\n";
#print "Number of 'N' seqs (discarded): $num_no\n";
print "Number of good seqs: $num_kept, num_short (less than $min_len nt): $num_short, num_contain_N: $num_contain_N\n";	
