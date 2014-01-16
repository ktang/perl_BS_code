#!/usr/bin/perl -w
#v1.0 July 15
# try to for each small region, narrow down first, then check gap

#method4:
#netDMC / basic mC >= cutoff
#basic mC = mentholated in both samples.


use strict;
use File::Spec;

my $debug_sub = 1;
my $debug = 0;
if($debug){
	print STDERR "debug:$debug\n";
	
}




my $usage = "$0 \n  <wig1>  <wig2> STDOUT \n\n";
die $usage unless (@ARGV == 2);

my $wig1 = shift or die "wig1";
my $wig2 = shift or die "wig2";

#my (%pos_hypers, %pos_hypos);
#my %pos_both;


my (%pos_wig1_only, %pos_wig2_only, %pos_both);

if($debug){
	print STDERR "read_wig...\n";
}
read_wig1($wig1, \%pos_wig1_only); 
read_wig2($wig2, \%pos_wig2_only, \%pos_wig1_only, \%pos_both);

my $num1 = cal_met(\%pos_wig1_only);
my $num2 = cal_met(\%pos_wig2_only);
my $num_both = cal_met(\%pos_both);

print "\n";
print $wig1, " only:\n";
print $num1, "\n";

print $wig2, " only:\n";
print $num2, "\n";

print"both:\n";
print $num_both, "\n";

exit;


sub read_wig1{
	my ($input, $ref) = @_;
	#print STDERR "read_wig reading $input...\t";
	open (IN, $input) or die "cannot open $input:$!";
	my $chr = "";
	while (<IN>){
		next if (/^track/);
		if (/variableStep\s+chrom=(\w+)/){
			$chr = lc $1;
			next;
		}
		chomp;
		my ($pos, $val) = split /\t/;
	#	$wigs{$chr} -> {$pos} = $val;
		$ref->{$chr}->{$pos} = 1;
	}
	close(IN);
#	print STDERR "DONE\n";
}

sub read_wig2{
	my ($input, $wig2_only_ref, $wig1_only_ref, $both_ref) = @_;
#	print STDERR "read_wig reading $input...\t";
	open (IN, $input) or die "cannot open $input:$!";
	my $chr = "";
	while (<IN>){
		next if (/^track/);
		if (/variableStep\s+chrom=(\w+)/){
			$chr = lc $1;
			next;
		}
		chomp;
		my ($pos, $val) = split /\t/;
		if(defined $wig1_only_ref->{$chr}->{$pos}){
			delete $wig1_only_ref->{$chr}->{$pos};
			$both_ref->{$chr}->{$pos} = 1;
		}else{
			$wig2_only_ref->{$chr}->{$pos} = 1;
		}
	}
	close(IN);
#	print STDERR "DONE\n";

}



sub cal_met{
	my ($ref) = @_ ;
	my $num = 0;
	foreach my $chr (sort keys %{$ref}){
		$num += scalar(keys %{${$ref}{$chr}})
	}
	return $num;
}
