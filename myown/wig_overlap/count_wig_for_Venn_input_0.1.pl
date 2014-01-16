#!/usr/bin/perl -w

#v0.1
# input dir and prefix, for mC and CG, CHG, CHH

use strict;
use File::Spec;

my $debug = 0;
if($debug){
	print STDERR "debug:$debug\n";
	
}

my $usage = "$0 \n <indir1> <wig_pre1> <indir2> <wig_pre2> STDOUT \n\n";
die $usage unless (@ARGV == 4);

#my $wig1 = shift or die "wig1";
#my $wig2 = shift or die "wig2";

my $dir1 = shift or die;
my $pre1 = shift or die;
my $dir2 = shift or die;
my $pre2 = shift or die;

die unless (-d $dir1);
die unless (-d $dir2);

check_file($dir1, $pre1, $dir2, $pre2, "mC");
check_file($dir1, $pre1, $dir2, $pre2, "mCG");
check_file($dir1, $pre1, $dir2, $pre2, "mCHG");
check_file($dir1, $pre1, $dir2, $pre2, "mCHH");

if($debug){
	print STDERR "file OK\n\n\n";
	exit;
}

Venn($dir1, $pre1, $dir2, $pre2, "mC");
Venn($dir1, $pre1, $dir2, $pre2, "mCG");
Venn($dir1, $pre1, $dir2, $pre2, "mCHG");
Venn($dir1, $pre1, $dir2, $pre2, "mCHH");


exit;

sub check_file{
	my ($sub_dir1, $sub_pre1, $sub_dir2, $sub_pre2, $label) = @_;
	my $file1 = File::Spec->catfile($sub_dir1, $sub_pre1 . "_" . $label . ".wig");
	my $file2 = File::Spec->catfile($sub_dir2, $sub_pre2 . "_" . $label . ".wig");
	
	die $file1 unless (-e $file1);
	die $file2 unless (-e $file2);
}

sub Venn{
	my ($sub_dir1, $sub_pre1, $sub_dir2, $sub_pre2, $label) = @_;
	my $wig1 = File::Spec->catfile($sub_dir1, $sub_pre1 . "_" . $label . ".wig");
	my $wig2 = File::Spec->catfile($sub_dir2, $sub_pre2 . "_" . $label . ".wig");
	
	my (%pos_wig1_only, %pos_wig2_only, %pos_both);

	read_wig1($wig1, \%pos_wig1_only); 
	read_wig2($wig2, \%pos_wig2_only, \%pos_wig1_only, \%pos_both);

	my $num1 = cal_met(\%pos_wig1_only);
	my $num2 = cal_met(\%pos_wig2_only);
	my $num_both = cal_met(\%pos_both);

	print "\n";
	print $sub_pre1, "_$label\t";
	print $num1, "\n";

	print $sub_pre2, "_$label\t";
	print $num2, "\n";

	print"both_$label\t";
	print $num_both, "\n";
}


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
