#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;

my $debug = 0;

#my ( $CG_num_cutoff, $CHG_num_cutoff, $CHH_num_cutoff ) = (5,5, 20);
#my ( $CG_diff_cutoff, $CHG_diff_cutoff, $CHH_diff_cutoff ) = (30, 10 ,10);

if($debug){
	print STDERR "debug = 1\n\n";
}
#my $usage = "$0 <input> <output>";
#die $usage unless(@ARGV == 2);

my $usage = "$0]\n <CG_num> <CG_diff> <CHG_num> <CHG_diff> <CHH_num> <CHH_diff>  <input> STDOUT \n\n";
die $usage unless(@ARGV == 7);

my $CG_num_cutoff    = shift or die ;
my $CG_diff_cutoff   = shift or die ;
my $CHG_num_cutoff   = shift or die ;
my $CHG_diff_cutoff  = shift or die ;
my $CHH_num_cutoff   = shift or die ;
my $CHH_diff_cutoff  = shift or die ;

my $input = shift or die;
#my $output = shift or die;

die unless (-e $input);
#die if( -e $output);

open(IN, $input) or die "cannot open $input: $!";

#die if(-e $output);
#open(OUT, ">$output") or die "cannot open $output: $!";

#my $CG_start = 3;# 4;
#my $per_start = 7;#8

my $head = <IN>;
print $head;

chomp $head;
my @hs = split "\t", $head;

my $CG_start = -1;
for my $i(0..$#hs){
	if($hs[$i] eq "CG_num" ){
		$CG_start  = $i;
		last;
	}
}
die if($CG_start == -1);

my $per_start = -1;
for my $i(0..$#hs){
	if($hs[$i] eq "mut_CG_per" or  $hs[$i] eq "wt_CG_per"){
		$per_start  = $i;
		last;
	}
}
die if($per_start == -1);

print STDERR " CG_start = $CG_start\t ";
print STDERR " per_start = $per_start\n\n ";

while(<IN>){
	chomp;
	my @a = split "\t";
	
#	my ( $cg_num, $chg_num, $chh_num ) = @a[26..28];
	my ( $cg_num, $chg_num, $chh_num ) = @a[$CG_start..($CG_start + 2)];
	
	my ( $cg_diff, $chg_diff, $chh_diff ) = (0) x 3;
	
#	$cg_diff  = $a[5] - $a[11]  if ($a[5]  ne "NA");
#	$chg_diff = $a[7] - $a[13]  if ($a[7]  ne "NA");
#	$chh_diff = $a[9] - $a[15]  if ($a[9]  ne "NA");
	
	$cg_diff  = $a[$per_start ] - $a[$per_start + 6]  if ($a[ $per_start ]  ne "NA");
	$chg_diff = $a[$per_start + 2] - $a[$per_start + 8]  if ($a[$per_start + 2]  ne "NA");
	$chh_diff = $a[$per_start + 4] - $a[$per_start + 10]  if ($a[$per_start + 4]  ne "NA");
	
	if($debug){
		print STDERR "number: ", join("\t",  ( $cg_num, $chg_num, $chh_num  )), "\n";
		print STDERR "diff: ",join("\t", ( $cg_diff, $chg_diff, $chh_diff ) ), "\n\n";
		exit;
	}
	
	if(  ($cg_num >= $CG_num_cutoff and $cg_diff >= $CG_diff_cutoff) 
		 or ($chg_num >= $CHG_num_cutoff and $chg_diff >=  $CHG_diff_cutoff ) 
		 or ($chh_num  >= $CHH_num_cutoff and $chh_diff >=  $CHH_diff_cutoff ) ){
		print $_, "\n";
	}
		
}

close(IN);
#close(OUT);

exit;

sub get_two_val{
	my ($x) = @_;
	$x =~ s/=//;
	my ($a, $b) = split "/", $x;
	return ($a, $b);
}