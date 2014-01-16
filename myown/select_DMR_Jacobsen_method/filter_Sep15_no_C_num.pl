#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;

my $debug = 0;

my ( $CG_num_cutoff, $CHG_num_cutoff, $CHH_num_cutoff ) = (5,5, 20);
my ( $CG_diff_cutoff, $CHG_diff_cutoff, $CHH_diff_cutoff ) = (30, 10 ,10);

if($debug){
	print STDERR "debug = 1\n\n";
}
#my $usage = "$0 <input> <output>";
#die $usage unless(@ARGV == 2);

my $usage = "$0]\n <input> STDOUT \n\n";
die $usage unless(@ARGV == 1);


my $input = shift or die;
#my $output = shift or die;

die unless (-e $input);
#die if( -e $output);

open(IN, $input) or die "cannot open $input: $!";

#die if(-e $output);
#open(OUT, ">$output") or die "cannot open $output: $!";

my $head = <IN>;
print $head;
#chomp $head;
#print  join("\t", ($head, "high_weigth_meth_level", "low_weigth_meth_level")), "\n";
#print  join("\t", ($head, "high_C_per", "low_C_per")), "\n";

#0			1	2	 3		  4			5		6		7				8		9	  	10		11			12		13				14		15
#chr	start	end	DMC_num	 wt_CG	wt_CG_per	wt_CHG	wt_CHG_per	wt_CHH	wt_CHH_per	mut_CG	mut_CG_per	mut_CHG	mut_CHG_per	 mut_CHH	mut_CHH_per

#chr	start	end	DMC_num	wt_CG	wt_CG_per	wt_CHG	wt_CHG_per	wt_CHH	wt_CHH_per	mut_CG	mut_CG_per	mut_CHG	mut_CHG_per	mut_CHH	mut_CHH_per |	16_avge_meth_level_wt
#avge_meth_level_mut	Gene	Strand	GeneType	GeneAnnot	TE	TEFamily	Intergenic	Promoter_1kb	26_CG_num	27_CHG_num	28_CHH_num	high_weigth_meth_level	low_weigth_meth_level

#0			1		2	 3		  4			5		6		7		8		9	  10		11		12		13		14		15
#chr1    42319   42670   20      5/9=    55.5556 2/35=   5.7143  68/451= 15.0776 0/8=    0.0000  0/39=   0.0000  4/324=  1.2346
while(<IN>){
	chomp;
	my @a = split "\t";
	
	my ( $cg_num, $chg_num, $chh_num ) = @a[26..28];
	
	my ( $cg_diff, $chg_diff, $chh_diff ) = (0) x 3;
	
	$cg_diff  = $a[5] - $a[11]  if ($a[5]  ne "NA");
	$chg_diff = $a[7] - $a[13]  if ($a[7]  ne "NA");
	$chh_diff = $a[9] - $a[15]  if ($a[9]  ne "NA");
	
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