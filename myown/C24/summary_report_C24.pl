#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;

my $debug = 0;

my $usage = "$0 <indir> <output>";
die $usage unless(@ARGV == 2 );

my $dir = shift or die "dir";
my $output = shift or die "output";

die unless (-d $dir);
die if (-e $output);

#my @pres = ( "037306005190A", "Num12A", "Num20A", "Num27A", 
#			 "Num77A", "akn1_2A", "ape_1A", "ape1l_2A",
#			 "arp_1A", "colA", "ibm1-4A", "Jmij-19A",
#			 "nrpd1-3A", "nrpel-11A", "P261A", "P313A",
#			 "P314A", "P31A", "pkl_1A", "shh1A"	);

my @pres = 	("C24_156A", "C24_204A", "C24_nrpd1aA", "C24_ros1nrpd1aA");		 

if($debug){
	@pres = ("colA");
}

#my $dir = "/Volumes/My_Book/20120427_ShangHai_data/handled_data/cytosine_coverage_info";
foreach my $pre(@pres){
	my $file = File::Spec->catfile($dir, $pre. "_cytosines_coverage_info.txt");
#	if($pre eq "ape1l_2A"){
#		$file = File::Spec->catfile($dir,  "ape1l_2A_not_col0_cytosines_coverage_info.txt");
#	}
	die unless (-e $file);
}

#my $output = File::Spec->catfile($dir, "basic_summary.txt");
die if (-e $output);
open(OUT, ">$output") or die;

print OUT join("\t", ("sample", "average_depth", "%_of_dep>=1",  "%_of_dep>=2", "CG_meth_level", "CHG_meth_level", "CHH_meth_level", "total_meth_level","error_rate")), "\n";
my $l;

foreach my $pre(@pres){
	my $file = File::Spec->catfile($dir, $pre. "_cytosines_coverage_info.txt");
#	if($pre eq "ape1l_2A"){
#		$file = File::Spec->catfile($dir,  "ape1l_2A_not_col0_cytosines_coverage_info.txt");
#	}
	die unless (-e $file);
	open(IN, $file) or die "cannot open $file:$!";
	$l = <IN>;# if($debug){print STDERR "1$l"}
	$l = <IN>; #if($debug){print STDERR "2$l"}
	my ($first, $second) = (0, 0);
	for my $i(1..5){
		$l = <IN>;
		chomp $l  ;
		my @a = split "\t", $l ;
		die unless ($a[0] eq "chr".$i);
		
		my ($t1, $t2) = get_two_value($a[-1]);
		$first += $t1;
		$second+= $t2;
		#if($debug){
		#	print STDERR "depth>=1: $t1 $t2\n";
		#}
	}
	if($debug){
		print STDERR "$first / $second", "\n\n";
	}

	my $dep1 = sprintf("%.3f", $first / $second * 100) . "%";
	
	for my $i(1..6){$l = <IN>}
	
	($first, $second) = (0, 0);
	for my $i(1..5){
		$l = <IN>;
		chomp $l ;
		my @a = split "\t", $l ;
		die unless ($a[0] eq "chr".$i);
		my ($t1, $t2) = get_two_value($a[-1]);
		$first += $t1;
		$second+= $t2;
		
	}
	if($debug){
		print STDERR "$first / $second", "\n\n";
	}
	my $dep2 = sprintf("%.3f", $first / $second * 100) . "%";
	
	for my $i(1..6){$l = <IN>}
	
	($first, $second) = (0, 0);
	for my $i(1..5){
		$l = <IN>;
		chomp $l ;
		my @a = split "\t", $l;
		die unless ($a[0] eq "chr".$i);
		my ($t1, $t2) = get_two_value($a[-1]);
		$first += $t1;
		$second+= $t2;
	}
	if($debug ){
		print STDERR "$first / $second", "\n\n";
	}
	my $avg_dep = sprintf("%.3f", $first / $second );
	
	
	for my $i(1..6){$l = <IN>}
	
	($first, $second) = (0, 0);
	
	
	my ($CG_first, $CG_second) = (0, 0);
	my ($CHG_first, $CHG_second) = (0, 0);
	my ($CHH_first, $CHH_second) = (0, 0);
	for my $i(1..5){
		$l = <IN>;
		chomp $l;
		my @a = split "\t",$l ;
		die unless ($a[0] eq "chr".$i);
		my ($t1, $t2) = get_two_value($a[-1]);
		$first += $t1;
		$second+= $t2;
		
		for my $j(1..2){
			my ($temp1, $temp2) = get_two_value($a[$j]);
			$CG_first += $temp1;
			$CG_second += $temp2;
		}
		for my $j(3..4){
			my ($temp1, $temp2) = get_two_value($a[$j]);
			$CHG_first += $temp1;
			$CHG_second += $temp2;
		}
		for my $j(5..6){
			my ($temp1, $temp2) = get_two_value($a[$j]);
			$CHH_first += $temp1;
			$CHH_second += $temp2;
		}
		
	}
	if($debug ){
		print STDERR "$first / $second", "\n";
		print STDERR "$CG_first / $CG_second ", "\n";
		print STDERR "$CHG_first / $CHG_second ", "\n";
		print STDERR "$CHH_first / $CHH_second ", "\n";
	}
	my $cg = sprintf("%.3f", $CG_first / $CG_second * 100) . "%";
	my $chg = sprintf("%.3f", $CHG_first / $CHG_second * 100) . "%";
	my $chh = sprintf("%.3f", $CHH_first / $CHH_second * 100) . "%";
	my $meth_level = sprintf("%.3f", $first / $second * 100) . "%";
	$l = <IN>;
	die unless ($l =~
	 /chrC/);
	chomp $l;
	my @a = split "\t", $l;
	my ($t1, $t2) = get_two_value($a[-1]);
	
	my $error = $t1 / $t2;
	close(IN);
	
	print OUT join("\t", ($pre,$avg_dep, $dep1, $dep2, $cg, $chg, $chh, $meth_level, $error )), "\n";
}
close(OUT);
#print OUT join("\t", ("sample", "average_depth", "%_of_dep>=1",  "%_of_dep>=2", "CG_meth_level", "CHG_meth_level", "CHH_meth_level", "total_meth_level","error_rate")), "\n";

exit;

sub get_two_value{
	my ($item) = @_;
	my @a = split "=", $item;
	return (split /\//, $a[0]);
}