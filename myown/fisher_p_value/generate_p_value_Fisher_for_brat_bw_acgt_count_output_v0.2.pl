#! /usr/bin/perl -w
#		   		  WT	    mut
#				  word2   ~word2
#mC	 	 word1    n11      n12 | n1p
#No.T	~word1    n21      n22 | n2p
#          --------------
#				  np1      np2   npp
use strict;
use File::Spec;
use Text::NSP::Measures::2D::Fisher2::twotailed;

# this script, 
# 1 only calculate depth $min_depth - $max_depth for both
# 2 only output p-value <= $p_cutoff;




#v1.0 use Dr. Liu's output as input to prove the correctness
#of the code.

#v1.1 use Col-0 high WT

#v1.2 use C24 input
#v1.2.1 as 1.2 is slow, so 1.2.1 only consider depth >=$min_depth and <= $max_depth.
#my ($npp, $n1p,    $np1,      $n11);
#   total  mC_sum   WT_depth  WT_mC;

my $debug = 0;
if($debug){
	print STDERR "debug = 1\n\n";
}

my $p_cutoff = 0.05;
my $min_depth = 4;
my $max_depth = 100;
my $const = 5000000;

print STDERR "depth between $min_depth and $max_depth\n\n\n";

#my $colA_forw = "/home/becky/SH_data_analysis/data_set/colA/colA_forw.txt";
#my $colA_rev  = "/home/becky/SH_data_analysis/data_set/colA/colA_rev.txt";

#my $colA_forw = "/Volumes/My_Book/20120427_ShangHai_data/handled_data/colA_trim/colA_brat_output/colA_forw.txt";
#my $colA_rev  = "/Volumes/My_Book/20120427_ShangHai_data/handled_data/colA_trim/colA_brat_output/colA_rev.txt";


my $usage = "$0 <dir> <WT_pre> <mut_pre> <outdir>";
die $usage unless (@ARGV == 4);

my $indir  = shift or die "indir";
my $WT_pre = shift or die "WT_pre";
my $pre	   = shift or die "pre";
my $outdir = shift or die "outdir";

my $colA_forw = File::Spec->catfile($indir, $WT_pre . "_forw.txt");
my $colA_rev  = File::Spec->catfile($indir, $WT_pre . "_rev.txt");

die unless (-e $colA_forw and -e $colA_rev);


die "wrong indir" unless (-d  $indir);
die "wrong outdir" unless (-d $outdir);

my $in_forw = File::Spec->catfile($indir, $pre . "_forw.txt");
my $in_rev  = File::Spec->catfile($indir, $pre . "_rev.txt");

if($debug){
	print  STDERR "input files:\n";
	print  STDERR join("\n", ($in_forw, $in_rev)), "\n\n";
	print STDERR "WT_file:\n";
	print STDERR join("\n", ($colA_forw, $colA_rev)) , "\n\n";
	
}

die "wrong input" unless (-e $in_forw and -e $in_rev);



my $output = File::Spec->catfile($outdir,$WT_pre ."_vs_". $pre . "_diff_bases_P" . $p_cutoff . "_depth_" . $min_depth . "_" . $max_depth . "_both.txt");
die "$output exists" if (-e $output);

if($debug){
	print STDERR "output:\t$output\n\n";
	exit;
}

open(OUT, ">$output") or die;

#chr1    0       0       CHH:0   0       +
#chr1    1       1       CHH:0   0       +
#chr1    2       2       CHH:0   0       +
#  0	 1		 2			3	 4		 5


my $i = 0;
my ($l_wt, $l_mut);

open(WT, $colA_forw) or die;
open(MUT, $in_forw) or die;


while ($l_wt = <WT>, $l_mut = <MUT>){
	$i++;
	if($i % $const == 0){
		print STDERR $i , " done\n";
	}
 	chomp $l_wt;
  	chomp $l_mut;

	my @a_wt  = split "\t", $l_wt;
 	my @a_mut = split "\t", $l_mut;
 	my ($type_wt, $wt_depth) = split ":" , $a_wt[3];
 	my ($type_mut,$mut_depth) = split ":", $a_mut[3];
 	
 	my ($chr, $pos, $strand, $type) = ($a_wt[0], $a_wt[1] + 1, $a_wt[5], $type_wt );
 	
 	next if($chr eq "chrC" or $chr eq "chrM" or $chr eq "pUC19");
 	
 	die $l_wt, $l_mut unless ( $a_mut[0] eq $a_wt[0]  and $a_mut[1] == $a_wt[1]  );
 	my ($npp, $n1p,    $np1,      $n11);

	if($wt_depth >= $min_depth and $mut_depth >= $min_depth and $wt_depth <= $max_depth  and $mut_depth <= $max_depth ){
		my $per_wt  = $a_wt[4];
		my $per_mut = $a_mut[4];
		
		my $wt_mC = round($wt_depth  * $per_wt );
		my $mut_mC= round($mut_depth * $per_mut);
	
		$npp = $wt_depth + $mut_depth;
		$n1p = $wt_mC + $mut_mC;
		$np1 = $wt_depth;
		$n11 = $wt_mC;
		next if($npp == 0);
		my  $p_value = calculateStatistic( n11=>$n11,
										   n1p=>$n1p,
						    			   np1=>$np1,
						    			   npp=>$npp);
						    			   
		#if($debug){
		#	if($i <= 20){
		#		print join("\t", ($chr, $pos, $strand, $type, $wt_mC, $wt_depth, $per_wt, $mut_mC, $mut_depth, $per_mut, $p_value)), "\n";
		#	}
		#}
						    			   
		my $errorCode;
		if( ($errorCode = getErrorCode())){
			print STDERR $errorCode." - ".getErrorMessage();
			die "$l_wt, $l_mut";
  		}
  	#	elsif($p_value <= $p_cutoff){
  		else{
  			if($p_value <= $p_cutoff){
				print OUT join("\t", ($chr, $pos, $strand, $type, $wt_mC, $wt_depth, $per_wt, $mut_mC, $mut_depth, $per_mut, $p_value)), "\n";
  			}
  		}
		
	}	
}

close(WT);
close(MUT);
 
$i = 0;
open(WT, $colA_rev) or die;
open(MUT, $in_rev) or die;

while ($l_wt = <WT>, $l_mut = <MUT>){
	$i++;
	if($i % $const == 0){
		print STDERR $i , " done\n";
	}
 	chomp $l_wt;
  	chomp $l_mut;

	my @a_wt  = split "\t", $l_wt;
 	my @a_mut = split "\t", $l_mut;
 	my ($type_wt, $wt_depth) = split ":" , $a_wt[3];
 	my ($type_mut,$mut_depth) = split ":", $a_mut[3];
 	
 	my ($chr, $pos, $strand, $type) = ($a_wt[0], $a_wt[1] + 1, $a_wt[5], $type_wt );
 	
 	next if($chr eq "chrC" or $chr eq "chrM" or $chr eq "pUC19");
 	
 	die $l_wt, $l_mut unless ( $a_mut[0] eq $a_wt[0]  and $a_mut[1] == $a_wt[1]  );
 	my ($npp, $n1p,    $np1,      $n11);

	if($wt_depth >= $min_depth and $mut_depth >= $min_depth and $wt_depth <= $max_depth  and $mut_depth <= $max_depth ){
		my $per_wt  = $a_wt[4];
		my $per_mut = $a_mut[4];
		
		my $wt_mC = round($wt_depth  * $per_wt );
		my $mut_mC= round($mut_depth * $per_mut);
	
		$npp = $wt_depth + $mut_depth;
		$n1p = $wt_mC + $mut_mC;
		$np1 = $wt_depth;
		$n11 = $wt_mC;
		next if($npp == 0);
		my  $p_value = calculateStatistic( n11=>$n11,
										   n1p=>$n1p,
						    			   np1=>$np1,
						    			   npp=>$npp);
						    			   
	#	if($debug){
	#		if($i <= 10){
	#			print join("\t", ($chr, $pos, $strand, $type, $wt_mC, $wt_depth, $per_wt, $mut_mC, $mut_depth, $per_mut, $p_value)), "\n";
	#		}
	#	}
						    			   
		my $errorCode;
		if( ($errorCode = getErrorCode())){
			print STDERR $errorCode." - ".getErrorMessage();
			die "$l_wt, $l_mut";
  		}
  	#	elsif($p_value <= $p_cutoff){
  		else{
  			if($p_value <= $p_cutoff){
				print OUT join("\t", ($chr, $pos, $strand, $type, $wt_mC, $wt_depth, $per_wt, $mut_mC, $mut_depth, $per_mut, $p_value)), "\n";
  			}
  		}
		
	}	
}


close(WT);
close(MUT);

close(OUT);


exit;

sub round{
    my($number) = shift;
    return int($number + .5 * ($number <=> 0)); # take care of negative numbers too
}
