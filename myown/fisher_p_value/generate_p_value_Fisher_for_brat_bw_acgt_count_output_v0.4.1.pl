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
#use Statistics::R;

print STDERR "\n\nuse perl cal p-value\n\n";


#v0.4.1
# modified for v0.4, use perl to cal p-val

#v0.4
# use R to calculate p-value
# restore depth cutoff 4-100 from v0.3 (same as v0.2)

#v0.3 
# no cutoff for depth
# command line paramter change

# this script, 
# 1 only calculate depth $min_depth - $max_depth for both
# 2 only output p-value <= $p_cutoff;


my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}

my $p_cutoff = 0.05;
my $min_depth = 4;
my $max_depth = 100;
#my $const = 5000000;

print STDERR "depth between $min_depth and $max_depth\n\n\n";

#my $usage = "$0 <WT_dir> <WT_pre> <mut_dir> <mut_pre> <outdir> <outpre>";
my $usage = "$0 \n<WT_dir> <WT_pre> <mut_dir> <mut_pre> <outdir>\n\n";
die $usage unless (@ARGV == 5);

my $WT_dir  = shift or die "WT_dir" ;
my $WT_pre  = shift or die "WT_pre" ;

my $mut_dir = shift or die "mut_dir";
my $mut_pre = shift or die "mut_pre";

my $outdir  = shift or die "outdir" ;
#my $outpre  = shift or die "outpre" ;


my $colA_forw = File::Spec->catfile($WT_dir, $WT_pre . "_forw.txt");
my $colA_rev  = File::Spec->catfile($WT_dir, $WT_pre . "_rev.txt");

die unless (-e $colA_forw and -e $colA_rev);


my $in_forw = File::Spec->catfile($mut_dir, $mut_pre . "_forw.txt");
my $in_rev  = File::Spec->catfile($mut_dir, $mut_pre . "_rev.txt");

#if($debug){
	print  STDERR "input files:\n";
	print  STDERR join("\n", ($in_forw, $in_rev)), "\n\n";
	print STDERR "WT_file:\n";
	print STDERR join("\n", ($colA_forw, $colA_rev)) , "\n\n";
	
#}

die "wrong input" unless (-e $in_forw and -e $in_rev);

#my $output = File::Spec->catfile($outdir,$WT_pre ."_vs_". $mut_pre . "_diff_bases_P" . $p_cutoff . "_NoDepthCutoff" . ".txt");
my $output = File::Spec->catfile($outdir,$WT_pre ."_vs_". $mut_pre . "_diff_bases_P" . $p_cutoff . "_depth_" . $min_depth . "_" . $max_depth . "_both.txt");
die "$output exists" if (-e $output);

print STDERR "output:\t$output\n\n";

if($debug){
	exit;
}

die "output exists" if(-e $output);
open(OUT, ">$output") or die;
close(OUT);

cal_p($colA_forw, $in_forw, $output, $min_depth , $max_depth, $p_cutoff );
cal_p($colA_rev,  $in_rev,  $output, $min_depth , $max_depth, $p_cutoff );

exit;

#chr1    0       0       CHH:0   0       +
#chr1    1       1       CHH:0   0       +
#chr1    2       2       CHH:0   0       +
#  0	 1		 2			3	 4		 5

sub cal_p{

	my ($ctrl_file, $mut_file, $outfile, $min_depth_sub , $max_depth_sub, $p_cutoff_sub) = @_;

	die unless (-e $ctrl_file );
	die unless (-e $mut_file  );
	die unless (-e $outfile   );

	my $const = 5000000;

	my $i = 0;
	my ($l_wt, $l_mut);

	open(WT,  $ctrl_file) or die;
	open(MUT, $mut_file ) or die;
	
	open(OUT, ">>$outfile") or die;

#	my $R = Statistics::R->new();

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
 	
 		die $l_wt, $l_mut unless ( $a_mut[1] == $a_wt[1]  );

 		my ($npp, $n1p,    $np1,      $n11);

		if($wt_depth >= $min_depth_sub and $mut_depth >= $min_depth_sub and $wt_depth <= $max_depth_sub  and $mut_depth <= $max_depth_sub ){
			my $per_wt  = $a_wt[4];
			my $per_mut = $a_mut[4];
		
			my $wt_mC = round($wt_depth  * $per_wt );
			my $mut_mC= round($mut_depth * $per_mut);
			
		#	my $p_value = 1;
		#	$R->set('a', $wt_mC);
    	#	$R->set('b', $wt_depth - $wt_mC);
   		#	$R->set('c', $mut_mC);
   		#	$R->set('d', $mut_depth - $mut_mC);
		#	$R->run(q`p = fisher.test(matrix( c(a,b,c,d ), ncol=2))$p.value`);
		#	$p_value = $R->get('p');
		
			$npp = $wt_depth + $mut_depth;
			$n1p = $wt_mC + $mut_mC;
			$np1 = $wt_depth;
			$n11 = $wt_mC;
			next if($npp == 0);
			my  $p_value = calculateStatistic( n11=>$n11,
											   n1p=>$n1p,
							    			   np1=>$np1,
							    			   npp=>$npp);
			my $errorCode;
			if( ($errorCode = getErrorCode())){
				print STDERR $errorCode." - ".getErrorMessage();
				die "$l_wt, $l_mut";
  			}
  	#	elsif($p_value <= $p_cutoff){
  			else{    			   
	 			if($p_value <= $p_cutoff_sub){
					print OUT join("\t", ($chr, $pos, $strand, $type, $wt_mC, $wt_depth, $per_wt, $mut_mC, $mut_depth, $per_mut, $p_value)), "\n";
  				}
  			}	
		}	
	}
#	$R->stop();

	close(WT);
	close(MUT);
 	close(OUT);
}



sub round{
    my($number) = shift;
    return int($number + .5 * ($number <=> 0)); # take care of negative numbers too
}
