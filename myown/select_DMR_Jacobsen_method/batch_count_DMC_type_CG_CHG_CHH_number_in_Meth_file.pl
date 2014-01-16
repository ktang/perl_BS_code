#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;



my $debug = 0;

#	print STDERR "\n\n nodupl tag\n\n";


my $script = "/Users/tang58/Kai_BS/myown/select_DMR_Jacobsen_method/count_DMC_type_CG_CHG_CHH_number_in_Meth_file.pl";
die unless (-e $script);

print STDERR "used script \n\n  $script \n\n";


#my $isMeth_dir = "/Users/tang58/misc/Huiming_Zhang/9_PolIV_PolV_ros1_rdd_ibm1/Oct16_Jacobsen_method_all_dep4/isMeth_dep4_6samples";
my $isMeth_dir = "../isMeth_dep4_6samples";
die unless (-d $isMeth_dir);

#my $wt_file = "/Users/tang58/misc/Huiming_Zhang/9_PolIV_PolV_ros1_rdd_ibm1/Oct16_Jacobsen_method_all_dep4/isMeth_dep4_6samples/colA_isMeth_depth4.txt";
my $wt_file = "../isMeth_dep4_6samples/colA_isMeth_depth4.txt";
 
my %files = ( "JKZ18_rdd_vs_colA_hyper_JacobsenMethod_Bin200_sliding50_gap100_BH_adj_p0.05_depth4_annotated_TAIR10.txt" => "JKZ18_rdd_isMeth_depth4.txt", 
			  "JKZ18_rdd_vs_colA_hypo_JacobsenMethod_Bin200_sliding50_gap100_BH_adj_p0.05_depth4_annotated_TAIR10.txt" => "JKZ18_rdd_isMeth_depth4.txt", 
			  "JKZ3_ros1-4_vs_colA_hyper_JacobsenMethod_Bin200_sliding50_gap100_BH_adj_p0.05_depth4_annotated_TAIR10.txt" => "JKZ3_ros1-4_isMeth_depth4.txt", 
			  "JKZ3_ros1-4_vs_colA_hypo_JacobsenMethod_Bin200_sliding50_gap100_BH_adj_p0.05_depth4_annotated_TAIR10.txt" => "JKZ3_ros1-4_isMeth_depth4.txt", 
			  "ibm1-4A_vs_colA_hyper_JacobsenMethod_Bin200_sliding50_gap100_BH_adj_p0.05_depth4_annotated_TAIR10.txt" => "ibm1-4A_isMeth_depth4.txt", 
			  "ibm1-4A_vs_colA_hypo_JacobsenMethod_Bin200_sliding50_gap100_BH_adj_p0.05_depth4_annotated_TAIR10.txt" => "ibm1-4A_isMeth_depth4.txt", 
			  "nrpd1-3A_vs_colA_hyper_JacobsenMethod_Bin200_sliding50_gap100_BH_adj_p0.05_depth4_annotated_TAIR10.txt" => "nrpd1-3A_isMeth_depth4.txt", 
			  "nrpd1-3A_vs_colA_hypo_JacobsenMethod_Bin200_sliding50_gap100_BH_adj_p0.05_depth4_annotated_TAIR10.txt" => "nrpd1-3A_isMeth_depth4.txt", 
			  "nrpe1-11A_vs_colA_hyper_JacobsenMethod_Bin200_sliding50_gap100_BH_adj_p0.05_depth4_annotated_TAIR10.txt" => "nrpe1-11A_isMeth_depth4.txt", 
			  "nrpe1-11A_vs_colA_hypo_JacobsenMethod_Bin200_sliding50_gap100_BH_adj_p0.05_depth4_annotated_TAIR10.txt" => "nrpe1-11A_isMeth_depth4.txt"
			);

my $usage = "$0 \n <indir> <outdir>  \n\n";
die $usage unless(@ARGV == 2);

my $indir = shift or die;
my $outdir = shift or die;
die unless (-d $indir);
die unless (-d $outdir);

#my $cutoff = shift or die;
#my $wt_file = shift or die;
die unless (-e $wt_file);

#my $wt_label = shift or die;

#opendir(DIR, $indir);
#my @files = grep /isMeth.+\.txt$/, readdir DIR;
#closedir (DIR);


foreach my $file(sort keys %files){
	my $input =  File::Spec->catfile($indir, $file);
	die unless (-e $input);
	
	my $isMeth_file_mut = File::Spec->catfile($isMeth_dir , $files{$file});
	die unless (-e $isMeth_file_mut);
	
	my $output;
	if($file =~ /(\S+)\.txt$/){
		$output = File::Spec->catfile($outdir, $1. "_DMC_type.txt");
	}
	
	#<bed_like_file> <isMeht_file_WT> <isMeth_file_mut> <output>
	my $cmd = "time perl $script $input $wt_file $isMeth_file_mut $output";

	
	print STDERR "\n", $cmd, "\n\n";
	if(!$debug){
		`$cmd`;
	}
}

exit;
