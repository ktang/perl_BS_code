#!/usr/bin/perl -w

use strict;
use File::Spec;

my $debug = 0;


my $usage = "$0 <list_label> <list> <outdir>";
die $usage unless(@ARGV == 3);

my $script_dir = "/Users/tang58/Kai_BS/myown/feature_methylation_level_calculation";
die unless (-d $script_dir);
my $script = File::Spec->catfile($script_dir, "1_cal_methylation_level_for_features_and_flanking_regions_for_two_replicates.pl");
die unless (-e $script);


my $list_label = shift or die;
my $input_list = shift or die;
my $outdir = shift or die "outdir";

die unless (-d $outdir);

#my @pres = ("pkl_1A", "shh1A", "nrpd1-3A", "nrpel-11A");
#<sample_label> <db_in1> <db_in2> <feature_bin_num> <flaking_bp> <flaking_bin_num> <input_list> <outdir> <TE_gene_label>



my %label_files = (
					"idm1-1" => [ "/Volumes/Macintosh_HD_2/idm1_new_met_data/UCR_server_data/Lane8_JKZ12_phd/s_8_isMeth.txt",
								  "/Volumes/Macintosh_HD_2/idm1_new_met_data/UCR_server_data/Lane7_JKZ13_phd/s_7_isMeth.txt"],
								  
					"ros1-4" => [ "/Volumes/Macintosh_HD_2/idm1_new_met_data/UCR_server_data/Lane4_JKZ3_ros1_4/s_4_isMeth.txt",
								  "/Volumes/Macintosh_HD_2/idm1_new_met_data/UCR_server_data/Lane5_JKZ4_ros1_4/s_5_isMeth.txt"],

					"rdd"	 => [ "/Volumes/Macintosh_HD_2/idm1_new_met_data/UCR_server_data/Lane3_JKZ17_rdd/s_3_isMeth.txt",
								  "/Volumes/Macintosh_HD_2/idm1_new_met_data/UCR_server_data/Lane5_JKZ18_rdd/s_5_isMeth.txt"]
				  );

my $feature_bin_num = 20;
my $flaking_bp = 2000;
my $flaking_bin_num = 20;


foreach my $pre(sort keys %label_files){
	my ($file_in1, $file_in2 )  = @{ $label_files{$pre} };
	die unless (-e $file_in1);
	die unless (-e $file_in2);
}

# <db_label> <input_db> <feature_bin_num> <flaking_bp> <flaking_bin_num> <input_list> <outdir> <outpre> 
foreach my $pre(sort keys %label_files){
	my ($file_in1, $file_in2 )  = @{ $label_files{$pre} };
	
	my $cmd = "time perl $script $pre $file_in1 $file_in2 $feature_bin_num $flaking_bp $flaking_bin_num $input_list $outdir $list_label";
	print STDERR $cmd, "\n\n";
	if(!$debug){
		`$cmd`;
	}
}

exit;