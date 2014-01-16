#!/usr/bin/perl -w

use strict;
use File::Spec;

my $debug = 0;


#my $usage = "$0 <list_label> <list> <outdir>";
#die $usage unless(@ARGV == 3);

my $usage = "$0 \n <list_file> <outdir> <list_label>  \n\n";
die $usage unless(@ARGV == 3);


# my $script_dir = "/Users/tang58/Kai_BS/myown/feature_methylation_level_calculation";
my $script_dir = "/Users/tang58/Kai_BS/myown/Huiming_downstream_analysis";
die unless (-d $script_dir);
my $script = File::Spec->catfile($script_dir, "cal_methylation_level_for_features_and_flanking_regions_for_SH_isMeth_file_format.pl");
die unless (-e $script);


my $list_file = shift or die;
my $outdir = shift or die "outdir";
my $list_label = shift or die;

die unless( -e $list_file );

die unless (-d $outdir);

#my @pres = ("pkl_1A", "shh1A", "nrpd1-3A", "nrpel-11A");
#<sample_label> <db_in1> <db_in2> <feature_bin_num> <flaking_bp> <flaking_bin_num> <input_list> <outdir> <TE_gene_label>



my %isMeth_files = (
#					"idm1-1" => [ "/Volumes/Macintosh_HD_2/idm1_new_met_data/UCR_server_data/Lane8_JKZ12_phd/s_8_isMeth.txt",
#								  "/Volumes/Macintosh_HD_2/idm1_new_met_data/UCR_server_data/Lane7_JKZ13_phd/s_7_isMeth.txt"],
								  
#					"ros1-4" => [ "/Volumes/Macintosh_HD_2/idm1_new_met_data/UCR_server_data/Lane4_JKZ3_ros1_4/s_4_isMeth.txt",
#								  "/Volumes/Macintosh_HD_2/idm1_new_met_data/UCR_server_data/Lane5_JKZ4_ros1_4/s_5_isMeth.txt"],
#
#					"rdd"	 => [ "/Volumes/Macintosh_HD_2/idm1_new_met_data/UCR_server_data/Lane3_JKZ17_rdd/s_3_isMeth.txt",
#								  "/Volumes/Macintosh_HD_2/idm1_new_met_data/UCR_server_data/Lane5_JKZ18_rdd/s_5_isMeth.txt"]
					"colA_nodupl2_dep5" =>
					"/Volumes/My_Book/20120723_HSZ12204_ARAyagH/nodupl/DMR_list/rrp6l/used/1_compare_with_PolIV_PolV/depth5_colA_nrpd_nrpe_rrp6L/isMeth_file/colA_nodupl2_isMeth_depth5.txt",
					 
					"nrpd1-3A_nodupl_dep5" => 
					"/Volumes/My_Book/20120723_HSZ12204_ARAyagH/nodupl/DMR_list/rrp6l/used/1_compare_with_PolIV_PolV/depth5_colA_nrpd_nrpe_rrp6L/isMeth_file/nrpd1-3A_nodupl_isMeth_depth5.txt",
					
					"nrpe1-11A_nodupl_dep5" => 
					"/Volumes/My_Book/20120723_HSZ12204_ARAyagH/nodupl/DMR_list/rrp6l/used/1_compare_with_PolIV_PolV/depth5_colA_nrpd_nrpe_rrp6L/isMeth_file/nrpe1-11A_nodupl_isMeth_depth5.txt",
					
					"rrp6L1-2A_nodupl_dep5" => 
					"/Volumes/My_Book/20120723_HSZ12204_ARAyagH/nodupl/DMR_list/rrp6l/used/1_compare_with_PolIV_PolV/depth5_colA_nrpd_nrpe_rrp6L/isMeth_file/rrp6L1-2A_nodupl_isMeth_depth5.txt"
				  );

my $feature_bin_num = 20;
my $flaking_bp = 2000;
my $flaking_bin_num = 20;


# <db_label> <input_db> <feature_bin_num> <flaking_bp> <flaking_bin_num> <input_list> <outdir> <outpre> 
foreach my $isMeth_label (sort keys %isMeth_files){
#	my ($file_in1, $file_in2 )  = @{ $label_files{$pre} };
	my $isMeth_file = $isMeth_files{$isMeth_label};
	die  $isMeth_label unless (-e $isMeth_file);
	#foreach my $list_label (sort keys %list_files){
	#	my $list_file = $list_files{ $list_label};
		die unless (-e $list_file);	
	
		my $pre = join("_in_", ($list_label, $isMeth_label));
		my $cmd = "time perl $script $feature_bin_num $flaking_bp $flaking_bin_num $isMeth_file $list_file $outdir $pre";
		print STDERR "\n",  $cmd, "\n\n";
		if(!$debug){
			`$cmd`;
		}
	#}
}

exit;
