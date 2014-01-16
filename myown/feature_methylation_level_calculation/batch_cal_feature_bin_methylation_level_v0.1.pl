#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);


#input outdir pre p-value
use strict;
use File::Spec;

my $debug = 0;

#my $indir = "/Volumes/My_Book/20120427_ShangHai_data/call_methylation/1_isMeth";
#die unless(-d $indir);
#my $p_value = 0.01;

my %lists = ( "all_TE" => "/Users/tang58/misc/Chengguo_Duan/2_nucleosome_analysis/results/TE_length/all_TE_coordinate_bed.txt",
			  "all_protein_coding_gene" => "/Users/tang58/misc/Chengguo_Duan/2_nucleosome_analysis/results/TE_length/all_protein_coding_gene_coordinate_bed.txt"
			);


my $usage = "$0 \n<list_label = all_TE / all_protein_coding_gene> <isMeth_dir> <outdir>\n\n";
die $usage unless(@ARGV == 3);

#my $script_dir = "/Users/tang58/misc/Chengguo_Duan/2_nucleosome_analysis/results/src/TE_bin_methylation_level";
#my $script = File::Spec->catfile($script_dir, "1_cal_methylation_level_for_features_and_flanking_regions_for_SH_isMeth_file_format.pl");
my $script = "/Users/tang58/Kai_BS/myown/feature_methylation_level_calculation/1_cal_methylation_level_for_features_and_flanking_regions_for_SH_isMeth_file_format.pl";
die unless (-e $script);


my $list_label = shift or die;

#my $indir = shift or die "indir";
#my $outdir = shift or die "outdir";
#my $p_cutoff = shift or die "p_value_cutoff";
my $indir = shift or die "isMeth_indir";
my $outdir = shift or die "outdir";

die unless (-d $indir and -d $outdir);
die unless (defined $lists{$list_label});
my $input_list = $lists{$list_label};


#opendir(DIR, $indir) or die;
#my @files = grep/colA_vs\S+diff_bases\S+\.txt/, readdir DIR;
#closedir(DIR);

#my @pres = ( "037306005190A", "Num12A", "Num27A", 
#			 "Num77A", "akn1_2A", "ape_1A", "ape1l_2A",
#			 "arp_1A", "ibm1-4A", "Jmij-19A",
#			 "nrpd1-3A", "nrpel-11A", "P261A", "P313A",
#			 "P314A", "P31A", "pkl_1A", "shh1A"	); # no colA
# my @pres = ("pkl_1A", "shh1A", "nrpd1-3A", "nrpel-11A");

my @pres = ("colA_nodupl2", "nrpd1-3A_nodupl", "nrpel-11A_nodupl", "rrp6L1-2A_nodupl");

my $feature_bin_num = 20;
my $flaking_bp = 2000;
my $flaking_bin_num = 20;


foreach my $pre(@pres){
	my $file = File::Spec->catfile($indir, $pre. "_isMeth_depth5.txt");
	die unless (-e $file);
}

# <db_label> <input_db> <feature_bin_num> <flaking_bp> <flaking_bin_num> <input_list> <outdir> <outpre> 
foreach my $pre(@pres){
	my $file = File::Spec->catfile($indir, $pre. "_isMeth_depth5.txt");
	die unless (-e $file);
	
	
	my $cmd = "perl $script $pre $file $feature_bin_num $flaking_bp $flaking_bin_num $input_list $outdir $list_label";
	print STDERR "\n", $cmd, "\n\n";
	if(!$debug){
		`$cmd`;
	}
}


exit;