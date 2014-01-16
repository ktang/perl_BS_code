#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;

my $debug = 0;

print STDERR "\n\n script may be needed modify\n\n grep \n\n";

if($debug){
	print STDERR "debug = 1\n\n";
}

my $winSize 		  = 200;
my $sliding_bp 		  = 50;
my $merged_gap_cutoff = 100;
my $netDMC_cutoff	  = 10;
my $percentage_cutoff = 100;


#my $script = "/Users/tang58/Kai_BS/myown/isMeth2wig_v0.2.pl";
my $script = "/Users/tang58/Kai_BS/myown/select_DMR_method5/select_DMR_output_list_for_individual_method5_output_list_v0.1.pl";
die unless (-e $script);

#my $control_wig = "/Volumes/My_Book/20120427_ShangHai_data/nodupl_July30/colA/colA_nodupl_chrC_wig/colA_nodupl2_mC.wig";
#my $ctrl_label  = "colA_nodupl2";

#my $control_wig = "/Volumes/My_Book/20120702_SH_extra_BS_Seq/nodupl/nodupl_chrC_error_rate_wig/col0B_nodupl_mC.wig";
#my $ctrl_label  = "col0B_nodupl";

#my $control_wig = "/Volumes/My_Book/20120723_HSZ12204_ARAyagH/nodupl/DMR_list/rrp6l/used/1_compare_with_PolIV_PolV/depth2_colA_nrpd_nrpe_rrp6L/wig/colA_nodupl2_depth2.wig";
#my $ctrl_label  = "colA_nodupl2_depth2";

my $control_wig = "/Volumes/My_Book/20120723_HSZ12204_ARAyagH/nodupl/DMR_list/rrp6l/used/1_compare_with_PolIV_PolV/depth5_colA_nrpd_nrpe_rrp6L/wig/colA_nodupl2_depth5.wig";
my $ctrl_label  = "colA_nodupl2_depth5";

die unless (-e $control_wig);


my $usage = "$0 <indir> <outdir>";
die $usage unless(@ARGV == 2);

my $indir = shift or die;
my $outdir = shift or die;
die unless (-d $indir);
die unless (-d $outdir);

opendir(DIR, $indir);
#my @files = grep /isMeth.+\.txt$/, readdir DIR;
my @files = grep /\.wig$/, readdir DIR;
closedir DIR;


my %pres;

foreach my $file(@files){
	if($file =~ /col/i){
		print STDERR "pass ", $file, "\n\n";
		next;
	}
	if($file =~ /(\S+)\.wig$/){
		my $pre = $1;
		$pres{$pre} = $file;
	}else{
		die $file;
	}
}


#my $dir = "/Volumes/My_Book/20120427_ShangHai_data/handled_data/cytosine_coverage_info";
# <control_pre> <control_input> <input> <outdir> <pre>
# <winSize> <sliding_bp> <merged_gap_cutoff> <netDMC_cutoff> <percentage_cutoff>

#<control_pre> <control_input> <input> <outdir> <pre> 
#<winSize> <sliding_bp> <merged_gap_cutoff> <big_score_cutoff> <increase_percentage_cutoff>
# 200-50-100-10-100

foreach my $pre(sort keys %pres){
	my $file = $pres{$pre};
#	my $cmd = "time perl $script $indir $outdir $pre $file";
#	my $cmd = "time perl $script $indir $outdir $pre" . "_nodupl". " $file";
	my $input = File::Spec->catfile($indir, $file);
	die unless (-e $input);
	my $cmd = "perl $script $ctrl_label $control_wig $input $outdir $pre $winSize $sliding_bp $merged_gap_cutoff $netDMC_cutoff $percentage_cutoff";
	print STDERR $cmd, "\n\n";
	if(!$debug){
		`$cmd`;
	}
}

exit;
