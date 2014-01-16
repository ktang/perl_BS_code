#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;

my $debug = 0;

print STDERR "\n\n script may be needed modify\n\n";

if($debug){
	print STDERR "debug = 1\n\n";
}

my $winSize 		  = 400;
my $sliding_bp 		  = 100;
my $merged_gap_cutoff = 200;
my $netDMC_cutoff	  = 20;
my $percentage_cutoff = 200;


#my $script = "/Users/tang58/Kai_BS/myown/isMeth2wig_v0.2.pl";
my $script = "/Users/tang58/Kai_BS/myown/select_DMR_method4/select_DMR_output_list_for_individual_method4_output_list_v1.1.pl";
die unless (-e $script);

#my $control_wig = "/Volumes/My_Book/20120427_ShangHai_data/nodupl_July30/colA/colA_nodupl_chrC_wig/colA_nodupl2_mC.wig";
#my $ctrl_label  = "colA_nodupl2";


my $control_wig = "/Volumes/My_Book/20120813_HSZ12204_ARAyagH/C24/wig_nodupl_chrC_error_rate/C24_WT_CF_nodupl_mC.wig";
my $ctrl_label  = "C24_WT_CF_nodupl";

die unless (-e $control_wig);


my $usage = "$0 <indir> <outdir>";
die $usage unless(@ARGV == 2);

my $indir = shift or die;
my $outdir = shift or die;
die unless (-d $indir);
die unless (-d $outdir);

opendir(DIR, $indir);
#my @files = grep /isMeth.+\.txt$/, readdir DIR;
my @files = grep /_mC\.wig$/, readdir DIR;
closedir DIR;


my %pres;

foreach my $file(@files){
	if($file =~ /C24/i){
		print STDERR "pass ", $file, "\n\n";
		next;
	}
	if($file =~ /(\S+)_mC/){
		my $pre = $1;
		$pres{$pre} = $file;
	}else{
		die $file;
	}
}


#my $dir = "/Volumes/My_Book/20120427_ShangHai_data/handled_data/cytosine_coverage_info";
# <control_pre> <control_input> <input> <outdir> <pre>
# <winSize> <sliding_bp> <merged_gap_cutoff> <netDMC_cutoff> <percentage_cutoff>

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
