#!/usr/bin/perl -w

use strict;
use File::Spec;

#my $usage = "$0 \n <WT_isMeth> <wt_label> <mut_isMeth> <mut_label> <outdir> <outpre>\n\n";
my $script = "/Users/tang58/Kai_BS/Jacobsen_Cell_method/step1_select_DMRs_Jacobsen_Cell_method_v0.0.pl";
die unless (-e $script);

my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 \n <indir_isMeth_mut_only> <WT_isMeth> <WT_label_used_in_file_name>  <outdir> \n\n";
die $usage unless(@ARGV == 4);

my $indir_isMeth = shift or die;
die unless (-d $indir_isMeth);

my $wt_isMet_file = shift or die;
die unless (-e $wt_isMet_file);
my $wt_label = shift or die;

my $outdir = shift or die;
die unless (-d $outdir);


opendir(METDIR, $indir_isMeth) or die;
my @isMet_files = grep /isMeth\S+\.txt$/, readdir METDIR;
closedir METDIR;



# JKZ18_rdd_isMeth
# JKZ18_rdd_vs

my @labels;

for my $i (0..$#isMet_files){
	my $isMet_file = $isMet_files[$i];
	if($isMet_file =~ /(\S+)_isMeth/){ $labels[$i] = $1;	}
}


print STDERR "isMet:\n";
print STDERR join("\n", @isMet_files ), "\n\n";
print STDERR join("\n", @labels ), "\n\n";


print STDERR "wt: $wt_isMet_file \n\n\n";

my $wt_file = $wt_isMet_file;


foreach my $i (0..$#isMet_files){
	
	my $isMet_file = File::Spec->catfile($indir_isMeth ,$isMet_files[$i]);
	my $mut_label = $labels[$i] ;
	my $pre = $mut_label . "_vs_" . $wt_label ;
	
	my $cmd = "perl $script $wt_file $wt_label $isMet_file  $mut_label $outdir $pre";
	#my $usage = "$0 \n <WT_isMeth> <wt_label> <mut_isMeth> <mut_label> <outdir> <outpre> \n\n";

	print STDERR $cmd, "\n\n";
	if(!$debug){
		`$cmd`;
	}

	
}

exit;

