#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;

# my $script = "/Users/tang58/Kai_BS/myown/select_DMR_Ecker_method/select_DMR_output_list_v1.0_sliding_window.pl";
my $script = "/Users/tang58/Kai_BS/IDM1_Science_DMR_one_step_Nov19_2013/2_select_DMR_output_list_IDM1_Science_method.pl";
#"/Users/tang58/Kai_BS/for_publish/select_DMR_output_list_sliding_window_for_publish_v0.0.pl";
die unless (-e $script);

print STDERR "script used:\n";
print STDERR $script, "\n\n\n";

my $debug = 0;

my $anchor_cutoff = 5;#2;
my $final_cutoff  = 10; #5;

#my $p_value = 0.01;

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 \n <indir> <outdir> <p_value> <dep> <postfix(edm2_paper)>\n\n";
die $usage unless(@ARGV == 5);

my $indir = shift or die;
die unless (-d $indir);

my $outdir = shift or die;
die unless (-d $outdir);

my $p_value = shift or die;
my $dep     = shift or die;

my $postfix = shift or die;

opendir(DIR, $indir) or die;
#my @files = grep /_both\.txt$/, readdir DIR;
my @files = grep /\.txt$/, readdir DIR;
closedir DIR;

#<input> <outdir> <pre> <p_value>
#col_0_pooled_vs_12_pooled_diff_bases_P0.01_depth_4_100_both.txt
foreach my $file (@files){
	if($file =~ /(\S+)_vs_(\S+)_diff_.+\.txt$/){
		my $input = File::Spec->catfile($indir, $file);
		
		my $pre = $2 . "_vs_" . $1;
		my $cmd = "time perl $script $input $outdir $pre $p_value $dep $anchor_cutoff $final_cutoff $postfix";
		print STDERR $cmd, "\n\n";
		if(!$debug){
			`$cmd`;
		}
	}
}

exit;