#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;

#my $script = "/Users/tang58/Kai_BS/myown/cytosines_coverage/cal_cytosines_coverage_for_acgt_count_output.pl";
my $script = "/Users/tang58/Kai_BS/for_publish/cal_meth_level_for_isMeth_v0.0.pl";

die unless (-e $script);

my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 \n <indir> <output> \n\n";
die $usage unless(@ARGV == 2);

my $indir = shift or die;
die unless (-e $indir);

#my $outdir = shift or die;
#die unless (-e $outdir);

my $output = shift or die;

opendir(DIR, $indir) or die;

my @files = grep /_isMeth_\S*\.txt$/, readdir DIR;
closedir DIR;

print STDERR join("\n", @files), "\n\n";

die if(-e $output);

if (!$debug){
	open(OUT, ">>$output") or die;

	print OUT join("\t", ( "sample",
			  "weighted_CG_level", "weighted_CHG_level", "weighted_CHH_level", "weighted_C_level",
			  "weighted_CG_div", "weighted_CHG_div", "weighted_CHH_div", "weighted_C_div",
			  
			  "frac_CG_level", "frac_CHG_level", "frac_CHH_level", "frac_C_level",
			  "mCG_num", "mCHG_num", "mCHH_num", "mC_num",
			  "CG_num", "CHG_num", "CHH_num", "C_num",
			  
			  "sum_CG_level", "sum_CHG_level", "sum_CHH_level", "sum_C_level",
			  "mean_CG_level", "mean_CHG_level", "mean_CHH_level", "mean_C_level"
			  
			  ) )  , "\n";

	close OUT;
}

# <isMeth_file> <sample_name> <head_yes_no> STDOUT

foreach my $file (@files){
	if($file =~ /(\S+)_isMeth_\S*\.txt$/){
		my $pre = $1;
		my $input = File::Spec->catfile($indir, $file);
		die unless (-e $input);
		my $cmd = "perl $script  $input  $pre no >> $output";
		print STDERR $cmd, "\n\n";
		if(!$debug){
			`$cmd`;
		}
	}
}

exit;
