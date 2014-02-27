#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;

#my $script = "/Users/tang58/Kai_BS/myown/cytosines_coverage/cal_cytosines_coverage_for_acgt_count_output.pl";
my $script = "/Users/tang58/Kai_BS/for_publish/9_beCalled_cal_methylation_level_for_features_and_flanking_regions_for_SH_isMeth_file_format_AllDepth_cutoff_v0.1_original_C_not_EckerSuggest.pl";

die unless (-e $script);

my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 \n <indir_isMeth> <bed_dir> <outdir>  <postfix(XXXPaper)> \n\n";
die $usage unless(@ARGV == 4);

my $indir = shift or die;
die unless (-e $indir);

my $bed_dir = shift or die;
die unless (-e $bed_dir);

my $outdir = shift or die;
die unless (-e $outdir);

my $postfix = shift or die "postfix";

#my $output = shift or die;

opendir(DIR, $indir) or die;
my @files = grep /_isMeth\S*\.txt$/, readdir DIR;
closedir DIR;

print STDERR join("\n", @files), "\n\n";

#die if(-e $output);

# <isMeth_file> <sample_label> <feature_bin_num(20)> <flaking_bp(2000)> <flaking_bin_num(20)> <input_bed_directory> <outdir>\n\n";

my ( $feature_bin_num , $flaking_bp, $flaking_bin_num ) = (20, 2000, 20);

foreach my $file (@files){
	if($file =~ /(\S+)_isMeth\S*\.txt$/){
		my $pre = $1;
		my $input = File::Spec->catfile($indir, $file);
		die unless (-e $input);
	#	my $cmd = "perl $script  $input  $pre no >> $output";
		my $cmd = "time perl $script $input $pre $feature_bin_num $flaking_bp $flaking_bin_num $bed_dir $outdir $postfix";
		print STDERR $cmd, "\n\n";
		if(!$debug){
			`$cmd`;
		}
	}
}

exit;
