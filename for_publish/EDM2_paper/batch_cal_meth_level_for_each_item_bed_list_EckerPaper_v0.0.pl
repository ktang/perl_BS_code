#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;

#my $script = "/Users/tang58/Kai_BS/myown/cytosines_coverage/cal_cytosines_coverage_for_acgt_count_output.pl";
my $script = "/Users/tang58/Kai_BS/for_publish/EDM2_paper/cal_meth_level_for_each_item_bed_list_EckerPaper_v0.0.pl";

die unless (-e $script);

my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 \n <indir_isMeth> <bed_list> <outpre> <outdir>  <postfix(XXXPaper)> \n\n";
die $usage unless(@ARGV == 5);

my $indir = shift or die;
die unless (-e $indir);

#my $bed_dir = shift or die;
#die unless (-e $bed_dir);
my $bed_file = shift or die;
die unless ( -e $bed_file );

my $outpre = shift or die;

my $outdir = shift or die;
die unless (-e $outdir);

my $postfix = shift or die "postfix";

#my $output = shift or die;

opendir(DIR, $indir) or die;
my @files = grep /_isMeth_\S*\.txt$/, readdir DIR;
closedir DIR;

print STDERR join("\n", @files), "\n\n";

foreach my $file (@files){
	if($file =~ /(\S+)_isMeth_\S*\.txt$/){
		my $pre = $1;
		my $input = File::Spec->catfile($indir, $file);
		die unless (-e $input);
		my $output = File::Spec->catfile($outdir, $outpre . "_meth_level_in_" . $pre . "_for_". $postfix. ".txt" );
		
	#	my $cmd = "perl $script  $input  $pre no >> $output";
	#	my $cmd = "time perl $script $input $pre $feature_bin_num $flaking_bp $flaking_bin_num $bed_dir $outdir $postfix";
	#m<isMeht_file> <sample_label> <bed_like_file> <output> \n\n";
		my $cmd = "perl $script $input $pre  $bed_file $output";
		print STDERR $cmd, "\n\n";
		if(!$debug){
			`$cmd`;
		}
	}
}

exit;
