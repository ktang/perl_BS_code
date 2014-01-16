#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;

#my $fold = 1.5;
#my $diff_cutoff = 5;


my $len_cutoff  = 100;
my $CG_cutoff  = 30;
my $CHG_cutoff  = 15;
my $CHH_cutoff  = 10;

my $script = "/Users/tang58/Kai_BS/myown/filter/v0.2/filter_len100_CG30_CHG15_CHH10_v0.2.pl";
die unless (-e $script);

print STDERR "script used:\n";
print STDERR $script, "\n\n\n";

my $debug = 0;

#my $p_value = 0.01;

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 \n <indir> <outdir> \n\n";
die $usage unless(@ARGV == 2);

my $indir = shift or die;
die unless (-d $indir);

my $outdir = shift or die;
die unless (-d $outdir);

opendir(DIR, $indir) or die;
#my @files = grep /_both\.txt$/, readdir DIR;
my @files = grep /\.txt$/, readdir DIR;
closedir DIR;


foreach my $file (@files){
	if($file =~ /(\S+)\.txt$/){
		my $input = File::Spec->catfile($indir, $file);
		my $output = File::Spec->catfile($outdir, $1. "_len". $len_cutoff . "_CG" . $CG_cutoff   ."_CHG". $CHG_cutoff ."_CHH". $CHH_cutoff .".txt");
		die if (-e $output);
	#	my $pre = $1;
		my $cmd = "perl $script $input  $len_cutoff $CG_cutoff $CHG_cutoff $CHH_cutoff > $output";
		print STDERR $cmd, "\n\n";
		if(!$debug){
			`$cmd`;
		}
	}
}

exit;