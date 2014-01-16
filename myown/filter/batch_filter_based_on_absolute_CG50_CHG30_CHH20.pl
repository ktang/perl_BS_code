#!/usr/bin/perl -w

#perl  batch_filter_based_on_absolute_CG50_CHG30_CHH20.pl /Users/tang58/misc/Huiming_Zhang/Sep12_rrp6L_list_Jacobsen_method/original/list/annotated list/

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;

#my $fold = 1.5;
#my $diff_cutoff = 5;


#my $len_cutoff  = 100;
my $len_cutoff  = 0.5;
my $CG_cutoff  = 50;
my $CHG_cutoff  = 30;
my $CHH_cutoff  = 20;

my $script = "/Users/tang58/Kai_BS/myown/filter/filter_absolute_CG50_CHG30_CHH20.pl";
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
		my $output = File::Spec->catfile($outdir, $1. "_absolute" . "_CG" . $CG_cutoff   ."_CHG". $CHG_cutoff ."_CHH". $CHH_cutoff .".txt");
		die if (-e $output);
	#	my $pre = $1;
		my $cmd = "time perl $script $input  $len_cutoff $CG_cutoff $CHG_cutoff $CHH_cutoff > $output";
		print STDERR $cmd, "\n\n";
		if(!$debug){
			`$cmd`;
		}
	}
}

exit;