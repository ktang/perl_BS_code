#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;
my $debug = 0;

my $trim_tool = "/Users/tang58/Software/BRAT/brat_bw-2.0.1/trim";
die unless (-e $trim_tool);

my $usage = "$0 <indir> <prefix>";
die $usage unless(@ARGV == 2);

my $indir = shift or die;
my $prefix = shift or die;
die unless (-d $indir);

my $q = 20;
my $L = 64;
my $m = 0;

my $in1 = File::Spec->catfile($indir, "1.fq.gz");
my $in2 = File::Spec->catfile($indir, "2.fq.gz");

die unless (-e $in1 and -e $in2);

my $raw_fq1 = $prefix . "_1_raw.fq";
my $raw_fq2 = $prefix . "_2_raw.fq";

my $cmd1 = "time zcat $in1 > $raw_fq1";
print STDERR $cmd1, "\n\n";
if(!$debug){
	die if (-e $raw_fq1);
	`$cmd1`;
}

my $cmd2 = "time zcat $in2 > $raw_fq2";
print STDERR $cmd2, "\n\n";
if(!$debug){
	die if(-e $raw_fq2);
	`$cmd2`;
}


my $trim_cmd = "time $trim_tool -1 $raw_fq1 -2 $raw_fq2 -P $prefix -q $q -L $L -m $m";
print STDERR $trim_cmd,"\n\n";

if(!$debug){
	`$trim_cmd`;
}
exit;