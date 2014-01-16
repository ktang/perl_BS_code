#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;

#my $script = "/Users/tang58/Kai_BS/myown/cytosines_coverage/cal_cytosines_coverage_for_acgt_count_output.pl";
my $script = "/Users/tang58/Kai_BS/myown/cytosines_coverage/cal_cytosines_coverage_for_acgt_count_output_v0.1.pl";

die unless (-e $script);

print STDERR "script used is: \n";

print STDERR $script, "\n\n";

my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 \n<indir> <outdir>\n\n";
die $usage unless(@ARGV == 2);

my $indir = shift or die;
die unless (-e $indir);

my $outdir = shift or die;
die unless (-e $outdir);

opendir(DIR, $indir) or die;

my @files = grep /_rev\.txt$/, readdir DIR;
closedir DIR;

print STDERR join("\n", @files), "\n\n";

foreach my $file (@files){
	if($file =~ /(\S+)_rev\.txt$/){
		my $pre = $1;
		my $cmd = "time perl $script $indir $pre $outdir";
		print STDERR $cmd, "\n\n";
		if(!$debug){
			`$cmd`;
		}
	}
}

exit;

sub round {
    my($number) = shift;
    #return int($number + .5);
    return int($number + .5 * ($number <=> 0)); # take care of negative numbers too

}
