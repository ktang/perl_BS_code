#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;

my $script = "/Users/tang58/Kai_BS/myown/fisher_p_value/generate_p_value_Fisher_for_brat_bw_acgt_count_output_v0.2.pl";

die unless (-e $script);

my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 <indir> <outdir>";
die $usage unless(@ARGV == 2);

my $indir = shift or die;
die unless (-d $indir);

my $outdir = shift or die;
die unless (-d $outdir);

opendir(DIR, $indir) or die;

my @files = grep /_rev\.txt$/, readdir DIR;
closedir DIR;

my $WT_pre = "";


foreach my $file (@files){
	if($file =~ /(\S+)_rev\.txt$/){
		my $pre = $1;
		$WT_pre = $pre if ($pre =~ /col/);
	}
}

foreach my $file (@files){
	if($file =~ /(\S+)_rev\.txt$/){
	
	
		my $pre = $1;
		next if ($pre =~ /col/);
		#my $usage = "$0 <dir> <WT_pre> <mut_pre> <outdir>";
		my $cmd = "time perl $script $indir $WT_pre $pre $outdir";
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
