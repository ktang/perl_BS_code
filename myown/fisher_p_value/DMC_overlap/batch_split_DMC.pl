#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;

my $script = "/Users/tang58/Kai_BS/myown/fisher_p_value/DMC_overlap/split_DMC.pl";
die unless (-e $script);

my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 \n <indir> <outdir>\n\n";
die $usage unless(@ARGV == 2);

my $indir = shift or die;
die unless (-d $indir);

my $outdir = shift or die;
die unless (-d $outdir);

opendir(DIR, $indir) or die;

my @files = grep /\.txt$/, readdir DIR;
closedir DIR;


foreach my $file (@files){
	if($file =~ /(\S+)\.txt$/){
		my $pre = $1;
		
		my $input = File::Spec->catfile($indir, $file);
		die unless (-e $input);
		
		#my $usage = "$0 <dir> <WT_pre> <mut_pre> <outdir>";
		my $cmd = " perl $script $input $outdir $pre ";
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
