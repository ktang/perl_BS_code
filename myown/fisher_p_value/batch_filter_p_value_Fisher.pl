#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;

#my $script = "/Users/tang58/Kai_BS/myown/fisher_p_value/generate_p_value_Fisher_for_brat_bw_acgt_count_output_v0.2.pl";
#die unless (-e $script);



my $debug = 0;

my $p_cutoff = 0.01;

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
my @files = grep /_vs_\S+\.txt$/, readdir DIR;
closedir DIR;

foreach my $file (@files){
	if($file =~ /(\S+)_P0\.05_(\S+)\.txt$/){
		my $input = File::Spec->catfile($indir, $file);
		my $pre1 = $1;
		my $pre2 = $2;
		
		my $output = File::Spec->catfile($outdir, $pre1. "_P0.01_" . $pre2 . ".txt");
		die if(-e $output);
		print STDERR "input : $file\n";
		print STDERR "output: $output\n\n";
		if(!$debug){
			open(IN, $input) or die;
			
			die if (-e $output);
			open(OUT, ">$output") or die;
			
			while(<IN>){
				chomp;
				my @a = split "\t";
				if($a[-1] <= $p_cutoff){
					print OUT join("\t", @a), "\n";
				}
			}
			close IN;
			close OUT;
		}
	}
}

exit;

sub round {
    my($number) = shift;
    #return int($number + .5);
    return int($number + .5 * ($number <=> 0)); # take care of negative numbers too

}
