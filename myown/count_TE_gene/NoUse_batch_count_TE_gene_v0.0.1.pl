#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;

my $debug = 0 ;

print STDERR "\n\n nodupl tag\n\n";


my $script = "/Users/tang58/Kai_BS/myown/count_TE_gene/count_num_in_annotation_gene_TE_intergenic_v0.0.1.pl";
die unless (-e $script);

my $usage = "$0 <indir> <output>";
die $usage unless(@ARGV == 2);

my $indir = shift or die;
my $output = shift or die;
die unless (-d $indir);

opendir(DIR, $indir);
my @files = grep /\.txt$/, readdir DIR;
closedir DIR;


#my $dir = "/Volumes/My_Book/20120427_ShangHai_data/handled_data/cytosine_coverage_info";

foreach my $file(@files){
#	my $cmd = "time perl $script $indir $outdir $pre $file";
	my $input = File::Spec->catfile($indir, $file);
	die unless (-e $input);
	my $cmd = "perl $script $input >> $output";
	print STDERR $cmd, "\n\n";
	if(!$debug){
		`$cmd`;
	}
}

exit;