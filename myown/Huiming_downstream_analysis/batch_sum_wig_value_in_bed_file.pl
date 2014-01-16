#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;

my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}

# <col_label> <in_wig> <in_bed_like_file> <output> 
my $script = "/Users/tang58/Kai_BS/myown/Huiming_downstream_analysis/sum_wig_value_in_bed_file.pl";
die unless (-e $script);

my $usage = "$0 \n <indir> <outdir> <wig_file> <col_label> \n\n";
die $usage unless(@ARGV == 4);

my $indir  = shift or die;
my $outdir = shift or die;

my $wig_file  = shift or die;

my $label  = shift or die;

die unless (-e $wig_file);

die unless (-d $indir);
die unless (-d $outdir);

opendir(DIR, $indir);
my @files = grep /\.txt$/, readdir DIR;
closedir DIR;

foreach my $infile(@files){
	if( $infile =~ /(\S+)\.txt$/){
		my $outfile = $1 . "_$label" . ".txt";
		
		my $input = File::Spec->catfile($indir,$infile);
		my $output = File::Spec->catfile($outdir, $outfile);
		
		die "wrong input" unless (-e $input);
		die "wrong output" if ( -e $output);

		my $cmd = "perl $script $label $wig_file $input $output";
		print STDERR $cmd, "\n\n";
		if(!$debug){
			`$cmd`;
		}
	}
}


exit;