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
my $script_0 = "/Users/tang58/Kai_BS/myown/Huiming_downstream_analysis/0_annotate_whether_in_promoter_v0.1.pl";
die unless (-e $script_0);
#my $usage = "$0 \n <1kb_1.5kb_2kb> <input_hyper_hypo_list> STDOUT \n\n";

my $script_1 = "/Users/tang58/Kai_BS/myown/Huiming_downstream_analysis/1_TE_near_or_not_gene_no_output_no_TE_limit100_v0.1.pl";
die unless (-e $script_1);

my $usage = "$0 \n <indir> <output> \n\n";
die $usage unless(@ARGV == 2);

my $indir  = shift or die;
#my $outdir = shift or die;
my $output = shift or die;

die unless (-d $indir);
die if(-e $output);

opendir(DIR, $indir);
my @files = grep /\.txt$/, readdir DIR;
closedir DIR;

open(OUT, ">>$output" ) or die;

my $gap = "1kb";

foreach my $infile(@files){
	if( $infile =~ /(\S+)\.txt$/){
		
		my $input = File::Spec->catfile($indir,$infile);
		
		die "wrong input" unless (-e $input);
		
		print OUT join("\n", ($infile, $gap)), "\n";
		
		my $cmd_0 = "$script_0 $gap $input >> $output";
		print STDERR $cmd_0, "\n\n";
		if(!$debug){
			`$cmd_0`;
		}
		
		print OUT "\n\n";
		
		my $cmd_1 = "$script_1 $input >> $output";
		print STDERR $cmd_1, "\n\n";
		if(!$debug){
			`$cmd_1`;
		}
		
		print OUT "\n\n";
	}
}

close OUT;

exit;