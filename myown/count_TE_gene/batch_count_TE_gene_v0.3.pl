#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;

my $debug = 0;

print STDERR "\n\n nodupl tag\n\n";


my $script = "/Users/tang58/Kai_BS/myown/count_TE_gene/count_num_in_annotation_gene_TE_intergenic_v0.3.pl";
die unless (-e $script);

print STDERR "\nscript used \n";
print STDERR $script, "\n\n";


my $usage = "$0 \n<indir> <output>\n\n";
die $usage unless(@ARGV == 2);

my $indir = shift or die;
my $output = shift or die;
die unless (-d $indir);

die if(-e $output);

if(!$debug){
	open(OUT, ">$output") or die;
	print OUT join("\t", ( "total", "TE", "gene","intergenic", "pseudogene", "other", "TE+protein","TE+inter", "all_intergenic_number", "sample" )), "\n";
	close OUT;
}

opendir(DIR, $indir);
my @files = grep /\.txt$/, readdir DIR;
closedir DIR;


#my $dir = "/Volumes/My_Book/20120427_ShangHai_data/handled_data/cytosine_coverage_info";

foreach my $file(@files){
#	my $cmd = "time perl $script $indir $outdir $pre $file";
	my $label = "";
	if($file =~ /(\S+)\.txt$/){
		$label = $1;
	}
	my $input = File::Spec->catfile($indir, $file);
	die unless (-e $input);
	my $cmd = "perl $script $input $label no >> $output";
	
	if(!$debug){
		`$cmd`;
	}else{
		print STDERR $cmd, "\n\n";
	}
}

exit;