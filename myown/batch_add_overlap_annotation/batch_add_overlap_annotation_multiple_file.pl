#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;

my $debug = 0 ;

my $script = "/Users/tang58/Kai_BS/myown/batch_add_overlap_annotation/add_overlap_annotation_multiple_file.pl";
die unless (-e $script);

print STDERR "\n\n script used is\n";
print STDERR $script, "\n\n";

my $usage = "$0 \n<indir> <outdir>\n\n";
die $usage unless(@ARGV == 2);

my $indir = shift or die;
my $outdir = shift or die;
die unless (-d $indir);
die unless (-d $outdir);

opendir(DIR, $indir);
my @files = grep /\.txt$/, readdir DIR;
closedir DIR;

#<file_name> <indir> <outdir> <output_name>
foreach my $file( @files ){
	my $out_name;
	
	if($file =~ /(\S+)\.txt$/){
		$out_name = $1 . "_all_overlap.txt";
	}else{
		die $file;
	}
	
	my $cmd = "perl $script $file $indir $outdir $out_name";
	print STDERR $cmd, "\n\n";
	if(!$debug){
		`$cmd`;
	}
}

exit;
