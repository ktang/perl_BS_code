#!/usr/bin/perl -w

#v0.2
# file name is not $file =~ /(\S+)_vs_colA\S*_(hyper|hypo)/
# not (\S+)\.txt
#but (\S+)_vs
use strict;
use File::Spec;

my $debug = 0 ;

my $script = "/Users/tang58/Kai_BS/myown/batch_add_overlap_annotation/add_overlap_annotation_multiple_file_v0.2.pl";
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
