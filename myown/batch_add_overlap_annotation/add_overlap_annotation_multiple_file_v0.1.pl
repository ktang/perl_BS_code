#!/usr/bin/perl -w

#v0.1
# file name is not $file =~ /(\S+)_vs_colA\S*_(hyper|hypo)/
# but (\S+)\.txt

use strict;
use File::Spec;

my $debug = 0 ;

my $script = "/Users/tang58/Kai_BS/myown/batch_add_overlap_annotation/add_overlap_annotation_SH.pl";
die unless (-e $script);

print STDERR "\n\n script used is\n";
print STDERR $script, "\n\n";

my $usage = "$0 \n <file_name> <indir> <outdir> <output_name>\n\n";
die $usage unless(@ARGV == 4);

my $orignal_file_name = shift or die;
my $indir = shift or die;
my $outdir = shift or die;
my $output_name = shift or die;


my $orignal_file = File::Spec->catfile($indir, $orignal_file_name);
die unless (-e $orignal_file);



die unless (-d $indir);
die unless (-d $outdir);

my $output = File::Spec->catfile($outdir, $output_name);
die if(-e $output);

opendir(DIR, $indir);
my @files = grep { $_ ne $orignal_file_name} (grep /\.txt$/, readdir DIR) ;
closedir DIR;

if($debug ){
	print STDERR join("\n", @files), "\n\n\n";
}

my $i = -1;

my $f0 = $orignal_file . ".temp0";

my $cp_cmd = "cp $orignal_file $f0";
print STDERR $cp_cmd, "\n\n";
if(!$debug){
	`$cp_cmd`;
}

foreach my $file (@files){
	$i++;
	
	my $bench_file = $orignal_file . ".temp" . $i;
	if(!$debug){
		die unless (-e $bench_file);
	}
	my $overlap_file =  File::Spec->catfile($indir, $file );
	die unless (-e $overlap_file);
	my $out_file =  $orignal_file . ".temp" . ($i+1);
	
	die if(-e $out_file);
	
	
	my $label ;
	#if($file =~ /(\S+)_vs_colA\S*_(hyper|hypo)/){
	if($file =~ /(\S+)\.txt$/){
		#$label = $1 . "_" . $2;
		$label = $1 ;
		
	}else{
		die $file;
	}
		
	my $cmd = "perl $script $bench_file $overlap_file $label $out_file";
	print STDERR $cmd, "\n\n";
	if(!$debug){
		`$cmd`;
	}
}

my $final = $orignal_file . ".temp" . ($i+1);

my $cp_out_cmd = "cp $final $output";
print STDERR $cp_out_cmd, "\n\n";
if(!$debug){
	`$cp_out_cmd`;
}

my $rm_cmd = "rm -f $indir/*txt.temp* ";
print STDERR $rm_cmd, "\n\n";
if(!$debug){
	`$rm_cmd`;
}


exit;
