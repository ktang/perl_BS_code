#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;

my $debug = 0;

#my $script = "/Users/tang58/Kai_BS/for_publish/batch_cal_meth_level/cal_meth_level_for_bed_list_EckerPaper_multiple_isMeth.pl";
my $script =  "/Users/tang58/Kai_BS/for_publish/batch_cal_meth_level_original_mC_NotEckerSuggest_v0.1_bed_may_overlap/cal_meth_level_for_bed_list_original_mC_NotEckerSuggest_multiple_isMeth_v0.1.pl";

die unless (-e $script);

print STDERR "\n\n script used is\n";
print STDERR $script, "\n\n";

my $usage = "$0 \n <dir_isMeth> <indir> <outdir>\n\n";
die $usage unless(@ARGV == 3);

my $dir_isMeth = shift or die;
my $indir = shift or die;
my $outdir = shift or die;
die unless (-d $indir);
die unless (-d $outdir);
die unless (-d $dir_isMeth);

opendir(DIR, $indir);
my @files = grep /\.txt$/, readdir DIR;
closedir DIR;

foreach my $file( @files ){
	my $out_name;
	
	if($file =~ /(\S+)\.txt$/){
		$out_name = $1 . "_with_meth_original_mC_NotE.txt"; #notE v0.1 may be overlap
	}else{
		die $file;
	}
	
	my $input_bed = File::Spec->catfile($indir, $file);
	#<isMeth_dir> <input_bed> <outdir> <output_name>
#	my $cmd = "perl $script $file $indir $outdir $out_name";
	my $cmd = "perl $script $dir_isMeth $input_bed $outdir $out_name";
	print STDERR $cmd, "\n\n";
	if(!$debug){
		`$cmd`;
	}
}

exit;
