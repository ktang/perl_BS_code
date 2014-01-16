#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;

#my $script = "/Users/tang58/Kai_BS/myown/select_DMR_Ecker_method/count_meth_level_in_Meth_file_for_bed_list.pl";
#my $script = "/Users/tang58/Kai_BS/myown/select_DMR_Ecker_method/count_meth_level_in_Meth_file_for_bed_list_v0.1.pl";
die unless (-e $script);


print STDERR "\n STOP STOP!!!\n";
print STDERR "Please check WT_isMeth_file and pre !!!\n\n";


print STDERR "script used:\n";
print STDERR $script, "\n\n\n";



my $debug = 0;

#my $p_value = 0.01;

if($debug){
	print STDERR "debug = 1\n\n";
}

my $WT_isMeth_file = "/Users/tang58/misc/Huiming_Zhang/9_PolIV_PolV_ros1_rdd_ibm1/isMeth_dep5_6samples/colA_isMeth_depth5.txt";
die unless (-e $WT_isMeth_file);

print STDERR "WT: $WT_isMeth_file\n\n";

my @pres = ( "JKZ134_ros1-4",  "ibm1-4A", "nrpd1-3A", "nrpe1-11A", "rdd_kaust_pool");

my $usage = "$0 \n <isMeth_file_dir> <bed_indir> <outdir> \n\n";
die $usage unless(@ARGV == 3);

my $isMeth_file_dir = shift or die;
die unless (-e $isMeth_file_dir);

my $indir = shift or die;
die unless (-d $indir);

my $outdir = shift or die;
die unless (-d $outdir);

#opendir(DIR, $indir) or die;
#my @files = grep /_both\.txt$/, readdir DIR;
#my @files = grep /\.txt$/, readdir DIR;
#closedir DIR;



my %posts = ("hyper" => "_dep5_hyper_P0.01_reduced_boundary_dep5_WinSize1000_gap1000_initialCutoff5_reportCutoff10.txt", 
			 "hypo"  => "_dep5_hypo_P0.01_reduced_boundary_dep5_WinSize1000_gap1000_initialCutoff5_reportCutoff10.txt");

#my $usage = "$0 \n <isMeht_file_wt> <isMeht_file_mut>  <bed_like_file> <hyper|hypo> <outdir> \n\n";

foreach my $pre (@pres){
	foreach my $label(keys %posts){
		my $mut_file = File::Spec->catfile( $isMeth_file_dir, $pre. "_isMeth_depth5.txt");
		
		if($pre eq "JKZ134_ros1-4"){
			$mut_file = File::Spec->catfile( $isMeth_file_dir, $pre. "_s1_isMeth_depth5.txt");
		}
		
		die $mut_file unless (-e $mut_file);
		
		my $bed_file = File::Spec->catfile( $indir, $pre. $posts{$label});
		
		die $bed_file unless (-e $bed_file);
		
		my $cmd = "time perl $script $WT_isMeth_file $mut_file $bed_file $label $outdir";
		print STDERR $cmd, "\n\n";
		if(!$debug){
			`$cmd`;
		}
		
		
	}
}

exit;