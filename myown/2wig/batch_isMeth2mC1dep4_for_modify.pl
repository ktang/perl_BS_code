#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;

my $script = "/Users/tang58/Kai_BS/myown/2wig/isMeth2wig_v0.5.pl";
die unless (-e $script);


print STDERR "\n STOP STOP!!!\n";
print STDERR "Please check WT_isMeth_file and pre !!!\n\n";


print STDERR "script used:\n";
print STDERR $script, "\n\n\n";



my $debug = 1;

#my $p_value = 0.01;

if($debug){
	print STDERR "debug = 1\n\n";
}



#my @pres = ( "JKZ134_ros1-4",  "ibm1-4A", "nrpd1-3A", "nrpe1-11A", "rdd_kaust_pool");

my %pres = ( 
			 "idm1idm2A" 		=> "/Volumes/My_Book/20120813_HSZ12204_ARAyagH/WQ_first/isMeth/idm1idm2A_isMeth_chrC_error_separately_called.txt",
			 "idm1nrpd1a-3A" 	=> "/Volumes/My_Book/20120813_HSZ12204_ARAyagH/WQ_first/isMeth/idm1nrpdla-3A_isMeth_chrC_error_separately_called.txt",
			 "idm2_JKZ12_hsp" 	=> "/Volumes/Macintosh_HD_2/idm1idm2_double_mutant_Sep1/re_analysis_of_Kaust_data/isMeth/JKZ12_hsp_s6_isMeth_chrC_error_separately_called.txt",
			 "idm2rdm3A" 		=> "/Volumes/My_Book/20120813_HSZ12204_ARAyagH/WQ_first/isMeth/idm2rdm3A_isMeth_chrC_error_separately_called.txt",
			 "idm2ros1-4A" 		=> "/Volumes/My_Book/20120813_HSZ12204_ARAyagH/WQ_first/isMeth/idm2ros1-4A_isMeth_chrC_error_separately_called.txt"
			);


#my $usage = "$0 \n <isMeth_file_dir> <bed_indir> <outdir> \n\n";
#die $usage unless(@ARGV == 3);

my $usage = "$0 \n  <outdir> \n\n";
die $usage unless(@ARGV == 1);


# <input_isMeth_name> <outdir> <out_pre> <mC_cutoff> <depth_cutoff>


my $outdir = shift or die;
die unless (-d $outdir);

foreach my $pre (sort keys %pres){
		
		my  $isMeth_file = $pres{$pre} ;
		die $isMeth_file unless (-e $isMeth_file);
		
		my $cmd = "time perl $script  $isMeth_file  $outdir $pre 1 4";
		print STDERR $cmd, "\n\n";
		if(!$debug){
			`$cmd`;
		}
	
}

exit;