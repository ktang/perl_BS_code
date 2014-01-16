#!/usr/bin/perl -w

use strict;

my $debug = 0;
print STDERR "debug:$debug\n\n";

my $usage = "$0 <do>";
die $usage unless (@ARGV == 1 and $ARGV[0] eq "do");

# my $script = "/Users/tang58/Kai_BS/myown/brat_bw_one_step/brat_bw_one_step_one_dir_nodupl_Col0.pl";
my $script = "/Volumes/Macintosh_HD_2/idm1idm2_double_mutant_Sep1/re_analysis_of_Kaust_data/brat_bw_one_step_one_dir_nodupl_Col0_Kaust.pl";


my $indir = "/Volumes/My_Book/Ath_BS_data/Run19_Illumina_GAIIx";

#<indir> <tar_file> <s1_file> <s2_file> <pre>

#"/Volumes/Macintosh_HD_2/BS_data_HiSeq/C24_background/JKZ19_C24_Luc_149" => ["JKZ19_C24_Luc_149_s2", 149],
my %files = (
			"Lane5_JKZ11_hsp.tar.gz" => [ "s_5_1_sequence.txt", "s_5_2_sequence.txt", "JKZ11_hsp_s5"],
			"Lane6_JKZ12_hsp.tar.gz" => [ "s_6_1_sequence.txt", "s_6_2_sequence.txt", "JKZ12_hsp_s6"],
			"Lane7_JKZ13_phd.tar.gz" => [ "s_7_1_sequence.txt", "s_7_2_sequence.txt", "JKZ13_phd_s7"],
			"Lane8_JKZ14_phd.tar.gz" => [ "s_8_1_sequence.txt", "s_8_2_sequence.txt", "JKZ14_phd_s8"]
);
			
			
foreach my $file (sort keys %files){
	my ($s1, $s2, $pre )= @{$files{$file}};
	
	my $cmd = "time perl $script $indir $file $s1 $s2 $pre";
	
	print STDERR $cmd,"\n\n";
	
	if(!$debug){
		`$cmd`;
	}
	
}