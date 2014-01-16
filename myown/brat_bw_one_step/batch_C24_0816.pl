#!/usr/bin/perl -w

use strict;

my $debug = 0;
print STDERR "debug:$debug\n\n";

my $usage = "$0 do";
die $usage unless (@ARGV == 1 and $ARGV[0] eq "do");

my $script = "/Users/tang58/Kai_BS/myown/brat_bw_one_step/brat_bw_one_step_one_dir_nodupl_C24.pl";

#"/Volumes/Macintosh_HD_2/BS_data_HiSeq/C24_background/JKZ19_C24_Luc_149" => ["JKZ19_C24_Luc_149_s2", 149],
# "/Volumes/My_Book/20120813_HSZ12204_ARAyagH/WQ_first/ape1L-1A/clean_data" => "ape1L-1A",
my %dirs = (
			 "/Volumes/My_Book/20120813_HSZ12204_ARAyagH/C24/c24A/clean_data" 		 => "C24_WT_CF",
			 "/Volumes/My_Book/20120813_HSZ12204_ARAyagH/C24/rdm16ros1-1/clean_data" => "rdm16ros1-1",
			 "/Volumes/My_Book/20120813_HSZ12204_ARAyagH/C24/ros1-1A/clean_data" 	 => "ros1-1_CF"
		   );


foreach my $dir (sort keys %dirs){
	die "no $dir" unless (-d $dir);
}
			
foreach my $dir (sort keys %dirs){
	my $pre = $dirs{$dir};
	
	my $cmd = "time perl $script $dir $pre ";
	
	print STDERR $cmd,"\n\n";
	
	if(!$debug){
		`$cmd`;
	}
	
}