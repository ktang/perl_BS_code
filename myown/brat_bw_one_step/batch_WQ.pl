#!/usr/bin/perl -w

use strict;

my $debug = 0;
print STDERR "debug:$debug\n\n";

my $script = "/Users/tang58/Kai_BS/myown/brat_bw_one_step/brat_bw_one_step_one_dir_nodupl_Col0.pl";

#"/Volumes/Macintosh_HD_2/BS_data_HiSeq/C24_background/JKZ19_C24_Luc_149" => ["JKZ19_C24_Luc_149_s2", 149],
my %dirs = ( "/Volumes/My_Book/20120813_HSZ12204_ARAyagH/WQ_first/ape1L-1A/clean_data" => "ape1L-1A",
			 "/Volumes/My_Book/20120813_HSZ12204_ARAyagH/WQ_first/ape2-2A/clean_data" => "ape2-2A",
			 "/Volumes/My_Book/20120813_HSZ12204_ARAyagH/WQ_first/idm1nrpdla-3A/clean_data" => "idm1nrpdla-3A",
			 "/Volumes/My_Book/20120813_HSZ12204_ARAyagH/WQ_first/idm2rdm3A//clean_data" => "idm2rdm3A",
			 "/Volumes/My_Book/20120813_HSZ12204_ARAyagH/WQ_first/idm2ros1-4A/clean_data" => "idm2ros1-4A",
			 "/Volumes/My_Book/20120813_HSZ12204_ARAyagH/WQ_first/idmlidm2A/clean_data" => "idm1idm2A"
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