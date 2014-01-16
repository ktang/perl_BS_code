#!/usr/bin/perl -w

use strict;

my $debug = 0;
print STDERR "debug:$debug\n\n";

my $usage = "$0 do";
die $usage unless (@ARGV == 1 and $ARGV[0] eq "do");

my $script = "/Users/tang58/Kai_BS/myown/brat_bw_one_step/brat_bw_one_step_one_dir_nodupl_Col0.pl";

#"/Volumes/Macintosh_HD_2/BS_data_HiSeq/C24_background/JKZ19_C24_Luc_149" => ["JKZ19_C24_Luc_149_s2", 149],
my %dirs = ( 
				"/Volumes/My_Book/20120813_HSZ12204_ARAyagH/col0/co1A/clean_data" => "Col0_WT_CF",
				"/Volumes/My_Book/20120813_HSZ12204_ARAyagH/col0/idm3A/clean_data" => "idm3A",
				"/Volumes/My_Book/20120813_HSZ12204_ARAyagH/col0/rdm16-2A/clean_data" => "rdm16-2A",
				"/Volumes/My_Book/20120813_HSZ12204_ARAyagH/col0/salk-070025A/clean_data" => "salk-070025A",
				"/Volumes/My_Book/20120813_HSZ12204_ARAyagH/col0/sk460A/clean_data" => "sk460A",
				"/Volumes/My_Book/20120813_HSZ12204_ARAyagH/col0/stalA/clean_data" => "sta1A",
				
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