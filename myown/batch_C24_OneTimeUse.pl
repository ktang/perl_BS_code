#!/usr/bin/perl -w

use strict;

my $debug = 0;
print STDERR "debug:$debug\n\n";

my $script = "/Users/tang58/Kai_BS/myown/Kai_batch_after_brat_before_cal_methy.pl";

#"/Volumes/Macintosh_HD_2/BS_data_HiSeq/C24_background/JKZ19_C24_Luc_149" => ["JKZ19_C24_Luc_149_s2", 149],
my %dirs = ( "/Volumes/Macintosh_HD_2/BS_data_HiSeq/C24_background/JKZ20_C24_Luc_150" => ["JKZ20_C24_Luc_150_s3", 150],
  		     "/Volumes/Macintosh_HD_2/BS_data_HiSeq/C24_background/JKZ21_ros1_1_151" => ["JKZ21_ros1_1_s5", 151],
  		     "/Volumes/Macintosh_HD_2/BS_data_HiSeq/C24_background/JKZ22_ros1_1_152" => ["JKZ22_ros1_1_152_s6", 152],
  		     "/Volumes/Macintosh_HD_2/BS_data_HiSeq/C24_background/JKZ23_ros3_1_153" => ["JKZ23_ros3_1_153_s7", 153],
  		     "/Volumes/Macintosh_HD_2/BS_data_HiSeq/C24_background/JKZ24_ros3_1_154" => ["JKZ24_ros3_1_154_s8", 154],
			);
			
foreach my $dir (sort keys %dirs){
	die "no $dir" unless (-d $dir);
}			
foreach my $dir (sort keys %dirs){
	my ($pre, $id ) = @{$dirs{$dir}};
	chdir($dir);
	print STDERR `pwd` , "\n";
	
	my $cmd = "date; time perl $script . $id $pre C24";
	print STDERR $cmd,"\n\n";
	
	if(!$debug){system("$cmd")}
	
}