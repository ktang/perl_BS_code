#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;

my $debug = 0;

my $debug_db = 0;

my %isMeth_db = (
					"JKZ131_Col0" 	=> "/Volumes/Macintosh_HD_2/idm1idm2_double_mutant_Sep1/re_map_HiSeq/isMeth/JKZ131_Col0_s8_isMeth_chrC_error_separately_called.txt",
					#"idm1_pool" 	=> "/Volumes/Macintosh_HD_2/idm1idm2_double_mutant_Sep1/Sep12/1_pool_Kaust_data/isMeth/idm1_kaust_pool_isMeth_chrC_error_separately_called.txt",
					"idm2_pool" 	=> "/Volumes/Macintosh_HD_2/idm1idm2_double_mutant_Sep1/Sep12/1_pool_Kaust_data/isMeth/idm2_kaust_pool_isMeth_chrC_error_separately_called.txt",
					#"ros1_4_pool" 	=> "/Volumes/Macintosh_HD_2/idm1idm2_double_mutant_Sep1/Sep12/1_pool_Kaust_data/isMeth/ros1_4_kaust_low_pool_isMeth_chrC_error_separately_called.txt",
					#"nrpd1_3B" 		=> "/Volumes/My_Book/20120702_SH_extra_BS_Seq/nodupl/isMeth_nodupl/nrpd1-3B_nodupl_isMeth_chrC_error_separately_called.txt",
					#"idm1nrpd1_3" 	=> "/Volumes/My_Book/20120813_HSZ12204_ARAyagH/WQ_first/isMeth/idm1nrpdla-3A_isMeth_chrC_error_separately_called.txt",
					#"idm2rdm3" 		=> "/Volumes/My_Book/20120813_HSZ12204_ARAyagH/WQ_first/isMeth/idm2rdm3A_isMeth_chrC_error_separately_called.txt",
					#"idm1idm2"		=> "/Volumes/My_Book/20120813_HSZ12204_ARAyagH/WQ_first/isMeth/idm1idm2A_isMeth_chrC_error_separately_called.txt",
					#"idm2ros1_4"	=> "/Volumes/My_Book/20120813_HSZ12204_ARAyagH/WQ_first/isMeth/idm2ros1-4A_isMeth_chrC_error_separately_called.txt"
					"rdd_pool"	=> "/Volumes/Macintosh_HD_2/idm1idm2_double_mutant_Sep1/Sep12/1_pool_Kaust_data/isMeth/rdd_kaust_pool_isMeth_chrC_error_separately_called.txt"	
				);
				
foreach my $db_pre(sort keys %isMeth_db){
	print STDERR $db_pre, "\t";
	die unless (-e $isMeth_db{$db_pre});
}

print STDERR "\n\n";

if($debug_db){
	%isMeth_db = ( "de1" => "/Users/tang58/Kai_BS/myown/for_boxplot/debug/isMeth_db/1_isMeth.txt",
				   "d2" => "/Users/tang58/Kai_BS/myown/for_boxplot/debug/isMeth_db/2_isMeth.txt"
				 );
				 
	foreach my $db_pre(sort keys %isMeth_db){
		print STDERR $db_pre, "\t";
		die unless (-e $isMeth_db{$db_pre});
	}

	print STDERR "\n\n";
}

my $script = "/Users/tang58/Kai_BS/myown/for_boxplot/count_meth_level_in_Meth_file_for_bed_list_single_isMeth.pl";
die unless (-e $script);

#my $dep_cutoff = 4;

print STDERR "\n\n script used is\n";
print STDERR $script, "\n\n";

my $usage = "$0 \n <dep_cutoff> <indir> <file_name>  <outdir> <output_name> <label_num> [ <lable> ...]\n\n";
die $usage unless(@ARGV >= 7);


my $dep_cutoff = shift or die;
my $indir = shift or die;
my $orignal_file_name = shift or die;
my $outdir = shift or die;
my $output_name = shift or die;

my $label_num = shift or die;
my $last_index = $label_num - 1;

die unless (@ARGV == $label_num);

my @labels;
#my @files;
#my @outputs;

for my $i(0..$last_index){
#	$files[$i]  = shift or die;
	$labels[$i] = shift or die;
	die $labels[$i] unless (defined $isMeth_db{$labels[$i]});
#	$outputs[$i] = File::Spec->catfile($outdir,  $labels[$i] . "_isMeth_depth" . $dep_cutoff . ".txt" );
#	die $files[$i] unless (-e $files[$i]);
#	die if (-e $outputs[$i]);
	
}




my $orignal_file = File::Spec->catfile($indir, $orignal_file_name);
die unless (-e $orignal_file);



die unless (-d $indir);
die unless (-d $outdir);

my $output = File::Spec->catfile($outdir, $output_name);
die if(-e $output);

my $i = -1;

my $f0 = $orignal_file . ".temp0";

my $cp_cmd = "cp $orignal_file $f0";
print STDERR $cp_cmd, "\n\n";
if(!$debug){
	`$cp_cmd`;
}

foreach my $label (@labels){
	$i++;
	
	my $isMeth_file = $isMeth_db{$label};
	my $in_file = $orignal_file . ".temp" . $i;
	if(!$debug){
		die unless (-e $in_file);
	}
	
	my $out_file =  $orignal_file . ".temp" . ($i+1);
	die if(-e $out_file);
	
	#my $cmd = "perl $script $bench_file $overlap_file $label $out_file";

	
	# <isMeht_file> <label> <bed_like_file> <output> <depth>
	my $cmd = "perl $script  $isMeth_file   $label  $in_file $out_file  $dep_cutoff ";
	
	print STDERR $cmd, "\n\n";
	if(!$debug){
		`$cmd`;
	}
}

my $final = $orignal_file . ".temp" . ($i+1);

if(!$debug){
	die unless (-e $final);
}

my $cp_out_cmd = "cp $final $output";
print STDERR $cp_out_cmd, "\n\n";
if(!$debug){
	`$cp_out_cmd`;
}

my $rm_cmd = "rm -f $indir/" . $orignal_file . ".temp*";
print STDERR $rm_cmd, "\n\n";
if(!$debug){
	`$rm_cmd`;
}


exit;
