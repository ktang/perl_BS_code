#!/usr/bin/perl -w

use strict;

my $debug = 0;
print STDERR "debug: $debug\n\n";

my $brat = "/Users/tang58/Software/BRAT/brat-1.1.18/brat";

my $cal_coverage_script = "/Users/tang58/Kai_BS/calculate_coverage_Kai.pl";

my $count_script = "/Users/tang58/Software/BRAT/brat-1.1.19/acgt-count";

my $cal_error_script = "/Users/tang58/Kai_BS/conversion_error.pl";

my $create_meth_db_input_script = "/Users/tang58/Kai_BS/create_meth_db_input_Kai.pl";

my $usage ="$0 <indir> <sample_id> <pre> <C24|Col0>";

die $usage unless (@ARGV == 4);

my ($indir, $sample_id, $pre, $flag) = @ARGV[0..3];

if($indir =~ /\/$/) {chop $indir}
#if($outdir =~ /\/$/) {chop $outdir}
my $outdir = $indir;

my $ref = "NONE";
my $pos_file = "";

if($flag eq "C24"){
	$ref = "/Users/tang58/DataBase/BRAT_reference_files/C24/C24_references_names.txt";
	$pos_file = "/Volumes/Macintosh_HD_2/BS_data_HiSeq/C24_background/ath_C24/ath_C24_position.txt";
}elsif($flag eq "Col0"){
	$ref = "/Users/tang58/DataBase/BRAT_reference_files/Col0/Col0_references_names.txt";
	$pos_file = "/Users/tang58/deep_seq_analysis/Bisulfite/work_with_Liu/data_UCR_server/ath_col0_position.txt";
}else{
	die $usage, "\nC24 or Col0, exactly";
}



die "wrong indir" unless (-d $indir);

die "wrong outdir" unless (-d $outdir);

#my $in_P1 = $pre."_1_pairs_brat_input.txt";
#my $in_P2 = $pre."_2_pairs_brat_input.txt";
#my $in_S1 = $pre."_1_single_brat_input.txt";
#my $in_S2 = $pre."_2_single_brat_input.txt";

#print STDERR join("\n", ($in_P1, $in_P2, $in_S1, $in_S2)), "\n\n";

#die "wrong input" unless (-e "$indir/$in_P1" and -e "$indir/$in_P2" and -e "$indir/$in_S1" and -e "$indir/$in_S2");

my $brat_out_p = $pre."_paired_brat_out.txt";
my $brat_out_S1 = $pre . "_1_brat_out.txt";
my $brat_out_S2 = $pre. "_2_brat_out.txt";

#print STDERR "output:\n";
#print STDERR join("\n", ($out_p, $out_S1, $out_S2)), "\n";
if(!$debug){
	die "brat output does not exists" unless ((-e "$outdir/$brat_out_p") and (-e "$outdir/$brat_out_S1") and (-e "$outdir/$brat_out_S2") );
}
my $log = $pre."_batch_log.txt";

##############################################
# calculate coverage
##############################################

my $covearge_file = $pre."_coverage_stats.txt";

die "$covearge_file exists" unless (! (-e "$outdir/$covearge_file"));

my $cmd_cal_coverage = "date; time perl $cal_coverage_script $outdir/$brat_out_p $outdir/$brat_out_S1 $outdir/$brat_out_S2 >> $outdir/$covearge_file";

print STDERR $cmd_cal_coverage, "\n\n";
if (!$debug) {system("$cmd_cal_coverage")}

##############################################
# count
##############################################
my ($single_list, $pair_list) = ($pre . "_single_res.txt", $pre . "_paired_res.txt");
die "list exists!" unless ( !(-e "$outdir/$single_list") and !(-e "$outdir/$pair_list"));
open(SOUT, ">$outdir/$single_list") or die "Can't open $single_list: $!";
open(POUT, ">$outdir/$pair_list") or die "Can't open $pair_list: $!";

print POUT  "$outdir/$brat_out_p", "\n";
print SOUT  "$outdir/$brat_out_S1", "\n";
print SOUT  "$outdir/$brat_out_S2", "\n";

close SOUT;
close POUT;

my $cmd_count = "date; time $count_script -r $ref -P $pre -p $outdir/$pair_list -s $outdir/$single_list -B";
print STDERR $cmd_count ,"\n\n";
if (!$debug) {system("$cmd_count")}

##############################################
# cal error
##############################################

my $forw_file = $pre."_forw.txt";
my $rev_file = $pre."_rev.txt";
if (!$debug) {die "rev/forw file does not exists" unless (-e "$outdir/$rev_file" and -e "$outdir/$forw_file")}

my $error_file = $pre."_conversion_error.txt";
die "error_file exists" unless (!(-e "$outdir/$error_file"));

my $cmd_cal_error_forw = "date; time perl $cal_error_script $forw_file >> $outdir/$error_file";
my $cmd_cal_error_rev  = "date; time perl $cal_error_script $rev_file  >> $outdir/$error_file";

print STDERR $cmd_cal_error_forw , "\n\n";
print STDERR $cmd_cal_error_rev , "\n\n";
if (!$debug) {
	system("$cmd_cal_error_forw");
	system("$cmd_cal_error_rev");
}


##############################################
# create_meth_db_input
##############################################
#date; time perl ~/Kai_BS/create_meth_db_input_Kai.pl /Users/tang58/deep_seq_analysis/Bisulfite/work_with_Liu/data_UCR_server/ath_col0_position.txt
# /Volumes/Macintosh_HD_2/Kai_BS/brat_output /Volumes/Macintosh_HD_2/Kai_BS/brat_output col0_s_8_JKZ131 131
# <pos_file> <indir> <outdir> <out_prefix> [<sample_id>] 
my $cmd_create_meth_db_input = "date; time perl $create_meth_db_input_script $pos_file $indir $outdir $pre $sample_id";

print STDERR $cmd_create_meth_db_input , "\n\n";

if (!$debug) {
	system("$cmd_create_meth_db_input");
}


exit;