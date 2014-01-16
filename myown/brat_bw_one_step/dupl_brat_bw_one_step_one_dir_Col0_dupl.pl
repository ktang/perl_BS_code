#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;

my $debug = 1;

if($debug){
	print STDERR "debug = 1\n\n";
}

my $brat_bw = "/Users/tang58/Software/BRAT/brat_bw-2.0.1/brat_bw";
die unless (-e $brat_bw);

my $trim_tool = "/Users/tang58/Software/BRAT/brat_bw-2.0.1/trim";
die unless (-e $trim_tool);

my $ref = "/Users/tang58/DataBase/TAIR_Col0_genome/index/BRAT_BW/brat_bw_arabdopsis_index";
die unless (-d $ref);

my $fas_ref = "/Users/tang58/DataBase/BRAT_reference_files/Col0/Col0_references_names.txt";
die unless (-e $fas_ref);

my $acgt_count = "/Users/tang58/Software/BRAT/brat_bw-2.0.1/acgt-count";
die unless (-e $acgt_count);

my $usage = "$0 <indir> <pre>";
die $usage unless(@ARGV == 2);

my $indir = shift or die "indir";
my $pre = shift or die "pre";
#my $outdir = shift or die "outdir";


opendir(DIR, $indir) or die "dir";

my @files = grep /fq\.gz$/, readdir DIR;

die unless (@files == 2);

die "1: $files[0]" unless ($files[0] =~ /_1\.fq\.gz$/);
die "2: $files[1]" unless ($files[1] =~ /_2\.fq\.gz$/);

my $gz1    = File::Spec->catfile($indir, $files[0]);
my $gz2    = File::Spec->catfile($indir, $files[1]);

my $fq_in1 = File::Spec->catfile($indir, $pre . "_1.fq");
my $fq_in2 = File::Spec->catfile($indir, $pre . "_2.fq");



if(!$debug){
	die unless (-e $gz1);
	die unless (-e $gz2);
	
	die if (-e $fq_in1);
	die if (-e $fq_in2);
}

######################
# 0 unzip
##########################

my $cmd_uzip1 = "time zcat $gz1 > $fq_in1";
print STDERR $cmd_uzip1, "\n\n";
if(!$debug){
	`$cmd_uzip1`;
}

if(!$debug){
	die unless (-e $gz1);
	die unless (-e $gz2);
	
	die unless (-e $fq_in1);
	die if (-e $fq_in2);
}

my $cmd_uzip2 = "time zcat $gz2 > $fq_in2";

print STDERR $cmd_uzip2, "\n\n";
if(!$debug){
	`$cmd_uzip2`;
}

if(!$debug){
	die unless (-e $gz1);
	die unless (-e $gz2);
	
	die unless (-e $fq_in1);
	die unless (-e $fq_in2);
}


#####################
# 1 brat_trim
#####################
my $q = 20;
my $L = 64;
my $m_trim = 0;

my $trim_cmd = "time $trim_tool -1 $fq_in1 -2 $fq_in2 -P $indir/$pre -q $q -L $L -m $m_trim";
print STDERR $trim_cmd,"\n\n";

if(!$debug){
	`$trim_cmd`;
}

#####################
# 2 brat_bw
#####################
my $i = 0;
my $a = 1000;
my $m_brat_bw = 2;

#PE reads 1,2
my $in_P1 = File::Spec->catfile($indir, $pre."_reads1.txt");
my $in_P2 = File::Spec->catfile($indir, $pre."_reads2.txt");

#SE mates 1,2
my $in_S1 = File::Spec->catfile($indir, $pre."_mates1.txt");
my $in_S2 = File::Spec->catfile($indir, $pre."_mates2.txt");


my $outdir = $indir;

my $out_p  = File::Spec->catfile($outdir, $pre."_paired_brat_bw_out.txt");
my $out_S1 = File::Spec->catfile($outdir, $pre . "_1_brat_bw_out.txt");
my $out_S2 = File::Spec->catfile($outdir, $pre. "_2_brat_bw_out.txt");
my $out_left = File::Spec->catfile($outdir, $pre . "_paired_brat_bw_out.txt.single_mates");

if(!$debug){
	die "wrong input" unless (-e $in_P1 and -e $in_P2 and -e $in_S1 and -e $in_S2);
	die "output exists" if(-e $out_p or -e $out_S1 or -e $out_S2);
}

my $log = File::Spec->catfile($outdir, $pre."_brat_bw_log.txt");


my $cmd_S1 = "date; time $brat_bw -s $in_S1 -P $ref -o $out_S1 -m $m_brat_bw  >> $log 2>&1";
print STDERR $cmd_S1,"\n\n";
open(OUT, ">>$log") or die;
print OUT $cmd_S1, "\n\n";
close OUT;

if(!$debug){
	`$cmd_S1`;
}


my $cmd_S2 = "date; time $brat_bw -s $in_S2 -P $ref -o $out_S2 -m $m_brat_bw -A >> $log 2>&1";
print STDERR $cmd_S2,"\n\n";
open(OUT, ">>$log") or die;
print OUT "\n\n", $cmd_S2, "\n\n";
close OUT;
if(!$debug){
	`$cmd_S2`;
}

my $cmd_P = "date; time $brat_bw -1 $in_P1 -2 $in_P2 -pe -o $out_p -P $ref -i $i -a $a -m $m_brat_bw  >> $log 2>&1";
print STDERR $cmd_P,"\n\n";
open(OUT, ">>$log") or die;
print OUT "\n\n",  $cmd_P, "\n\n";
close OUT;
if(!$debug){
	`$cmd_P`;
}


#####################
# 3 acgt_count
#####################

if(!$debug){
	die unless (-e $out_p);
	die unless (-e $out_S1);
	die unless (-e $out_S2);
}


my $single_list = File::Spec->catfile($indir, $pre."_single_list.txt");
die if (-e $single_list);
my $pair_list   = File::Spec->catfile($indir, $pre."_pair_list.txt");
die if (-e $pair_list );

die if (-e $single_list);
open(SOUT, ">$single_list") or die "Can't open $single_list: $!";

die if (-e $pair_list );
open(POUT, ">$pair_list") or die "Can't open $pair_list: $!";

print SOUT "$out_S1","\n";
print SOUT "$out_S2","\n";
print SOUT $out_left. "\t";
close SOUT;

print POUT "$out_p", "\n";
close POUT;

my $cmd_acgt_count = "date; time $acgt_count -r $fas_ref -P $indir/$pre -p $pair_list -s $single_list -B";

print STDERR $cmd_acgt_count, "\n\n";
if(!$debug){
	`$cmd_acgt_count`;
}

print STDERR "\n\n OK \n\n";


exit;