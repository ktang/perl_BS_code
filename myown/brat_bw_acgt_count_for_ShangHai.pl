#!/usr/bin/perl -w

#/Users/tang58/Kai_BS/myown
use strict;
use File::Spec;

my $debug = 1;

#my $acgt_count = "/Users/tang58/Software/BRAT/brat-1.1.19/acgt-count";

my $acgt_count = "/Users/tang58/Software/BRAT/brat_bw-2.0.1/acgt-count";
die unless (-e $acgt_count);

my $ref = "/Users/tang58/DataBase/BRAT_reference_files/Col0/Col0_references_names.txt";
die unless (-e $ref);

my $usage ="$0 <dir> <pre>";

die $usage unless (@ARGV == 2);

#my ($indir, $outdir, $pre, $flag) = @ARGV[0..3];

my $dir = shift or die "dir";
my $pre = shift or die "pre";

die "wrong dir" unless (-d $dir);

my $in_P    = File::Spec->catfile($dir, $pre ."_paired_brat_bw_out.txt");

my $in_S1   = File::Spec->catfile($dir, $pre . "_1_brat_bw_out.txt"   );
my $in_S2   = File::Spec->catfile($dir, $pre . "_2_brat_bw_out.txt"    );
my $in_left = File::Spec->catfile($dir, $pre . "_paired_brat_bw_out.txt.single_mates");

print STDERR "brat_file:\n";
print STDERR join("\n", ($in_P, $in_S1, $in_S2, $in_left)), "\n\n";

die "wrong input" unless (-e $in_P and -e $in_S1 and -e $in_S2 and -e $in_left);

#my ($single_list, $pair_list) = ($pre."_single_list.txt", $pre."_pair_list.txt");
my $pair_list   = File::Spec->catfile($dir, $pre . "_pair_list.txt");
my $single_list = File::Spec->catfile($dir, $pre . "_single_list.txt");

my $log_file = File::Spec->catfile($dir, $pre . "acgt_count_log.txt");

die "wrong list file" if (-e $pair_list or -e $single_list);

open(POUT, ">$pair_list") or die "Can't open $pair_list: $!";
print POUT $in_P, "\n";
close (POUT);


open(SOUT, ">$single_list") or die "Can't open $single_list: $!";

print SOUT $in_S1,"\n";
print SOUT $in_S2,"\n";
print SOUT $in_left, "\n";
close SOUT;

my $cmd = "time $acgt_count -r $ref -P $pre -p $pair_list -s $single_list -B >> $log_file 2>&1";

print STDERR $cmd, "\n\n";
if(!$debug){
	`$cmd`;
}
exit;