#!/usr/bin/perl -w

#outuse: for Becky's srver or linux;
# the ref file are in diff dir

use strict;
use File::Spec;

#my $brat = "/Users/tang58/Software/BRAT/brat-1.1.18/brat";
my $brat_bw = "/Users/tang58/Software/BRAT/brat_bw-2.0.1/brat_bw";
die unless (-e $brat_bw);
my $debug = 1;

my $usage ="$0 <indir> <outdir> <pre>";

die $usage unless (@ARGV == 3);

#my ($indir, $outdir, $pre, $flag) = @ARGV[0..3];
my $indir = shift or die "indir";
my $outdir = shift or die "outdir";
my $pre = shift or die "prefix";

#
my $ref = "/Users/tang58/DataBase/TAIR_Col0_genome/index/BRAT_BW/brat_bw_arabdopsis_index";
die unless (-d $ref);
my $i = 0;
my $a = 1000;
my $m = 2;
#

die "wrong indir" unless (-d $indir);
die "wrong outdir" unless (-d $outdir);

#my $in_P1 = $pre."_1_pairs_brat_input.txt";
#my $in_P2 = $pre."_2_pairs_brat_input.txt";
#my $in_S1 = $pre."_1_single_brat_input.txt";
#my $in_S2 = $pre."_2_single_brat_input.txt";

#PE reads 1,2
my $in_P1 = File::Spec->catfile($indir, $pre."_reads1.txt");
my $in_P2 = File::Spec->catfile($indir, $pre."_reads2.txt");

#SE mates 1,2
my $in_S1 = File::Spec->catfile($indir, $pre."_mates1.txt");
my $in_S2 = File::Spec->catfile($indir, $pre."_mates2.txt");

print STDERR "paired-end files:\n";
print STDERR join("\n", ($in_P1, $in_P2)), "\n\n";

print STDERR "left files:\n";
print STDERR join("\n", ( $in_S1, $in_S2)), "\n\n";

die "wrong input" unless (-e $in_P1 and -e $in_P2 and -e $in_S1 and -e $in_S2);

my $out_p  = File::Spec->catfile($outdir, $pre."_paired_brat_bw_out.txt");
my $out_S1 = File::Spec->catfile($outdir, $pre . "_1_brat_bw_out.txt");
my $out_S2 = File::Spec->catfile($outdir, $pre. "_2_brat_bw_out.txt");

print STDERR "output:\n";
print STDERR join("\n", ($out_p, $out_S1, $out_S2)), "\n\n";

#die "output exists" unless (!(-e "$outdir/$out_p") and !(-e "$outdir/$out_S1") and !(-e "$outdir/$out_S2") );
die "output exists" if(-e $out_p or -e $out_S1 or -e $out_S2);
my $log = File::Spec->catfile($outdir, $pre."_brat_bw_log.txt");


my $cmd_S1 = "date; time $brat_bw -s $in_S1 -P $ref -o $out_S1 -m $m  >> $log 2>&1";
print STDERR $cmd_S1,"\n\n";
open(OUT, ">>$log") or die;
print OUT $cmd_S1, "\n\n";
close OUT;

if(!$debug){
	`$cmd_S1`;
}


my $cmd_S2 = "date; time $brat_bw -s $in_S2 -P $ref -o $out_S2 -m $m -A >> $log 2>&1";
print STDERR $cmd_S2,"\n\n";
open(OUT, ">>$log") or die;
print OUT "\n\n", $cmd_S2, "\n\n";
close OUT;
if(!$debug){
	`$cmd_S2`;
}

my $cmd_P = "date; time $brat_bw -1 $in_P1 -2 $in_P2 -pe -o $out_p -P $ref -i $i -a $a -m $m  >> $log 2>&1";
print STDERR $cmd_P,"\n\n";
open(OUT, ">>$log") or die;
print OUT "\n\n",  $cmd_P, "\n\n";
close OUT;
if(!$debug){
	`$cmd_P`;
}

exit;