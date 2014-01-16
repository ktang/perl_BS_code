#!/usr/bin/perl -w

use strict;

my $acgt_count = "/Users/tang58/Software/BRAT/brat-1.1.19/acgt-count";

my $usage ="$0 <indir> <outdir> <pre> <C24|Col0>";

die $usage unless (@ARGV == 4);

my ($indir, $outdir, $pre, $flag) = @ARGV[0..3];

my $ref = "NONE";

if($flag eq "C24"){
	$ref = "";
}elsif($flag eq "Col0"){
	$ref = "/Users/tang58/DataBase/BRAT_reference_files/Col0/Col0_references_names.txt";
}else{
	die $usage, "\nC24 or Col0, exactly";
}

die "wrong indir" unless (-d $indir);

die "wrong outdir" unless (-d $outdir);


my $in_p = $pre."_paired_brat_out.txt";
my $in_S1 = $pre . "_1_brat_out.txt";
my $in_S2 = $pre. "_2_brat_out.txt";

print STDERR "brat_file:\n";
print STDERR join("\n", ($in_p, $in_S1, $in_S2)), "\n";

die "wrong input" unless ((-e "$indir/$in_p") and (-e "$indir/$in_S1") and (-e "$indir/$in_S2") );

my ($single_list, $pair_list) = ($pre."_single_list.txt", $pre."_pair_list.txt");

die "wrong list file" unless (!(-e "$indir/$single_list") and !(-e "$indir/$pair_list"));

open(SOUT, ">$indir/$single_list") or die "Can't open $single_list: $!";
open(POUT, ">$indir/$pair_list") or die "Can't open $pair_list: $!";

print SOUT "$indir/$in_S1","\n";
print SOUT "$indir/$in_S2","\n";

print POUT "$indir/$in_p", "\n";

close SOUT;
close POUT;

my $cmd = "date; time $acgt_count -r $ref -P $pre -p $indir/$pair_list -s $indir/$single_list -B";

print STDERR "$cmd\n";
system("$cmd");


exit;