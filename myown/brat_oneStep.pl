#!/usr/bin/perl -w

use strict;

my $brat = "/Users/tang58/Software/BRAT/brat-1.1.18/brat";

my $usage ="$0 <indir> <outdir> <pre> <C24|Col0>";

die $usage unless (@ARGV == 4);

my ($indir, $outdir, $pre, $flag) = @ARGV[0..3];

my $ref = "NONE";

if($flag eq "C24"){
	$ref = "/Users/tang58/DataBase/BRAT_reference_files/C24/C24_references_names.txt";
}elsif($flag eq "Col0"){
	$ref = "/Users/tang58/DataBase/BRAT_reference_files/Col0/Col0_references_names.txt";
}else{
	die $usage, "\nC24 or Col0, exactly";
}

die "wrong indir" unless (-d $indir);

die "wrong outdir" unless (-d $outdir);

my $in_P1 = $pre."_1_pairs_brat_input.txt";
my $in_P2 = $pre."_2_pairs_brat_input.txt";
my $in_S1 = $pre."_1_single_brat_input.txt";
my $in_S2 = $pre."_2_single_brat_input.txt";

print STDERR join("\n", ($in_P1, $in_P2, $in_S1, $in_S2)), "\n\n";

die "wrong input" unless (-e "$indir/$in_P1" and -e "$indir/$in_P2" and -e "$indir/$in_S1" and -e "$indir/$in_S2");

my $out_p = $pre."_paired_brat_out.txt";
my $out_S1 = $pre . "_1_brat_out.txt";
my $out_S2 = $pre. "_2_brat_out.txt";

print STDERR "output:\n";
print STDERR join("\n", ($out_p, $out_S1, $out_S2)), "\n";

die "output exists" unless (!(-e "$outdir/$out_p") and !(-e "$outdir/$out_S1") and !(-e "$outdir/$out_S2") );

my $log = $pre."_brat_log.txt";

my $cmd_S1 = "date; time $brat -r $ref -s $indir/$in_S1 -m 2 -bs -o $outdir/$out_S1  >> $outdir/$log 2>&1";
my $cmd_S2 = "date; time $brat -r $ref -s $indir/$in_S2 -m 2 -bs -o $outdir/$out_S2 -A >> $outdir/$log 2>&1";
my $cmd_P = "date; time $brat -r $ref -1 $indir/$in_P1 -2 $indir/$in_P2 -u -i 0 -a 1000 -m 2 -pe -bs -o $outdir/$out_p >> $outdir/$log 2>&1 ";

print STDERR "$cmd_P\n";
system("$cmd_P");

print STDERR "$cmd_S1\n";
system ("$cmd_S1");

print STDERR "$cmd_S2\n";
system("$cmd_S2");

exit;