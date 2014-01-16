#!/usr/bin/perl -w

# v0.1:
# only consider 5 chrs

use strict;
use File::Spec;

#my %chr_len = ("chr1"=>30427671, "chr2"=>19698289, "chr3"=>23459830, "chr4"=>18585056, "chr5"=>26975502, "chrM"=>366924, "chrC"=>154478, "pUC19"=>2686);

my $debug = 0;
my $constant = 2000000;
if ($debug) {$constant = 300 ;}

if($debug){
	print  STDERR "debug = $debug\n\n";
}

my $usage = "$0 \n <indir>  <pre> <outdir> \n\n";
die $usage unless (@ARGV == 3);

my $indir  = shift or die "indir";
my $pre    = shift or die "pre";
my $outdir = shift or die "outdir";

#die "wrong dir" unless (-d $dir);
die "wrong indir" unless  (-d $indir);
die "wrong outdir" unless (-d $outdir);

my $in_forw = File::Spec->catfile($indir, $pre . "_forw.txt");
my $in_rev  = File::Spec->catfile($indir, $pre . "_rev.txt");

print  STDERR "input files:\n";
print  STDERR join("\n", ($in_forw, $in_rev)), "\n\n";

die "wrong input" unless (-e $in_forw and -e $in_rev );

my $output = File::Spec->catfile($outdir,$pre. "_cytosines_depth.txt");
die "$output exists" if (-e $output);
print STDERR "output:\t$output\n\n";

die if($debug);


open (OUT, ">>$output") or die;

my %depths;

#chr1    79      79      CHH:27  0.148148        -
# 0		  1		 2			3		4			 5
open(FORW, $in_forw) or die;
open(REV, $in_rev) or die;

while(<FORW>){
	chomp;
	my @a = split "\t";
	last if ($a[0] eq "chrC");
	my ($type, $depth) = split ":", $a[3];
	$depths{$depth}++;	
}
close(FORW);

while(<REV>){
	chomp;
	my @a = split "\t";
	last if ($a[0] eq "chrC");
	my ($type, $depth) = split ":", $a[3];
	$depths{$depth}++;
}
close(REV);

print OUT join("\t", ("depth", "count")), "\n";
foreach my $key (sort {$a<=>$b} keys %depths){
	print OUT join ("\t", ($key, $depths{$key})), "\n";
}

close(OUT);

exit;
