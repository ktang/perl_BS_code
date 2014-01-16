#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;

my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 \n<dir1> <pre1> <dir2> <pre2> <outdir> <outpre>\n\n";
die $usage unless(@ARGV == 6);

my $dir1 = shift or die;
my $pre1 = shift or die;
my $dir2 = shift or die;
my $pre2 = shift or die;
my $outdir = shift or die;
my $outpre = shift or die;

my $forw_in1 = File::Spec->catfile ($dir1 ,$pre1 . "_forw.txt");
my $forw_in2 = File::Spec->catfile ($dir2 ,$pre2 . "_forw.txt");
my $forw_out = File::Spec->catfile ($outdir, $outpre . "_forw.txt");

die unless (-e $forw_in1);
die unless (-e $forw_in2);
die if (-e $forw_out);


my $rev_in1 = File::Spec->catfile ($dir1 ,$pre1 . "_rev.txt");
my $rev_in2 = File::Spec->catfile ($dir2 ,$pre2 . "_rev.txt");
my $rev_out = File::Spec->catfile ($outdir, $outpre . "_rev.txt");

die unless (-e $rev_in1);
die unless (-e $rev_in2);
die if (-e $rev_out);

if($debug){
	
	print STDERR join("\n", ($forw_out, $rev_out)), "\n\n";
	
	print STDERR "OK\n\n";
	exit;
}

combine($forw_in1, $forw_in2, $forw_out);
combine($rev_in1,  $rev_in2,  $rev_out);

exit;
#chr1    33      33      CHG:13  0.153846        -
#chr1    79      79      CHH:27  0.148148        -
#chr1    109     109     CG:17   0.882353        -

sub combine{
	my ($in1, $in2, $output) = @_;
	die unless (-e $in1);
	die unless (-e $in2);
	die if (-e $output);
	open(IN1, $in1) or die;
	open(IN2, $in2) or die;
	open(OUT, ">$output") or die;
	
	while(my $l1 = <IN1>, my $l2 = <IN2>){
		chomp $l1;
		chomp $l2;
		
		my @a1 = split "\t", $l1;
		my @a2 = split "\t", $l2;
		die unless ($a1[1] == $a2[1]);
		
		my ($t1, $dep1) = split ":", $a1[3];
		my ($t2, $dep2) = split ":", $a2[3];
		
		my $mC1 = round( $dep1 * $a1[4]);
		my $mC2 = round( $dep2 * $a2[4]);
		my $dep_sum = $dep1 + $dep2;
		my $per = 0;
		if($dep_sum >0){
			$per = sprintf ( '%.6g', ($mC1 + $mC2) / $dep_sum);
		}
		$a1[3] = join(":", ($t1, $dep_sum));
		$a1[4] = $per;
		print OUT join("\t", @a1),"\n"; 
		
	}
	
	close IN1;
	close IN2;
	close OUT;
}

sub round {
    my($number) = shift;
    #return int($number + .5);
    return int($number + .5 * ($number <=> 0)); # take care of negative numbers too
}
