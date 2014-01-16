#!/usr/bin/perl -w
#v1.0 July 15
# try to for each small region, narrow down first, then check gap

#method4:
#netDMC / basic mC >= cutoff
#basic mC = mentholated in both samples.


use strict;
use File::Spec;

my $debug_sub = 1;
my $debug = 0;
if($debug){
	print STDERR "debug:$debug\n";
	
}




#my $usage = "$0 \n  <wig1>  <wig2> STDOUT \n\n";
my $usage = "$0 \n  <isMet1> <label1>  <isMet2> <label2> <min_cutoff> <max_cutoff> STDOUT \n\n";
die $usage unless (@ARGV == 6);

my $file1      = shift or die  "file1" ;
my $label1 	   = shift or die  "label1";
my $file2 	   = shift or die  "file2" ;
my $label2	   = shift or die  "label2";
my $min_cutoff = shift or die;
my $max_cutoff = shift or die;

die unless ($max_cutoff > $min_cutoff and $min_cutoff > 0);

die unless (-e $file1);
die unless (-e $file2);

my (%num1_only, %num2_only, %both);
#0		 1			2		3		4		5		6				7
#chr     pos     strand  type    num_C   depth   percentage      isMeth
#chr1    1       +       CHH     0       0       0       0
#chr1    2       +       CHH     0       0       0       0

open(IN1, $file1) or die;
open(IN2, $file2) or die;

my $head1 = <IN1> ;
my $head2 = <IN2> ;
chomp $head1;
chomp $head2;


my @h1 = split "\t", $head1;
my @h2 = split "\t", $head2;

die unless ($h1[3] eq "type" and $h1[5] eq "depth" and $h1[7] eq "isMeth" );
die unless ($h2[3] eq "type" and $h2[5] eq "depth" and $h2[7] eq "isMeth" );

my ($l1, $l2);

while($l1 = <IN1>, $l2 = <IN2>){
	chomp $l1;
	chomp $l2;
	my @a1 = split "\t", $l1;
	my @a2 = split "\t", $l2;
	die unless ($a1[1] == $a2[1]);
	
	my $dep1 = $a1[5];
	my $dep2 = $a2[5];
	
	if ($dep1 >= $min_cutoff and $dep2 >= $min_cutoff and $dep1 <= $max_cutoff and $dep2 <= $max_cutoff){
		my $type = $a1[3];
		my $isMet1 = $a1[7];
		my $isMet2 = $a2[7];
		if($isMet1 == 1 and $isMet2 == 1) {
			$both{$type} ++;
			$both{C}++;
		}elsif($isMet1 == 1 and $isMet2 == 0) {
			$num1_only{$type} ++ ; 
			$num1_only{C}++ ; 
		}elsif($isMet1 == 0 and $isMet2 == 1) {
			$num2_only{$type} ++ ; 
			$num2_only{C}++ ; 
		}elsif($isMet1 == 0 and $isMet2 == 0) {
			next;
		}else{
			print STDERR $l1 , "\n", $l2, "\n\n";
			exit;
		}
	}
}
close(IN1);
close(IN2);

my @types = ("C", "CG", "CHG", "CHH");

print "depth_cutoff:", $min_cutoff, "-", $max_cutoff, "\n";
print join("\t", ("type", $label1, $label2, "both")), "\n";



foreach my $type (@types){
	my  ($num1, $num2, $num_both) = (0) x 3;
	if(defined $num1_only{$type}) {$num1     = $num1_only{$type}}
	if(defined $num2_only{$type}) {$num2     = $num2_only{$type}}
	if(defined $both{$type} )     {$num_both = $both{$type} }
	print join("\t", ($type,$num1, $num2, $num_both )), "\n";
}



exit;