#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;

my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 \n <in1> <label1> <in2> <label2> STDOUT";
die $usage unless(@ARGV == 4);

my $first_file = shift or die;
die unless (-e $first_file);

my $label1 = shift or die;

my $second_file = shift or die;
die unless (-e $second_file);

my $label2 = shift or die;

my %DMCs;

my $num1 = read_first_file($first_file, \%DMCs);

my ($num2, $both) = read_second_file($second_file, \%DMCs);
print STDERR "\n";

print STDERR "$label1\t", $num1, "\n";
print STDERR "$label2\t", $num2, "\n\n";





print STDERR "overlap between $label1 and $label2:\n", $both, "\n";
print STDERR "$label1 only:\t", $num1  - $both, "\n";
print STDERR "$label2 only:\t", $num2 - $both, "\n\n";



exit;

sub read_first_file {
	my ($file, $ref ) = @_;
	die unless (-e $file);
	open(IN, $file) or die;
	my $i = 0;
	while(<IN>){
		$i++;
		chomp;
		my @a = split "\t";
		my ($chr, $pos) = @a[0..1];
		$ref->{$chr}->[$pos] = 1;
	}	
	close(IN);
	return $i;
}

sub read_second_file{
	my ($file, $ref ) = @_;
	die unless (-e $file);
	open(IN, $file) or die;
	my $i = 0;
	my $both = 0;
	while(<IN>){
		$i++;
		chomp;
		my @a = split "\t";
		my ($chr, $pos) = @a[0..1];
		if(defined 	$ref->{$chr}->[$pos] ){
			$both ++;
		}
	}	
	close(IN);
	return ($i, $both);
}