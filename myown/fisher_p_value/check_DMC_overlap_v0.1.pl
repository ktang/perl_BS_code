#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

#v0.1
#used for three file
# 7 numbers

use strict;
use File::Spec;

my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 \n <in1> <label1> <in2> <label2> <in3> <label3> STDOUT";
die $usage unless(@ARGV == 6);

my $first_file = shift or die;
die unless (-e $first_file);

my $label1 = shift or die;

my $second_file = shift or die;
die unless (-e $second_file);

my $label2 = shift or die;


my $third_file = shift or die;
die unless (-e $third_file);

my $label3 = shift or die;



my (%DMC1 , %DMC2 , %DMC3 , %total_DMC);

my $num1 = read_file($first_file, \%DMC1, \%total_DMC ,  1);
my $num2 = read_file($second_file, \%DMC2, \%total_DMC , 2);
my $num3 = read_file($third_file, \%DMC3, \%total_DMC ,  4);

# my $num1 = read_first_file($first_file, \%DMCs);

my @numbers = (0) x 8;

foreach my $chr(sort keys %total_DMC){
	foreach my $pos(keys %{$total_DMC{$chr}}){
		my $index = $total_DMC{$chr}->{$pos};
		$numbers[$index]++;
	}
}

#1:label1 only
#2:label2 only
#3:label1 and label2
#4:label3 only
#5:label1 and label3
#6:label2 and label3
#7:all three

print STDERR "\n";

print STDERR "$label1\t", $num1, "\n";
print STDERR "$label2\t", $num2, "\n";
print STDERR "$label3\t", $num3, "\n\n";


print STDERR "$label1 only:", $numbers[1], "\n";
print STDERR "$label1 and $label2 only:", $numbers[3], "\n";
print STDERR "$label2 only:", $numbers[2], "\n";
print STDERR "$label1 and $label3 only:", $numbers[5], "\n";
print STDERR "all three:", $numbers[7],  "\n";
print STDERR "$label2 and $label3 only:", $numbers[6], "\n";
print STDERR "$label3 only:", $numbers[4], "\n\n";

exit;

#read_file($first_file, \%DMC, \%total_DMC , );
sub read_file{
	my ($file, $ref , $total_ref, $num) = @_;
	die unless (-e $file);
	open(IN, $file) or die;
	my $i = 0;

	while(<IN>){
		$i++;
		chomp;
		my @a = split "\t";
		my ($chr, $pos) = @a[0..1];
		$ref->{$chr}->[$pos] = 1;
		$total_ref->{$chr}->{$pos} += $num;
	}	
	close(IN);
	return $i;
}

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