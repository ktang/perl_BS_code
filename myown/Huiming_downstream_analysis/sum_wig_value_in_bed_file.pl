#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;

my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = " \n $0 <col_label> <in_wig> <in_bed_like_file> <output> \n\n ";
die $usage unless(@ARGV == 4);

my $label    = shift or die;
my $wig_file = shift or die;
my $bed_file = shift or die;
my $output   = shift or die;

die unless (-e $wig_file);
die unless (-e $bed_file);
die if (-e $output);

my @lists;
my %wigs;

if($debug){
	print STDERR "\n\nOK\n\n";
	exit;
}

read_bed($bed_file, \@lists);
read_wig($wig_file, \%wigs);
cal_sum_in_bed( \@lists, \%wigs );

output($output, \@lists, $label);

exit;

#read_bed($bed_file, \@lists);
sub read_bed{
	my ($file, $ref_a) = @_;
	
	die unless (-e $file);
	
	open(IN, $file) or die "$file : $!";
	
	@{$ref_a} = <IN>;
	
	close(IN);	
}

#read_wig($wig_file, \%wigs);
sub read_wig{
	my ( $file, $ref ) = @_;
	die unless (-e $file);
	open(IN, $file) or die "cannot open $file: $!";
	
	my $chr = "";
	while (<IN>){
		next if (/^track/);
		if (/variableStep\s+chrom=(\w+)/){
			$chr = lc $1;
			next;
		}
		chomp;
		my ($pos, $val) = split /\t/;
	#	$wigs{$chr} -> {$pos} = $val;
		$ref->{$chr}->{$pos} = abs($val);
	}
	close(IN);
}

#cal_sum_in_bed( \@lists, \%wigs );
sub cal_sum_in_bed{
	my ($ref_a, $ref_h) = @_;
	my $last_index = scalar(@{$ref_a}) - 1;
	for my $i(1..$last_index){
		my $l = $ref_a ->[$i];
		chomp $l;
		my @a = split "\t", $l;
		my ($chr, $start, $end) = @a[0..2];
		my $sum = 0;
		for my $j ($start .. $end){
			if(defined $ref_h->{$chr}->{$j} ){
				$sum += $ref_h->{$chr}->{$j} ;
			}
		}
		$ref_a ->[$i] = join("\t", (@a, $sum)); 
	}
}

#output($output, \@lists);
sub output{
	my ($file, $ref_a, $label_sub) = @_;
	die if (-e $file);
	open(OUT, ">>$file") or die;
	my $head = $ref_a->[0];
	
	chomp $head;
	
	print OUT join("\t", ($head, $label_sub)), "\n";
	
	my $last_index = scalar(@{$ref_a}) - 1;
	for my $i(1..$last_index){
		print OUT $ref_a->[$i], "\n";
	}
	close(OUT);
}