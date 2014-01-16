#!/usr/bin/perl -w

# get_p_value_for_list.pl
# get original p-value for the list

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;

my $debug = 1;

my $const = 100000;

my $p_cutoff = 0.00001;

if($debug){
	print STDERR "debug = 1\n\n";
}

my $usage = "$0 \n <input_list> <p-value_list> STDOUT \n\n";
die $usage unless(@ARGV == 2);


my $input = shift or die;

my $p_value_file = shift or die;

die unless (-e $input);
die unless (-e $p_value_file);

my %records;
my @lists;

read_list($input, \%records, \@lists);

open(IN, $p_value_file)  or die;
my $h = <IN>;
#print $h;

my $i = 0;

while(<IN>){
	$i++;
	if($i % $const == 0){
		print STDERR $i , " DONE\n";
	}
	chomp;
	my @a = split "\t";
	my ($chr, $start) = @a[0..1];
	if(defined $records{$chr}->[$start]){
	#	print $_, "\n";
		if($a[3] <= $p_cutoff or $a[5] <= $p_cutoff or $a[7] <= $p_cutoff){
			my $index = $records{$chr}->[$start];
			$lists[$index] = 1;
		}
	}
}

close(IN);

foreach my $j (0..$#lists){
	if($lists[$j] =~ /chr/ ){
		#print  $lists[$j] , "\n";
		my @a = split "_", $lists[$j];
		print $a[0]  .":". $a[1] ."-". $a[2], "\n"; 
	}
}

exit;

#read_list($input, \%records);

sub read_list{
	my ($file, $ref, $ref_array) = @_;
	open(IN, $file) or die;
	my $head = <IN>;
	
	my $l = -1;
	
	while(<IN>){
		$l++;
		chomp;
		my @a = split "\t";
		my ($chr, $start, $end) = @a[0..2];
		$ref_array ->[$l] = join("_", ($chr, $start, $end));

		my $index_s = int (($start - 1) / 50);
		my $index_e = int(($end - 1) / 50 );
		for my $i( $index_s - 4..$index_e +4 ){
			if($i>=0){
				my $pos = $i * 50 + 1;
				$ref->{$chr}->[$pos] = $l;
			}
		}
	}	
	close(IN);
	print STDERR "Done\n\n";
}