#!/usr/bin/perl -w

use strict;
use File::Spec;

my $debug = 0;

my $bin_size     = 500000;
my $sliding_size = 500000;
my $covered_num = int ($bin_size / $sliding_size);

my %chr_len = ("chr1"=>30427671, "chr2"=>19698289, "chr3"=>23459830, "chr4"=>18585056, "chr5"=>26975502);

my %chr_len_cumulative = ("chr1"=> 0 , "chr2"=>30427671, "chr3"=>50125960, "chr4"=>73585790, "chr5"=>92170846);

my %bin_last_index;
# bin start from 0
for my $chr (sort keys %chr_len){
	my $len = $chr_len{$chr};
	my $last = int( ($len - 1) / $sliding_size);
	$bin_last_index{$chr} = $last;
}
if($debug){
	foreach my $chr (sort keys %bin_last_index){
		print STDERR join("\t", ($chr, $bin_last_index{$chr})), "\n";
	}
}




my $usage = "$0 \n <DMR_list> <outdir> <outpre> \n\n";
#label can be ros1,ros2,ros3,ros4,rdd,hsp,phd,zdp
die $usage unless (@ARGV == 3);

my $input = shift or die;
my $outdir = shift or die;
my $outpre = shift or die;

die unless (-d $outdir);
die unless (-e $input);

my $output = File::Spec->catfile($outdir, $outpre . "_DMR_number_WinSize" . $bin_size . "_sliding" . $sliding_size . ".txt");

die if( -e $output);


if($debug){
	print STDERR  $output, "\nOK\n\n";
	exit;
}


my %DMR_numbers;

read_DMR_file( $input, \%DMR_numbers  );

output( $output, \%DMR_numbers, \%bin_last_index, \%chr_len);



exit;

#read_DMR_file( $input, \%DMR_numbers  );

sub read_DMR_file{
	my ($file, $ref) = @_;
	die unless (-e $file);
	
	open(IN, $file) or die "cannot open $file : $!";
	my $head = <IN>;
	while (<IN>){
		chomp;
		my @a = split "\t";
		my ($chr, $start, $end) = @a[0..2];
		my $pos = int( ( $start + $end) / 2 );
		
		
		my $base_index = int( ($pos - 1) / $sliding_size);
		
		for my $i (($base_index - $covered_num)..($base_index + $covered_num)){
			if($i >= 0){
				my $interval_start = $i * $sliding_size + 1;
				my $interval_end   = $interval_start + $bin_size - 1;
				if($pos >= $interval_start and $pos <= $interval_end){
					
					$ref->{$chr}->[$i] ++;
				}
			}
		}
		
	}
	close(IN);
}

# open (OUT, ">$output") or die "cannot open $output: $!";
# print OUT join("\t", ("chr", "start", "end", "num_CG_plus", "num_CG_minus", "num_CHG_plus", "num_CHG_minus", "num_CHH_plus", "num_CHH_minus")), "\n";

#output( $output, \%DMR_numbers, \%bin_last_index );

sub output{
	my ($file, $ref_num, $ref_index, $ref_len) = @_;
	die if(-e $file);
	open (OUT, ">$file") or die "cannot open $file: $!";
	print OUT join("\t", ("chr", "start", "end", "DMR_number", "cumulative_pos" ) ) , "\n";
	
	foreach my $chr (sort keys %{$ref_index}){
		my $last_index = $ref_index->{$chr};
		
		for my $i (0..( $last_index - 1 )){
			my $start = $i * $sliding_size + 1;
			my $end   = $start + $bin_size - 1;
			my $num = 0;
			if (defined $ref_num->{$chr}->[$i]){
				$num = $ref_num->{$chr}->[$i];
			}
			
			print OUT join("\t", ($chr, $start, $end, $num, $chr_len_cumulative{$chr} + int( ( $start + $end) / 2 ) ) ) , "\n";
		}
		
		my $start = $last_index * $sliding_size + 1;
		my $end   = $ref_len->{$chr};
		
		my $num = 0;
		if (defined $ref_num->{$chr}->[$last_index]){
			$num = $ref_num->{$chr}->[$last_index];
		}
		print OUT join("\t", ($chr, $start, $end, $num ,$chr_len_cumulative{$chr} + int( ( $start + $end) / 2 ) ) ) , "\n";
		
	}
	
	close(OUT);
}
