#!/usr/bin/perl -w
 
use strict;
use File::Spec;

my $debug = 0;

my $usage = "$0 \n<bed_like_file> <isMeht_file> <output>\n\n";
die $usage unless (@ARGV == 3);

my $bed_file = shift or die;
my $isMeth_file = shift or die;
my $output = shift or die;


die unless (-e $bed_file);
die unless (-e $isMeth_file);
die if (-e $output);

if($debug){
	print STDERR "\n\nOK\n\n";
	exit;
}



my @bed_list;
my %pos;
read_bed_file($bed_file, \@bed_list);
record_region_pos(\@bed_list, \%pos);

my %nums;

read_isMeth_file($isMeth_file, \%pos , \%nums);
%pos = ();

output($output, \@bed_list, \%nums);

exit;

sub output{
	my ($file, $ref_list, $num_ref) = @_;
	die if ( -e $file);
	open(OUT, ">>$file") or die;
	my $last_index = scalar(@{$ref_list}) - 1;
	my $head = $ref_list->[0];
	print OUT join("\t", ($head, "CG_num", "CHG_num", "CHH_num")), "\n";
	
	foreach my $i (1..$last_index){
		my ($CG_num, $CHG_num, $CHH_num ) = (0) x 3;
		if(defined $num_ref->{$i}->{"CG"}) {$CG_num = $num_ref->{$i}->{"CG"}}
		if(defined $num_ref->{$i}->{"CHG"}) {$CHG_num = $num_ref->{$i}->{"CHG"}}
		if(defined $num_ref->{$i}->{"CHH"}) {$CHH_num = $num_ref->{$i}->{"CHH"}}
		print OUT join("\t", ($ref_list->[$i], $CG_num, $CHG_num, $CHH_num)), "\n";
	}	
	close(OUT);	
}

sub read_bed_file{
	my ($file, $ref) = @_;
	
	die unless (-e $file);
	
	open(IN, $file) or die;
	
	my $i = -1;
	while(<IN>){
		$i++;
		chomp;
		$ref->[$i] = $_;
	}
	
	close(IN);
}

sub record_region_pos{
	my ($ref_list, $ref_h) = @_;
	my $last_index = scalar(@{$ref_list}) - 1;
	
	foreach my $i (1..$last_index){
		my @a = split "\t", $ref_list->[$i];
		my ($chr, $start, $end ) = @a[0..2];
		
		for my $j($start..$end){
			$ref_h->{$chr}->[$j] = $i;
		}
	}
}

#chr	pos	strand	type	num_C	depth	percentage	isMeth
#chr1	100	+	CHH	0	11	0	0
#chr1	101	+	CHH	0	12	0	0
#chr1	102	+	CHH	0	12	0	1
sub read_isMeth_file{
	my ($file, $pos_ref, $num_ref) = @_;
	
	die unless (-e $file);
	open(IN, $file) or die;
	my $head = <IN>;
	while(<IN>){
		chomp;
		my @a = split "\t";
		my ($chr, $pos, $type) = ($a[0], $a[1], $a[3] );
		if(defined $pos_ref->{$chr}->[$pos]){
			my $index = $pos_ref->{$chr}->[$pos];
			$num_ref->{$index}->{$type}++;
		}
	}
	
	close(IN);
}