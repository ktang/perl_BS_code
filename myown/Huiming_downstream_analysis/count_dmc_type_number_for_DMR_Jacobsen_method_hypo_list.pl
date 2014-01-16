#!/usr/bin/perl -w
 
use strict;
use File::Spec;

my $debug = 0;
if($debug){
	print STDERR "debug = 1\n\n";
}

my @types = ("CG", "CHG", "CHH");
my $const = 2000000;
my $depth_cutoff = 5;
#my $p_value_cutoff = 0.05;
my %chr_len = ("chr1"=>30427671, "chr2"=>19698289, "chr3"=>23459830, "chr4"=>18585056, "chr5"=>26975502);

my @chrs = ("chr1", "chr2", "chr3", "chr4", "chr5");


my $usage = "$0 \n <WT_isMeth> <mut_isMeth> <in_bed_list> <output>\n\n";

die $usage unless (@ARGV == 4);

my $wt_in  = shift or die;
my $mut_in = shift or die;
my $input = shift or die;
my $output = shift or die;
die $wt_in  unless (-e $wt_in);
die $mut_in unless (-e $mut_in);

print STDERR "output: $output\n\n";
die unless (-e $input);
die if (-e $output);

open (OUT, ">$output") or die;


my %regions_in_list; 
my @border_list;

read_list($input, \@border_list);

my @heads = split "\t", $border_list[0];
print OUT join("\t", (@heads, "DMC_CG", "DMC_CHG", "DMC_CHH" )), "\n";

foreach my $i (1..$#border_list){
#	my ($chr, $start, $end, $num) = split "_", $border_list[$i];
	my @a = split "\t", $border_list[$i];
	my ($chr, $start, $end) = @a[0..2];
	for my $j($start..$end){
		$regions_in_list{$chr}->[$j] = $i;
	}
}


my ($l_wt, $l_mut);
#chr     pos     strand  type    num_C   depth   percentage      isMeth
#chr1    1       +       CHH      0       0       0	           0
# 0      1       2         3      4       5       6                7

open (WT,  $wt_in )  or die "cannot open  $wt_in: $!";
open (MUT, $mut_in ) or die "cannot open $mut_in: $!";
my $h_wt  = <WT>;
my $h_mut = <MUT>;

my $line = 0;


#my (%mC_list, %dep_list);
#my %per_list;

my @numbers;

while($l_wt = <WT>, $l_mut = <MUT>){
	$line++;
	if($line % $const == 0){
		print STDERR $line, "\tdone\n";
	}
	chomp $l_wt;
	chomp $l_mut;
	my @a_wt  = split "\t", $l_wt;
	my @a_mut = split "\t", $l_mut;
	my $pos = $a_wt[1];

#	die unless ($a_wt[1] == $a_mut[1]);
	my $chr = $a_wt[0];
	my $type = $a_wt[3];
	
	next unless (defined $regions_in_list{$chr}->[$pos]);
	my $order = $regions_in_list{$chr}->[$pos];	
	
	my $dep_wt = $a_wt[5];
	my $dep_mut = $a_mut[5];

	if($dep_wt >= $depth_cutoff and $dep_mut >= $depth_cutoff){
		my $mC_wt = $a_wt[4];
	#	my $mC_mut = $a_mut[4];
		my $per_wt  =  $a_wt[6];
		my $per_mut =  $a_mut[6];
		
		my $type = $a_wt[3];
		
	#	$mC_list{"WT"}->{$type}->[$order] += $mC_wt;
	#	$mC_list{"mut"}->{$type}->[$order] += $mC_mut;
	#	push @{ $per_list{"WT"}->[$order] } , $per_wt;
	#	push @{ $per_list{"mut"}->[$order] } , $per_mut;
	
		#if($per_wt == 0){
		#	if($mC_mut >= 2){
		#		$hyper_DMCs{$chr}->[$pos] = 1;
		#	}
		#}elsif($per_mut / $per_wt >= 2){
		#	$hyper_DMCs{$chr}->[$pos] = 1;
		#}
	
		if($per_mut == 0){
			if($mC_wt >= 2){
				#$hypo_DMCs{$chr}->[$pos] = 1;
				$numbers[$order]->{$type} ++;
			}
		}elsif($per_wt / $per_mut >= 2){
			#$hypo_DMCs{$chr}->[$pos] = 1;
			$numbers[$order]->{$type} ++;
		}
	}
}
close(WT);
close(MUT);


foreach my $i (1..$#border_list){
	#my ($chr, $start, $end, $num) = split "_", $border_list[$i];
	my @a = split "\t", $border_list[$i];
	
	my ($cg, $chg, $chh) = (0) x 3;
	
	if(defined $numbers[$i]->{"CG"}) {  $cg = $numbers[$i]->{"CG"}}
	if(defined $numbers[$i]->{"CHG"}) { $chg = $numbers[$i]->{"CHG"}}
	if(defined $numbers[$i]->{"CHH"}) { $chh = $numbers[$i]->{"CHH"}}
	
	print OUT join("\t", (@a, $cg, $chg, $chh)), "\n";
	
}


close(OUT);

exit;

#read_list($input, \@list_array);
sub read_list{
	my ($file, $ref) = @_;
	die unless (-e $file);
	my $i = -1;
	open(IN, $file) or die;
	while(<IN>){
		$i++;
		chomp;
	#	my @a = split "\t";
	#	$ref->[$i] = join("_", @a);
		$ref->[$i] = $_;
	}
}


sub cal_mean{
	my ($ref) = @_;
	my $sum = 0;
	my $num = 0;
	
	foreach my $item (@{$ref}){
		$sum+=$item;
		$num++;
	}
	if($num == 0){
		return "NONE";
	}else{
		return ($sum / $num);
	}
}