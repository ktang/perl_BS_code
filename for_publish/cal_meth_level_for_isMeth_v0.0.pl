#!/usr/bin/perl -w

use strict;
use File::Spec;

my $debug = 0;
my $constant = 2000000;
if ($debug) {
	$constant = 300 ;
}

if($debug){
	print  STDERR "debug = $debug\n\n";
}

my %labels = ("yes"=>1, "no"=>1);

my $usage = "$0 \n <isMeth_file> <sample_name> <head_yes_no> STDOUT\n\n";
die $usage unless (@ARGV == 3);

#my $indir = shift or die "indir";
#my $infile = shift or die "infile";
#my $outdir = shift or die "outdir";
#my $pre = shift or die "pre";

my $isMeth = shift or die;
my $sample = shift or die;
my $print_head = shift or die;

die $usage unless ( defined $labels{$print_head} );

my @types = ("total", "CG", "CHG", "CHH");

my (%called_mC_nums,%C_type_nums); # Fraction of methylated cytosines
my %mean_meth_level_sum; #Mean methylation level
my ( %seq_mC_sum,  %seq_depth_sum ); #Weighted methylation level

Initialization(\%called_mC_nums);
Initialization(\%C_type_nums);
Initialization(\%mean_meth_level_sum);
Initialization(\%seq_mC_sum);
Initialization(\%seq_depth_sum);


my $max_dep = 0;

open(IN, $isMeth) or die "cannot open $isMeth: $!";
# 0		 1			2	  3			4		5		6				7
#chr     pos     strand  type    num_C   depth   percentage      isMeth
#chr1    334     +       CHH     0       20      0       0
#chr1    336     +       CHH     0       20      0       0
#0	 1	2	 3	 4	 5	 6	 7

my $head = <IN>;
chomp $head;
my @h = split "\t", $head;

die unless ($h[-1] eq "isMeth");

if($debug){
	close IN;
	print STDERR "\n\nOK\n\n";
	exit;
}

while(<IN>){
	chomp;
	my @a = split "\t";
			
	#my $per = $a[6];
	#my ($type, $depth) = ( $a[3], $a[5]);
	#my $C_num = $a[4];
	
	my ($type, $num_C, $depth,   $per,  $isMeth) = @a[3..7];
	
	$C_type_nums{$type}++;
	$C_type_nums{"total"}++;
	
	
	if($depth > $max_dep){
		$max_dep = $depth;
	}
	if($isMeth == 0){
		$num_C = 0;
		$per   = 0;
	}elsif( $isMeth == 1){
		$called_mC_nums{$type}++;
		$called_mC_nums{ "total" }++;
	}else{
		die $_, "\n\n";
	}
	
	$mean_meth_level_sum{$type} += $per;
	$mean_meth_level_sum{ "total" } += $per;
	
	$seq_mC_sum{$type} += $num_C;
	$seq_mC_sum{"total" } += $num_C;
	
	
	$seq_depth_sum{$type} += $depth;
	$seq_depth_sum{"total" } += $depth;
	
	
	
#my (%called_mC_nums,%C_type_nums); # Fraction of methylated cytosines
#my %mean_meth_level_sum; #Mean methylation level
#my ( %seq_mC_sum,  %seq_depth_sum ); #Weighted methylation level
	
	
}

close (IN);

print STDERR "max_depth = $max_dep \n\n";

#weighted; fraction;  mean

check(\%called_mC_nums);
check(\%C_type_nums);
check(\%mean_meth_level_sum);
check(\%seq_mC_sum);
check(\%seq_depth_sum);

my @weighted_level = ("NA") x 4;
my @weighted_div = ("0/0=") x 4;
my @frac_level = ("NA") x 4;
my @mean_level = ("NA") x 4;

my @called_mC_nums_a = (0) x 4;
my @C_type_nums_a    = (0) x 4;
my @mean_meth_level_sum_a = (0) x 4;

weight_meth(\%seq_mC_sum, \%seq_depth_sum, \@weighted_level , \@weighted_div);

call_meth(\%called_mC_nums, \%C_type_nums, \@frac_level);
call_meth(\%mean_meth_level_sum , \%C_type_nums, \@mean_level);



assign(\%called_mC_nums, \@called_mC_nums_a);
assign(\%C_type_nums, \@C_type_nums_a);
assign(\%mean_meth_level_sum , \@mean_meth_level_sum_a);

if(  $print_head eq "yes" ){
	print join("\t", ( "sample",
			  "weighted_CG_level", "weighted_CHG_level", "weighted_CHH_level", "weighted_C_level",
			  "weighted_CG_div", "weighted_CHG_div", "weighted_CHH_div", "weighted_C_div",
			  
			  "frac_CG_level", "frac_CHG_level", "frac_CHH_level", "frac_C_level",
			  "mCG_num", "mCHG_num", "mCHH_num", "mC_num",
			  "CG_num", "CHG_num", "CHH_num", "C_num",
			  
			  "sum_CG_level", "sum_CHG_level", "sum_CHH_level", "sum_C_level",
			  "mean_CG_level", "mean_CHG_level", "mean_CHH_level", "mean_C_level"
			  
			  ) )  , "\n";
}

print join("\t", ($sample,
		  @weighted_level,
		  @weighted_div,
		  @frac_level,
		  @called_mC_nums_a,
		  @C_type_nums_a,
		  @mean_meth_level_sum_a,
		  @mean_level		  
		  )), "\n";


exit;

sub assign{
	my ($first_ref, $ref_a) = @_;
	
	my @types_sub = ("CG", "CHG", "CHH", "total");
	my $i = -1;
	foreach my $type_sub (@types_sub){
		$i++;
		$ref_a->[$i] = $first_ref->{$type_sub};
	}
}

sub call_meth{
	my ($first_ref, $second_ref, $ref_a) = @_;
	
	my @types_sub = ("CG", "CHG", "CHH", "total");
	my $i = -1;
	foreach my $type_sub (@types_sub){
		$i++;
		my ($first, $second) = ($first_ref->{$type_sub}  , $second_ref->{$type_sub} ); 
		if( $second != 0){
			$ref_a->[$i] = sprintf( "%.3f", 100 * $first / $second );
		}
	}
	
}

# weight_meth (\%seq_mC_sum, \%seq_depth_sum, \@weighted_level , \@weighted_div);
sub weight_meth{
	my ($mC_ref_h, $dep_ref_h, $level_ref_a, $div_ref_a) = @_;
	my @types_sub = ("CG", "CHG", "CHH", "total");
	my $i = -1;
	foreach my $type_sub (@types_sub){
		$i++;
		my ($mC_sub, $dep_sub) = ($mC_ref_h->{$type_sub}  , $dep_ref_h->{$type_sub} ); 
		if( $dep_sub != 0){
			$level_ref_a->[$i] = sprintf( "%.3f", 100 * $mC_sub /$dep_sub );
			$div_ref_a->[$i] = "$mC_sub/$dep_sub=";
		}
	}
}

sub Initialization{
	my ($ref) = @_;
	my @types_sub = ("total", "CG", "CHG", "CHH");
	foreach my $type_sub (@types_sub){
		$ref-> {$type_sub} = 0;
	}
}

sub check{
	my ($ref) = @_;
	my @types_sub = ("CG", "CHG", "CHH");
	my $sum = 0;
	foreach my $type_sub (@types_sub){
		$sum += $ref-> {$type_sub} ;
	}
	
	my $t = $ref->{"total"};
	
	if ( abs( $sum - $t ) > 0.01 ){
		print STDERR "maybe wrong \n";
		print STDERR "sum = $sum != $t \n\n";
	}
}