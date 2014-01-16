#!/usr/bin/perl -w

#chr     start   end     hypo_DMC_num    hypo_DMC_type   hyper_DMC_num   hyper_DMC_type  informative_C_num       CG      CHG     CHH     wt_CG   wt_CG_per       wt_CHG  wt_CHG_per      wt_CHH  wt_CHH_per      mut_CG  mut_CG_per      mut_CHG mut_CHG_per     mut_CHH mut_CHH_per     wt_C_per        mut_C_per       wt_MethLevel_perC       mut_MethLevel_perC
#chr1    880854  881979  11      11+0+0  2       2+0+0   477     57      116     304     977/1727=       56.5721 3/3669= 0.0818  13/9362=        0.1389  571/1202=       47.5042 5/2499= 0.2001  5/5767= 0.0867  6.7286  6.1365  6.9620  5.9933
	#																						11						13				15						17						19				21

use strict;
#use Text::NSP::Measures::2D::Fisher2::twotailed;

use Text::NSP::Measures::2D::Fisher2::right;
#high >= low


#my $usage = "$0 \n <input> <index4/7> STDOUT\n\n";
#die $usage unless (@ARGV == 2);

my $usage = "$0 \n <input>  STDOUT\n\n";
die $usage unless (@ARGV == 1);

my $input = shift or die "input";

die unless (-e $input);

my $index = 11;  #shift or die; #7;
open(IN, $input) or die;

my $head = <IN>;

chomp $head;
my @pts = split "\t", $head;
die unless ($pts[$index] =~ /_CG$/);

print join("\t", (@pts, "p_CG", "p_CHG", "p_CHH")), "\n";

while(<IN>){
	chomp;
	my @a = split "\t";
	
	my @p = (1) x 3;
	
	for my $i(0..2){
		my ($mC_high, $dep_high) = get_two($a[$index + 2*$i]);
		my ($mC_low, $dep_low)   = get_two($a[$index + 6 + 2*$i]);
		
		# my $npp = $mC1 + $non_mC1 + $mC2 + $non_mC2; #37;
		# my $n1p = $mC1 + $mC2; #27;
		# my $np1 = $mC1 + $non_mC1; #20; 
		# my $n11 = $mC1; #14;
		
		my $npp =  $dep_high + $dep_low;
		my $n1p =  $mC_high  + $mC_low;
		my $np1 =  $dep_high;
		my $n11 =  $mC_high;
		
		if( $npp != 0){
			my  $p_value = calculateStatistic( n11=>$n11,
              		                       n1p=>$n1p,
                    		               np1=>$np1,
                            		       npp=>$npp);
			$p[$i] = $p_value;
		}
	}
	print join("\t", (@a, @p)), "\n";
}

close(IN);


exit;

sub get_two{
	my ($temp) = @_;
	
	die unless ($temp =~/=/);
	$temp =~ s/=//;
	return split "\/", $temp;
}