#!/usr/bin/perl -w

#v0.0_filter_first_then_merge
# get raw 200-bp from p-adj file, filter with 2fold reduce 
# then combine them
# then filter 10 DMC
use strict;
use File::Spec;


my $print_debug = 1;

my $debug = 0;
if($debug){
	print STDERR "debug = 1\n\n";
}
my $bin_size     = 200;
my $sliding_size = 50;
my $const = 4000000;
my $allowed_gap  = 200;
my $covered_num = int ($bin_size / $sliding_size);

my $DMC_cutoff  = 10;
my $depth_cutoff = 4;
my $p_value_cutoff = 0.05;
my %chr_len = ("chr1"=>30427671, "chr2"=>19698289, "chr3"=>23459830, "chr4"=>18585056, "chr5"=>26975502);
my @types = ( "CG", "CHG", "CHH");
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


my @chrs = ("chr1", "chr2", "chr3", "chr4", "chr5");

#my $usage = "$0 <WT_isMeth> <mut_isMeth> <outdir> <outpre>";
my $usage = "$0 <WT_isMeth> <mut_isMeth> <in_adj_p_file> <outpre>";
die $usage unless (@ARGV == 4);

my $wt_in  			= shift or die;
my $mut_in 			= shift or die;
my $temp_p_adj_file = shift or die;
#my $output 			= shift or die;
my $outpre 			= shift or die;
my $output =  $outpre."_BinSize". $bin_size . "_sliding" . $sliding_size . "_BH_adj_p_value". $p_value_cutoff . "_depth". $depth_cutoff. "_filter2fold_first.txt";



die $wt_in  unless (-e $wt_in);
die $mut_in unless (-e $mut_in);

die $temp_p_adj_file unless (-e $temp_p_adj_file);
die "$output exists" if(-e $output);
print STDERR "output: $output\n\n";

if($debug){
	exit;
}

open(OUT, ">>$output") or die;


#my (%CG_num, %CHG_num, %CHH_num); 		# number of CXX that meet the requirement that both depth >= $depth_cutoff
#my (%CG_depth, %CHG_depth, %CHH_depth);	#depth sum
#my (%mCG, %mCHG, %mCHH);               	#mCXX sum



my @raw_list;
my %records_index; #$recores_index{$chr}->[i] = 1;

#############
#	read list
############################
open(ADJ, $temp_p_adj_file) or die "cannot open $temp_p_adj_file: $!";
my $h = <ADJ>;
while(<ADJ>){
	chomp;
	my @a = split "\t";
	my ($chr, $start, $end) = @a[0..2];
	if(  $a[3] <= $p_value_cutoff or  $a[4] <= $p_value_cutoff or $a[5] <= $p_value_cutoff ){
		push @raw_list, join("\t", ($chr, $start, $end));
		my $base_index = int( ($start - 1) / $sliding_size);
		$records_index{$chr}->[$base_index] = 1;
	} 
}
close(ADJ);

print STDERR "raw_list:" , scalar(@raw_list), "\n";
if ( $print_debug ){
	for my $i(0..19){
		print STDERR $raw_list[$i], "\n";
	}
}

open (WT,  $wt_in  )  or die "cannot open  $wt_in: $!";
open (MUT, $mut_in )  or die "cannot open $mut_in: $!";

my $h_wt  = <WT>;
my $h_mut = <MUT>;
my ( $l_wt, $l_mut );
my $line = 0;

my %small_window_per_list;
my %DMCs;

while($l_wt = <WT>, $l_mut = <MUT>){
	$line++;
	if($line % $const == 0){
		print STDERR $line, "\tdone\n";
	}
	chomp $l_wt;
	chomp $l_mut;
	my @a_wt  = split "\t", $l_wt;
	my @a_mut = split "\t", $l_mut;
	
	die unless ($a_wt[1] == $a_mut[1]);
	my $chr = $a_wt[0];

	
	my $dep_wt = $a_wt[5];
	my $dep_mut = $a_mut[5];
	
	next if($chr =~ /chrC|chrM|pUC19/);
		
	if($dep_wt >= $depth_cutoff and $dep_mut >= $depth_cutoff){
		my $pos = $a_wt[1];
		my $type = $a_wt[3];

		#my $mC_wt = $a_wt[4];
		#my $mC_mut = $a_mut[4];
	
		my $per_wt  = $a_wt[6];
		my $per_mut = $a_mut[6];
		
		if($per_mut != 0){
			if($per_wt / $per_mut <= 0.5){
				$DMCs{$chr}->[$pos] = $type;
			}
		}

		
			
		my $base_index = int( ($pos - 1) / $sliding_size);
		
		for my $i (($base_index - $covered_num)..($base_index + $covered_num)){
			if($i >= 0 and defined $records_index{$chr}->[$i]) {
				my $interval_start = $i * $sliding_size + 1;
				my $interval_end   = $interval_start + $bin_size - 1;
				if($pos >= $interval_start and $pos <= $interval_end){
					push @{ $small_window_per_list{"WT"}->[$i] }, $per_wt;
					push @{ $small_window_per_list{"mut"}->[$i] }, $per_mut;
				}
			}
		}
	}
}
close(WT);
close(MUT);

#####################
#	filter
##################################
my @filter_raw_list = ();


foreach my $i (0..$#raw_list){
	#my ($chr, $start, $end, $num) = split "\t", $border_list[$i];
	my ( $chr, $start, $end ) = split "\t", $raw_list[$i];
	
	my $base_index = int( ($start - 1) / $sliding_size);
	
	die "( $chr, $start, $end ) " unless (defined $records_index{$chr}->[$base_index] );
	
	die unless (defined $small_window_per_list{"WT"}->[$base_index]);
	die unless (defined $small_window_per_list{"WT"}->[$base_index]);
	
	my $avge_wt  =  cal_mean(\@{$small_window_per_list{"WT"}->[$base_index]});
	my $avge_mut =  cal_mean(\@{$small_window_per_list{"mut"}->[$base_index]});
	
	if($avge_wt eq "NONE" or $avge_mut eq "NONE" ){
		print STDERR join("\t", ( $chr, $start, $end , $avge_wt, $avge_mut)), "\n";
	}elsif($avge_mut > 0){
		if($avge_wt / $avge_mut <= 0.5){
			#print OUT join("\t", (@a,  $avge_mut ,$avge_wt)), "\n";
			push @filter_raw_list, join("\t", ( $chr, $start, $end,  $avge_mut ,$avge_wt ));
		}
	}
}

#####################
#	generate merged_filter_raw_list
##################################

my @merged_filter_raw_list;
my %regions_in_list;

print STDERR "filter_raw_list:" , scalar(@filter_raw_list), "\n";
if ( $print_debug ){
	for my $i(0..19){
		print STDERR $filter_raw_list[$i], "\n";
	}
}



my ($last_chr, $last_start, $last_end ) = ("chr0", 0, 0);


foreach my $i(0..$#filter_raw_list){
	my ( $chr, $start, $end,  $avge_mut ,$avge_wt ) = split "\t", $filter_raw_list[$i];
	if($chr ne $last_chr || $start > $last_end + $allowed_gap + 1 ){
			if($last_chr ne "chr0") {
				push @merged_filter_raw_list ,    join("\t", ( $last_chr, $last_start, $last_end)) ;
			}
			($last_chr, $last_start, $last_end ) =  ($chr, $start, $end);
		}else{
			$last_end = $end;
	}
}


###########
#	generate final_border_list
#####################
my @final_border_list;

print STDERR "merged_filter_raw_list:" , scalar(@merged_filter_raw_list), "\n";

foreach my $i ( 0..$#merged_filter_raw_list ){
	my ($chr, $start, $end) = split "\t", $merged_filter_raw_list[$i];
	my $DMC_num = 0;
	for my $j ($start..$end){
		if(defined $DMCs{$chr}->[$j]){
			$DMC_num++;
		}
	}
	if($DMC_num >= $DMC_cutoff ){
		my ($real_start, $real_end) = ($start, $end) ;
		for my $j($start..$end){
			if(defined $DMCs{$chr}->[$j]){
				$real_start = $j;
				last;
			}
		}
		for (my $j = $end; $j >= $real_start; $j--){
			if(defined $DMCs{$chr}->[$j]){
				$real_end = $j;
				last;
			}
		}
		push @final_border_list, join("\t" , ($chr, $real_start, $real_end, $DMC_num));
	}
}

print STDERR "final_border_list:" , scalar(@final_border_list), "\n";


foreach my $i (0..$#final_border_list){
	my ($chr, $start, $end) = split "\t", $final_border_list[$i];
	for my $j($start..$end){
		$regions_in_list{$chr}->[$j] = $i;
	}
}


#####################
#	read again
##################################

open (WT,  $wt_in )  or die "cannot open  $wt_in: $!";
open (MUT, $mut_in ) or die "cannot open $mut_in: $!";

$h_wt  = <WT>;
$h_mut = <MUT>;

 $line = 0;
my (%mC_list, %dep_list);
my %per_list;

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
		my $mC_mut = $a_mut[4];
		
		my $per_wt  =  $a_wt[6];
		my $per_mut =  $a_mut[6];
		
		$mC_list{"WT"}->{$type}->[$order] += $mC_wt;
		$mC_list{"mut"}->{$type}->[$order] += $mC_mut;

		$dep_list{"WT"}->{$type}->[$order] += $dep_wt;
		$dep_list{"mut"}->{$type}->[$order] += $dep_mut;
	
		push @{ $per_list{"WT"}->[$order] } , $per_wt;
		push @{ $per_list{"mut"}->[$order] } , $per_mut;
	}
}
close(WT);
close(MUT);

foreach my $i (0..$#final_border_list){
	my ($chr, $start, $end, $num) = split "\t", $final_border_list[$i];
	my (%divide_mut, %divide_wt);
	my (%per_mut, %per_wt);

	foreach my $type (@types){
		my ($mC_wt , $dep_wt ) = (0) x 2;
		my ($mC_mut, $dep_mut) = (0) x 2;
		
		if(defined $dep_list{"WT"} ->{$type}->[$i] ) { $dep_wt  = $dep_list{"WT"} ->{$type}->[$i] }
		if(defined $dep_list{"mut"}->{$type}->[$i] ) { $dep_mut = $dep_list{"mut"}->{$type}->[$i] }
		if(defined $mC_list {"WT"} ->{$type}->[$i] ) { $mC_wt   = $mC_list {"WT"} ->{$type}->[$i] }
		if(defined $mC_list {"mut"}->{$type}->[$i] ) { $mC_mut  = $mC_list {"mut"}->{$type}->[$i] }
		$divide_mut{$type} = "$mC_mut/$dep_mut=";
		$divide_wt {$type} = "$mC_wt/$dep_wt=";
		if($dep_mut != 0) {
			$per_mut{$type} = sprintf("%.4f",  100 * $mC_mut/$dep_mut);
		}else{
			$per_mut{$type} = "NA";
		}	
		if($dep_wt != 0) {
			$per_wt{$type} = sprintf("%.4f",  100 * $mC_wt/$dep_wt);
		}else{
			$per_wt{$type} = "NA";
		}	
	}
	my $avge_wt  =  cal_mean(\@{$per_list{"WT"}->[$i]});
	my $avge_mut =  cal_mean(\@{$per_list{"mut"}->[$i]});
	
	if($avge_wt eq "NONE" or $avge_mut eq "NONE" ){
		print STDERR join("\t", ($chr, $start, $end, $num, 
							 $divide_mut{"CG"}, $per_mut{"CG"},  $divide_mut{"CHG"}, $per_mut{"CHG"}, $divide_mut{"CHH"}, $per_mut{"CHH"},
					 		 $divide_wt{"CG"}, $per_wt{"CG"},    $divide_wt{"CHG"}, $per_wt{"CHG"},   $divide_wt{"CHH"}, $per_wt{"CHH"} , $avge_mut ,$avge_wt )), "\n";

	}elsif($avge_mut > 0){
		if($avge_wt / $avge_mut <= 0.5){
			print OUT join("\t", ($chr, $start, $end, $num, 
							 $divide_mut{"CG"}, $per_mut{"CG"},  $divide_mut{"CHG"}, $per_mut{"CHG"}, $divide_mut{"CHH"}, $per_mut{"CHH"},
					 		 $divide_wt{"CG"}, $per_wt{"CG"},    $divide_wt{"CHG"}, $per_wt{"CHG"},   $divide_wt{"CHH"}, $per_wt{"CHH"} , $avge_mut ,$avge_wt )), "\n";
		}
	}
}


close( OUT );
exit;

#if($per_mut != 0){
#			if($per_wt / $per_mut <= 0.5){
#				$DMCs{$chr}->[$pos] = $type;
#			}
#		}


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
		#return ($sum / $num);
		return ( sprintf ( "%.4f", 100 * $sum / $num));
	}
}
