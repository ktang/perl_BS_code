#!/usr/bin/perl -w

# v2.0
# use (hyper-hypo) instead hyper
# use (hypo-hyper) instead hypo

#v1.1
#add two col in the output total_bp_hyper total_bp_hypo

#oringinal
#modified from 
#/Volumes/Macintosh_HD_2/idm1_new_met_data/script/Liu_algorithm_Oct1_data_setreproduce_Liu_Algorithm_hyper_and_hypo_v1.2.2.pl

#v1.0
# just use wig file as input, hyper_DMC defined as mC in mut but not in colA

use strict;
use File::Spec;
#$ | = 1;

my $debug = 0;
if($debug){
	print STDERR "debug:$debug\n";
	print STDERR "input is special,read the script\n";
}


my $window_size = 100;
my $allowed_gap = 100;
my $cutoff_for_small_region =10;
my $cutoff_for_merged_region = 25;

my %chr_len = ("chr1"=>30427671, "chr2"=>19698289, "chr3"=>23459830, "chr4"=>18585056, "chr5"=>26975502);


my $usage = "$0 <input> <outdir> <pre>";
die $usage unless (@ARGV == 3);

my $input = shift or die "input";
my $outdir = shift or die "outdir";
my $pre = shift or die "pre";
#my $p_flag = shift or die "p-value";

die unless (-e $input and -d $outdir);

my $colA_wig = "/Volumes/My_Book/20120427_ShangHai_data/call_methylation/1_isMeth2wig/colA_mC.wig";
die "colA_wig" unless (-e $colA_wig);



my $postfix = "wig"."_netDMC_WinSize" . $window_size . "_gap" . $allowed_gap . "_initialCutoff" . $cutoff_for_small_region . "_reportCutoff"
              . $cutoff_for_merged_region . ".txt";

my $hyper_output = File::Spec->catfile($outdir, $pre. "_hyper_" . $postfix);
my $hypo_output  = File::Spec->catfile($outdir, $pre. "_hypo_" . $postfix);

print STDERR "output:\n";
print STDERR join("\n", ($hyper_output, $hypo_output)), "\n\n";

die "output exists" if(-e  $hyper_output or -e $hypo_output);


my (%pos_hypers, %pos_hypos);


read_wig($colA_wig, \%pos_hypos); #mC in colA only is hypo DMC
read_mut_wig($input, \%pos_hypers, \%pos_hypos);


	
	
#cal_density(\%records, \%density);
#cal_density{#from 0, $window_size

my (%den_hypers, %den_hypos);

cal_density(\%pos_hypers, \%den_hypers, $window_size );
cal_density(\%pos_hypos,  \%den_hypos,  $window_size );


# select 
#my ($ref_den, $ref_pos, $ref_line) = @_;

	my (@line_hyper, @line_hypo);


#	select_region(\%pos_hypers, \%den_hypers, \@line_hyper, $window_size, $allowed_gap, $cutoff_for_small_region, $cutoff_for_merged_region); #old
#	select_region(\%pos_hypos,  \%den_hypos,  \@line_hypo , $window_size, $allowed_gap, $cutoff_for_small_region, $cutoff_for_merged_region); #old

	select_region(\%pos_hypers, \%den_hypers, \@line_hyper, $window_size, $allowed_gap, $cutoff_for_small_region, $cutoff_for_merged_region, \%pos_hypos,  \%den_hypos); 
	select_region(\%pos_hypos,  \%den_hypos,  \@line_hypo , $window_size, $allowed_gap, $cutoff_for_small_region, $cutoff_for_merged_region, \%pos_hypers, \%den_hypers);

	my ($num_hyper, $num_hypo);

	$num_hyper = scalar (@line_hyper);
	$num_hypo  = scalar (@line_hypo);
	
#	my $bp_hyper = sum_length(\@line_hyper);
#	my $bp_hypo  = sum_length(\@line_hypo);
	
	#my $bp_per_hyper = sprintf("%.1f", $bp_hyper / $num_hyper) ;
	#my $bp_per_hypo  = sprintf("%.1f", $bp_hypo / $num_hypo) ;

#	my ($bp_per_hyper , $bp_per_hypo) = (0, 0);
	
#	if($num_hyper > 0){ $bp_per_hyper = sprintf("%.1f", $bp_hyper / $num_hyper) ;}
#	if($num_hypo > 0) { $bp_per_hypo  = sprintf("%.1f", $bp_hypo / $num_hypo) ;}

#	print STDERR "hyper_DMR: $num_hyper\n";
#	print STDERR "hypo_DMR:  $num_hypo\n";
	
output_hyper(\@line_hyper, $hyper_output);
output_hypo(\@line_hypo,  $hypo_output);


exit;

#read_wig($colA_wig, \%pos_hypos); #mC in colA only is hypo DMC
#read_mut_wig($input, \%pos_hypers, \%pos_hypos);

sub output_hyper{
	my ($ref, $output) = @_;
	die "$output exists" if (-e $output);
	
	open (OUT, ">$output") or die "cannot open $output: $!";
	print OUT join("\t", ("chr", "start", "end", "hyper_mC_num")), "\n";
	for my $ref_region(@{$ref}){
		#my ($chr, $start, $end, $num) = @{$ref_region2};
		print OUT join("\t",  @{$ref_region}), "\n";
	}
	close(OUT);	
}

sub output_hypo{
	my ($ref, $output) = @_;
	die "$output exists" if (-e $output);
	
	open (OUT, ">$output") or die "cannot open $output: $!";
	print OUT join("\t", ("chr", "start", "end", "hypo_mC_num")), "\n";
	for my $ref_region(@{$ref}){
		#my ($chr, $start, $end, $num) = @{$ref_region2};
		print OUT join("\t",  @{$ref_region}), "\n";
	}
	close(OUT);	
}


sub sum_length{
	my ($ref) = @_;
	my $len = 0;
	foreach my $ref_region(@{$ref}){
		$len += ($ref_region->[2] - $ref_region->[1] + 1);
	}
	return $len;
}

sub read_wig{
	my ($input,$ref) = @_;
	print STDERR "read_wig reading $input...\t";
	open (IN, $input) or die "cannot open $input:$!";
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
		$ref->{$chr}->{$pos} = 1;
	}
	close(IN);
	print STDERR "DONE\n";
}

sub read_mut_wig{
	my ($input, $hyper_ref, $hypo_ref) = @_;
	print STDERR "read_wig reading $input...\t";
	open (IN, $input) or die "cannot open $input:$!";
	my $chr = "";
	while (<IN>){
		next if (/^track/);
		if (/variableStep\s+chrom=(\w+)/){
			$chr = lc $1;
			next;
		}
		chomp;
		my ($pos, $val) = split /\t/;
		if(defined $hypo_ref->{$chr}->{$pos}){
			delete $hypo_ref->{$chr}->{$pos};
		}else{
			$hyper_ref->{$chr}->{$pos} = 1;
		}
	}
	close(IN);
	print STDERR "DONE\n";

}



sub overlap{# first in second ; return percentage
	my ($ref_first, $ref_second) = @_;
	
	my %beds;
	
	my $num_overlap = 0;
	for my $ref_region2(@{$ref_second}){
		my ($chr, $start, $end, $num) = @{$ref_region2};
		for my $i ($start..$end){
			$beds{$chr}->[$i] = 1;
		}		
	}

	for my $ref_region1(@{$ref_first}){
		my ($chr, $start, $end, $num) = @{$ref_region1};
		
		my $flag = 0;
		
		for my $i ($start..$end){
			if(defined $beds{$chr}->[$i]){
				$flag = 1;
				last;
			}
		}
		
		if ($flag == 1){$num_overlap ++}

	}
	
	my $total = scalar (@{$ref_first});
	if($total == 0){
		return (0, "0.00");
	}
	else{
	#	return (sprintf ("%.1f",100 * $num_overlap / $total)."%" )
		return ($num_overlap, sprintf ("%.2f",100 * $num_overlap / $total));
	}
}

#	select_region(\%pos_hypers, \%den_hypers, \@line_hyper, $window_size, $allowed_gap, $cutoff_for_small_region, $cutoff_for_merged_region, \%pos_hypos,  \%den_hypos); 
#	select_region(\%pos_hypos,  \%den_hypos,  \@line_hypo , $window_size, $allowed_gap, $cutoff_for_small_region, $cutoff_for_merged_region, \%pos_hypers, \%den_hypers);


sub select_region{# $window_size, $allowed_gap, $cutoff_for_small_region, $cutoff_for_merged_region
	my ($ref_pos, $ref_den, $ref_line, $window_size_sub, $allowed_gap_sub, $cutoff_for_small_region_sub, $cutoff_for_merged_region_sub, $ref_pos_opp, $ref_den_opp) = @_;
	my ($last_start, $last_end)  = (0, 0);
	my $last_chr = "NONE";
	
	foreach my $chr(sort keys %chr_len){
		my $len = $chr_len{$chr};
	
		foreach my $i(0..(int($len/$window_size_sub)+1)){
			my $temp_den = 0;
			
			if(defined ${$ref_den}{$chr}->[$i] ){
				$temp_den = ${$ref_den}{$chr}->[$i];
				
				if(defined $ref_den_opp->{$chr}->[$i]){
					$temp_den -= ($ref_den_opp->{$chr}->[$i]);
				}
			}
		
			if ( $temp_den >= $cutoff_for_small_region_sub ){
				
				#my ($new_start, $new_end) = ($i*$window_size+1, ($i+1) * $window_size);#regular_start, end
				
				my ($regular_start, $regular_end) = ($i*$window_size_sub + 1, ($i+1) * $window_size_sub);#regular_start, end from win_size
			 
				my ($new_start, $new_end) = (-1, -1);
				
				for my $k($regular_start..$regular_end){
					if (defined ${$ref_pos}{$chr}->{$k}){
						$new_start = $k;
						last;
					}
				}
				
				for(my $k = $regular_end; $k >= $regular_start; $k --){
					if (defined ${$ref_pos}{$chr}->{$k}){
						$new_end = $k;
						last;
					}
				}
				
				die "$regular_start, $regular_end" if ($new_start == -1 or $new_end == -1);
			
				if($chr ne $last_chr || $new_start > $last_end + $allowed_gap_sub){ # not overlap
					my $num = 0;
					my $num_opp = 0;
				
					foreach my $k ($last_start..$last_end){
						if(defined ${$ref_pos}{$last_chr}->{$k}){ $num++; }
						if(defined $ref_pos_opp->{$last_chr}->{$k}) {$num_opp++;}
					}
				
					if(  ($num - $num_opp) >= $cutoff_for_merged_region_sub ){
						push (@{$ref_line}, [$last_chr, $last_start, $last_end, "$num-$num_opp"]);
						
					}
					($last_chr, $last_start, $last_end) = ($chr, $new_start, $new_end);
				}
				
				else{
					#$last_end = ($i+1) * $window_size;
					my $flag = 0;
					for(my $k = ($i+1) * $window_size_sub; $k >= $i*$window_size_sub; $k --){
						if (defined ${$ref_pos}{$chr}->{$k}){
							$last_end = $k;
							$flag = 1;
							last;
						}
					}
					
					die "$i" unless ($flag);
					
					
				}
			}#if ( $temp_den >= $cutoff_for_small_region )
		}
	}
}
#										WT									mut
#chr1    108     +       CHG     9       11      0.818182        1       8       0.125   0.00547749464158133
# 0		  1		 2		  3		 4		  5			6			 7		 8			9			10

sub read_file{
	my ($file, $ref_hyper, $ref_hypo) = @_;
#	print STDERR "reading $file...\t";
	open (IN, $file) or die "cannot open $file: $!";

	while(<IN>){
		chomp;
		my @a = split "\t";
		if($a[9] > $a[6]){
			${$ref_hyper}{$a[0]}->{$a[1]} = 1;
		}			
		else{
			${$ref_hypo}{$a[0]}->{$a[1]} = 1;
		}
	}
	close(IN);
}

sub cal_density{#from 0, $window_size
	my ($ref_file, $ref_den, $window_size_sub) = @_;
	
	#print STDERR "cal_density...\t";
#	foreach my $chr (sort keys %{$ref_file} ){
#		foreach my $pos (sort {$a <=> $b} keys %{ ${$ref_file}{$chr}} ){
	foreach my $chr (keys %{$ref_file} ){
		foreach my $pos ( keys %{ ${$ref_file}{$chr}} ){
			
			
			my $base_index = int( ($pos - 1) / $window_size_sub);
			${$ref_den}{$chr}->[$base_index]++;
		}
	}
	#print STDERR "DONE\n";
}


sub cal_met{
	my ($ref) = @_;
	my $num = 0;
	foreach my $chr (sort keys %{$ref}){
		$num += scalar(keys %{${$ref}{$chr}})
	}
	return $num;
}