#!/usr/bin/perl -w
#v1.2
# just change input order

#v1.0 July 15
# try to for each small region, narrow down first, then check gap

#method4:
#netDMC / basic mC >= cutoff
#basic mC = mentholated in both samples.


use strict;
use File::Spec;
$| = 1;


my $debug = 0;


my $debug_sub = 0;
if($debug){
	print STDERR "debug:$debug\n";
	print STDERR "input is special,read the script\n";
}




my %chr_len = ("chr1"=>30427671, "chr2"=>19698289, "chr3"=>23459830, "chr4"=>18585056, "chr5"=>26975502);


#my $usage = "$0 <input> <outdir> <pre> <win> <gap> <small> <final>";
#die $usage unless (@ARGV == 7);

#my $usage = "$0 <control_pre> <input> <outdir> <pre> <winSize> <sliding_bp> <merged_gap_cutoff> <netDMC_cutoff> <percentage_cutoff>";
#die $usage unless (@ARGV == 9);

my $usage = "$0 \n <winSize> <sliding_bp> <merged_gap_cutoff> <netDMC_cutoff> <percentage_cutoff> <control_input> <control_pre> <outdir> <input>  <pre> \n\n";
die $usage unless (@ARGV == 10);

my $window_size 	  = shift or die "window_size";
my $sliding_bp		  = shift or die "sliding_bp";
my $merged_gap_cutoff = shift or die "merged_gap_cutoff";
my $netDMC_cutoff     = shift or die "netDMC_cutoff";
my $percentage_cutoff = shift or die "percentage_cutoff";

my $colA_wig          = shift or die "control_input";
my $control_pre		  = shift or die "control_pre";

my $outdir 			  = shift or die "outdir";

my $input  			  = shift or die "input";
my $pre    			  = shift or die "pre";

die unless (-e $input and -d $outdir);

#my $colA_wig = "/Volumes/My_Book/20120427_ShangHai_data/call_methylation/1_isMeth2wig/colA_mC.wig";
#my $colA_wig = "/Volumes/Macintosh_HD_2/BS_data_HiSeq/C24_background/cal_methy/wig/JKZ19_C24_Luc.wig";
#my $colA_wig = "/Volumes/Macintosh_HD_2/BS_data_HiSeq/C24_background/brat_bw_new_analysis/nodupl/nodupl_pooled_wig/C24_WT_nodupl_mC.wig";
#my $colA_wig = "/Volumes/My_Book/20120702_SH_extra_BS_Seq/wig/chrM_error_rate_less/col_0_pooled_mC.wig";
die "colA_wig" unless (-e $colA_wig);



my $postfix = "_method4_win" . $window_size . "_sliding" . $sliding_bp . "_gap" . $merged_gap_cutoff . "_netDMC" . $netDMC_cutoff . "_per" . $percentage_cutoff .".txt";

my $hyper_output = File::Spec->catfile($outdir, $pre. "_hyper_vs_" . $control_pre	 . $postfix);
my $hypo_output  = File::Spec->catfile($outdir, $pre. "_hypo_vs_" . $control_pre	 . $postfix);
#my $output = File::Spec->catfile($outdir, $pre. $postfix);

print STDERR "output:\n";
print STDERR join("\n", ($hyper_output, $hypo_output)), "\n\n";

die "output exists" if(-e  $hyper_output or -e $hypo_output);
#die "output exists" if(-e  $output);

#open (OUT, ">$output") or die;


my (%pos_hypers, %pos_hypos);
my %pos_both;


if($debug){
	print STDERR "read_wig...\n";
}
read_wig($colA_wig, \%pos_hypos); #mC in colA only is hypo DMC
#read_mut_wig($input, \%pos_hypers, \%pos_hypos);
read_mut_wig($input, \%pos_hypers, \%pos_hypos, \%pos_both);


	
	
#cal_density(\%records, \%density);
#cal_density{#from 0, $window_size

my (%den_hypers, %den_hypos);
my %den_both;

if($debug){
	print STDERR "cal_density...\n";
}
cal_density(\%pos_hypers, \%den_hypers, $window_size, $sliding_bp );
cal_density(\%pos_hypos,  \%den_hypos,  $window_size, $sliding_bp );
cal_density(\%pos_both,   \%den_both,   $window_size, $sliding_bp );

#print OUT join("\t", ("chr", "start", "curr_end" , "hyperDMC_num", "hypoDMC_num", "netDMC_num", "mC_num", "percentage(netDMC/mC)")), "\n";


####################################################
# head
####################################################
=head
foreach my $chr (sort keys %chr_len){
	my $length = $chr_len{$chr};
#	my $num_of_interval = int($length / $sliding_bp) + 1;
	for my $i (0..int($length / $sliding_bp)){
		my $curr_start = $i * $sliding_bp + 1;
		my $curr_end = $curr_start + $window_size - 1;
		if($curr_end > $length ){
			$curr_end = $length;
		}
		
		my ($hyperDMC_num, $hypoDMC_num, $netDMC_num, $mC_num, $per) = (0) x 5;
		if(defined $den_hypers{$chr}->[$i]) {$hyperDMC_num = $den_hypers{$chr}->[$i] ; }
		if(defined $den_hypos{$chr}->[$i] ) {$hypoDMC_num  = $den_hypos{$chr}->[$i]  ; }
		if(defined $den_both{$chr}->[$i]  ) {$mC_num       = $den_both{$chr}->[$i]   ; }
	
		$netDMC_num = $hyperDMC_num - $hypoDMC_num;
		if($mC_num != 0){
			$per = sprintf("%.2f", $netDMC_num / $mC_num * 100);
		}
		else{
			if($netDMC_num != 0){
				$per = "$netDMC_num / 0";
			}
		}
	
		print OUT join("\t", ($chr, $curr_start, $curr_end, $hyperDMC_num, $hypoDMC_num, $netDMC_num, $mC_num, $per )), "\n";
	
	}
	
}

=cut
####################################################
#   cut
####################################################

# select 
#my ($ref_den, $ref_pos, $ref_line) = @_;







	my (@line_hyper, @line_hypo);


#	select_region(\%pos_hypers, \%den_hypers, \@line_hyper, $window_size, $allowed_gap, $cutoff_for_small_region, $cutoff_for_merged_region); #old
#	select_region(\%pos_hypos,  \%den_hypos,  \@line_hypo , $window_size, $allowed_gap, $cutoff_for_small_region, $cutoff_for_merged_region); #old

#	select_region(\%pos_hypers, \%den_hypers, \@line_hyper, $window_size, $allowed_gap, $cutoff_for_small_region, $cutoff_for_merged_region, \%pos_hypos,  \%den_hypos); 
#	select_region(\%pos_hypos,  \%den_hypos,  \@line_hypo , $window_size, $allowed_gap, $cutoff_for_small_region, $cutoff_for_merged_region, \%pos_hypers, \%den_hypers);
if($debug){
	print STDERR "select_region...\n";
}
	select_region(\%pos_hypers, \%den_hypers, \@line_hyper, $window_size, $sliding_bp, $merged_gap_cutoff, $netDMC_cutoff, $percentage_cutoff, 
				  \%pos_hypos,  \%den_hypos, \%pos_both, \%den_both);
				  
	select_region(\%pos_hypos,  \%den_hypos, \@line_hypo,   $window_size, $sliding_bp, $merged_gap_cutoff, $netDMC_cutoff, $percentage_cutoff, 
			      \%pos_hypers, \%den_hypers, \%pos_both, \%den_both);
				  
				  

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

	print STDERR "hyper_DMR: $num_hyper\n";
	print STDERR "hypo_DMR:  $num_hypo\n\n";

if($debug){
	print STDERR "output...\n";
}
output_hyper(\@line_hyper, $hyper_output);
output_hypo(\@line_hypo,  $hypo_output);


exit;

#read_wig($colA_wig, \%pos_hypos); #mC in colA only is hypo DMC
#read_mut_wig($input, \%pos_hypers, \%pos_hypos);

sub output_hyper{
	my ($ref, $output) = @_;
	die "$output exists" if (-e $output);
	
	open (OUT, ">$output") or die "cannot open $output: $!";
	print OUT join("\t", ("chr", "start", "end", "hyper_mC", "hypo_mC", "netmC", "basic_mC", "percentage")), "\n";
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
	print OUT join("\t", ("chr", "start", "end", "hypo_mC", "hyper_mC","netmC", "basic_mC", "percentage")), "\n";
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
	#print STDERR "read_wig reading $input...\t";
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
#	print STDERR "DONE\n";
}

sub read_mut_wig{
	my ($input, $hyper_ref, $hypo_ref, $both_ref) = @_;
#	print STDERR "read_wig reading $input...\t";
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
			$both_ref->{$chr}->{$pos} = 1;
		}else{
			$hyper_ref->{$chr}->{$pos} = 1;
		}
	}
	close(IN);
#	print STDERR "DONE\n";

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


sub select_region{
	#select_region(\%pos_hypers, \%den_hypers, \@line_hyper, $window_size, $sliding_bp, $merged_gap_cutoff, $netDMC_cutoff, $percentage_cutoff, 
#				  \%pos_hypos,  \%den_hypos, \%pos_both,   \%den_both);

	my ($ref_pos, $ref_den, $ref_line, $window_size_sub, $sliding_bp_sub, $merged_gap_cutoff_sub, $netDMC_cutoff_sub, $percentage_cutoff_sub, 
		$ref_pos_opp, $ref_den_opp, $ref_pos_both, $ref_den_both) = @_;
	my ($last_chr, $last_start, $last_end)  = ("chr0", 0, 0);

	#	my @chr_array = ("chr2");
	
	foreach my $chr(sort keys %chr_len){
	#foreach my $chr (@chr_array){
		my $len = $chr_len{$chr};
	
		foreach my $i(0..(int( ($len - 1) / $sliding_bp_sub ) ) ){
		
			my ($posDMC_num, $oppDMC_num, $mC_num) = (0) x 3; #pos = positive; opp = opponent
			
					
			if(defined $ref_den->{$chr}->[$i]     ) {$posDMC_num = $ref_den->{$chr}->[$i]     ; }
			if(defined $ref_den_opp->{$chr}->[$i] ) {$oppDMC_num = $ref_den_opp->{$chr}->[$i] ; }
			if(defined $ref_den_both->{$chr}->[$i]) {$mC_num	 = $ref_den_both->{$chr}->[$i]; }
			
			my $netDMC_num = $posDMC_num - $oppDMC_num;
			my $unit_percentage = $percentage_cutoff_sub;
			
			if($mC_num > 0) {
				$unit_percentage = sprintf("%.2f", $netDMC_num / $mC_num * 100);
			}
					
			if ( $netDMC_num >= $netDMC_cutoff_sub and $unit_percentage >= $percentage_cutoff_sub){
				
				my $win_start = $i * $sliding_bp_sub + 1;
				my $win_end   = $win_start + $window_size_sub - 1;
				my $narrow_donw_win_start = $win_start;
				my $narrow_donw_win_end   = $win_end;
				
				#$narrow_donw_win_start and end
				for my $k($win_start..$win_end){
					if (defined $ref_pos->{$chr}->{$k}){
						$narrow_donw_win_start = $k;
						last;
					}
				}
				for(my $k = $win_end; $k >= $win_start; $k--){
					if (defined $ref_pos->{$chr}->{$k}){
						$narrow_donw_win_end = $k;
						last;
					}
				}
				
				if($chr ne $last_chr || $narrow_donw_win_start > $last_end + $merged_gap_cutoff_sub + 1){ # not overlap
					my $num = 0;
					my $num_opp = 0;
					my $num_mC  = 0;
					
			
					
					if($last_start == 0){
						
												
						($last_chr, $last_start, $last_end) = ($chr, $narrow_donw_win_start, $narrow_donw_win_end);
						
						if($debug_sub){
							print STDERR "last=0: ", join("\t", ($last_chr, $last_start, $last_end )) , "\n";
						}
						
						next;
					}
				
					foreach my $k ($last_start..$last_end){
						if(defined $ref_pos->{$last_chr}->{$k}      ) { $num++;    }
						if(defined $ref_pos_opp->{$last_chr}->{$k}  ) { $num_opp++;}
						if(defined $ref_pos_both->{$last_chr}->{$k} ) { $num_mC++; }
					}

					my $final_percentage = $percentage_cutoff_sub;
					
					
					my $final_netDMC_num = $num - $num_opp;
					if($num_mC > 0){
						$final_percentage = sprintf("%.2f", $final_netDMC_num / $num_mC * 100);
					}
				
				
					if($final_netDMC_num >= $netDMC_cutoff_sub and $final_percentage >=  $percentage_cutoff_sub){
							push (@{$ref_line}, [$last_chr, $last_start, $last_end, $num, $num_opp, $final_netDMC_num, $num_mC, $final_percentage]);
						#	($last_chr, $last_start, $last_end) = ($chr, $real_start, $real_end);

					}else{
						#print STDERR "final: $last_chr, $last_start, $last_end, $num, $num_opp, $final_netDMC_num, $num_mC, $final_percentage", "\n";
					}
				
					#($last_chr, $last_start, $last_end) = ($chr, $real_start, $real_end);
					#($last_chr, $last_start, $last_end) = ($chr, $win_start, $win_end);
					($last_chr, $last_start, $last_end) = ($chr, $narrow_donw_win_start, $narrow_donw_win_end);
				}
				
				else{
					#$last_end = ($i+1) * $sliding_bp_sub  + $window_size_sub;
					$last_end = $narrow_donw_win_end;
					
				}
			}#if ( $temp_den >= $cutoff_for_small_region )
		}
	}
}

sub cal_density{#from 0, $window_size
	my ($ref_file, $ref_den, $window_size_sub, $sliding_bp_sub) = @_;
	
	my $covered_num = int ($window_size_sub / $sliding_bp_sub);
	
	#print STDERR "cal_density...\t";
#	foreach my $chr (sort keys %{$ref_file} ){
#		foreach my $pos (sort {$a <=> $b} keys %{ ${$ref_file}{$chr}} ){
	foreach my $chr (keys %{$ref_file} ){
		foreach my $pos ( keys %{ ${$ref_file}{$chr}} ){
		#	my $base_index = int( ($pos - 1) / $window_size_sub);
			my $base_index = int( ($pos - 1) / $sliding_bp_sub);
			
			for my $i (($base_index - $covered_num)..($base_index + $covered_num)){
				if($i >= 0){
					my $interval_start = $i * $sliding_bp_sub + 1;
					my $interval_end   = $interval_start + $window_size_sub - 1;
					if($pos >= $interval_start and $pos <= $interval_end){
						$ref_den->{$chr}->[$i]++;
					}
				}
			}
			
		#	${$ref_den}{$chr}->[$base_index]++;
		}
	}
	#print STDERR "DONE\n";
}


sub cal_met{
	my ($ref) = @_ ;
	my $num = 0;
	foreach my $chr (sort keys %{$ref}){
		$num += scalar(keys %{${$ref}{$chr}})
	}
	return $num;
}
