#!/usr/bin/perl -w
#method 5
# for each window ,cal scores(sum of methylation level) in the window
#

#v0.1 for small score =0 , big score >=3 is OK



use strict;
use File::Spec;
$| = 1;

my $debug_print = 1;

my $debug_sub = 1;
my $debug = 1;
if($debug){
	print STDERR "debug:$debug\n";
#	print STDERR "input is special,read the script\n";
}




my %chr_len = ("chr1"=>30427671, "chr2"=>19698289, "chr3"=>23459830, "chr4"=>18585056, "chr5"=>26975502);


#my $usage = "$0 <input> <outdir> <pre> <win> <gap> <small> <final>";
#die $usage unless (@ARGV == 7);

#my $usage = "$0 <control_pre> <input> <outdir> <pre> <winSize> <sliding_bp> <merged_gap_cutoff> <netDMC_cutoff> <percentage_cutoff>";
#die $usage unless (@ARGV == 9);

my $usage = "$0 \n <control_pre> <control_input> <input> <outdir> <pre> <winSize> <sliding_bp> <merged_gap_cutoff> <big_score_cutoff> <increase_percentage_cutoff>\n\n";
die $usage unless (@ARGV == 10);


my $control_pre		  = shift or die "control_pre";
my $colA_wig          = shift or die "control_input";
my $input  			  = shift or die "input";
my $outdir 			  = shift or die "output";
my $pre    			  = shift or die "pre";
my $window_size 	  = shift or die "window_size";
my $sliding_bp		  = shift or die "sliding_bp";
my $merged_gap_cutoff = shift or die "merged_gap_cutoff";
#my $netDMC_cutoff     = shift or die "netDMC_cutoff";
my $big_score_cutoff     = shift or die "big_score_cutoff";
my $percentage_cutoff = shift or die "percentage_cutoff";

my $zero_score_cutoff = 3;

die unless (-e $input and -d $outdir);

#my $colA_wig = "/Volumes/My_Book/20120427_ShangHai_data/call_methylation/1_isMeth2wig/colA_mC.wig";
#my $colA_wig = "/Volumes/Macintosh_HD_2/BS_data_HiSeq/C24_background/cal_methy/wig/JKZ19_C24_Luc.wig";
#my $colA_wig = "/Volumes/Macintosh_HD_2/BS_data_HiSeq/C24_background/brat_bw_new_analysis/nodupl/nodupl_pooled_wig/C24_WT_nodupl_mC.wig";
#my $colA_wig = "/Volumes/My_Book/20120702_SH_extra_BS_Seq/wig/chrM_error_rate_less/col_0_pooled_mC.wig";
die "colA_wig" unless (-e $colA_wig);



my $postfix = "_method5_win" . $window_size . "_sliding" . $sliding_bp . "_gap" . $merged_gap_cutoff . "_netDMC" . $big_score_cutoff . "_per" . $percentage_cutoff . "_zero_score_cutoff3.txt";

my $hyper_output = File::Spec->catfile($outdir, $pre. "_hyper_vs_" . $control_pre	 . $postfix);
my $hypo_output  = File::Spec->catfile($outdir, $pre. "_hypo_vs_" . $control_pre	 . $postfix);
#my $output = File::Spec->catfile($outdir, $pre. $postfix);

print STDERR "output:\n";
print STDERR join("\n", ($hyper_output, $hypo_output)), "\n\n";

die "output exists" if(-e  $hyper_output or -e $hypo_output);
#die "output exists" if(-e  $output);

#open (OUT, ">$output") or die;


#my (%pos_hypers, %pos_hypos);
#my %pos_both;


if($debug){
	print STDERR "read_wig...\n";
}

#read_wig($colA_wig, \%pos_hypos); #mC in colA only is hypo DMC
#read_mut_wig($input, \%pos_hypers, \%pos_hypos);
#read_mut_wig($input, \%pos_hypers, \%pos_hypos, \%pos_both);
my (%pos_wt, %pos_mut);

read_wig5($colA_wig, \%pos_wt);
read_wig5($input, \%pos_mut);
	
	
#cal_density(\%records, \%density);
#cal_density{#from 0, $window_size

#my (%den_hypers, %den_hypos);
#my %den_both;

if($debug){
	print STDERR "cal_density...\n";
}
#cal_density(\%pos_hypers, \%den_hypers, $window_size, $sliding_bp );
#cal_density(\%pos_hypos,  \%den_hypos,  $window_size, $sliding_bp );
#cal_density(\%pos_both,   \%den_both,   $window_size, $sliding_bp );

my (%den_wt, %den_mut);

cal_density5( \%pos_wt,   \%den_wt,   $window_size, $sliding_bp );
cal_density5( \%pos_mut,  \%den_mut,  $window_size, $sliding_bp );




my (@line_hyper, @line_hypo);


if($debug){
	print STDERR "select_region...\n";
}

#hyper in mut
select_region5( $window_size, $sliding_bp, $merged_gap_cutoff, $big_score_cutoff, $percentage_cutoff, \@line_hyper,  \%pos_mut,  \%den_mut,  \%pos_wt,   \%den_wt) ;
#hypo in wt
select_region5( $window_size, $sliding_bp, $merged_gap_cutoff, $big_score_cutoff, $percentage_cutoff, \@line_hypo,   \%pos_wt,   \%den_wt,   \%pos_mut,  \%den_mut);			  

my ($num_hyper, $num_hypo);

$num_hyper = scalar (@line_hyper);
$num_hypo  = scalar (@line_hypo);

print STDERR "hyper_DMR: $num_hyper\n";
print STDERR "hypo_DMR:  $num_hypo\n";

if($debug){
	print STDERR "output...\n";
}

output_hyper5 (\@line_hyper, $hyper_output);
output_hypo5  (\@line_hypo,  $hypo_output);
exit;

#read_wig5($colA_wig, \%pos_wt);
#read_wig5($input, \%pos_mut);
sub read_wig5{
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
#	print STDERR "DONE\n";
}

#read_wig($colA_wig, \%pos_hypos); #mC in colA only is hypo DMC
#read_mut_wig($input, \%pos_hypers, \%pos_hypos);

sub output_hyper5{
	my ($ref, $output) = @_;
	die "$output exists" if (-e $output);
	
	open (OUT, ">$output") or die "cannot open $output: $!";
	#print OUT join("\t", ("chr", "start", "end", "hyper_mC", "hypo_mC", "netmC", "basic_mC", "percentage")), "\n";
	print OUT join("\t", ("chr", "start", "end", "mut_score", "wt_score", "increased_percentage")), "\n";
	for my $ref_region(@{$ref}){
		#my ($chr, $start, $end, $num) = @{$ref_region2};
		print OUT join("\t",  @{$ref_region}), "\n";
	}
	close(OUT);	
}

sub output_hypo5{
	my ($ref, $output) = @_;
	die "$output exists" if (-e $output);
	
	open (OUT, ">$output") or die "cannot open $output: $!";
	#print OUT join("\t", ("chr", "start", "end", "hypo_mC", "hyper_mC","netmC", "basic_mC", "percentage")), "\n";
	print OUT join("\t", ("chr", "start", "end", "wt_score", "mut_score", "increased_percentage")), "\n";

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


#hyper in mut
#select_region5( $window_size, $sliding_bp, $merged_gap_cutoff, $big_score_cutoff, $percentage_cutoff, \@line_hyper,  \%pos_mut,  \%den_mut,  \%pos_wt,   \%den_wt) ;
#hypo in mut
#select_region5( $window_size, $sliding_bp, $merged_gap_cutoff, $big_score_cutoff, $percentage_cutoff, \@line_hypo,   \%pos_wt,   \%den_wt,   \%pos_mut,  \%den_mut);	

sub select_region5{
	my ( $window_size_sub, $sliding_bp_sub, $merged_gap_cutoff_sub, $big_score_cutoff_sub, $percentage_cutoff_sub, 
		 $ref_line, $ref_pos_hyper, $ref_den_hyper, $ref_pos_hypo, $ref_den_hypo) = @_;
	
	my $hypo_zero_label = -99999999;
	
	my ($last_chr, $last_start, $last_end)  = ("chr0", 0, 0);
	
	foreach my $chr(sort keys %chr_len){
	#foreach my $chr (@chr_array)
		my $len = $chr_len{$chr};
	
		foreach my $i(0..(int( ($len - 1) / $sliding_bp_sub ) ) ){
			my ($score_hyper, $score_hypo) = (0) x 2;
			
			if(defined $ref_den_hyper->{$chr}->[$i]) { $score_hyper  = $ref_den_hyper ->{$chr}->[$i]}
			if(defined $ref_den_hypo ->{$chr}->[$i]) { $score_hypo   = $ref_den_hypo ->{$chr}->[$i]}
			
		#	my $net_score = $score_hyper - $score_hypo;
			#my $unit_percentage = $percentage_cutoff_sub;
			my $unit_percentage = $hypo_zero_label ;
			
			if($score_hypo > 0) {
				$unit_percentage = sprintf("%.2f", 100 * ( $score_hyper - $score_hypo  ) / $score_hypo );
			}
			
			if( $debug_print == 1 ){
				if( $chr eq "chr5" and ( $i >= 170215 and $i <= 170227) ){
					print STDERR join("\t", ($i * $sliding_bp_sub + 1,$i * $sliding_bp_sub + $window_size_sub  , $score_hyper, $score_hypo )), "\n";
				}
			}
			
			if ( ($score_hyper >= $big_score_cutoff_sub and $unit_percentage >= $percentage_cutoff_sub) or 
				 ($unit_percentage == $hypo_zero_label and $score_hyper >= $zero_score_cutoff)){
				my $win_start = $i * $sliding_bp_sub + 1;
				my $win_end   = $win_start + $window_size_sub - 1;
				my $narrow_donw_win_start = $win_start;
				my $narrow_donw_win_end   = $win_end;
				#$narrow_donw_win_start and end
				for my $k($win_start..$win_end){
					if (defined $ref_pos_hyper->{$chr}->{$k}){
						$narrow_donw_win_start = $k;
						last;
					}
				}#for
				for(my $k = $win_end; $k >= $win_start; $k--){
					if (defined $ref_pos_hyper->{$chr}->{$k}){
						$narrow_donw_win_end = $k;
						last;
					}
				}#for
				
				
				
				
				if($chr ne $last_chr || $narrow_donw_win_start > $last_end + $merged_gap_cutoff_sub + 1){ # not overlap
					my $score_hyper_total = 0;
					my $score_hypo_total = 0;
					if($last_start == 0){
						($last_chr, $last_start, $last_end) = ($chr, $narrow_donw_win_start, $narrow_donw_win_end);
						if($debug_sub){
							print STDERR "last=0: ", join("\t", ($last_chr, $last_start, $last_end )) , "\n";
						}
						next;
					}

					foreach my $k ($last_start..$last_end){
						if(defined $ref_pos_hyper->{$last_chr}->{$k}    ) { $score_hyper_total += $ref_pos_hyper->{$last_chr}->{$k};    }
						if(defined $ref_pos_hypo->{$last_chr}->{$k}     ) { $score_hypo_total  += $ref_pos_hypo->{$last_chr}->{$k} ;}
					}
					
					#my $final_percentage = $percentage_cutoff_sub;
					my $final_percentage = $hypo_zero_label;
					if($score_hypo_total > 0 ){
						$final_percentage = sprintf("%.2f", 100 * ($score_hyper_total - $score_hypo_total) / $score_hypo_total);
					}
					
					if( ( $score_hyper_total >= $big_score_cutoff_sub and $final_percentage >=  $percentage_cutoff_sub) 
							or ( $final_percentage == $hypo_zero_label and  $score_hyper_total >= $zero_score_cutoff) ){
								
							if($final_percentage == $hypo_zero_label) {$final_percentage = 10000}
							push (@{$ref_line}, [$last_chr, $last_start, $last_end, $score_hyper_total, $score_hypo_total, $final_percentage]);
						#	($last_chr, $last_start, $last_end) = ($chr, $real_start, $real_end);

					}				
					($last_chr, $last_start, $last_end) = ($chr, $narrow_donw_win_start, $narrow_donw_win_end);
				
				}#if ( $score_hyper >= $big_score_cutoff_sub and $unit_percentage >= $percentage_cutoff_sub)
			
				else{
					#$last_end = ($i+1) * $sliding_bp_sub  + $window_size_sub;
					$last_end = $narrow_donw_win_end;
					
				}
				
			}
		}#foreach my $i(0..(int( ($len - 1) / $sliding_bp_sub ) ) )
	}
}

#cal_density5(\%pos_wt,   \%den_wt,   $window_size, $sliding_bp );
#cal_density5(\%pos_mut,  \%den_mut,  $window_size, $sliding_bp );

sub cal_density5{
	my ($ref_pos, $ref_den, $window_size_sub, $sliding_bp_sub) = @_;
	my $covered_num = int ($window_size_sub / $sliding_bp_sub);
	
	foreach my $chr (keys %{$ref_pos} ){
		foreach my $pos ( keys %{ ${$ref_pos}{$chr}} ){
			my $base_index = int( ($pos - 1) / $sliding_bp_sub);
			for my $i (($base_index - $covered_num)..($base_index + $covered_num)){
				if($i >= 0){
					my $interval_start = $i * $sliding_bp_sub + 1;
					my $interval_end   = $interval_start + $window_size_sub - 1;
					if($pos >= $interval_start and $pos <= $interval_end){
						$ref_den->{$chr}->[$i] += ($ref_pos->{$chr}->{$pos});
					}
				}
			}
		}
	}
}

sub cal_met{
	my ($ref) = @_ ;
	my $num = 0;
	foreach my $chr (sort keys %{$ref}){
		$num += scalar(keys %{${$ref}{$chr}})
	}
	return $num;
}