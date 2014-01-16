#!/usr/bin/perl -w



# v0.1
# output format is more detail
# add DMC type number and opposite DMC number

#modified from 
#/Volumes/Macintosh_HD_2/idm1_new_met_data/script/Liu_algorithm_Oct1_data_setreproduce_Liu_Algorithm_hyper_and_hypo_v1.2.2.pl
use strict;
use File::Spec;

my $debug = 0;
if($debug){
	print STDERR "debug:$debug\n";
	print STDERR "input is special,read the script\n";
}

my $window_size = 1000;
my $allowed_gap = 1000;
my $cutoff_for_small_region = 5;
my $cutoff_for_merged_region = 10;

my %chr_len = ("chr1"=>30427671, "chr2"=>19698289, "chr3"=>23459830, "chr4"=>18585056, "chr5"=>26975502);


my $usage = "$0\n <input> <outdir> <pre> <p_value> <dep>\n\n";
die $usage unless (@ARGV == 5);

my $input = shift or die "input";
my $outdir = shift or die "outdir";
my $pre = shift or die "pre";
my $p_flag = shift or die "p-value";
my $dep = shift or die "dep";
die unless (-e $input and -d $outdir);

#my ($window_size, $allowed_gap, $cutoff_for_small_region, $cutoff_for_merged_region) = @ARGV[0..3];


my $postfix = "P". $p_flag ."_reduced_boundary_dep" . $dep . "_WinSize" . $window_size . "_gap" . $allowed_gap . "_initialCutoff" . $cutoff_for_small_region . "_reportCutoff"
              . $cutoff_for_merged_region . ".txt";

my $hyper_output = File::Spec->catfile($outdir, $pre. "_hyper_" . $postfix);
my $hypo_output  = File::Spec->catfile($outdir, $pre. "_hypo_" . $postfix);

print STDERR "output:\n";
print STDERR join("\n", ($hyper_output, $hypo_output)), "\n\n";

die "output exists" if(-e  $hyper_output or -e $hypo_output);

if ($debug){
	print STDERR "\nOK\n\n";
	exit;
}
#my (%pos_rdd_hypers, %pos_rdd_hypos, %pos_ros1_hypers, %pos_ros1_hypos, %pos_idm1_hypers, %pos_idm1_hypos);
# read_file	my ($file, $ref_hyper, $ref_hypo) = @_;

my (%pos_hypers, %pos_hypos);

read_file($input, \%pos_hypers, \%pos_hypos);

	print STDERR "hyper and hypo num: ", join ("\t", ( cal_met(\%pos_hypers), cal_met(\%pos_hypos) ) ), "\n";

#my %density;
#my (%den_rdd_hypers, %den_rdd_hypos, %den_ros1_hypers, %den_ros1_hypos, %den_idm1_hypers, %den_idm1_hypos);
#my (%list_rdd_hypers, %list_rdd_hypos, %list_ros1_hypers, %list_ros1_hypos, %list_idm1_hypers, %list_idm1_hypos);

my (%den_hypers, %den_hypos);
my (%list_hypers, %list_hypos);



#cal_density(\%records, \%density);

cal_density(\%pos_hypers, \%den_hypers);
cal_density(\%pos_hypos,  \%den_hypos );


# select 
#my ($ref_den, $ref_pos, $ref_line) = @_;

my (@line_hyper, @line_hypo);


#select_region(\%pos_hypers, \%den_hypers, \@line_hyper );
#select_region(\%pos_hypos,  \%den_hypos,  \@line_hypo );

select_region2(\%pos_hypers, \%den_hypers, \@line_hyper , \%pos_hypos );
select_region2(\%pos_hypos,  \%den_hypos,  \@line_hypo,   \%pos_hypers );



#my ($num_rdd_hyper, $num_rdd_hypo, $num_ros1_hyper, $num_ros1_hypo, $num_idm1_hyper, $num_idm1_hypo);

my ($num_hyper, $num_hypo);

$num_hyper = scalar (@line_hyper);
$num_hypo  = scalar (@line_hypo);

print STDERR "hyper_DMR: $num_hyper\n";
print STDERR "hypo_DMR:  $num_hypo\n";
	
output_hyper(\@line_hyper, $hyper_output);
output_hypo(\@line_hypo,  $hypo_output);

exit;


#										WT									mut
#chr1    108     +       CHG     9       11      0.818182        1       8       0.125   0.00547749464158133
# 0		  1		 2		  3		 4		  5			6			 7		 8			9			10
#chr1    224     +       CHH     1       29      0.0344828       12      35      0.342857        0.00366488289766808
##############
#	1 read_file
#########################
sub read_file{
	my ($file, $ref_hyper, $ref_hypo) = @_;
#	print STDERR "reading $file...\t";
	open (IN, $file) or die "cannot open $file: $!";

	while(<IN>){
		chomp;
		my @a = split "\t";
		if($a[9] > $a[6]){
			${$ref_hyper}{$a[0]}->{$a[1]} = $a[3];
		}			
		else{
			${$ref_hypo}{$a[0]}->{$a[1]} = $a[3];
		}
	}
	close(IN);
}


##############
#	2 cal_met
#########################
sub cal_met{
	my ($ref) = @_;
	my $num = 0;
	foreach my $chr (sort keys %{$ref}){
		$num += scalar(keys %{${$ref}{$chr}})
	}
	return $num;
}


##############
#	3 cal_density
#########################
sub cal_density{#from 0
	my ($ref_file, $ref_den) = @_;
	
	#print STDERR "cal_density...\t";
	foreach my $chr (sort keys %{$ref_file} ){
		foreach my $pos (keys %{ ${$ref_file}{$chr}} ){
			my $base_index = int( ($pos - 1) / $window_size);
			${$ref_den}{$chr}->[$base_index]++;
		}
	}
	#print STDERR "DONE\n";
}


##############
#	4 select_region
#########################
sub select_region2{
	my ($ref_pos, $ref_den, $ref_line, $opp_pos_ref) = @_;
	my ($last_start, $last_end)  = (0, 0);
	my $last_chr = "NONE";
	
	foreach my $chr(sort keys %chr_len){
		my $len = $chr_len{$chr};
		foreach my $i(0..(int($len/$window_size)+1)){
			my $temp_den = 0;
			if(defined ${$ref_den}{$chr}->[$i] ){ $temp_den = ${$ref_den}{$chr}->[$i]; }
		
			if ( $temp_den >= $cutoff_for_small_region ){
				my ($regular_start, $regular_end) = ($i*$window_size+1, ($i+1) * $window_size);#regular_start, end from win_size
				my ($new_start, $new_end) = (-1, -1);
				
				for my $k($regular_start..$regular_end){
					if (defined ${$ref_pos}{$chr}->{$k}){
						$new_start = $k; 	last;
					}
				}				
				for(my $k = $regular_end; $k >= $regular_start; $k --){
					if (defined ${$ref_pos}{$chr}->{$k}){
						$new_end = $k; 	last;
					}
				}
				die "$regular_start, $regular_end" if ($new_start == -1 or $new_end == -1);
				if($chr ne $last_chr || $new_start > $last_end + $allowed_gap + 1){ # not overlap
					my ($positive_num, $opposite_num) = (0) x 2;
					my (%pos_nums, %oppo_nums);
					foreach my $k ($last_start..$last_end){
						if(defined ${$ref_pos}{$last_chr}->{$k}){ 
							#$num++;
							$positive_num++;
							my $type = ${$ref_pos}{$last_chr}->{$k};
							$pos_nums{$type}++;
						}
						if(defined $opp_pos_ref->{$last_chr}->{$k}){
							my $type = $opp_pos_ref->{$last_chr}->{$k};
							$opposite_num++;
							$oppo_nums{$type}++;
						}
					}
					my @types = ("CG", "CHG", "CHH");
					foreach my $t(@types){
						unless(defined $pos_nums{$t} ) { $pos_nums{$t}  = 0 }
						unless(defined $oppo_nums{$t}) { $oppo_nums{$t} = 0 }
					}
				
					if(  $positive_num >= $cutoff_for_merged_region ){
						push @{$ref_line}, [$last_chr, $last_start, $last_end, 
											$positive_num, join("+", ( $pos_nums{"CG"}, $pos_nums{"CHG"}, $pos_nums{"CHH"} )) , 
											$opposite_num, join("+", ( $oppo_nums{"CG"}, $oppo_nums{"CHG"}, $oppo_nums{"CHH"} )) 
										   ];
					}
					($last_chr, $last_start, $last_end) = ($chr, $new_start, $new_end);
				}else{
					my $flag = 0;
					for(my $k = ($i+1) * $window_size; $k >= $i*$window_size; $k --){
						if (defined ${$ref_pos}{$chr}->{$k}){
							$last_end = $k;
							$flag = 1;
							last;
						}
					}					
					die "$i" unless ($flag);
				}
			}#if ( $temp_den >= $cutoff_for_small_region )
		}#		foreach my $i(0..(int($len/$window_size)+1))
	}
#	print STDERR "last one: " , join("\t", ($last_chr, $last_start, $last_end)), "\n\n";
	my ($positive_num, $opposite_num) = (0) x 2;
	my (%pos_nums, %oppo_nums);
	foreach my $k ($last_start..$last_end){
		if(defined ${$ref_pos}{$last_chr}->{$k}){ 
			$positive_num++;
			my $type = ${$ref_pos}{$last_chr}->{$k};
			$pos_nums{$type}++;
		}
		if(defined $opp_pos_ref->{$last_chr}->{$k}){
			my $type = $opp_pos_ref->{$last_chr}->{$k};
			$opposite_num++;
			$oppo_nums{$type}++;
		}
	}
	my @types = ("CG", "CHG", "CHH");
	foreach my $t(@types){
		unless(defined $pos_nums{$t} ) { $pos_nums{$t}  = 0 }
		unless(defined $oppo_nums{$t}) { $oppo_nums{$t} = 0 }
	}
	if(  $positive_num >= $cutoff_for_merged_region ){
		push @{$ref_line}, [$last_chr, $last_start, $last_end, 
							$positive_num, join("+", ( $pos_nums{"CG"}, $pos_nums{"CHG"}, $pos_nums{"CHH"} )) , 
							$opposite_num, join("+", ( $oppo_nums{"CG"}, $oppo_nums{"CHG"}, $oppo_nums{"CHH"} )) 
						   ];
	}

}



=head
#sub select_region{
	my ($ref_pos, $ref_den, $ref_line) = @_;
	my ($last_start, $last_end)  = (0, 0);
	my $last_chr = "NONE";
	
	foreach my $chr(sort keys %chr_len){
		my $len = $chr_len{$chr};
	
		foreach my $i(0..(int($len/$window_size)+1)){
			my $temp_den = 0;
			
			if(defined ${$ref_den}{$chr}->[$i] ){
				$temp_den = ${$ref_den}{$chr}->[$i] ;
			}
		
			if ( $temp_den >= $cutoff_for_small_region ){
				
				#my ($new_start, $new_end) = ($i*$window_size+1, ($i+1) * $window_size);#regular_start, end
				
				my ($regular_start, $regular_end) = ($i*$window_size+1, ($i+1) * $window_size);#regular_start, end from win_size
			
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
			
				if($chr ne $last_chr || $new_start > $last_end + $allowed_gap + 1){ # not overlap
					my $num = 0;
				
					foreach my $k ($last_start..$last_end){
						if(defined ${$ref_pos}{$last_chr}->{$k}){ $num++; }
					}
				
					if(  $num >= $cutoff_for_merged_region ){
						push (@{$ref_line}, [$last_chr, $last_start, $last_end, $num]);
						
					}
					($last_chr, $last_start, $last_end) = ($chr, $new_start, $new_end);
				}
				
				else{
					#$last_end = ($i+1) * $window_size;
					my $flag = 0;
					for(my $k = ($i+1) * $window_size; $k >= $i*$window_size; $k --){
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
=cut


#############
#	5 last output
##########################
sub output_hyper{
	my ($ref, $output) = @_;
	die "$output exists" if (-e $output);
	
	open (OUT, ">$output") or die "cannot open $output: $!";
	print OUT join("\t", ("chr", "start", "end", "hyper_DMC_num", "hyper_DMC_type", "hypo_DMC_num", "hypo_DMC_type" )), "\n";
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
	print OUT join("\t", ("chr", "start", "end", "hypo_DMC_num", "hypo_DMC_type", "hyper_DMC_num", "hyper_DMC_type")), "\n";
	for my $ref_region(@{$ref}){
		#my ($chr, $start, $end, $num) = @{$ref_region2};
		print OUT join("\t",  @{$ref_region}), "\n";
	}
	close(OUT);	
}
