#!/usr/bin/perl -w

#meth_level_patterns_Zilberman_single_isMeth_original_mC_NotEckerSuggest_v0.0.pl
# modified from
# /Users/tang58/Kai_BS/for_publish/9_beCalled_cal_methylation_level_for_features_and_flanking_regions_for_SH_isMeth_file_format_AllDepth_cutoff_v0.1_original_C_not_EckerSuggest.pl

#Gene and TE meta analysis (ends analysis). A. thaliana TAIR-annotated genes or transposons 
#were aligned at the 5' end or the 3' end. For genes and TEs, we discarded from the analysis 1500 
#bp or 250 bp, respectively, from the end opposite to the one used for alignment to avoid 
#averaging the edges of shorter genes and TEs with the middles of longer sequences.

#one cal_pattern function is enough, just need two input list, one is for five_p, the other is for three_p
# and when read input for three_p, just reverse the strand
use strict;
use File::Spec;

my $dep_cutoff = 4;
my $bin_size = 100;
my $bp_outside_gene_body = 4000;
my $bp_inside_gene_body = 4000;
my %discarded_bp = (
                    "Gene"      => 1500,
                    "TE"        => 250
                    );

print STDERR "\n\n\n Please change varialbe detail's value !!!\n\n";

my $debug = 0;
if($debug){
	print STDERR "debug = 1\n\n";
}

my $detail = 0;

if($debug == 0){
	$detail = 0;
}

#my  $usage = "$0 \n <isMeth_file> <sample_label> <feature_bin_num(20)> <flaking_bp(2000)> <flaking_bin_num(20)> <input_bed_directory> <outdir>  <postfix(XXXPaper)> \n\n";
#die $usage unless (@ARGV == 8);

my  $usage = "$0 \n <isMeth_file> <sample_label> <feature(TE/Gene)> <input_bed_directory> <outdir>  <postfix(XXXPaper)> \n\n";
die $usage unless (@ARGV == 6);

my $isMeth_file		= shift or die "isMeth_file";
my $sample_label	= shift or die;
my $feature 		= shift or die;

die unless (defined $discarded_bp{$feature});
#die $usage, "\n\nWT\n\n" unless ($WT_label eq "WT");

#my $feature_bin_num = shift or die "feature_bin_num";
#my $flaking_bp = shift or die "flaking_bp";
#my $flaking_bin_num = shift or die "flaking_bin_num";

my $indir = shift or die "indir";
my $outdir   = shift or die "outdir";

my $postfix = shift or die "postfix";


die unless (-d $indir);
die unless (-d $outdir);
die unless (-e $isMeth_file );


opendir(DIR, $indir) or die "cannot open indir";
my @bed_lists = grep /txt$/, readdir DIR;
close DIR;

print STDERR "\n", join("\n", @bed_lists), "\n\n";

my @outputs;
my @feature_pres;

my $i = -1;

foreach my $file( @bed_lists ){
	$i++;
	if($file =~ /(\S+)_coordinate_bed.txt$/){
		$feature_pres[$i] = $1;
	#	$outputs[$i] = File::Spec->catfile($outdir,
	#		$feature_pres[$i] . "_in_" . $sample_label . "_binNum" . $feature_bin_num . "_flaking" . $flaking_bp . "bp_binNum" . $flaking_bin_num . "_" . $postfix . "_notEckerSuggest.txt");
	$outputs[$i] = File::Spec->catfile($outdir,
					   $feature_pres[$i] . "_in_" . $sample_label . "_Zilberman_method_meth_pattern_" . $postfix . "_original_mC_notEckerSuggest.txt");
		die $outputs[$i]  if (-e $outputs[$i] ) ;
	}
	else{
		die $file;
	}
}
print STDERR join("\n", @outputs), "\n\n";
print STDERR join("\n", @feature_pres), "\n\n";

die "wrong number " unless (@feature_pres == @bed_lists and @outputs == @bed_lists);


if ($debug){
	
	print join("\n", @feature_pres), "\n\n";
	print join("\n", @bed_lists), "\n\n";
	print join("\n", @outputs), "\n\n";
	
	print STDERR "\n\nOK\n\n";
	exit;
}

#my (@target_list , @left_list);
my @input_list;# each member is also a array record the gene/TE 



my $last_index = $#bed_lists; #file number - 1
for $i (0..$last_index){
	my $input = File::Spec->catfile($indir, $bed_lists[$i]);
	die unless (-e $input);
	print STDERR "reading list ", $bed_lists[$i] , "\t";
	read_list2($input, \@{$input_list[$i]}  );
	print STDERR "DONE\n";
}


#my (@total_depths, @total_mCs) ;# $array[i]->{CG/CHG/CHH}->[0..(2*flaking_bin_num + feature_bin_num -1 )] += 
#chr     pos     strand  type    num_C   depth   percentage      isMeth
#chr1    109     +       CG      10      11      0.909091        1
#chr1    114     +       CHG     2       12      0.166667        1
my (@total_depths_5p, @total_mCs_5p) ;# $array[i]->{CG/CHG/CHH}->[0..(2*flaking_bin_num + feature_bin_num -1 )] += 
my (@total_depths_3p, @total_mCs_3p) ;# $array[i]->{CG/CHG/CHH}->[0..(2*flaking_bin_num + feature_bin_num -1 )] += 



for my $i_chr (1..5){
	print STDERR  "reading chr$i_chr\n";
	my (@positions, @depths_temp, @mCs_temp);
	open(DB, $isMeth_file) or die "db";
	my $curr_chr = "chr" . $i_chr;
	
#	my $line = 0;
#	my $const = 3000000;
	my $db_head = <DB>;
	die unless ($db_head =~ /isMeth/);
	
	while(<DB>){

		next unless(/$curr_chr/i);
		my @a = split "\t";
		next unless ( $a[5] >= $dep_cutoff );
	#	die unless (lc($a[0]) eq $curr_chr);
		die unless ( $a[0] eq $curr_chr);
		
		my $pos = $a[1];
		my $type = $a[3];
		$positions[$pos]   = $type;
		$depths_temp[$pos] = $a[5];
		$mCs_temp[$pos] = $a[4];
	}
	close(DB);
	print STDERR "\n\n";

	
	for $i (0..$last_index){
		print STDERR "cal $i \n";
#my $bin_size = 100;
#my $bp_outside_gene_body = 4000;
#my $bp_inside_gene_body = 4000;
#$discarded_bp{$feature}
		
		#one cal_pattern function is enough, just need two input list, one is for five_p, the other is for three_p
		# and when read input for three_p, just reverse the strand (+/- => -/+)
		cal_pattern( $curr_chr,  \@{$input_list[$i]}, \@positions, \@depths_temp, \@mCs_temp,  \%{$total_depths_5p[$i]}, \%{$total_mCs_5p[$i]},
			    $discarded_bp{$feature}, $bin_size, $bp_outside_gene_body, $bp_inside_gene_body, "five_prime"  );
		
		cal_pattern( $curr_chr,  \@{$input_list[$i]}, \@positions, \@depths_temp, \@mCs_temp,  \%{$total_depths_3p[$i]}, \%{$total_mCs_3p[$i]},
			    $discarded_bp{$feature}, $bin_size, $bp_outside_gene_body, $bp_inside_gene_body, "three_prime" );
		
	}
	
}

my $outside_gene_interval_number =  int ( ( $bp_outside_gene_body - 1 )/ $bin_size) + 1;
my $inside_gene_interval_number  =  int ( ($bp_inside_gene_body - 1) /$bin_size  ) + 1;


for $i (0..$last_index){
	my $output = $outputs[$i];
#	output2($output, \%{$total_depths[$i]}, \%{$total_mCs[$i]}, $total_intervel_num);
	output_Zilberman( $output, \%{$total_depths_5p[$i]}, \%{$total_mCs_5p[$i]},\%{$total_depths_3p[$i]}, \%{$total_mCs_3p[$i]}, $outside_gene_interval_number, $inside_gene_interval_number);
}

print STDERR "\nOK\n\n";

exit;

# cal_pattern( $curr_chr,  \@{$input_list[$i]}, \@positions, \@depths_temp, \@mCs_temp,  \%{$total_depths[$i]}, \%{$total_mCs[$i]},
#			    $discarded_bp{$feature}, $bin_size, $bp_outside_gene_body, $bp_inside_gene_body  );

sub cal_pattern{
	#my ($chr_sub,         $list_refa, $pos_refa,   $depths_temp_refa, $mCs_temp_refa, $total_depth_refh,  $total_mC_refh, $feature_bin_num_sub, $flaking_bp_sub, $flaking_bin_num_sub )  = @_;

	my ($chr_sub, $list_refa, $pos_refa, $depths_temp_refa, $mCs_temp_refa, $total_depth_refh, $total_mC_refh,
	    $discarded_bp_sub, $bin_size_sub, $bp_outside_gene_body_sub, $bp_inside_gene_body_sub,
	    $prime_end
	    ) = @_;
	die unless ( $prime_end =~ /five_prime|three_prime/  );
	
	my %reverse_strand = ( "+" => "-", "-" => "+" );
#	|======>
#-N..-1 | 1,2,3..
# <======|
#N..3,2,1|-1,-2..
#	|

#start always > end
	foreach my $i(0..( scalar(@{$list_refa}) -1 )){
		my $this = $list_refa->[$i];
		next unless ($this =~ /$chr_sub/);
	#	my ($chr, $start, $end, $length, $strand, $id) = split "_", $this;
	
		my ($chr, $start, $end, $strand) = split "_", $this;
		die unless($chr eq $chr_sub);
		if ($prime_end eq "three_prime"){
			$strand = $reverse_strand{$strand};
		}
		
		if( $strand eq "+"){
		#######
		# outside
			{
			my $outside_start = (  ( $start - $bp_outside_gene_body_sub) >= 1  )?($start - $bp_outside_gene_body_sub):(1);
			my $outside_interval_num = int (  ($start - $outside_start - 1)/$bin_size_sub ) + 1;
			
			if ($detail) { print STDERR join("\t", ("+", "outside_start", $outside_start, "outside_interval_num",$outside_interval_num)) ,"\n" }
			
			for my $j (1..($outside_interval_num -1) ){
				my $interval_start =  $start - $j * $bin_size_sub ;
				my $interval_end  = $interval_start + $bin_size_sub - 1 ;
				if ($detail) { print STDERR join("\t", ("+", $j, "interval_start", $interval_start, "interval_end",$interval_end)) ,"\n" }
				for my $k ($interval_start..$interval_end){
					#$ref->{-$j} ->{"mC"} = ;
					if(defined $pos_refa->[$k]) {
						my $type = $pos_refa->[$k];
						$total_depth_refh->{$type}->{-$j} +=	$depths_temp_refa->[$k];
						$total_mC_refh->{$type}->{-$j}	  +=	$mCs_temp_refa->[$k];
						
						$total_depth_refh->{C}->{-$j}	+=	$depths_temp_refa->[$k];
						$total_mC_refh->{C}->{-$j}	+=	$mCs_temp_refa->[$k];
						
					}
				}
				
			}
			#################
			#last one interval
			my $j = $outside_interval_num;
			my $interval_start = $outside_start;
			my $interval_end   =$start - ( $j - 1 ) * $bin_size_sub - 1;
			if ($detail) { print STDERR join("\t", ("+last_out", $j, "interval_start", $interval_start, "interval_end",$interval_end)) ,"\n" }
			for my $k ($interval_start..$interval_end){
				#$ref->{-$j} ->{"mC"} = ;
				if(defined $pos_refa->[$k]) {
					my $type = $pos_refa->[$k];
					$total_depth_refh->{$type}->{-$j} +=	$depths_temp_refa->[$k];
					$total_mC_refh->{$type}->{-$j}	  +=	$mCs_temp_refa->[$k];
					
					$total_depth_refh->{C}->{-$j}	+=	$depths_temp_refa->[$k];
					$total_mC_refh->{C}->{-$j}	+=	$mCs_temp_refa->[$k];
				}
			}
			#end of last one interval
			############
			}
		# end of outside
		############
		
		###########
		#	inside
			{
				my $inside_end = ( ($start + $bp_inside_gene_body_sub) <= ($end - $discarded_bp_sub) )?( $start + $bp_inside_gene_body_sub ):($end - $discarded_bp_sub);
				my $inside_interval_num  = int ( ( $inside_end - $start - 1) /$bin_size_sub  ) + 1;
				if ($detail) { print STDERR join("\t", ("+", "inside_end", $inside_end, "inside_interval_num",$inside_interval_num)) ,"\n" }
				for my $j (1..($inside_interval_num -1) ){
					my $interval_start = $start + ($j - 1) * $bin_size_sub;
					my $interval_end   = $interval_start + $bin_size_sub - 1 ;
					if ($detail) { print STDERR join("\t", ("+in", $j, "interval_start", $interval_start, "interval_end",$interval_end)) ,"\n" }
					for my $k ($interval_start..$interval_end){
						if(defined $pos_refa->[$k]) {
							my $type = $pos_refa->[$k];
							$total_depth_refh->{$type}->{$j} +=	$depths_temp_refa->[$k];
							$total_mC_refh->{$type}->{$j}	  +=	$mCs_temp_refa->[$k];
					
							$total_depth_refh->{C}->{$j}	+=	$depths_temp_refa->[$k];
							$total_mC_refh->{C}->{$j}	+=	$mCs_temp_refa->[$k];
						}	
					}
				}
					
				#################
				#last one interval
				my $j = $inside_interval_num;
				my $interval_start = $start + ($j - 1) * $bin_size_sub;
				my $interval_end = $inside_end;
				if ($detail) { print STDERR join("\t", ("+last_in", $j, "interval_start", $interval_start, "interval_end",$interval_end)) ,"\n" }

				for my $k ($interval_start..$interval_end){
					if(defined $pos_refa->[$k]) {
						my $type = $pos_refa->[$k];
						$total_depth_refh->{$type}->{$j} +=	$depths_temp_refa->[$k];
						$total_mC_refh->{$type}->{$j}	  +=	$mCs_temp_refa->[$k];
				
						$total_depth_refh->{C}->{$j}	+=	$depths_temp_refa->[$k];
						$total_mC_refh->{C}->{$j}	+=	$mCs_temp_refa->[$k];
					}	
				}
				#end of last one interval
				############
			}
		#end of inside
		#########
		
		}

# inside |  outside	
# <======|
#N..3,2,1|-1,-2..
#	|
		
		elsif( $strand eq "-"){
		#######
		# outside
			{
				my $outside_end = $end + $bp_outside_gene_body_sub;
				my $outside_interval_num = int ( ($outside_end - $end - 1 ) / $bin_size_sub ) + 1;
				if ($detail) { print STDERR join("\t", ("-", "outside_end", $outside_end, "outside_interval_num",$outside_interval_num)) ,"\n" }
				for my $j (1..($outside_interval_num -1) ){
					my $interval_start = $end + 1 + ($j - 1) * $bin_size_sub;
					my $interval_end   = $interval_start + $bin_size_sub - 1;
				
					if ($detail) { print STDERR join("\t", ("-", $j, "interval_start", $interval_start, "interval_end",$interval_end)) ,"\n" }

					
					for my $k ($interval_start..$interval_end){
						#$ref->{-$j} ->{"mC"} = ;
						if(defined $pos_refa->[$k]) {
							my $type = $pos_refa->[$k];
							$total_depth_refh->{$type}->{-$j} +=	$depths_temp_refa->[$k];
							$total_mC_refh->{$type}->{-$j}	  +=	$mCs_temp_refa->[$k];
						
							$total_depth_refh->{C}->{-$j}	+=	$depths_temp_refa->[$k];
							$total_mC_refh->{C}->{-$j}	+=	$mCs_temp_refa->[$k];
						}
					}				
				}
					
				#################
				#last one interval
				my $j 		   = $outside_interval_num;
				my $interval_start = $end + 1 + ($j - 1) * $bin_size_sub;
				my $interval_end   = $outside_end;
				if ($detail) { print STDERR join("\t", ("-last_out", $j, "interval_start", $interval_start, "interval_end",$interval_end)) ,"\n" }

				for my $k ($interval_start..$interval_end){
					#$ref->{-$j} ->{"mC"} = ;
					if(defined $pos_refa->[$k]) {
						my $type = $pos_refa->[$k];
						$total_depth_refh->{$type}->{-$j} +=	$depths_temp_refa->[$k];
						$total_mC_refh->{$type}->{-$j}	  +=	$mCs_temp_refa->[$k];
					
						$total_depth_refh->{C}->{-$j}	+=	$depths_temp_refa->[$k];
						$total_mC_refh->{C}->{-$j}	+=	$mCs_temp_refa->[$k];
					}
				}
				# end of outside last one interval
				############
			}
		# end of outside
		############
		
		###########
		#	inside
			{
				my $inside_start = ( ($start + $discarded_bp_sub) > ($end - $bp_inside_gene_body_sub))?( $start + $discarded_bp_sub):( $end - $bp_inside_gene_body_sub);
				my $inside_interval_num  = int ( ( $end - $inside_start - 1) / $bin_size_sub ) + 1;
				if ($detail) { print STDERR join("\t", ("-", "inside_start", $inside_start, "inside_interval_num",$inside_interval_num)) ,"\n" }

				for my $j (1..($inside_interval_num -1) ){
					my $interval_end   = $end - ($j - 1) *  $bin_size_sub;
					my $interval_start = $interval_end - $bin_size_sub + 1;
					if ($detail) { print STDERR join("\t", ("-in", $j, "interval_start", $interval_start, "interval_end",$interval_end)) ,"\n" }
					
					for my $k ($interval_start..$interval_end){
						if(defined $pos_refa->[$k]) {
							my $type = $pos_refa->[$k];
							$total_depth_refh->{$type}->{$j} +=	$depths_temp_refa->[$k];
							$total_mC_refh->{$type}->{$j}	  +=	$mCs_temp_refa->[$k];
				
							$total_depth_refh->{C}->{$j}	+=	$depths_temp_refa->[$k];
							$total_mC_refh->{C}->{$j}	+=	$mCs_temp_refa->[$k];
						}	
					}
				}
				
				my $j = $inside_interval_num;
				my $interval_start = $inside_start;
				my $interval_end  = $end - ($j - 1) *  $bin_size_sub;
				if ($detail) { print STDERR join("\t", ("-last_in", $j, "interval_start", $interval_start, "interval_end",$interval_end)) ,"\n" }
				
				for my $k ($interval_start..$interval_end){
					if(defined $pos_refa->[$k]) {
						my $type = $pos_refa->[$k];
						$total_depth_refh->{$type}->{$j} +=	$depths_temp_refa->[$k];
						$total_mC_refh->{$type}->{$j}	 +=	$mCs_temp_refa->[$k];
						$total_depth_refh->{C}->{$j}	 +=	$depths_temp_refa->[$k];
						$total_mC_refh->{C}->{$j}	 +=	$mCs_temp_refa->[$k];
					}	
				}
			}
		#end of inside
		#########
		
		}
		else{
			die "strand: ", $this;
		}
	
	}
	
}


#	output_Zilberman( $output, \%{$total_depths_5p[$i]}, \%{$total_mCs_5p[$i]},\%{$total_depths_3p[$i]}, \%{$total_mCs_3p[$i]}, $outside_gene_interval_number, $inside_gene_interval_number);
sub output_Zilberman{
	my ($output_sub, $total_depths_5p_sub, $total_mCs_5p_sub, $total_depths_3p_sub, $total_mCs_3p_sub, $outside_gene_interval_number_sub, $inside_gene_interval_number_sub ) = @_;
	die if(-e $output_sub);

	
	my @types_3  = ("CG", "CHG", "CHH");
	my @types_4  = ("CG", "CHG", "CHH", "C");
	
	my @indexs = (  (-$outside_gene_interval_number_sub)..-1, 1..$inside_gene_interval_number_sub);
	
	for my $i ( @indexs ){
		for my $type (@types_4){
			unless (defined $total_depths_5p_sub->{$type}->{$i})	{ 	$total_depths_5p_sub->{$type}->{$i}= 0;}
			unless (defined $total_mCs_5p_sub->{$type}->{$i})	{ 	$total_mCs_5p_sub->{$type}->{$i}= 0;}
			unless (defined $total_depths_3p_sub->{$type}->{$i})	{ 	$total_depths_3p_sub->{$type}->{$i}= 0;}
			unless (defined $total_mCs_3p_sub->{$type}->{$i})	{ 	$total_mCs_3p_sub->{$type}->{$i}= 0;}
		}
	}
	
	check_sum( $total_depths_5p_sub, \@indexs );
	check_sum( $total_mCs_5p_sub, \@indexs );
	check_sum( $total_depths_3p_sub, \@indexs );
	check_sum( $total_mCs_3p_sub, \@indexs );
	
	
	open(OUT, ">>$output_sub") or die "cannot open $output_sub:$!";
	
#	print OUT join("\t", ("order", "mCG_num", "CG_total_num", "CG_per",
#			"mCHG_num", "CHG_total_num", "CHG_per",
#			"mCHH_num", "CHH_total_num", "CHH_per",
#			"mC_num", "C_total_num", "C_per")), "\n";

	print OUT join("\t", ("order",
			"mCG_num_5", "CG_total_num_5", "CG_per_5",
			"mCHG_num_5", "CHG_total_num_5", "CHG_per_5",
			"mCHH_num_5", "CHH_total_num_5", "CHH_per_5",
			"mC_num_5", "C_total_num_5", "C_per_5",
			
			"mCG_num_3", "CG_total_num_3", "CG_per_3",
			"mCHG_num_3", "CHG_total_num_3", "CHG_per_3",
			"mCHH_num_3", "CHH_total_num_3", "CHH_per_3",
			"mC_num_3", "C_total_num_3", "C_per_3"
			)), "\n";
	for my $i ( @indexs ){
		my (%depths_sub_5, %mC_sub_5, %pers_sub_5);
		my (%depths_sub_3, %mC_sub_3, %pers_sub_3);
		for my $type (@types_4){
			$depths_sub_5 {$type} = $total_depths_5p_sub ->{$type}->{$i};
			$mC_sub_5 {$type} = $total_mCs_5p_sub ->{$type}->{$i};
			
			$depths_sub_3 {$type} = $total_depths_3p_sub ->{$type}->{$i};
			$mC_sub_3 {$type} = $total_mCs_3p_sub ->{$type}->{$i};
			
			$pers_sub_5{$type} = 0;
			$pers_sub_3{$type} = 0;
			if ( $depths_sub_5 {$type} > 0 ) { $pers_sub_5{$type} =  sprintf ("%.4f", 100 * $mC_sub_5{$type} / $depths_sub_5 {$type}) }
			if ( $depths_sub_3 {$type} > 0 ) { $pers_sub_3{$type} =  sprintf ("%.4f", 100 * $mC_sub_3{$type} / $depths_sub_3 {$type}) }
		}
		print OUT $i, "\t";
		for my $type (@types_4){
			print OUT join("\t", ($mC_sub_5{$type}, $depths_sub_5{$type}, $pers_sub_5{$type})), "\t";
		}
		
		for my $type (@types_3){
			print OUT join("\t", ($mC_sub_3{$type}, $depths_sub_3{$type}, $pers_sub_3{$type})), "\t";
		}
		my $type = "C";
		print OUT join("\t", ($mC_sub_3{$type}, $depths_sub_3{$type}, $pers_sub_3{$type})), "\n";

	}
	close OUT;
}

#	check_sum( $total_depths_5p_sub, $outside_gene_interval_number_sub, $inside_gene_interval_number_sub );
sub check_sum{
	my($ref, $ref_a ) = @_;
	
	my @index = @{$ref_a};
	
	my @types_3  = ("CG", "CHG", "CHH");
	for my $i ( @index){
		my $sum = 0;
		for my $type (@types_3){
			$sum += $ref->{$type}->{$i};
		}
		if( $sum != $ref->{C}->{$i}){
			print STDERR join("\t", ("wrong", $i) ),"\n";
#			print STDERR join("\t", ($i)), "\n";
		}
	}
}



sub read_list2{
	my ($file, $ref) = @_;
	die unless (-e $file);
	open(IN, $file)  or die "cannot open $file:$!\n\n";
	my $i = -1;
	my $head = <IN>;
	chomp $head;
	my @temp = split "\t", $head;
	unless(lc($temp[0]) eq "chr" or $temp[0] =~/^#/){
		$i++;
		$ref->[$i] = join("_", @temp[0..3]);
	}
	while(<IN>){
		$i++;
		chomp;
		my @a = split "\t";
		$a[0] = lc($a[0]);
		$ref->[$i] = join("_", @a[0..3]);# chr	start	end strand	ID; 0..3: no ID
	}
	close(IN);
	
	print STDERR "read_list2: $file index(file_line_nunber - 2): $i \n\n";
}
