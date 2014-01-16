#!/usr/bin/perl -w

#cal_meth_level_for_each_item_bed_list_EckerPaper_v0.0
# main purpose of this script is calculate methylation level for
# each of the gene in the genome;
# to generate scatterplot and corelation r

 
use strict;
use File::Spec;

my $debug = 0;

my $depth_cutoff = 4;

#my $usage = "$0 \n <isMeht_file_wt> <isMeht_file_mut>  <bed_like_file> <hyper|hypo> <outdir> \n\n";
#die $usage unless (@ARGV == 5);

my $usage = "$0 \n <isMeht_file> <sample_label> <bed_like_file> <output> \n\n";
die $usage unless (@ARGV == 4);

my $isMeth_file 	= shift or die;
my $sample_label 	= shift or die;

#my $isMeth_file_WT	= shift or die;
#my $isMeth_file_mut	= shift or die;
my $bed_file		= shift or die;
#my $hyper_hypo_label	= shift or die;
#my $outdir		= shift or die;
# my $depth_cutoff	= shift or die;

#my $output = shift or die;

my $output = shift or die;
#die $usage , "\n", "hyper|hypo" unless ( $hyper_hypo_label =~ /hyper|hypo/);

die unless (-e $bed_file);
die unless (-e $isMeth_file);

die if (-e $output);

if($debug){
	print STDERR "\n\n";
#	print STDERR join("\n", ($isMeth_file_WT, $isMeth_file_mut, $bed_file, $output )), "\n";
	print STDERR join("\n", ($isMeth_file, $bed_file, $output )), "\n";
	print STDERR "\n\nOK\n\n";
	exit;
}

my @bed_list;

my $last_one_index = read_bed_file($bed_file, \@bed_list);

print STDERR "last_one_index(head = 0): line_num: ", $last_one_index, "\n\n";


my ( %C_type_nums,  %called_mC_nums, %meth_level_sum,  %seq_mC_sum,  %seq_depth_sum ); #Weighted methylation level

# %{wt/mut}->{1/2}->{CG/CHG/CHH/total}

Initialization_type_num(\%C_type_nums, $last_one_index);

Initialization(\%called_mC_nums, $last_one_index);
Initialization(\%meth_level_sum, $last_one_index);
Initialization(\%seq_mC_sum, $last_one_index);
Initialization(\%seq_depth_sum, $last_one_index);


for my $i_chr (1..5){
	print STDERR  "reading chr$i_chr\n";
	
	my ( @positions, @depths_temp, @mCs_temp );
	my ( @isMeth,  @per  );
	
	open(DB, $isMeth_file) or die "db";

	my $curr_chr = "chr" . $i_chr;
	
	#my $line = 0;
	#my $const = 3000000;
	
	my $db_head = <DB>;
	die "db_head" unless ($db_head =~ /isMeth/);
	
	while(<DB>){
	#	$line++;
	#	if($line % $const == 0){
			#print STDERR $line, "\t";
	#	}
		next unless(/$curr_chr/i);
		my @a = split "\t";
		next unless ( $a[5] >= $depth_cutoff );
		die unless (lc($a[0]) eq $curr_chr);
		
		my $pos = $a[1];
		my $type = $a[3];
		$positions[$pos]   = $type;
		$depths_temp[$pos] = $a[5];
	#chr     pos     strand  type    num_C   depth   percentage      isMeth
	#chr1    647     +       CG      17      18      0.944444        1
#	0	 1	2	 3	 4	5	6		7
		if( $a[-1] == 1 ){
			$mCs_temp[$pos] = $a[4];
			$isMeth[$pos]   = 1;
			$per[$pos]      = $a[6];
		}elsif( $a[-1] == 0 ){
			$mCs_temp[$pos] = 0;
			$isMeth[$pos]   = 0;
			$per[$pos]      = 0;
		}else{
			die "isMeth: $_\n\n";
		}
				
	}
	close(DB);
#	print STDERR "\n\n";
	
	get_info_each_chr(
			  \@bed_list,	$curr_chr ,
			  \@positions, \@depths_temp, \@mCs_temp, \@isMeth,  \@per,
			  \%C_type_nums,  \%called_mC_nums, \%meth_level_sum,  \%seq_mC_sum,  \%seq_depth_sum,
			  $depth_cutoff	
			  
			  );
}




check_type_num(\%C_type_nums, $last_one_index );


check(\%called_mC_nums, $last_one_index );
check(\%meth_level_sum, $last_one_index );
check(\%seq_mC_sum, $last_one_index );
check(\%seq_depth_sum, $last_one_index );

output_list ($output, \@bed_list, 
			  \%C_type_nums,  \%called_mC_nums, \%meth_level_sum,  \%seq_mC_sum,  \%seq_depth_sum,
			  $sample_label
			  );


exit;

sub cal_meth{
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


sub assign{
	my ($first_ref, $ref_a) = @_;
	
	my @types_sub = ("CG", "CHG", "CHH", "total");
	my $i = -1;
	foreach my $type_sub (@types_sub){
		$i++;
		$ref_a->[$i] = $first_ref->{$type_sub};
	}
}



#done
sub check_type_num{
	my ($ref, $last_index_sub) = @_;
	#my @labs_sub = ("wt", "mut");
	my @types_sub = ("CG", "CHG", "CHH");
	
		for my $i (1..$last_index_sub){
			
			my $sum = 0;
			foreach my $type_sub (@types_sub){
				$sum  +=  $ref->{$i}->{$type_sub} ;
			}
			
			my $t = $ref->{$i}->{"total"};
			if ( abs( $sum - $t ) > 0.01 ){
				print STDERR "maybe wrong  \n", $i, "\n\n";
				print STDERR "sum = $sum != $t \n\n";
			}
			
			
		}
	
}


#modified
sub check{
	my ($ref, $last_index_sub) = @_;
	#my @labs_sub = ("wt", "mut");
	my @labs_sub = ( "mut");

	my @types_sub = ("CG", "CHG", "CHH");
	
	foreach my $lab_sub (@labs_sub){
		for my $i (1..$last_index_sub){
			
			my $sum = 0;
			foreach my $type_sub (@types_sub){
				$sum  +=  $ref->{$lab_sub}->{$i}->{$type_sub} ;
			}
			
			my $t = $ref->{$lab_sub}->{$i}->{"total"};
			if ( abs( $sum - $t ) > 0.01 ){
				print STDERR "maybe wrong  \n", join("\t", ( $lab_sub, $i )), "\n\n";
				print STDERR "sum = $sum != $t \n\n";
			}
			
			
		}
	}
}


#modified
sub Initialization_type_num{
	my ($ref, $last_index_sub) = @_;
#	my @labs_sub = ("wt", "mut");
	my @types_sub = ("CG", "CHG", "CHH", "total");
	for my $i (1..$last_index_sub){
		foreach my $type_sub (@types_sub){
			$ref->{$i}->{$type_sub} = 0;
			
		}
	}
}


# Initialization(\%seq_depth_sum, $last_one_index);

# modified
sub Initialization{
	my ($ref, $last_index_sub) = @_;
	#my @labs_sub = ("wt", "mut");
	my @labs_sub = ( "mut");
	
	my @types_sub = ("CG", "CHG", "CHH", "total");
	foreach my $i (1..$last_index_sub){
		foreach my $type_sub (@types_sub){
			foreach my $lab_sub (@labs_sub){
				$ref->{$lab_sub}->{$i}->{$type_sub} = 0;
			}
		}
	}
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
	
	return $i;
}



sub output_list{
	my ($file, $ref_list, 	
	    $C_type_nums_ref,  $called_mC_nums_ref, $meth_level_sum_ref,  $seq_mC_sum_ref,  $seq_depth_sum_ref,
	    $sample_label_sub
	) = @_;
	
	my @types = ("CG", "CHG", "CHH", "total");

	die if ( -e $file);
	open(OUT, ">>$file") or die;
	my $last_index = scalar(@{$ref_list}) - 1;
	my $head = $ref_list->[0];
	
#	print OUT join("\t", ($head, @types,
#			      "wmCG_mut", "wmCHG_mut", "wmCHH_mut", "wmC_mut",
#			      "wmCG_wt", "wmCHG_wt", "wmCHH_wt", "wmC_wt",
#			      "fmCG_mut", "fmCHG_mut", "fmCHH_mut", "fmC_mut",
#			      "fmCG_wt", "fmCHG_wt", "fmCHH_wt", "fmC_wt",
#			      "mmCG_mut", "mmCHG_mut", "mmCHH_mut", "mmC_mut",
#			      "mmCG_wt", "mmCHG_wt", "mmCHH_wt", "mmC_wt"
#	)), "\n";
	
	print OUT join("\t", ($head, "CG_CHG_CHH", "mCG_mCHG_mCHH",
			      "wmCG_" . $sample_label_sub, "wmCHG_" . $sample_label_sub, "wmCHH_" . $sample_label_sub, "wmC_" . $sample_label_sub,
			      "fmCG_" . $sample_label_sub, "fmCHG_" . $sample_label_sub, "fmCHH_" . $sample_label_sub, "fmC_" . $sample_label_sub,
			      "mmCG_" . $sample_label_sub, "mmCHG_" . $sample_label_sub, "mmCHH_" . $sample_label_sub, "mmC_" . $sample_label_sub,
	
	)), "\n";
	
	foreach my $i (1..$last_index){
	#	my @weighted_level_wt = ("NA") x 4;
	#	my @frac_level_wt = ("NA") x 4;
	#	my @mean_level_wt = ("NA") x 4;
		
		my @weighted_level_mut = ("NA") x 4;
		my @frac_level_mut = ("NA") x 4;
		my @mean_level_mut = ("NA") x 4;

		my @C_type_nums_a    = (0) x 4;
		
		assign(\%{$C_type_nums_ref -> {$i}}, \@C_type_nums_a);
		
		# cal_meth(\%called_mC_nums, \%C_type_nums, \@frac_level);

	#	cal_meth( \% {$seq_mC_sum_ref->{"wt"}->{$i} }, \% {$seq_depth_sum_ref->{"wt"}->{$i} }, \@weighted_level_wt );
	#	cal_meth( \% {$called_mC_nums_ref->{"wt"}->{$i} }, \% { $C_type_nums_ref->{$i} }, \@frac_level_wt );
	#	cal_meth( \% {$meth_level_sum_ref->{"wt"}->{$i} }, \% { $C_type_nums_ref->{$i} }, \@mean_level_wt );

		cal_meth( \% {$seq_mC_sum_ref->{"mut"}->{$i} }, \% {$seq_depth_sum_ref->{"mut"}->{$i} }, \@weighted_level_mut );
		cal_meth( \% {$called_mC_nums_ref->{"mut"}->{$i} }, \% { $C_type_nums_ref->{$i} }, \@frac_level_mut );
		cal_meth( \% {$meth_level_sum_ref->{"mut"}->{$i} }, \% { $C_type_nums_ref->{$i} }, \@mean_level_mut );
		
		my $C_num_output = join("+", ( $C_type_nums_ref-> {$i}->{CG}, $C_type_nums_ref-> {$i}->{CHG}, $C_type_nums_ref-> {$i}->{CHH} ) );
		
		my $mC_num_output = join("+", ( $called_mC_nums_ref->{ "mut" } -> {$i}->{CG}, $called_mC_nums_ref->{ "mut" } -> {$i}->{CHG}, $called_mC_nums_ref->{ "mut" } -> {$i}->{CHH} ) );

		
		print OUT join("\t", ( $ref_list->[$i], $C_num_output, $mC_num_output, 
				#       @C_type_nums_a,
				      @weighted_level_mut,
				 #     @weighted_level_wt,
				      
				      @frac_level_mut,
				  #    @frac_level_wt,
				      
				      @mean_level_mut
				   #   @mean_level_wt
				)), "\n";		
	}	
	close(OUT);	
}




#get_info_each_chr(	$curr_chr ,  \@positions, \@depths_temp, \@mCs_temp,
#				\%C_type_nums,  \%called_mC_nums, \%meth_level_sum,  \%seq_mC_sum,  \%seq_depth_sum,
#				$depth_cutoff	
#			  );

sub  get_info_each_chr{
	my ($ref_list, $chr_sub,
	    $pos_refa, $depths_temp_refa, $mCs_temp_refa,  $isMeth_temp_refa,  $per_temp_refa,
		$C_type_nums_ref,  $called_mC_nums_ref, $meth_level_sum_ref,  $seq_mC_sum_ref,  $seq_depth_sum_ref,
		$depth_cutoff_sub
	    ) = @_;
	
	my @list_sub = @{$ref_list};
	for my $index(1..$#list_sub){
		my $this = $list_sub[$index];
		next unless( $this =~ /$chr_sub/i);
		
		my @a_sub = split "\t", $this;
		die unless (lc($a_sub[0] ) eq $chr_sub);
		
		my ($start,  $end) = @a_sub[1..2];
		
		for my $j($start..$end){
			if(defined $pos_refa->[$j]){
				my $type = $pos_refa->[$j];
				
				my $per_mut = $per_temp_refa->[$j];
				my $mC_mut  = $mCs_temp_refa->[$j];
				my $dep_mut = $depths_temp_refa->[$j];
				
				$C_type_nums_ref-> {$index}->{$type}++;
				$C_type_nums_ref-> {$index}->{"total"}++;
				
				$meth_level_sum_ref->{"mut"}->{$index}->{$type} += $per_mut;
				$meth_level_sum_ref->{"mut"}->{$index}->{ "total" } += $per_mut;
				$seq_mC_sum_ref->{"mut"}->{$index}->{$type} += $mC_mut;
				$seq_mC_sum_ref->{"mut"}->{$index}->{"total" } += $mC_mut;
				$seq_depth_sum_ref->{"mut"}->{$index}->{$type} += $dep_mut;
				$seq_depth_sum_ref->{"mut"}->{$index}->{"total" } += $dep_mut;
				if( $isMeth_temp_refa ->[$j] == 1){
					$called_mC_nums_ref->{ "mut" } -> {$index}-> {$type}++;
					$called_mC_nums_ref->{"mut"} -> {$index}->{ "total" }++;
				}
			}
		}
		
		
	}
}