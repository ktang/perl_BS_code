#!/usr/bin/perl -w

#v0.3
# from v0.1
#differece: use Ecker's suggestion, if isMet = 0, set mC = 0;
 
use strict;
use File::Spec;

my $debug = 0;

my $depth_cutoff = 4;

my $usage = "$0 \n <isMeht_file_wt> <isMeht_file_mut>  <bed_like_file> <hyper|hypo> <outdir> \n\n";
die $usage unless (@ARGV == 5);

my $isMeth_file_WT	= shift or die;
my $isMeth_file_mut	= shift or die;
my $bed_file		= shift or die;
my $hyper_hypo_label	= shift or die;
my $outdir		= shift or die;
# my $depth_cutoff	= shift or die;

#my $output = shift or die;

die $usage , "\n", "hyper|hypo" unless ( $hyper_hypo_label =~ /hyper|hypo/);

die unless (-e $bed_file);
#die unless (-e $isMeth_file);
die unless (-e $isMeth_file_WT);
die unless (-e $isMeth_file_mut);

die unless ($bed_file =~ /$hyper_hypo_label/);

my  ($volume,$directories,$file) = File::Spec->splitpath( $bed_file );

my $outfile ;

if($file =~ /(\S+)\.txt$/){
	$outfile = $1 . "_detail_EckerPaper.txt";
}

my $output = File::Spec->catfile($outdir, $outfile);

die if (-e $output);

if($debug){
	print STDERR "\n\n";
	print STDERR join("\n", ($isMeth_file_WT, $isMeth_file_mut, $bed_file, $output )), "\n";
	print STDERR "\n\nOK\n\n";
	exit;
}


my @bed_list;
my %pos;


my $last_one_index = read_bed_file($bed_file, \@bed_list);
# head is index 0 and "\n" is chompped.

print STDERR "last: line_num: ", $last_one_index, "\n\n";

record_region_pos(\@bed_list, \%pos);

# weight; frac; mean 3 groups


#my %C_type_nums;
#my %called_mC_nums; # Fraction of methylated cytosines
#my %meth_level_sum; #Mean methylation level
my ( %C_type_nums,  %called_mC_nums, %meth_level_sum,  %seq_mC_sum,  %seq_depth_sum ); #Weighted methylation level

# %{wt/mut}->{1/2}->{CG/CHG/CHH/total}

Initialization_type_num(\%C_type_nums, $last_one_index);

Initialization(\%called_mC_nums, $last_one_index);
Initialization(\%meth_level_sum, $last_one_index);
Initialization(\%seq_mC_sum, $last_one_index);
Initialization(\%seq_depth_sum, $last_one_index);

read_isMeth_get_info( $isMeth_file_WT,  $isMeth_file_mut,  \%pos,
			\%C_type_nums,  \%called_mC_nums, \%meth_level_sum,  \%seq_mC_sum,  \%seq_depth_sum,
			$depth_cutoff		
);


check_type_num(\%C_type_nums, $last_one_index );


check(\%called_mC_nums, $last_one_index );
check(\%meth_level_sum, $last_one_index );
check(\%seq_mC_sum, $last_one_index );
check(\%seq_depth_sum, $last_one_index );



if($hyper_hypo_label eq "hyper"){

	output_hyper_list($output, \@bed_list, 
			  \%C_type_nums,  \%called_mC_nums, \%meth_level_sum,  \%seq_mC_sum,  \%seq_depth_sum
			  );
	
}

elsif( $hyper_hypo_label eq "hypo"){
	output_hypo_list($output, \@bed_list, 
			\%C_type_nums,  \%called_mC_nums, \%meth_level_sum,  \%seq_mC_sum,  \%seq_depth_sum 
			);

}else{
	die $usage;
}

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

sub check{
	my ($ref, $last_index_sub) = @_;
	my @labs_sub = ("wt", "mut");
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

sub Initialization_type_num{
	my ($ref, $last_index_sub) = @_;
	my @labs_sub = ("wt", "mut");
	my @types_sub = ("CG", "CHG", "CHH", "total");
	for my $i (1..$last_index_sub){
		foreach my $type_sub (@types_sub){
			$ref->{$i}->{$type_sub} = 0;
			
		}
	}
}


# Initialization(\%seq_depth_sum, $last_one_index);
sub Initialization{
	my ($ref, $last_index_sub) = @_;
	my @labs_sub = ("wt", "mut");
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


sub record_region_pos{
	my ($ref_list, $ref_h) = @_;
	my $last_index = scalar(@{$ref_list}) - 1;
	
	foreach my $i (1..$last_index){
		my @a = split "\t", $ref_list->[$i];
		my ($chr, $start, $end ) = @a[0..2];
		
		for my $j($start..$end){
			$ref_h->{$chr}->[$j] = $i; #record index
		}
	}
}


sub output_hyper_list{
	my ($file, $ref_list, 	
	    $C_type_nums_ref,  $called_mC_nums_ref, $meth_level_sum_ref,  $seq_mC_sum_ref,  $seq_depth_sum_ref
	) = @_;
	
	my @types = ("CG", "CHG", "CHH", "total");

	die if ( -e $file);
	open(OUT, ">>$file") or die;
	my $last_index = scalar(@{$ref_list}) - 1;
	my $head = $ref_list->[0];
	
#	print OUT join("\t", ($head,  "CG", "CHG", "CHH", "informative_C_num",
#				"mut_CG", "mut_CG_per", "mut_CHG", "mut_CHG_per", "mut_CHH", "mut_CHH_per",
#				 "wt_CG",  "wt_CG_per",  "wt_CHG",  "wt_CHG_per",  "wt_CHH", "wt_CHH_per",
#				 "mut_C_per", "wt_C_per",
#				 "mut_MethLevel_perC", "wt_MethLevel_perC"
#				 
#				 )), "\n";

	print OUT join("\t", ($head, @types,
			      "wmCG_mut", "wmCHG_mut", "wmCHH_mut", "wmC_mut",
			      "wmCG_wt", "wmCHG_wt", "wmCHH_wt", "wmC_wt",
			      "fmCG_mut", "fmCHG_mut", "fmCHH_mut", "fmC_mut",
			      "fmCG_wt", "fmCHG_wt", "fmCHH_wt", "fmC_wt",
			      "mmCG_mut", "mmCHG_mut", "mmCHH_mut", "mmC_mut",
			      "mmCG_wt", "mmCHG_wt", "mmCHH_wt", "mmC_wt"
	
	)), "\n";
	
	foreach my $i (1..$last_index){
		my @weighted_level_wt = ("NA") x 4;
		my @frac_level_wt = ("NA") x 4;
		my @mean_level_wt = ("NA") x 4;
		
		my @weighted_level_mut = ("NA") x 4;
		my @frac_level_mut = ("NA") x 4;
		my @mean_level_mut = ("NA") x 4;

		my @C_type_nums_a    = (0) x 4;
		
		assign(\%{$C_type_nums_ref -> {$i}}, \@C_type_nums_a);
		
		# cal_meth(\%called_mC_nums, \%C_type_nums, \@frac_level);

		cal_meth( \% {$seq_mC_sum_ref->{"wt"}->{$i} }, \% {$seq_depth_sum_ref->{"wt"}->{$i} }, \@weighted_level_wt );
		cal_meth( \% {$called_mC_nums_ref->{"wt"}->{$i} }, \% { $C_type_nums_ref->{$i} }, \@frac_level_wt );
		cal_meth( \% {$meth_level_sum_ref->{"wt"}->{$i} }, \% { $C_type_nums_ref->{$i} }, \@mean_level_wt );

		cal_meth( \% {$seq_mC_sum_ref->{"mut"}->{$i} }, \% {$seq_depth_sum_ref->{"mut"}->{$i} }, \@weighted_level_mut );
		cal_meth( \% {$called_mC_nums_ref->{"mut"}->{$i} }, \% { $C_type_nums_ref->{$i} }, \@frac_level_mut );
		cal_meth( \% {$meth_level_sum_ref->{"mut"}->{$i} }, \% { $C_type_nums_ref->{$i} }, \@mean_level_mut );
		
		
		print OUT join("\t", ( $ref_list->[$i],
				       @C_type_nums_a,
				      @weighted_level_mut,
				      @weighted_level_wt,
				      
				      @frac_level_mut,
				      @frac_level_wt,
				      
				      @mean_level_mut,
				      @mean_level_wt
				)), "\n";		
	}	
	close(OUT);	
}

sub output_hypo_list{
	my ($file, $ref_list, 	
	    $C_type_nums_ref,  $called_mC_nums_ref, $meth_level_sum_ref,  $seq_mC_sum_ref,  $seq_depth_sum_ref
	) = @_;
	
	my @types = ("CG", "CHG", "CHH", "total");

	die if ( -e $file);
	open(OUT, ">>$file") or die;
	my $last_index = scalar(@{$ref_list}) - 1;
	my $head = $ref_list->[0];
	
	print OUT join("\t", ($head, @types,
			      "wmCG_wt", "wmCHG_wt", "wmCHH_wt", "wmC_wt",
			      "wmCG_mut", "wmCHG_mut", "wmCHH_mut", "wmC_mut",
			      
			      "fmCG_wt", "fmCHG_wt", "fmCHH_wt", "fmC_wt",
			      "fmCG_mut", "fmCHG_mut", "fmCHH_mut", "fmC_mut",

			      "mmCG_wt", "mmCHG_wt", "mmCHH_wt", "mmC_wt",
			      "mmCG_mut", "mmCHG_mut", "mmCHH_mut", "mmC_mut"
	)), "\n";
	
	foreach my $i (1..$last_index){
		my @weighted_level_wt = ("NA") x 4;
		my @frac_level_wt = ("NA") x 4;
		my @mean_level_wt = ("NA") x 4;
		
		my @weighted_level_mut = ("NA") x 4;
		my @frac_level_mut = ("NA") x 4;
		my @mean_level_mut = ("NA") x 4;

		my @C_type_nums_a    = (0) x 4;
		
		assign(\%{$C_type_nums_ref -> {$i}}, \@C_type_nums_a);
		
		# cal_meth(\%called_mC_nums, \%C_type_nums, \@frac_level);

		cal_meth( \% {$seq_mC_sum_ref->{"wt"}->{$i} }, \% {$seq_depth_sum_ref->{"wt"}->{$i} }, \@weighted_level_wt );
		cal_meth( \% {$called_mC_nums_ref->{"wt"}->{$i} }, \% { $C_type_nums_ref->{$i} }, \@frac_level_wt );
		cal_meth( \% {$meth_level_sum_ref->{"wt"}->{$i} }, \% { $C_type_nums_ref->{$i} }, \@mean_level_wt );

		cal_meth( \% {$seq_mC_sum_ref->{"mut"}->{$i} }, \% {$seq_depth_sum_ref->{"mut"}->{$i} }, \@weighted_level_mut );
		cal_meth( \% {$called_mC_nums_ref->{"mut"}->{$i} }, \% { $C_type_nums_ref->{$i} }, \@frac_level_mut );
		cal_meth( \% {$meth_level_sum_ref->{"mut"}->{$i} }, \% { $C_type_nums_ref->{$i} }, \@mean_level_mut );
		
		
		print OUT join("\t", ( $ref_list->[$i],
				       @C_type_nums_a,
				       
				      @weighted_level_wt,
				      @weighted_level_mut,
				      
				      @frac_level_wt,
				      @frac_level_mut,
				      
				      @mean_level_wt,
				      @mean_level_mut
				)), "\n";		
	}	
	close(OUT);	
}




# read_isMeth_get_info ( $isMeth_file_WT,  $isMeth_file_mut,  \%pos,
#			\%C_type_nums,  \%called_mC_nums, \%meth_level_sum,  \%seq_mC_sum,  \%seq_depth_sum,
#			$depth_cutoff		
#);

sub read_isMeth_get_info{
	my ($isMeth_file_WT_sub,  $isMeth_file_mut_sub,
		$pos_ref,
		$C_type_nums_ref,  $called_mC_nums_ref, $meth_level_sum_ref,  $seq_mC_sum_ref,  $seq_depth_sum_ref,
		$depth_cutoff_sub) = @_;
	
	die unless (-e $isMeth_file_WT_sub);
	die unless (-e $isMeth_file_mut_sub) ;
	
	open(WT, $isMeth_file_WT_sub) or die;
	open(MUT, $isMeth_file_mut_sub) or die;
	
	my $h_wt = <WT>;
	my $h_mut  = <MUT>;
	
	my ($l_wt, $l_mut);
	
	while( $l_wt = <WT>, $l_mut = <MUT>){
		chomp $l_wt;
		my @a_wt = split "\t", $l_wt;
		my ($chr, $pos) = @a_wt[0..1];
		next unless (defined $pos_ref->{$chr}->[$pos] );
		
		chomp $l_mut;
		my @a_mut  = split "\t", $l_mut;
		die $l_wt, "\n", $l_mut unless ($pos == $a_mut[1]);
		
		my ($mC_wt , $dep_wt,  $per_wt ) = @a_wt[4..6]; 
		my ($mC_mut, $dep_mut, $per_mut ) = @a_mut[4..6]; 
		
		if($dep_wt >=  $depth_cutoff_sub and $dep_mut >= $depth_cutoff_sub ){
			
		# %{wt/mut}->{1/2}->{CG/CHG/CHH/total}
			my $index = $pos_ref->{$chr}->[$pos];
			my $type = $a_mut[3];
			
			if( $a_wt[-1]  == 0){
				$mC_wt =  0;
				$per_wt = 0;
			}elsif(  $a_wt[-1]  ==  1 ){
				$called_mC_nums_ref->{ "wt" } -> {$index}-> {$type}++;
				$called_mC_nums_ref->{"wt"} -> {$index}->{ "total" }++;
			}else{
				die $l_wt, "\n", $l_mut;
			}
		
			if( $a_mut[-1] == 0){
				$mC_mut = 0;
				$per_mut = 0;
			}elsif(  $a_mut[-1]  ==  1 ){
				$called_mC_nums_ref->{ "mut" } -> {$index}-> {$type}++;
				$called_mC_nums_ref->{"mut"} -> {$index}->{ "total" }++;
			}else{
				die $l_wt, "\n", $l_mut;
			}
			
			$C_type_nums_ref-> {$index}->{$type}++;
			$C_type_nums_ref-> {$index}->{"total"}++;
			
			
			$meth_level_sum_ref->{"wt"}->{$index}->{$type} += $per_wt;
			$meth_level_sum_ref->{"wt"}->{$index}->{ "total" } += $per_wt;
			$seq_mC_sum_ref->{"wt"}->{$index}->{$type} += $mC_wt;
			$seq_mC_sum_ref->{"wt"}->{$index}->{"total" } += $mC_wt;
			$seq_depth_sum_ref->{"wt"}->{$index}->{$type} += $dep_wt;
			$seq_depth_sum_ref->{"wt"}->{$index}->{"total" } += $dep_wt;
	
			$meth_level_sum_ref->{"mut"}->{$index}->{$type} += $per_mut;
			$meth_level_sum_ref->{"mut"}->{$index}->{ "total" } += $per_mut;
			$seq_mC_sum_ref->{"mut"}->{$index}->{$type} += $mC_mut;
			$seq_mC_sum_ref->{"mut"}->{$index}->{"total" } += $mC_mut;
			$seq_depth_sum_ref->{"mut"}->{$index}->{$type} += $dep_mut;
			$seq_depth_sum_ref->{"mut"}->{$index}->{"total" } += $dep_mut;
		}
	}
	
	close(WT);
	close(MUT);
}