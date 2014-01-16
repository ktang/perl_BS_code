#!/usr/bin/perl -w
 
use strict;
use File::Spec;

my $debug = 0;

my $usage = "$0 \n <isMeht_file_wt> <isMeht_file_mut>  <bed_like_file> <hyper|hypo> <outdir> \n\n";
die $usage unless (@ARGV == 5);

my $isMeth_file_WT   = shift or die;
my $isMeth_file_mut  = shift or die;
my $bed_file 		 = shift or die;
my $hyper_hypo_label = shift or die;
my $outdir			 = shift or die;

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
	$outfile = $1 . "_detail.txt";
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
read_bed_file($bed_file, \@bed_list);
record_region_pos(\@bed_list, \%pos);


if($hyper_hypo_label eq "hyper"){
	my ( %mC_list_high, %dep_list_high, @percentage_list_high );
	my ( %mC_list_low,  %dep_list_low,  @percentage_list_low );
	
	my %C_num_stat;
	
	#read_isMeth_get_mC_dep high first, then low
	read_isMeth_get_mC_dep2( $isMeth_file_mut,  $isMeth_file_WT, \%pos, 
							\%mC_list_high, \%dep_list_high, \@percentage_list_high,
							\%mC_list_low,  \%dep_list_low,  \@percentage_list_low, \%C_num_stat);
							
	output_hyper_list($output, \@bed_list, 
					  \%mC_list_high, \%dep_list_high, \@percentage_list_high,
					  \%mC_list_low,  \%dep_list_low,  \@percentage_list_low, \%C_num_stat);
	
}elsif( $hyper_hypo_label eq "hypo"){
	my ( %mC_list_high, %dep_list_high, @percentage_list_high );
	my ( %mC_list_low,  %dep_list_low,  @percentage_list_low );
	
	my %C_num_stat;

	
	read_isMeth_get_mC_dep2( $isMeth_file_WT, $isMeth_file_mut,\%pos, 
							\%mC_list_high, \%dep_list_high, \@percentage_list_high,
							\%mC_list_low,  \%dep_list_low,  \@percentage_list_low, \%C_num_stat);
	
	output_hypo_list($output, \@bed_list, 
					 \%mC_list_high, \%dep_list_high, \@percentage_list_high,
					 \%mC_list_low,  \%dep_list_low,  \@percentage_list_low, \%C_num_stat);

}else{
	die $usage;
}

exit;

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
}


sub record_region_pos{
	my ($ref_list, $ref_h) = @_;
	my $last_index = scalar(@{$ref_list}) - 1;
	
	foreach my $i (1..$last_index){
		my @a = split "\t", $ref_list->[$i];
		my ($chr, $start, $end ) = @a[0..2];
		
		for my $j($start..$end){
			$ref_h->{$chr}->[$j] = $i;
		}
	}
}


sub	output_hyper_list{
	my ($file, $ref_list, 	
		$high_mC_list_ref, $high_dep_list_ref, $high_percentage_list_ref,
		$low_mC_list_ref,  $low_dep_list_ref,  $low_percentage_list_ref,
		$C_num_ref) = @_;
	
	my @types = ("CG", "CHG", "CHH");

	die if ( -e $file);
	open(OUT, ">>$file") or die;
	my $last_index = scalar(@{$ref_list}) - 1;
	my $head = $ref_list->[0];
	
	print OUT join("\t", ($head, "informative_C_num", "CG", "CHG", "CHH",
						  "mut_CG", "mut_CG_per", "mut_CHG", "mut_CHG_per", "mut_CHH", "mut_CHH_per",
								 "wt_CG",  "wt_CG_per",  "wt_CHG",  "wt_CHG_per",  "wt_CHH", "wt_CHH_per",
								 "mut_C_per", "wt_C_per",
								 "mut_MethLevel_perC", "wt_MethLevel_perC" )), "\n";
	
	foreach my $i (1..$last_index){
		
		my %high_divide; # record like "$mC_mut/$dep_mut="
		my %high_per;
		
		my %low_divide; # record like "$mC_mut/$dep_mut="
		my %low_per;
		
		my ($high_mC_total, $high_dep_total) = (0) x 2;
		my ($low_mC_total , $low_dep_total ) = (0) x 2;
		
		foreach my $type(@types){
			my ($high_mC , $high_dep ) = (0) x 2;
			my ($low_mC , $low_dep ) = (0) x 2;
			
			if(defined $high_mC_list_ref ->{$type}->[$i]) {  $high_mC  =  $high_mC_list_ref->{$type}->[$i] }
			if(defined $high_dep_list_ref->{$type}->[$i]) { $high_dep  =  $high_dep_list_ref->{$type}->[$i]}
			
			if(defined $low_mC_list_ref ->{$type}->[$i]) {  $low_mC  =  $low_mC_list_ref->{$type}->[$i] }
			if(defined $low_dep_list_ref->{$type}->[$i]) { $low_dep  =  $low_dep_list_ref->{$type}->[$i]}
			

			$high_divide{$type} = "$high_mC/$high_dep=";
			$low_divide{$type} = "$low_mC/$low_dep=";
			
			
			$high_mC_total	+= $high_mC;
			$high_dep_total	+= $high_dep ;
			$low_mC_total	+= $low_mC ;
			$low_dep_total	+= $low_dep;
			
			
			if($high_dep != 0) { $high_per{$type} = sprintf("%.4f",  100 * $high_mC/$high_dep);		}
			else          { $high_per{$type} = "NA";  }	
			
			if($low_dep != 0) { $low_per{$type} = sprintf("%.4f",  100 * $low_mC/$low_dep);		}
			else          { $low_per{$type} = "NA";  }	
		}
		
		my $high_avge  =  cal_mean(\@{$high_percentage_list_ref->[$i]});
		my $low_avge  =  cal_mean(\@{$low_percentage_list_ref->[$i]});
		
		
		my $high_C_per = sprintf("%.4f",  100 * $high_mC_total/$high_dep_total);	
		my $low_C_per =  sprintf("%.4f",  100 * $low_mC_total/$low_dep_total);	
		
		#print OUT join("\t", ($ref_list->[$i], $CG_num, $CHG_num, $CHH_num)), "\n";
		
		my ($C, $CG, $CHG, $CHH) = (0) x 4;
		
		if( defined $C_num_ref->{"C"}->[$i] ){ $C = $C_num_ref->{"C"}->[$i] }
		if( defined $C_num_ref->{"CG"}->[$i] ){ $CG = $C_num_ref->{"CG"}->[$i] }
		if( defined $C_num_ref->{"CHG"}->[$i] ){ $CHG = $C_num_ref->{"CHG"}->[$i] }
		if( defined $C_num_ref->{"CHH"}->[$i] ){ $CHH = $C_num_ref->{"CHH"}->[$i] }
		
				
		print OUT join("\t", ($ref_list->[$i], $C, $CG, $CHG, $CHH, 
							 $high_divide{"CG"}, $high_per{"CG"},    $high_divide{"CHG"}, $high_per{"CHG"},   $high_divide{"CHH"}, $high_per{"CHH"} , 
							 $low_divide{"CG"}, $low_per{"CG"},    $low_divide{"CHG"}, $low_per{"CHG"},   $low_divide{"CHH"}, $low_per{"CHH"} , 
							 $high_C_per, $low_C_per, $high_avge , $low_avge )), "\n";
	}	
	close(OUT);	

	
}

sub output_hypo_list{
	my ($file, $ref_list, 	
		$high_mC_list_ref, $high_dep_list_ref, $high_percentage_list_ref,
		$low_mC_list_ref,  $low_dep_list_ref,  $low_percentage_list_ref,
		$C_num_ref) = @_;
	
	my @types = ("CG", "CHG", "CHH");

	die if ( -e $file);
	open(OUT, ">>$file") or die;
	my $last_index = scalar(@{$ref_list}) - 1;
	my $head = $ref_list->[0];
	
	print OUT join("\t", ($head,  "informative_C_num", "CG", "CHG", "CHH", 
								"wt_CG",  "wt_CG_per",  "wt_CHG",  "wt_CHG_per",  "wt_CHH", "wt_CHH_per", 
								"mut_CG", "mut_CG_per", "mut_CHG", "mut_CHG_per", "mut_CHH", "mut_CHH_per",
								 
								  "wt_C_per", "mut_C_per",
								  "wt_MethLevel_perC","mut_MethLevel_perC")), "\n";
	
	foreach my $i (1..$last_index){
		
		my %high_divide; # record like "$mC_mut/$dep_mut="
		my %high_per;
		
		my %low_divide; # record like "$mC_mut/$dep_mut="
		my %low_per;
		
		my ($high_mC_total, $high_dep_total) = (0) x 2;
		my ($low_mC_total , $low_dep_total ) = (0) x 2;
		
		foreach my $type(@types){
			my ($high_mC , $high_dep ) = (0) x 2;
			my ($low_mC , $low_dep ) = (0) x 2;
			
			if(defined $high_mC_list_ref ->{$type}->[$i]) {  $high_mC  =  $high_mC_list_ref->{$type}->[$i] }
			if(defined $high_dep_list_ref->{$type}->[$i]) { $high_dep  =  $high_dep_list_ref->{$type}->[$i]}
			
			if(defined $low_mC_list_ref ->{$type}->[$i]) {  $low_mC  =  $low_mC_list_ref->{$type}->[$i] }
			if(defined $low_dep_list_ref->{$type}->[$i]) { $low_dep  =  $low_dep_list_ref->{$type}->[$i]}
			

			$high_divide{$type} = "$high_mC/$high_dep=";
			$low_divide{$type} = "$low_mC/$low_dep=";
			
			
			$high_mC_total	+= $high_mC;
			$high_dep_total	+= $high_dep ;
			$low_mC_total	+= $low_mC ;
			$low_dep_total	+= $low_dep;
			
			
			if($high_dep != 0) { $high_per{$type} = sprintf("%.4f",  100 * $high_mC/$high_dep);		}
			else          { $high_per{$type} = "NA";  }	
			
			if($low_dep != 0) { $low_per{$type} = sprintf("%.4f",  100 * $low_mC/$low_dep);		}
			else          { $low_per{$type} = "NA";  }	
		}
		
		my $high_avge  =  cal_mean(\@{$high_percentage_list_ref->[$i]});
		my $low_avge  =  cal_mean(\@{$low_percentage_list_ref->[$i]});
		
		
		my $high_C_per = sprintf("%.4f",  100 * $high_mC_total/$high_dep_total);	
		my $low_C_per =  sprintf("%.4f",  100 * $low_mC_total/$low_dep_total);	
		
		#print OUT join("\t", ($ref_list->[$i], $CG_num, $CHG_num, $CHH_num)), "\n";
		
		my ($C, $CG, $CHG, $CHH) = (0) x 4;
		
		if( defined $C_num_ref->{"C"}->[$i] ){ $C = $C_num_ref->{"C"}->[$i] }
		if( defined $C_num_ref->{"CG"}->[$i] ){ $CG = $C_num_ref->{"CG"}->[$i] }
		if( defined $C_num_ref->{"CHG"}->[$i] ){ $CHG = $C_num_ref->{"CHG"}->[$i] }
		if( defined $C_num_ref->{"CHH"}->[$i] ){ $CHH = $C_num_ref->{"CHH"}->[$i] }
		
		print OUT join("\t", ($ref_list->[$i],  $C, $CG, $CHG, $CHH,
							 $high_divide{"CG"}, $high_per{"CG"},    $high_divide{"CHG"}, $high_per{"CHG"},   $high_divide{"CHH"}, $high_per{"CHH"} , 
							 $low_divide{"CG"}, $low_per{"CG"},    $low_divide{"CHG"}, $low_per{"CHG"},   $low_divide{"CHH"}, $low_per{"CHH"} , 
							 $high_C_per, $low_C_per, $high_avge , $low_avge )), "\n";
	}	
	close(OUT);	


}


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
		return ( sprintf ( "%.4f", 100 * $sum / $num));
	}
}

#read_isMeth_get_mC_dep2(   $isMeth_file_WT, $isMeth_file_mut,\%pos, 
#							\%mC_list_high, \%dep_list_high, \@percentage_list_high,
#							\%mC_list_low,  \%dep_list_low,  \@percentage_list_low );


sub read_isMeth_get_mC_dep2{
	my ($high_file, $low_file, $region_in_list_ref, 
		$high_mC_list_ref, $high_dep_list_ref, $high_percentage_list_ref,
		$low_mC_list_ref,  $low_dep_list_ref,  $low_percentage_list_ref ,
		$C_num_ref ) = @_;
		
	die unless (-e $high_file);
	die unless (-e $low_file) ;
	
	open(HIGH, $high_file) or die;
	open(LOW,  $low_file ) or die;
	
	my $h_high = <HIGH>;
	my $h_low  = <LOW>;
	
	my ($l_high, $l_low);
	while( $l_high = <HIGH>, $l_low = <LOW>){
		chomp $l_high;
		my @a_high = split "\t", $l_high;
		my ($chr, $pos) = @a_high[0..1];
		next unless (defined $region_in_list_ref->{$chr}->[$pos] );
		
		chomp $l_low;
		my @a_low  = split "\t", $l_low;
		die $l_high, "\n", $l_low unless ($pos == $a_low[1]);
		
		my $type = $a_low[3];
		my $index = $region_in_list_ref->{$chr}->[$pos];

		$C_num_ref->{$type}->[$index] ++;
		$C_num_ref->{"C"}->[$index] ++;
		
		my ($mC_high , $dep_high, $per_high ) = @a_high[4..6]; 
		$high_mC_list_ref->{$type}->[$index] += $mC_high;
		$high_dep_list_ref->{$type}->[$index] += $dep_high;
		push @{$high_percentage_list_ref->[$index]} , $per_high;
		
		my ($mC_low , $dep_low, $per_low ) = @a_low[4..6]; 
		$low_mC_list_ref->{$type}->[$index] += $mC_low;
		$low_dep_list_ref->{$type}->[$index] += $dep_low;
		push @{$low_percentage_list_ref->[$index]} , $per_low;
		
	}
	
	close(HIGH);
	close(LOW);	
}