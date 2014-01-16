#!/usr/bin/perl -w
 
use strict;
use File::Spec;

my $debug = 0;

my $usage = "$0 \n<bed_like_file> <isMeht_file> <output> <label>\n\n";
die $usage unless (@ARGV == 4);

my $bed_file = shift or die;
my $isMeth_file = shift or die;
my $output = shift or die;
my $label = shift or die;




die unless (-e $bed_file);
die unless (-e $isMeth_file);
die if (-e $output);

if($debug){
	print STDERR "\n\nOK\n\n";
	exit;
}



my @bed_list;
my %pos;
read_bed_file($bed_file, \@bed_list);
record_region_pos(\@bed_list, \%pos);

#my %nums;

my ( %mC_list, %dep_list, @percentage_list );

#read_isMeth_file($isMeth_file, \%pos , \%nums);
read_isMeth_get_mC_dep( $isMeth_file, \%pos, \%mC_list, \%dep_list, \@percentage_list );

%pos = ();

#output($output, \@bed_list, \%nums);
output2($output, \@bed_list, \%mC_list, \%dep_list, \@percentage_list , $label) ;

exit;

sub output2{
	my ($file, $ref_list, 	$mC_list_ref, $dep_list_ref, $percentage_list_ref, $label_sub) = @_;

	my @types = ("CG", "CHG", "CHH");

	die if ( -e $file);
	open(OUT, ">>$file") or die;
	my $last_index = scalar(@{$ref_list}) - 1;
	my $head = $ref_list->[0];
	print OUT join("\t", ($head,  $label_sub."_CG", $label_sub."_CG_per",$label_sub."_CHG", $label_sub."_CHG_per",$label_sub."_CHH", $label_sub. "_CHH_per",
			  			   "avge_meth_level_" . $label_sub )), "\n";
	
	foreach my $i (1..$last_index){
		
		my %divide; # record like "$mC_mut/$dep_mut="
		my %per;
		
		foreach my $type(@types){
			my ($mC , $dep ) = (0) x 2;
			
			if(defined $mC_list_ref ->{$type}->[$i]) {  $mC  =  $mC_list_ref->{$type}->[$i] }
			if(defined $dep_list_ref->{$type}->[$i]) { $dep  =  $dep_list_ref->{$type}->[$i]}
			
			$divide{$type} = "$mC/$dep=";
			
			if($dep != 0) { $per{$type} = sprintf("%.4f",  100 * $mC/$dep);		}
			else          { $per{$type} = "NA";  }	

			
		}
		
		
		my $avge  =  cal_mean(\@{$percentage_list_ref->[$i]});
	
	
		#print OUT join("\t", ($ref_list->[$i], $CG_num, $CHG_num, $CHH_num)), "\n";
		print OUT join("\t", ($ref_list->[$i],  $divide{"CG"}, $per{"CG"},    $divide{"CHG"}, $per{"CHG"},   $divide{"CHH"}, $per{"CHH"} , $avge )), "\n";
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

sub read_isMeth_get_mC_dep{
	my ($file, $region_in_list_ref, 
	$mC_list_ref, $dep_list_ref, $percentage_list_ref ) = @_;
	
	die unless (-e $file);
	
	open(IN, $file) or die;
	
	my $h = <IN>;
	
	
	while (<IN>){
		chomp;
		
		my @a = split "\t";
		my ($chr, $pos) = @a[0..1];
		next unless (defined $region_in_list_ref->{$chr}->[$pos] );
		my $type = $a[3];
		my $index = $region_in_list_ref->{$chr}->[$pos];
		
		######
		#
		############
	#	$C_num_ref->{$index}->{$type} ++;
		#######################
		
		my ($mC , $dep, $per ) = @a[4..6]; 
		$mC_list_ref->{$type}->[$index] += $mC;
		$dep_list_ref->{$type}->[$index] += $dep;
		push @{$percentage_list_ref->[$index]} , $per;
	}
	close(IN);
}


sub output{
	my ($file, $ref_list, $num_ref) = @_;
	die if ( -e $file);
	open(OUT, ">>$file") or die;
	my $last_index = scalar(@{$ref_list}) - 1;
	my $head = $ref_list->[0];
	print OUT join("\t", ($head, "CG_num", "CHG_num", "CHH_num")), "\n";
	
	foreach my $i (1..$last_index){
		my ($CG_num, $CHG_num, $CHH_num ) = (0) x 3;
		if(defined $num_ref->{$i}->{"CG"}) {$CG_num = $num_ref->{$i}->{"CG"}}
		if(defined $num_ref->{$i}->{"CHG"}) {$CHG_num = $num_ref->{$i}->{"CHG"}}
		if(defined $num_ref->{$i}->{"CHH"}) {$CHH_num = $num_ref->{$i}->{"CHH"}}
		print OUT join("\t", ($ref_list->[$i], $CG_num, $CHG_num, $CHH_num)), "\n";
	}	
	close(OUT);	
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

#chr	pos	strand	type	num_C	depth	percentage	isMeth
#chr1	100	+	CHH	0	11	0	0
#chr1	101	+	CHH	0	12	0	0
#chr1	102	+	CHH	0	12	0	1
sub read_isMeth_file{
	my ($file, $pos_ref, $num_ref) = @_;
	
	die unless (-e $file);
	open(IN, $file) or die;
	my $head = <IN>;
	while(<IN>){
		chomp;
		my @a = split "\t";
		my ($chr, $pos, $type) = ($a[0], $a[1], $a[3] );
		if(defined $pos_ref->{$chr}->[$pos]){
			my $index = $pos_ref->{$chr}->[$pos];
			$num_ref->{$index}->{$type}++;
		}
	}
	
	close(IN);
}