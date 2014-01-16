#!/usr/bin/perl -w
 
use strict;
use File::Spec;

my $debug = 0;

my $usage = "$0 \n<bed_like_file> <isMeht_file_WT> <isMeth_file_mut> <output>\n\n";
die $usage unless (@ARGV == 4);

my $bed_file = shift or die;
#my $isMeth_file = shift or die;
my $wt_in  = shift or die;
my $mut_in = shift or die;
my $output = shift or die;

die $wt_in  unless (-e $wt_in);
die $mut_in unless (-e $mut_in);

die unless (-e $bed_file);
#die unless (-e $isMeth_file);
die if (-e $output);

if($debug){
	print STDERR "\n\nOK\n\n";
	exit;
}

open (WT,  $wt_in )  or die "cannot open  $wt_in: $!";
open (MUT, $mut_in ) or die "cannot open $mut_in: $!";

# 0      1       2         3      4       5       6                7
#chr     pos     strand  type    num_C   depth   percentage      isMeth
#chr1    1       +       CHH      0       0       0	               0
# 0      1       2         3      4       5       6                7




my @bed_list;
my %pos_h;
read_bed_file($bed_file, \@bed_list);
record_region_pos(\@bed_list, \%pos_h);


my (%hyper_DMCs, %hypo_DMCs);
my ($l_wt, $l_mut);

my $h_wt  = <WT>;
my $h_mut = <MUT>;

#my $line = 0;

while($l_wt = <WT>, $l_mut = <MUT>){

	chomp $l_wt;
	chomp $l_mut;
	my @a_wt  = split "\t", $l_wt;
	my @a_mut = split "\t", $l_mut;
	
	die unless ($a_wt[1] == $a_mut[1]);
	my $chr = $a_wt[0];
	my $pos = $a_wt[1];
	next unless (defined $pos_h{$chr}->[$pos]);
	my $type = $a_wt[3];
	
	my $dep_wt = $a_wt[5];
	my $dep_mut = $a_mut[5];
	
	my $mC_wt = $a_wt[4];
	my $mC_mut = $a_mut[4];
	
	my $per_wt  = $a_wt[6];
	my $per_mut = $a_mut[6];

	if($per_wt == 0){
		if($mC_mut >= 2){
	#		$hyper_DMCs{$chr}->[$pos] = 1;
			$hyper_DMCs{$chr}->[$pos] = $type;
		}
	}elsif($per_mut / $per_wt >= 2){
	#	$hyper_DMCs{$chr}->[$pos] = 1;
		$hyper_DMCs{$chr}->[$pos] = $type;
	}
	
	if($per_mut == 0){
		if($mC_wt >= 2){
		#	$hypo_DMCs{$chr}->[$pos] = 1;
			$hypo_DMCs{$chr}->[$pos] = $type;
		}
	}elsif($per_wt / $per_mut >= 2){
	#	$hypo_DMCs{$chr}->[$pos] = 1;
		$hypo_DMCs{$chr}->[$pos] = $type;
	}
}

close(WT);
close(MUT);

%pos_h=(); # empty the hash



#output($output, \@bed_list, \%nums);

output_DMC_type($output, \@bed_list, \%hyper_DMCs, \%hypo_DMCs);

exit;
#output_DMC_type($output, \@bed_list, \%hyper_DMCs, \%hypo_DMCs);

sub output_DMC_type{
	my ($file, $ref_list, $hyper_ref, $hypo_ref) = @_;
	die if ( -e $file);
	open(OUT, ">>$file") or die;
	my $last_index = scalar(@{$ref_list}) - 1;
	my $head = $ref_list->[0];
	print OUT join("\t", ($head, "hyper_DMC", "hyper_DMC_type", "hypo_DMC", "hypo_DMC_type")), "\n";
	
	foreach my $index (1..$last_index){
		
		my %hyper_num = ( "C" =>0, "CG" => 0, "CHG" => 0, "CHH" => 0 );
		my %hypo_num = ( "C" =>0, "CG" => 0, "CHG" => 0, "CHH" => 0 );
		
		my @a = split "\t", $ref_list->[$index];
		my ($chr, $start, $end) = @a[0..2];
		for my $j ($start..$end){
			if(defined $hyper_ref->{$chr}->[$j]){
				my $type = $hyper_ref->{$chr}->[$j];
				$hyper_num{C}++;
				$hyper_num{$type}++;
			}
			if(defined $hypo_ref->{$chr}->[$j]){
				my $type = $hypo_ref->{$chr}->[$j];
				$hypo_num{C}++;
				$hypo_num{$type}++;
			}
		}
		my $hyper_label = join("+", ( $hyper_num{CG}, $hyper_num{CHG}, $hyper_num{CHH} ));
		my $hypo_label  = join("+", ( $hypo_num{CG}, $hypo_num{CHG}, $hypo_num{CHH} ));
		print OUT join("\t", (@a, $hyper_num{C}, $hyper_label, $hypo_num{C}, $hypo_label  )), "\n";
		
	}
	
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