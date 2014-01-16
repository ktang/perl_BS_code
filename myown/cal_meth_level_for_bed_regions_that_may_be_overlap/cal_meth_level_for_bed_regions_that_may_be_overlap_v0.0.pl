#!/usr/bin/perl -w

#input db is 
# /Volumes/Macintosh_HD_2/idm1_new_met_data/UCR_server_data/JKZ131_high_WT_right/Lane8_JKZ_Col0_131_s_8_isMeth.txt

#input_liet format
# chr start end strand

use strict;
use File::Spec;

my $debug = 0;
if($debug){
	print STDERR "debug = 1\n\n";
}


print STDERR "\ninput bed must have head", "\n\n";

my $detail = 0;

if($debug == 0){
	$detail = 0;
}


#my  $usage = "$0 \n<dep_cutoff> <sample_label> <isMeth_file> <feature_bin_num> <flaking_bp> <flaking_bin_num> <input_bed_list> <outdir> <feature_pre>\n\n";
#die $usage unless(@ARGV == 9);

my  $usage = "$0 \n<dep_cutoff> <sample_label> <isMeth_file>  <input_bed_list> <outdir> <feature_pre | like 1kb>\n\n";
die $usage unless(@ARGV == 6);



my $dep_cutoff = shift or die;
my $sample_label = shift or die;
my $isMeth_file = shift or die "isMeth_file";


my $input_list_file = shift or die "input_list_file";
my $outdir   = shift or die "outdir";
my $feature_pre   = shift or die "outpre";

die unless (-e $input_list_file);
die unless (-d $outdir);

die unless (-e $isMeth_file and -e $input_list_file);


#my $output = File::Spec->catfile($outdir, $outpre . "_in_" . $label . "_binNum" . $feature_bin_num . "_flaking" . $flaking_bp . "bp_binNum" . $flaking_bin_num . ".txt");
#my $output = File::Spec->catfile($outdir,
#			 $feature_pre . "_in_" . $sample_label . "_binNum" . $feature_bin_num . "_flaking" . $flaking_bp . "bp_binNum" . $flaking_bin_num . "_dep" . $dep_cutoff .".txt");

my $output = File::Spec->catfile($outdir,
			 $feature_pre . "_in_" . $sample_label .  "_dep" . $dep_cutoff ."_meth_level.txt");

print STDERR "output:\n$output\n\n";
die if (-e $output);

if ($debug){
	#print STDERR $output, "\n\nOK\n\n";
	exit;
}

#my (@target_list , @left_list);
my @input_list;

#read_list($target_file, \@target_list);
#read_list($left_file,   \@left_list);

read_list($input_list_file, \@input_list);


my (%total_depths, %total_mCs); #$hash{C/CG/CHG/CHH}->[index] += 

my %numbers; # $hash{C/CG/CHG/CHH}->[index]

#chr     pos     strand  type    num_C   depth   percentage      isMeth
#chr1    109     +       CG      10      11      0.909091        1
#chr1    114     +       CHG     2       12      0.166667        1


for my $i_chr (1..5){
	print STDERR  "reading chr$i_chr\n";
#	my @positions;
#	my (@nums_wt, @nums_mut);
#	my (@mCs_wt, @mCs_mut);
	#my (@positions, @nums_wt, @nums_mut, @mCs_wt, @mCs_mut);	

#	my (@positions, @depths_temp, @mCs_temp);
	my (@types_temp, @depths_temp, @mCs_temp);
	
	open(DB, $isMeth_file) or die "db";

	my $curr_chr = "chr" . $i_chr;
	
	my $line = 0;
	my $const = 3000000;
	
	my $db_head = <DB>;
	
	while(<DB>){
		$line++;
		if($line % $const == 0){
			print STDERR $line, "\t";
		}
		next unless(/$curr_chr/i);
		my @a = split "\t";
		next unless ( $a[5] >= $dep_cutoff );
		die unless (lc($a[0]) eq $curr_chr);

		my $pos = $a[1];
		my $type = $a[3];
		
		$types_temp[$pos]   = $type;
		$depths_temp[$pos] = $a[5];
		$mCs_temp[$pos]    = $a[4];		
	}
	close(DB);
	print STDERR "\n\n";

	cal_meth_level($curr_chr, \@input_list, \@types_temp, \@depths_temp, \@mCs_temp,  \%total_depths, \%total_mCs, \%numbers );
	
}

output_list($output, \@input_list, \%total_depths, \%total_mCs, \%numbers );

exit;

#output_list($output, \@input_list, \%total_depths, \%total_mCs, \%numbers );

sub output_list{
	my ($file, $ref_list_a, $ref_total_depths_h, $ref_total_mCs_h, $ref_number_h) = @_;
	die if(-e $file);
	open(OUT, ">>$file") or die;
	
	my $head = $ref_list_a->[0];
	
	my $last_index = scalar(@{$ref_list_a}) - 1;
	
	print OUT join("\t", ($head, "informative_C_num", "CG_num", "CHG_num", "CHH_num",
							 "CG", "CG_per", "CHG", "CHG_per", "CHH", "CHH_per", "all_C_percentage"  )), "\n";
	
	my @types = ("C", "CG", "CHG", "CHH");

	
	foreach my $i (1..$last_index){
	
		my %divide; # record like "$mC_mut/$dep_mut="
		my %per;
		my %nums;

		
		foreach my $type(@types){
			my ($mC , $dep, $number ) = (0) x 3;
			
			if(defined $ref_total_mCs_h ->{$type}->[$i]) {  $mC  =  $ref_total_mCs_h->{$type}->[$i] }
			if(defined $ref_total_depths_h->{$type}->[$i]) { $dep  =  $ref_total_depths_h->{$type}->[$i]}
			if(defined $ref_number_h->{$type}->[$i])     {$number  = $ref_number_h->{$type}->[$i] }
			
			$nums{$type} = $number ;
			
			$divide{$type} = "$mC/$dep=";
			if($dep != 0) { $per{$type} = sprintf("%.4f",  100 * $mC/$dep);		}
			else          { $per{$type} = "NA";  }	
		}
		
		print OUT join("\t", ($ref_list_a->[$i], $nums{C}, $nums{CG}, $nums{CHG}, $nums{CHH}, 
											   $divide{CG}, $per{CG}, 
											   $divide{CHG}, $per{CHG}, 
											   $divide{CHH}, $per{CHH}, 
											   $per{C} ) ) , "\n";
	}
	
	close(OUT);
}

# cal_meth_level($curr_chr, \@input_list, \@types_temp, \@depths_temp, \@mCs_temp,  \%total_depths, \%total_mCs, \%numbers );

sub cal_meth_level{
	my ($chr_sub, $ref_list_a, $ref_type_a, $ref_temp_dep_a, $ref_temp_mC_a, $ref_total_depths_h, $ref_total_mCs_h, $ref_number_h) = @_;
	
	my $last_index = scalar(@{$ref_list_a}) - 1;
	foreach my $index (1..$last_index){
		my $line = $ref_list_a->[$index];
		my @a =  split "\t", $line;
		my ($chr, $start, $end ) = @a[0..2];
		
		next unless ($chr eq $chr_sub);
		
		for my $i ($start..$end){
			if(defined $ref_type_a->[$i]){
				my $type = $ref_type_a->[$i];
				my $dep  = $ref_temp_dep_a->[$i];
				my $mC   = $ref_temp_mC_a->[$i];
				
				$ref_number_h->{C}->[$index] ++;
				$ref_number_h->{$type}->[$index]++;
				
				$ref_total_depths_h->{C}->[$index] += $dep;
				$ref_total_mCs_h   ->{C}->[$index] += $mC;
				
				$ref_total_depths_h->{$type}->[$index] += $dep;
				$ref_total_mCs_h->{$type}->[$index]    += $mC;
			}#if(defined $ref_type_a->[$i])
		}#for my $i ($start..$end)
	}#foreach my $index
}



#sub read_list2{
#	my ($file, $ref) = @_;
#	die unless (-e $file);
#	open(IN, $file)  or die "cannot open $file:$!\n\n";
#	my $i = -1;
#	my $head = <IN>;
#	chomp $head;
#	my @temp = split "\t", $head;
#	unless(lc($temp[0]) eq "chr"){
#		$i++;
#		$ref->[$i] = join("_", @temp[0..3]);
#	}
#	while(<IN>){
#		$i++;
#		chomp;
#		my @a = split "\t";
#		$a[0] = lc($a[0]);
#		$ref->[$i] = join("_", @a[0..3]);
#	}
#	close(IN);
#}#

sub read_list{
	my ($file, $ref) = @_;
	die unless (-e $file);
	open(IN, $file)  or die ;
#	my $head = <IN>;
	my $i = -1;
	while(<IN>){
		$i++;
		chomp;
		my @a = split "\t";
		#$ref->[$i] = join("_", @a[0..5]);
		$ref->[$i] = $_;
	}
	close(IN);
}

sub round {
    my($number) = shift;
    #return int($number + .5);
    return int($number + 0.5 * ($number <=> 0)); # take care of negative numbers too

}