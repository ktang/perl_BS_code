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

my $detail = 0;

if($debug == 0){
	$detail = 0;
}


my  $usage = "$0 \n<dep_cutoff> <sample_label> <isMeth_file> <feature_bin_num> <flaking_bp> <flaking_bin_num> <input_bed_list> <outdir> <feature_pre>\n\n";
die $usage unless(@ARGV == 9);

my $dep_cutoff = shift or die;
my $sample_label = shift or die;
my $isMeth_file = shift or die "isMeth_file";
#die $usage, "\n\nWT\n\n" unless ($WT_label eq "WT");

my $feature_bin_num = shift or die "feature_bin_num";
my $flaking_bp = shift or die "flaking_bp";
my $flaking_bin_num = shift or die "flaking_bin_num";
#( $feature_bin_num, $flaking_bp, $flaking_bin_num )

my $input_list_file = shift or die "input_list_file";
my $outdir   = shift or die "outdir";
my $feature_pre   = shift or die "outpre";

die unless (-e $input_list_file);
die unless (-d $outdir);

die unless (-e $isMeth_file and -e $input_list_file);


#my $output = File::Spec->catfile($outdir, $outpre . "_in_" . $label . "_binNum" . $feature_bin_num . "_flaking" . $flaking_bp . "bp_binNum" . $flaking_bin_num . ".txt");
my $output = File::Spec->catfile($outdir,
			 $feature_pre . "_in_" . $sample_label . "_binNum" . $feature_bin_num . "_flaking" . $flaking_bp . "bp_binNum" . $flaking_bin_num . "_dep" . $dep_cutoff .".txt");

print STDERR "output:\n$output\n\n";
die if (-e $output);

if ($debug){
	print STDERR "\n\nOK\n\n";
	exit;
}

#my (@target_list , @left_list);
my @input_list;

#read_list($target_file, \@target_list);
#read_list($left_file,   \@left_list);

read_list2($input_list_file, \@input_list);


my (%total_depths, %total_mCs); #$hash{CG/CHG/CHH}->[0..(2*flaking_bin_num + feature_bin_num -1 )] += 

#chr     pos     strand  type    num_C   depth   percentage      isMeth
#chr1    109     +       CG      10      11      0.909091        1
#chr1    114     +       CHG     2       12      0.166667        1


for my $i_chr (1..5){
	print STDERR  "reading chr$i_chr\n";
#	my @positions;
#	my (@nums_wt, @nums_mut);
#	my (@mCs_wt, @mCs_mut);
	#my (@positions, @nums_wt, @nums_mut, @mCs_wt, @mCs_mut);	
	my (@positions, @depths_temp, @mCs_temp);
	
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
#		my $type = $a[5];
		my $type = $a[3];
		
		$positions[$pos]   = $type;
	#	$depths_temp[$pos] = $a[2];
	#	$mCs_temp[$pos]    = round($a[2] * $a[3]);
		$depths_temp[$pos] = $a[5];
		$mCs_temp[$pos]    = $a[4];		
	}
	close(DB);
	print STDERR "\n\n";

#	cal_list($curr_chr, \@target_list, \@positions, \@nums_wt, \@nums_mut,    \@mCs_wt,     \@mCs_mut,    \%total_wt_target, \%total_mut_target, \%mCXX_wt_target, \%mCXX_mut_target);
	cal_list($curr_chr, \@input_list, \@positions, \@depths_temp, \@mCs_temp,  \%total_depths, \%total_mCs, $feature_bin_num, $flaking_bp, $flaking_bin_num );
	
	#if(!$debug){
	#	cal_list($curr_chr, \@left_list,   \@positions, \@nums_wt, \@nums_mut,    \@mCs_wt,     \@mCs_mut,    \%total_wt_left,   \%total_mut_left,   \%mCXX_wt_left,   \%mCXX_mut_left);
	#}
#	       my($chr_sub, $list_refa,     $pos_refa,$num_wt_refa, $num_mut_refa, $mC_wt_refa, $mC_mut_refa, $total_wt_refh,     $total_mut_refh,    $mCXX_wt_refh,     $mCXX_mut_refh  ) 
}




#output($out_target, \%total_wt_target, \%total_mut_target, \%mCXX_wt_target, \%mCXX_mut_target);
#if(!$debug){
#	output($out_left  , \%total_wt_left,   \%total_mut_left,   \%mCXX_wt_left,   \%mCXX_mut_left);
#}
my $total_intervel_num = 2 * $flaking_bin_num + $feature_bin_num; 		
output2($output, \%total_depths, \%total_mCs, $total_intervel_num);

exit;

sub output2{
	my ($file, $total_depth_refh,  $total_mC_refh, $total_intervel_num_sub) = @_;
	my @types  = ("CG", "CHG", "CHH");
	die if(-e $file);
	
	open(OUT, ">>$file") or die "cannot open $file:$!";
	
	print OUT join("\t", ("order", "mCG_num", "CG_total_num", "CG_per",
									"mCHG_num", "CHG_total_num", "CHG_per",
									"mCHH_num", "CHH_total_num", "CHH_per",
									"mC_num", "C_total_num", "C_per")), "\n";
	
	for my $i(0..($total_intervel_num - 1)){
		my (%depths_sub, %mC_sub, %pers_sub);
		
		for my $type (@types){
			$depths_sub{$type} = 0;
			$mC_sub{$type} = 0;
			$pers_sub{$type} = 0;
						
			if(defined $total_depth_refh ->{$type}->[$i]) { $depths_sub {$type} = $total_depth_refh ->{$type}->[$i]}
			if(defined $total_mC_refh ->{$type}->[$i]) { $mC_sub {$type} = $total_mC_refh ->{$type}->[$i]}
						
			if($depths_sub {$type}  > 0){$pers_sub{$type} = sprintf ("%.4f", 100 * $mC_sub{$type} / $depths_sub {$type}) }
			
		}
		
		my $all_mC_sub = $mC_sub{CG} + $mC_sub{CHG} + $mC_sub{CHH};
		my $all_depth_sub = $depths_sub{CG} + $depths_sub{CHG} + $depths_sub{CHH};
		my $per_C = 0;
		if($all_depth_sub > 0){
			$per_C = sprintf ("%.4f", 100 * $all_mC_sub / $all_depth_sub);
		}
		
				
		print OUT $i + 1, "\t";
		for my $type (@types){
			print OUT join("\t", ($mC_sub{$type}, $depths_sub{$type}, $pers_sub{$type})), "\t";
		}
		
		print OUT join("\t", ($all_mC_sub, $all_depth_sub, $per_C)),"\n";
	}
	close(OUT);	 
}


#up:   0-9;
#TE:   10-19;
#down: 20-29;

sub cal_list{
#	my ($chr_sub, $list_refa,     $pos_refa,$num_wt_refa, $num_mut_refa, $mC_wt_refa, $mC_mut_refa, $total_wt_refh,     $total_mut_refh,    $mCXX_wt_refh,     $mCXX_mut_refh  )  = @_;
#	cal_list($curr_chr, \@input_list, \@positions, \@depths_temp,    \@mCs_temp,     \%total_depths,     \%total_mCs, $feature_bin_num, $flaking_bp, $flaking_bin_num );
	my ($chr_sub,         $list_refa, $pos_refa,   $depths_temp_refa, $mCs_temp_refa, $total_depth_refh,  $total_mC_refh, $feature_bin_num_sub, $flaking_bp_sub, $flaking_bin_num_sub )  = @_;
#	my $interval_num = 10;
#	my $flaking_bp = 1000;
#	my $flaking_interval = int( $flaking_bp  / $interval_num  );

	my $flaking_interval = int( $flaking_bp_sub / $flaking_bin_num_sub);
	
	#total intervel = 2* $flaking_bin_num + $feature_bin_num
	
#	my $total_interval_num = 2* $flaking_bin_num + $feature_bin_num;
	
	foreach my $j(0..( scalar(@{$list_refa}) -1 )){
		my $this = $list_refa->[$j];
		next unless ($this =~ /$chr_sub/);
	#	my ($chr, $start, $end, $length, $strand, $id) = split "_", $this;
	
		my ($chr, $start, $end, $strand) = split "_", $this;
		die unless($chr eq $chr_sub);
		my $part_len = int( ($end - $start + 1) / $feature_bin_num_sub);# interval length
		
		if($strand eq "+"){
			
			##############
			#	feature body
			##############
			
			for my $k(0..($feature_bin_num_sub - 2)){
				my $temp_start = $start + $k * $part_len;
				my $temp_end   = $temp_start + $part_len - 1;
				
				if($detail){
					print STDERR "+ TE:\t";
					print STDERR join("\t", ($temp_start, $temp_end) ), "\n";
				}
				
				for my $l( $temp_start..$temp_end ){
					if(defined $pos_refa->[$l]){
						my $type = $pos_refa->[$l];
						$total_depth_refh->{$type}-> [$k + $flaking_bin_num_sub] += $depths_temp_refa->[$l] ;# notice : change []
						$total_mC_refh ->{$type}-> [$k + $flaking_bin_num_sub] += $mCs_temp_refa->[$l] ;     # notice : change []
					
					}#if defined
				}#for l
			}#for k
			
			my $temp_start = $start + ($feature_bin_num_sub - 1) * $part_len;
			my $temp_end   = $end;
			if($detail){
				print STDERR "+ TE last:\t";
				print STDERR join("\t", ($temp_start, $temp_end) ), "\n";
			}
			#for my $l( ($start + ($interval_num - 1) * $part_len)..$end){
			for my $l( $temp_start..$temp_end ){
				if(defined $pos_refa->[$l]){
					my $type = $pos_refa->[$l];
					
					$total_depth_refh->{$type}-> [($feature_bin_num_sub - 1) + $flaking_bin_num_sub] += $depths_temp_refa->[$l] ;# notice : change []
					$total_mC_refh ->{$type}-> [ ($feature_bin_num_sub - 1) + $flaking_bin_num_sub] += $mCs_temp_refa->[$l] ;     # notice : change []
				}
			}
			
			###########################
			# for up stream  $flaking_bp
			##########################
			my $up_start = $start - $flaking_bp;
			if($up_start <= 0){
				next; 
				print STDERR $this, "\n";
			}
			
			for my $k(0..($flaking_bin_num_sub - 1)){
				my $temp_start = $up_start + $k * $flaking_interval;
				my $temp_end   = $temp_start + $flaking_interval - 1;
				
				if($detail){
					print STDERR "+ up:\t";
					print STDERR join("\t", ($temp_start, $temp_end) ), "\n";
				}
				
				for my $l( $temp_start..$temp_end ){
					if(defined $pos_refa->[$l]){
						my $type = $pos_refa->[$l];
						
						$total_depth_refh->{$type}-> [$k ] += $depths_temp_refa->[$l] ;# notice : change []
						$total_mC_refh ->{$type}-> [$k ] += $mCs_temp_refa->[$l] ;
					}#if defined
				}
				
			}
			
			###########################
			# for down stream  $flaking_bp
			##########################
			my $down_start = $end + 1;
			for my $k(0..($flaking_bin_num_sub - 1)){
				my $temp_start = $down_start + $k * $flaking_interval;
				my $temp_end   = $temp_start + $flaking_interval - 1;
				
				if($detail){
					print STDERR "+ down:\t";
					print STDERR join("\t", ($temp_start, $temp_end) ), "\n";
				}
					
				for my $l( $temp_start..$temp_end ){
					if(defined $pos_refa->[$l]){
						my $type = $pos_refa->[$l];
						$total_depth_refh->{$type}-> [$k + $flaking_bin_num_sub + $feature_bin_num_sub] += $depths_temp_refa->[$l] ;# notice : change []
						$total_mC_refh ->{$type}-> [$k + $flaking_bin_num_sub + $feature_bin_num_sub] += $mCs_temp_refa->[$l] ;
					}#if defined
				}
			}
			if($detail){
						print STDERR "\n\n\n";
			}
		}
		elsif($strand eq "-"){
			
			###########################
			#  feature body
			##########################
			
			for my $k(0..($feature_bin_num_sub - 2)){
				my $temp_end   = $end - $k * $part_len;
				my $temp_start = $temp_end  + 1 - $part_len;
				
				
				if($detail){
					print STDERR "- TE:\t";
					print STDERR join("\t", ($temp_start, $temp_end) ), "\n";
				}
				
				
				for my $l( $temp_start..$temp_end ){
					if(defined $pos_refa->[$l]){
						my $type = $pos_refa->[$l];
						$total_depth_refh->{$type}-> [$k + $flaking_bin_num_sub ] += $depths_temp_refa->[$l] ;# notice : change []
						$total_mC_refh ->{$type}-> [$k + $flaking_bin_num_sub] += $mCs_temp_refa->[$l] ;
					}#if defined
				}#for l
			}#for k
			
			#########last interval
			my $temp_start = $start;
			my $temp_end   = $end - ($feature_bin_num_sub - 1) * $part_len;
			
			if($detail){
				print STDERR "- TE last:\t";
				print STDERR join("\t", ($temp_start, $temp_end) ), "\n";
			}
			
			for my $l( $temp_start..$temp_end ){
				if(defined $pos_refa->[$l]){
					my $type = $pos_refa->[$l];
					$total_depth_refh->{$type}-> [($feature_bin_num_sub -  1) + $flaking_bin_num_sub ] += $depths_temp_refa->[$l] ;# notice : change []
					$total_mC_refh ->{$type}-> [ ($feature_bin_num_sub - 1) + $flaking_bin_num_sub] += $mCs_temp_refa->[$l] ;
				}#if defined
			}#for l
			
			###########################
			# for up stream  $flaking_bp
			##########################
			my $up_end = $end + $flaking_bp;
			for my $k(0..($flaking_bin_num_sub - 1)){
				my $temp_end   = $up_end - $k * $flaking_interval ;
				my $temp_start = $temp_end  + 1 - $flaking_interval;
				
				if($detail){
					print STDERR "- up:\t";
					print STDERR join("\t", ($temp_start, $temp_end) ), "\n";
				}
				
				for my $l( $temp_start..$temp_end ){
					if(defined $pos_refa->[$l]){
						my $type = $pos_refa->[$l];
						$total_depth_refh->{$type}-> [$k] += $depths_temp_refa->[$l] ;# notice : change []
						$total_mC_refh ->{$type}-> [$k ] += $mCs_temp_refa->[$l] ;     # notice : change []
					}#if defined
				}#for l
			}#for k
			
			###########################
			# for down stream  $flaking_bp
			##########################
			if($start -$flaking_bp <=0 ){
				next; 
				print STDERR $this, "\n";
			}
			
			my $down_end = $start - 1;
			for my $k(0..($flaking_bin_num_sub - 1)){
				my $temp_end   = $down_end - $k * $flaking_interval ;
				my $temp_start = $temp_end  + 1 - $flaking_interval;
				
				if($detail){
					print STDERR "- down:\t";
					print STDERR join("\t", ($temp_start, $temp_end) ), "\n";
				}
				
				for my $l( $temp_start..$temp_end ){
					if(defined $pos_refa->[$l]){
						my $type = $pos_refa->[$l];
						$total_depth_refh->{$type}-> [$k + $flaking_bin_num_sub + $feature_bin_num_sub] += $depths_temp_refa->[$l] ;# notice : change []
						$total_mC_refh ->{$type}-> [$k + $flaking_bin_num_sub + $feature_bin_num_sub] += $mCs_temp_refa->[$l] ;     # notice : change []
					}#if defined
				}#for l
				
			}#for k
			
		}else{
			die $this;
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
	unless(lc($temp[0]) eq "chr"){
		$i++;
		$ref->[$i] = join("_", @temp[0..3]);
	}
	while(<IN>){
		$i++;
		chomp;
		my @a = split "\t";
		$a[0] = lc($a[0]);
		$ref->[$i] = join("_", @a[0..3]);
	}
	close(IN);
}

#sub read_list{
#	my ($file, $ref) = @_;
#	die unless (-e $file);
#	open(IN, $file)  or die ;
#	my $head = <IN>;
#	my $i = -1;
#	while(<IN>){
#		$i++;
#		chomp;
#		my @a = split "\t";
#		$ref->[$i] = join("_", @a[0..5]);
#	}
#	close(IN);
#}

sub round {
    my($number) = shift;
    #return int($number + .5);
    return int($number + 0.5 * ($number <=> 0)); # take care of negative numbers too

}