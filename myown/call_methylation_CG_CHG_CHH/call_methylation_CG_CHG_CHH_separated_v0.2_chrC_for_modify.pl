#!/usr/bin/perl -w

# for modify
# debug to see why cutoff for depth = 2 is NA

#v0.2 from v0.0 directly not from v0.1
# as v0.1 output a lot of 
#print STDERR $_, "\n";
#v0.2 mask them

# calcualte methylation status using Lister et al. 2009 method
# do  separate different position types (CG, CHG, CHH)

# v0.1 
# add strand col

use strict;
use File::Spec;
use Math::CDF qw(:all);
my $FDR = 0.01;
my $print_values = 0;

my $debug = 1;

if($debug){
	print STDERR "\n\ndebug= 1\n\n";
}

my $usage = "$0 <indir> <pre> <outdir> <meth_max> <error_rate>";
die $usage unless (@ARGV == 5);

my $indir = shift or die "indir";
my $pre = shift or die "pre";
my $outdir = shift or die "outdir";
my $meth_max = shift or die "meth_max"; # artificial max depth for nodupl 200;
my $rate = shift or die "rate";         #error rate;

die "wrong indir" unless (-d $indir);
die "wrong outdir" unless (-d $outdir);

my $in_forw = File::Spec->catfile($indir, $pre . "_forw.txt");
my $in_rev  = File::Spec->catfile($indir, $pre . "_rev.txt");

print  STDERR "input files:\n";
print  STDERR join("\n", ($in_forw, $in_rev)), "\n\n";

die "wrong input" unless (-e $in_forw and -e $in_rev );

my $cutoff_file = File::Spec->catfile($outdir,$pre. "_chrC_error_cutoff_separately_called2.txt");
my $output = 	  File::Spec->catfile($outdir,$pre. "_isMeth_chrC_error_separately_called2.txt");

my $depth_file = File::Spec->catfile($outdir,$pre. "_depth_mC.txt");

die "depth_file" if(-e $depth_file);

die "$output exists" if (-e $output);
die "cutoff file exists" if(-e $cutoff_file);
if($debug){
	print STDERR "output:\t$output\n\n";
	exit;
}

#die if($debug);

open(OUT, ">>$output") or die;

#my $real_meth_max = 0;
#my @accu_meth;
#my @total_meth;
#my @meth_level;
#$meth_level[$depth]->[$num_C]++;
#$total_meth[$depth]++;

my ( $real_meth_max_CG, $real_meth_max_CHG, $real_meth_max_CHH )= (0, 0, 0); #not used in nodupl call methylation
my ( @total_meth_CG, @total_meth_CHG, @total_meth_CHH  );					 #		$total_meth_CG[$depth]++;
my ( @meth_level_CG, @meth_level_CHG, @meth_level_CHH  );					 #	$meth_level_CG[$depth]->[$num_C]++;



#chr1    0       0       CHH:0   0       +
#chr1    1       1       CHH:0   0       +
#chr1    2       2       CHH:0   0       +

#last if /chrC/
#print STDERR "max_pos_id is: $max_pos_id\n";


print STDERR "reading $in_forw..\t";
open(MFF, $in_forw) or die "Can't open $in_forw: $!";
while(<MFF>){
#	next if (/SampleID/);
	chomp;
	
	my @a = split "\t";
	my $chr = $a[0];
	my $per = $a[4];
	my ($type, $depth) = split ":", $a[3];
	my $num_C = round ($depth * $per); 
	
	last if ($chr eq "chrC");
	
#	$meth_level[$depth]->[$num_C]++;
#	$total_meth[$depth]++;
#		
#	if($real_meth_max < $depth){
#		$real_meth_max = $depth;
#	}
	if($type eq "CG"){
		$meth_level_CG[$depth]->[$num_C]++;
		$total_meth_CG[$depth]++;
		if($real_meth_max_CG < $depth){
			$real_meth_max_CG = $depth;
		}
	}
	
	elsif( $type eq "CHG" ){
		$meth_level_CHG[$depth]->[$num_C]++;
		$total_meth_CHG[$depth]++;
		if($real_meth_max_CHG < $depth){
			$real_meth_max_CHG = $depth;
		}
	}
	
	elsif( $type eq "CHH" ){
		$meth_level_CHH[$depth]->[$num_C]++;
		$total_meth_CHH[$depth]++;
		if($real_meth_max_CHH < $depth){
			$real_meth_max_CHH = $depth;
		}
	}else{
		die $_;
	}
	
}
close MFF;
print STDERR "Done\n";


print STDERR "reading $in_rev..\t";

open(MFR, $in_rev) or die;
while(<MFR>){
	chomp;

	my @a = split "\t";
	my $chr = $a[0];
	my $per = $a[4];
	my ($type, $depth) = split ":", $a[3];
	my $num_C = round ($depth * $per); 
	
	last if ($chr eq "chrC");
	
	if($type eq "CG"){
		$meth_level_CG[$depth]->[$num_C]++;
		$total_meth_CG[$depth]++;
		if($real_meth_max_CG < $depth){
			$real_meth_max_CG = $depth;
		}
	}
	
	elsif( $type eq "CHG" ){
		$meth_level_CHG[$depth]->[$num_C]++;
		$total_meth_CHG[$depth]++;
		if($real_meth_max_CHG < $depth){
			$real_meth_max_CHG = $depth;
		}
	}
	
	elsif( $type eq "CHH" ){
		$meth_level_CHH[$depth]->[$num_C]++;
		$total_meth_CHH[$depth]++;
		if($real_meth_max_CHH < $depth){
			$real_meth_max_CHH = $depth;
		}
	}else{
		die $_;
	}
	
}
close(MFR);
print STDERR "Done\n\n";

print STDERR "File: $pre, Conversion Error: $rate\n";			
print STDERR "Max depth: $meth_max\n";


my $real_max = max($real_meth_max_CG, $real_meth_max_CHG, $real_meth_max_CHH );


output_depth($depth_file, \@meth_level_CG, \@meth_level_CHG, \@meth_level_CHH , $real_max );

my ( @accu_meth_CG,  @accu_meth_CHG,  @accu_meth_CHH    );

cal_accumulated_number(\@total_meth_CG ,  \@meth_level_CG, \@accu_meth_CG);
#cal_accumulated_number(\@total_meth_CG , \@meth_level_CHG, \@accu_meth_CHG);
cal_accumulated_number(\@total_meth_CHG , \@meth_level_CHG, \@accu_meth_CHG);
cal_accumulated_number(\@total_meth_CHH , \@meth_level_CHH, \@accu_meth_CHH);



undef @meth_level_CG;# empty memory
undef @meth_level_CHG;# empty memory
undef @meth_level_CHH;# empty memory

#	my( $total_meth_ref, $accu_meth_ref, $cutoff_ref ) = @_;
my ( @cutoff_CG, @cutoff_CHG, @cutoff_CHH );
cal_cutoff(\@total_meth_CG,  \@accu_meth_CG,  \@cutoff_CG);
cal_cutoff(\@total_meth_CHG, \@accu_meth_CHG, \@cutoff_CHG);
cal_cutoff(\@total_meth_CHH, \@accu_meth_CHH, \@cutoff_CHH);

undef @total_meth_CG;
undef @total_meth_CHG;
undef @total_meth_CHH;



open (CUT, ">>$cutoff_file") or die "cannot open $cutoff_file: $!";

#print CUT "cutoff for $pre:\n";
print CUT join("\t", ($pre."_depth", "CG_cutoff", "CHG_cutoff", "CHH_cutoff")),"\n";

for my $index (0..$meth_max){
#	if (!defined $cutoff[$index]){
#		$cutoff[$index] = 0;
#	}

	my ( $cut_CG, $cut_CHG, $cut_CHH ) = ("NA", "NA", "NA");
	
	if(defined $cutoff_CG[$index] ) { $cut_CG  =  $cutoff_CG[$index]  }
	if(defined $cutoff_CHG[$index] ) {$cut_CHG =  $cutoff_CHG[$index]  }
	if(defined $cutoff_CHH[$index] ) {$cut_CHH =  $cutoff_CHH[$index]  }
	print CUT join("\t", ($index , $cut_CG, $cut_CHG, $cut_CHH )), "\n";
}

close(CUT);


my %num_pos; #
#chr1    0       0       CHH:0   0       +
#chr1    1       1       CHH:0   0       +
#chr1    2       2       CHH:0   0       +
print OUT join("\t", ("chr", "pos", "strand", "type", "num_C", "depth", "percentage", "isMeth")),"\n";


print STDERR "reading $in_forw..\t";

open(MFF, $in_forw) or die "Can't open $in_forw: $!";
while(<MFF>){
	chomp;
	
	my @a = split "\t";
	my $chr = $a[0];
	my $pos = $a[1] + 1;
	my $per = $a[4];
	my $strand = $a[5];
	my ($type, $depth) = split ":", $a[3];
	my $num_C = round ($depth * $per); 
	
	last if ($chr eq "chrC");
	
	my $isMeth = 0;
	
	if($depth > $meth_max){
		print STDERR "depth > max: $_\n";
		if( $per > $rate +  0.05 ){
			$isMeth = 1;
			$num_pos{$type}++;
		}
	}
	else{
				
		if( $type eq "CG" ){
			
			if($cutoff_CG[$depth] eq "NA"){
				print STDERR $_, "\n";
				if($per > $rate +  0.05){
					$isMeth = 1;
					$num_pos{$type}++;
				}
			}elsif( $num_C >= $cutoff_CG[$depth] ){
				$isMeth = 1;
				$num_pos{$type}++;
			}
		}
		elsif( $type eq "CHG"){
			if($cutoff_CHG[$depth] eq "NA"){
				print STDERR $_, "\n";
				if($per > $rate +  0.05){
					$isMeth = 1;
					$num_pos{$type}++;
				}
			}elsif( $num_C >= $cutoff_CHG[$depth] ){
				$isMeth = 1;
				$num_pos{$type}++;
			}
		}
		elsif($type eq "CHH" ){
			if($cutoff_CHH[$depth] eq "NA"){
				print STDERR $_, "\n";
				if($per > $rate +  0.05){
					$isMeth = 1;
					$num_pos{$type}++;
				}
			}elsif( $num_C >= $cutoff_CHH[$depth] ){
				$isMeth = 1;
				$num_pos{$type}++;
			}
		}
		else{
			die $_;	
		}
	}
	
	print OUT join("\t", ($chr, $pos, $strand, $type, $num_C, $depth, $per, $isMeth)), "\n";
	
}
close (MFF);
print STDERR "Done\n";

print STDERR "reading $in_rev..\t";

open(MFR, $in_rev) or die;
while(<MFR>){
	chomp;
	
	my @a = split "\t";
	my $chr = $a[0];
	my $pos = $a[1] + 1;
	my $per = $a[4];
	my $strand = $a[5];
	my ($type, $depth) = split ":", $a[3];
	my $num_C = round ($depth * $per); 
	
	last if ($chr eq "chrC");
	my $isMeth = 0;
	if($depth > $meth_max){
		if( $per > $rate +  0.05 ){
			$isMeth = 1;
			$num_pos{$type}++;
		}
	}
	else{
				
		if( $type eq "CG" ){
			
			if($cutoff_CG[$depth] eq "NA"){
	#			print STDERR $_, "\n";
				if($per > $rate +  0.05){
					$isMeth = 1;
					$num_pos{$type}++;
				}
			}elsif( $num_C >= $cutoff_CG[$depth] ){
				$isMeth = 1;
				$num_pos{$type}++;
			}
		}
		elsif( $type eq "CHG"){
			if($cutoff_CHG[$depth] eq "NA"){
		#		print STDERR $_, "\n";
				if($per > $rate +  0.05){
					$isMeth = 1;
					$num_pos{$type}++;
				}
			}elsif( $num_C >= $cutoff_CHG[$depth] ){
				$isMeth = 1;
				$num_pos{$type}++;
			}
		}
		elsif($type eq "CHH" ){
			if($cutoff_CHH[$depth] eq "NA"){
#				print STDERR $_, "\n";
				if($per > $rate +  0.05){
					$isMeth = 1;
					$num_pos{$type}++;
				}
			}elsif( $num_C >= $cutoff_CHH[$depth] ){
				$isMeth = 1;
				$num_pos{$type}++;
			}
		}
		else{
			die $_;	
		}
	}
	print OUT join("\t", ($chr, $pos,$strand, $type, $num_C, $depth, $per, $isMeth)), "\n";

}
close(MFR);
print STDERR "Done\n\n";


print STDERR "Total mC position in ", $pre, "\n";
my $total = 0;
foreach my $t(sort keys %num_pos){
    print STDERR join("\t",( $t,  $num_pos{$t}) ), "\n";
	$total += $num_pos{$t};
}
print STDERR "Total\t$total\n";

exit;

sub max{
	my @a = @_;
	my @b = sort {$b <=> $a} @a;
	return $b[0];
}

#output_depth($depth_file, \@total_meth_CG, \@total_meth_CG, \@total_meth_CHH );
#$meth_level_CHH[$depth]->[$num_C]
sub output_depth{
	my ($file, $CG_ref, $CHG_ref, $CHH_ref, $max) = @_;
	
	die if(-e $file);
	open(DEP, ">>$file") or die;
	print DEP join("\t", ("depth", "mC", "CG", "CHG", "CHH")), "\n";
	for my $i(0..$max){
		for my $j(0..$i){
			my ($num_cg, $num_chg, $num_chh) = (0, 0, 0);
			if(defined $CG_ref ->[$i]->[$j])  {$num_cg = $CG_ref->[$i]->[$j]}
			if(defined $CHG_ref->[$i]->[$j])  {$num_chg = $CHG_ref->[$i]->[$j]}
			if(defined $CHH_ref->[$i]->[$j])  {$num_chh = $CHH_ref->[$i]->[$j]}
			if( $num_cg + $num_chg + $num_chh  > 0){
				print DEP join("\t", ( $i, $j, $num_cg, $num_chg, $num_chh )), "\n";
			}
		}
	}
	
	close(DEP);
}


#2/2 YES or NO??
sub	cal_cutoff{
	my( $total_meth_ref, $accu_meth_ref, $cutoff_ref ) = @_;
	print STDERR "find cutoff at each depth level..\t";
	# find cutoff at each depth level
	foreach my $i(0..$meth_max){
		next unless(defined $total_meth_ref->[$i]);
		if($i < 2){
			$cutoff_ref->[$i] = $i+1;
			#next;
		}else{
			$cutoff_ref->[$i] = "NA";
			my $start = int ($i * $rate + 0.999999);
			my $j = $start;
			
		#	if($i == 2){
			#	print STDERR "start..i = $start, $i\n\n";
		#	}
			
		#	foreach $j($start..$i){
			while($j <= $i){	
				#if($i == 2){
				#	print STDERR "in body : j = $j\n\n";
				#}
				
				my $prob = pbinom($j, $i, $rate);
				$prob = $prob - pbinom($j-1, $i, $rate);
				my ($num_unC, $num_mC) = (0,0);
				$num_unC = $accu_meth_ref->[$i]->[$j-1];
				if(!defined $num_unC){
					$num_unC = 0;
				}
				$num_mC = $total_meth_ref->[$i] - $num_unC;
				my $left = $prob * $num_unC;
				my $right = $FDR * $num_mC;
				
				#if($i == 2){
				#	print STDERR "left:	$left = $prob * $num_unC; \n\n";
				#	print STDERR "right:	$right = $FDR * $num_mC\n\n";
				#}
				
				if($left < $right){
				#	if($i == 2){
				#		print STDERR "assign: $j\n\n";
				#	}
					$cutoff_ref->[$i] = $j;
					last;
				}
				
				$j++;
			}
			
		#	if($i == 2){
		#		print STDERR "last: j = $j\n\n";
		#	}
			
			
			#if($j > $i){
			if($cutoff_ref->[$i] eq "NA"){
				$cutoff_ref->[$i] = $i + 1;
		#		if($i == 2){
		#			print STDERR "assign here as last: $j\n\n";
		#		}
			}
			#}
			
		}
	}
	print STDERR "Done\n";
}

sub round {
    my ($number) = shift;
    #return int($number + .5);
    return int($number + 0.5 * ($number <=> 0)); # take care of negative numbers too
}



sub cal_accumulated_number{#cal_accumulated_number(\@total_meth_CHH , \@meth_level_CHH, \@accu_meth_CHH);
	my($total_meth_ref, $meth_level_ref, $accu_meth_ref) = @_;
	
	print STDERR "calculate accumulated number..\t";
	# calculate accumulated number of positions with <= k unconverted Cs
		foreach my $i(1..$meth_max){
			next unless (defined $total_meth_ref->[$i]);

			if(!defined $meth_level_ref->[$i]->[0]){
				$meth_level_ref->[$i]->[0] = 0;
			}
			$accu_meth_ref->[$i]->[0] = $meth_level_ref->[$i]->[0];
			foreach my $j(1..$i){
				if(!defined $meth_level_ref->[$i]->[$j]){
						$accu_meth_ref->[$i]->[$j] = $accu_meth_ref->[$i]->[$j-1]			
				}
				else{
						$accu_meth_ref->[$i]->[$j] = $accu_meth_ref->[$i]->[$j-1] +
													$meth_level_ref->[$i]->[$j];
				}
			}
		}
	print STDERR "Done\n\n";
	
	print STDERR "depth = 2 : $total_meth_ref->[2] \n";
	print STDERR "0: ",$meth_level_ref->[2]->[0] ,"\t1: ",$meth_level_ref->[2]->[1], "\t2: ", $meth_level_ref->[2]->[2], "\n";
	print STDERR "accum:0: ",$accu_meth_ref->[2]->[0] ,"\t1: ",$accu_meth_ref->[2]->[1], "\t2: ", $accu_meth_ref->[2]->[2], "\n\n";

}