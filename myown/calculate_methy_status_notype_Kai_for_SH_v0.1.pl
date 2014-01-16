#!/usr/bin/perl -w
# calcualte methylation status using Lister et al. 2009 method
# do not separate different position types (CG, CHG, CHH)

# v0.1 
# add strand col

use strict;
use File::Spec;
use Math::CDF qw(:all);
my $FDR = 0.01;
my $print_values = 0;

my $debug = 0;

if($debug){
	print STDERR "\n\ndebug= 1\n\n\n";
}

#my $usage = "$0 <meth file> <error rate> <cutoff_file> ";
#die $usage unless(@ARGV >= 3);
#my ($meth_file, $rate) = @ARGV[0..1];
#my $cutoff_file = $ARGV[2];

my $usage = "$0 <dir> <pre> <meth_max> <error_rate>";
die $usage unless (@ARGV == 4);

my $dir = shift or die "dir";
my $pre = shift or die "pre";
#my $outfile = shift or die "outfile";
#my $cutoff_file = shift or die "cutoff_file";
my $meth_max = shift or die "meth_max";
my $rate = shift or die "rate";

die "wrong dir" unless (-d $dir);

my $in_forw = File::Spec->catfile($dir, $pre . "_forw.txt");
my $in_rev  = File::Spec->catfile($dir, $pre . "_rev.txt");



print  STDERR "input files:\n";
print  STDERR join("\n", ($in_forw, $in_rev)), "\n\n";

die "wrong input" unless (-e $in_forw and -e $in_rev );

my $cutoff_file = File::Spec->catfile($dir,$pre. "_cutoff.txt");

#my $outdir = "/Volumes/My_Book/20120427_ShangHai_data/call_methylation/1_isMeth";
#die unless (-e $outdir);
#my $output = File::Spec->catfile($outdir,$pre. "_isMeth.txt");
my $output = File::Spec->catfile($dir,$pre. "_isMeth.txt");

if($debug){
	$output = File::Spec->catfile($dir,$pre. "_isMeth.txt");
}


die "$output exists" if (-e $output);
die "cutoff file exists" if(-e $cutoff_file);
print STDERR "output:\t$output\n\n";

open(OUT, ">$output") or die;

my $real_meth_max = 0;
my @accu_meth;
my @total_meth;
my @meth_level;
#my $max_pos_id = 42859753; #don't want to include chrC, chrM, pUPC19
#if(@ARGV > 3){
#	$max_pos_id = $ARGV[3];
#}

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
#	my ($sample_id, $pos_id, $depth, $num_C, $percent, $type) = split /\t/;
	#next if($pos_id > $max_pos_id);
	
	my @a = split "\t";
	my $chr = $a[0];
	my $per = $a[4];
	my ($type, $depth) = split ":", $a[3];
	my $num_C = round ($depth * $per); 
	
	last if ($chr eq "chrC");
	
	$meth_level[$depth]->[$num_C]++;
	$total_meth[$depth]++;
		
	if($real_meth_max < $depth){
		$real_meth_max = $depth;
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
	
	$meth_level[$depth]->[$num_C]++;
	$total_meth[$depth]++;
		
	if($real_meth_max < $depth){
		$real_meth_max = $depth;
	}

}
close(MFR);
print STDERR "Done\n\n";

print STDERR "File: $pre, Conversion Error: $rate\n";			
print STDERR "Max depth: $meth_max\n";

print STDERR "calculate accumulated number..\t";
# calculate accumulated number of positions with <= k unconverted Cs
	foreach my $i(0..$meth_max){
		next if($i == 0 or !defined $total_meth[$i]);

		if(!defined $meth_level[$i]->[0]){
			$meth_level[$i]->[0] = 0;
		}
		$accu_meth[$i]->[0] = $meth_level[$i]->[0];
		foreach my $j(1..$i){
			if(!defined $meth_level[$i]->[$j]){
					$accu_meth[$i]->[$j] = $accu_meth[$i]->[$j-1]			
			}
			else{
					$accu_meth[$i]->[$j] = $accu_meth[$i]->[$j-1] +
											$meth_level[$i]->[$j];
			}
		}
	}
print STDERR "Done\n\n";

print STDERR "depth = 2 : $total_meth[2]\n";
print STDERR "0: ",$meth_level[2]->[0] ,"\t1: ",$meth_level[2]->[1], "\t2: ", $meth_level[2]->[2], "\n";
print STDERR "accum:0: ",$accu_meth[2]->[0] ,"\t1: ",$accu_meth[2]->[1], "\t2: ", $accu_meth[2]->[2], "\n\n";

undef @meth_level;# empty memory

print STDERR "Empth meth_leverl\n\nfind cutoff at each depth level..\t";

# find cutoff at each depth level
my @cutoff;
	foreach my $i(0..$meth_max){
		next if(!defined $total_meth[$i]);
		if($i < 2){
			$cutoff[$i] = $i+1;
			#next;
		}else{
			$cutoff[$i] = 4;
			
			my $start = int ($i * $rate + 0.999999);
			foreach my $j($start..$i){
			
				my $prob = pbinom($j, $i, $rate);
				$prob = $prob - pbinom($j-1, $i, $rate);
				my ($num_unC, $num_mC) = (0,0);
				$num_unC = $accu_meth[$i]->[$j-1];
				if(!defined $num_unC){
					$num_unC = 0;
				}
				$num_mC = $total_meth[$i] - $num_unC;
				my $left = $prob * $num_unC;
				my $right = $FDR * $num_mC;
				if($left < $right){
					$cutoff[$i] = $j;
					last;
				}
			}
		}
	}
print STDERR "Done\n";
undef @total_meth;
print STDERR "Empty total_meth\n\n";

open (CUT, ">>$cutoff_file") or die "cannot open $cutoff_file: $!";

print CUT "cutoff for $pre:\n";
print CUT join("\t", ("depth", "cutoff")),"\n";

for my $index (0..$meth_max){
#	if (!defined $cutoff[$index]){
#		$cutoff[$index] = 0;
#	}
	print CUT join("\t", ($index , $cutoff[$index])), "\n" if (defined $cutoff[$index]);
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
		if( $type eq "CG" and $per >= 0.5 ){
			$isMeth = 1;
			$num_pos{$type}++;
		}elsif( $type eq "CHG" and $per >= 0.25){
			$isMeth = 1;
			$num_pos{$type}++;
		}elsif($type eq "CHH" and $per >= 0.15){
			$isMeth = 1;
			$num_pos{$type}++;
		}
	}
	else{
		if($depth >= 2 && $num_C >= 1 && $num_C >= $cutoff[$depth]){
			$isMeth = 1;
			$num_pos{$type}++;
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
		if( $type eq "CG" and $per >= 0.5 ){
			$isMeth = 1;
			$num_pos{$type}++;
		}elsif( $type eq "CHG" and $per >= 0.25){
			$isMeth = 1;
			$num_pos{$type}++;
		}elsif($type eq "CHH" and $per >= 0.15){
			$isMeth = 1;
			$num_pos{$type}++;
		}
	}
	else{
		if($depth >= 2 && $num_C >= 1 && $num_C >= $cutoff[$depth]){
			$isMeth = 1;
			$num_pos{$type}++;
		}
	}
	print OUT join("\t", ($chr, $pos,$strand, $type, $num_C, $depth, $per, $isMeth)), "\n";

}
close(MFR);
print STDERR "Done\n\n";


print STDERR "Total mC position in ", $pre, "\n";
my $total = 0;
foreach my $t(sort keys %num_pos){
    print STDERR "$t:", $num_pos{$t}, "\n";
	$total += $num_pos{$t};
}
print STDERR "Total: $total\n";

exit;

sub round {
    my ($number) = shift;
    #return int($number + .5);
    return int($number + 0.5 * ($number <=> 0)); # take care of negative numbers too
}