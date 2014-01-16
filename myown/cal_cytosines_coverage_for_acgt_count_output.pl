#!/usr/bin/perl -w
# given acgt-count output, calculate the depth of coverage at each cytosine position
# and percentage of total cytosines that was covered by at least one or two reads.

#chr1    33      33      CHG:13  0.153846        -
#chr1    79      79      CHH:27  0.148148        -
# 0		  1		 2			3		4			 5

use strict;
use File::Spec;

#my %chr_len = ("chr1"=>30427671, "chr2"=>19698289, "chr3"=>23459830, "chr4"=>18585056, "chr5"=>26975502, "chrM"=>366924, "chrC"=>154478, "pUC19"=>2686);

my $debug = 0;
my $constant = 2000000;
if ($debug) {$constant = 300 ;}

if($debug){
	print  STDERR "debug = $debug\n\n";
}

my $usage = "$0 <dir>  <pre>";
die $usage unless (@ARGV == 2);

my $dir = shift or die "dir";
my $pre = shift or die "pre";

die "wrong dir" unless (-d $dir);

my $in_forw = File::Spec->catfile($dir, $pre . "_forw.txt");
my $in_rev  = File::Spec->catfile($dir, $pre . "_rev.txt");

print  STDERR "input files:\n";
print  STDERR join("\n", ($in_forw, $in_rev)), "\n\n";

die "wrong input" unless (-e $in_forw and -e $in_rev );

my $output = File::Spec->catfile($dir,$pre. "_cytosines_coverage_info.txt");
die "$output exists" if (-e $output);
print STDERR "output:\t$output\n\n";

die if($debug);

#die "die debug\n\n" if ($debug);

open (OUT, ">>$output") or die;

my (%num_forw, %num_rev); #record cytosine number in reference %num{$chr}->{CG}
my (%cover_1W, %cover_2W, %cover_1C, %cover_2C, %cover_depth_forw, %cover_depth_rev);
my (%cytos_forw, %cytos_rev); #record read cytosine depth in each cytosine refrence

#chr1    79      79      CHH:27  0.148148        -
# 0		  1		 2			3		4			 5
open(FORW, $in_forw) or die;
open(REV, $in_rev) or die;

while(<FORW>){
	chomp;
	my @a = split "\t";
	my $chr = $a[0];
	my $per = $a[4];
	my ($type, $depth) = split ":", $a[3];
	my $C_num = round ($depth * $per);
	
#	if ($debug){
#		print STDERR join("\t", ($type, $C_num, $depth)), "\n";
#	}
	
	$num_forw{$chr}->{$type}++;
	$cytos_forw{$chr}->{$type} += $C_num;
	
	$cover_depth_forw{$chr}->{$type} += $depth;
	
	if($depth >= 1){
		$cover_1W{$chr}->{$type}++;
	}
	if($depth >= 2){
		$cover_2W{$chr}->{$type}++;
	}
	
	
}
close(FORW);

while(<REV>){
	chomp;
	my @a = split "\t";
	my $chr = $a[0];
	my $per = $a[4];
	my ($type, $depth) = split ":", $a[3];
	my $C_num = round ($depth * $per);

	$num_rev{$chr}->{$type}++;
	$cytos_rev{$chr}->{$type} += $C_num;
	
	$cover_depth_rev{$chr}->{$type} += $depth;
	
	if($depth >= 1){
		$cover_1C{$chr}->{$type}++;
	}
	if($depth >= 2){
		$cover_2C{$chr}->{$type}++;
	}
	
}
close(REV);

my @types = ("CG", "CHG", "CHH");

print OUT "Coverage >= 1:\n";
#print OUT "On Watson strand:\n";

print OUT join("\t", ("chr", "CG+", "CG-", "CHG+", "CHG-", "CHH+", "CHH-",  "forward", "reverse","total")), "\n";
#my @chrs = ("chr1", "chr2", "chr3", "chr4", "chr5", "chrC", "chrM", "pUC19");
foreach my $chr(sort keys %num_rev){
	print OUT $chr;
#	numerator / denominator

	#as numerator
	my $total = 0;
	my $total_f = 0;
	my $total_r = 0;
	
	# server as denominator
	my $total_num = 0;
	my $total_num_f = 0;
	my $total_num_r = 0;
	
	
	foreach my $type(@types){
		my ($f, $r) = (0, 0);
		
		if(defined $cover_1W{$chr}->{$type}){$f = $cover_1W{$chr}->{$type}}
		if(defined $cover_1C{$chr}->{$type}){$r = $cover_1C{$chr}->{$type}}
		
		$total += ($f + $r);
		$total_f += $f;
		$total_r += $r;
		
		$total_num   += ($num_forw{$chr}->{$type} + $num_rev{$chr}->{$type} );
		$total_num_f += $num_forw{$chr}->{$type} ;
		$total_num_r += $num_rev{$chr}->{$type};
		
		
		print OUT "\t", $f, "/", $num_forw{$chr}->{$type}, "=", sprintf("%.3f",  $f / $num_forw{$chr}->{$type} * 100), "%";
		print OUT "\t", $r, "/", $num_rev{$chr}->{$type}, "=",  sprintf("%.3f", $r / $num_rev{$chr}->{$type} * 100), "%";
		
	}
	print OUT "\t";
	
	print OUT $total_f, "/", $total_num_f, "=", sprintf("%.3f", $total_f  / $total_num_f * 100), "%\t";
	print OUT $total_r, "/", $total_num_r, "=", sprintf("%.3f", $total_r  / $total_num_r * 100), "%\t";
	print OUT $total, "/", $total_num, "=", sprintf("%.3f", $total / $total_num * 100), "%\n";
	
	#print OUT join("\t", ($chr,  , $cover_1W{$chr}/$chr_len{$chr} * 100, "%\n";
}

print OUT "\n";

print OUT "Coverage >= 2:\n";
print OUT join("\t", ("chr", "CG+", "CG-", "CHG+", "CHG-", "CHH+", "CHH-",  "forward", "reverse","total")), "\n";

foreach my $chr(sort keys %num_rev){
	print OUT $chr;
	# server as denominator
	my $total_num = 0;
	my $total_num_f = 0;
	my $total_num_r = 0;
	
	#as numerator
	my $total = 0;
	my $total_f = 0;
	my $total_r = 0;	
	foreach my $type(@types){
		my ($f, $r) = (0, 0);
		
		if(defined $cover_2W{$chr}->{$type}){$f = $cover_2W{$chr}->{$type}}
		if(defined $cover_2C{$chr}->{$type}){$r = $cover_2C{$chr}->{$type}}
		
		$total += ($f + $r);
		$total_f += $f;
		$total_r += $r;
		
		$total_num   += ($num_forw{$chr}->{$type} + $num_rev{$chr}->{$type} );
		$total_num_f += $num_forw{$chr}->{$type} ;
		$total_num_r += $num_rev{$chr}->{$type};
		
		
		print OUT "\t", $f, "/", $num_forw{$chr}->{$type}, "=", sprintf("%.3f",  $f / $num_forw{$chr}->{$type} * 100), "%";
		print OUT "\t", $r, "/", $num_rev{$chr}->{$type}, "=",  sprintf("%.3f", $r / $num_rev{$chr}->{$type} * 100), "%";
		
	}
	print OUT "\t";
	
	print OUT $total_f, "/", $total_num_f, "=", sprintf("%.3f", $total_f  / $total_num_f * 100), "%\t";
	print OUT $total_r, "/", $total_num_r, "=", sprintf("%.3f", $total_r  / $total_num_r * 100), "%\t";
	print OUT $total, "/", $total_num, "=", sprintf("%.3f", $total / $total_num * 100), "%\n";
	

	
	#print OUT join("\t", ($chr,  , $cover_1W{$chr}/$chr_len{$chr} * 100, "%\n";
}

print OUT "\n";

print OUT "Depth_coverage:\n";
print OUT join("\t", ("chr", "CG+", "CG-", "CHG+", "CHG-", "CHH+", "CHH-", "forward", "reverse", "total")), "\n";

foreach my $chr(sort keys %num_rev){
	print OUT $chr;
	# server as denominator
	my $total_num = 0;
	my $total_num_f = 0;
	my $total_num_r = 0;
	
	#as numerator
	my $total = 0;
	my $total_f = 0;
	my $total_r = 0;		
	foreach my $type(@types){
		my ($f, $r) = (0, 0);
		
		if(defined $cover_depth_forw{$chr}->{$type}){$f = $cover_depth_forw{$chr}->{$type}}
		if(defined $cover_depth_rev{$chr}->{$type}){$r = $cover_depth_rev{$chr}->{$type}}
		
		$total += ($f + $r);
		$total_f += $f;
		$total_r += $r;
		
		$total_num   += ($num_forw{$chr}->{$type} + $num_rev{$chr}->{$type} );
		$total_num_f += $num_forw{$chr}->{$type} ;
		$total_num_r += $num_rev{$chr}->{$type};
		
		
		print OUT "\t", $f, "/", $num_forw{$chr}->{$type}, "=", sprintf("%.3f", $f / $num_forw{$chr}->{$type} );
		print OUT "\t", $r, "/", $num_rev{$chr}->{$type}, "=",  sprintf("%.3f", $r / $num_rev{$chr}->{$type} );
		
	}
	print OUT "\t";
	
	print OUT $total_f, "/", $total_num_f, "=", sprintf("%.3f", $total_f / $total_num_f ), "\t";
	print OUT $total_r, "/", $total_num_r, "=", sprintf("%.3f", $total_r / $total_num_r), "\t";
	print OUT $total,   "/", $total_num,   "=", sprintf("%.3f", $total   / $total_num ), "\n";
	
	
	#print OUT join("\t", ($chr,  , $cover_1W{$chr}/$chr_len{$chr} * 100, "%\n";
}

print OUT "\n";

print OUT "Methylation_level:\n";
print OUT join("\t", ("chr", "CG+", "CG-", "CHG+", "CHG-", "CHH+", "CHH-", "forward", "reverse", "total")), "\n";
foreach my $chr(sort keys %num_rev){
	print OUT $chr;
	# server as denominator
	my $total_num = 0;
	my $total_num_f = 0;
	my $total_num_r = 0;
	
	#as numerator
	my $total = 0;
	my $total_f = 0;
	my $total_r = 0;	
	
	foreach my $type(@types){
		my ($f, $r) = (0, 0);
		
		if(defined $cytos_forw{$chr}->{$type}){$f = $cytos_forw{$chr}->{$type}}
		if(defined $cytos_rev{$chr}->{$type}){$r = $cytos_rev{$chr}->{$type}}
		
		$total += ($f + $r);
		$total_f += $f;
		$total_r += $r;
		
		$total_num   += ($cover_depth_forw{$chr}->{$type} + $cover_depth_rev{$chr}->{$type} );
		$total_num_f += $cover_depth_forw{$chr}->{$type} ;
		$total_num_r += $cover_depth_rev{$chr}->{$type};
		
		
		print OUT "\t", $f, "/", $cover_depth_forw{$chr}->{$type}, "=", sprintf("%.3f",  $f / $cover_depth_forw{$chr}->{$type} * 100), "%";
		print OUT "\t", $r, "/", $cover_depth_rev{$chr}->{$type}, "=",  sprintf("%.3f", $r / $cover_depth_rev{$chr}->{$type} * 100), "%";
		
	}
	print OUT "\t";
	
	print OUT $total_f, "/", $total_num_f, "=", sprintf("%.3f", $total_f  / $total_num_f * 100), "%\t";
	print OUT $total_r, "/", $total_num_r, "=", sprintf("%.3f", $total_r  / $total_num_r * 100), "%\t";
	print OUT $total, "/", $total_num, "=", sprintf("%.3f", $total / $total_num * 100), "%\n";
	
	
	#print OUT join("\t", ($chr,  , $cover_1W{$chr}/$chr_len{$chr} * 100, "%\n";
}



close(OUT);

exit;

sub round {
    my ($number) = shift;
    #return int($number + .5);
    return int($number + 0.5 * ($number <=> 0)); # take care of negative numbers too
}