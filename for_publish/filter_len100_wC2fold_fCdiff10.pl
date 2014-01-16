#!/usr/bin/perl -w

use strict;


#my $usage = "$0 \n <input> <length_cutoff> <CG_cutoff> <CHG_cutoff> <CHH_cutoff> STDOUT\n\n";
#die $usage unless (@ARGV == 5);


my $usage = "$0 \n <input> <length_cutoff> <fold> <fC_diff> STDOUT\n\n";
die $usage unless (@ARGV == 4);



my $input = shift or die "input";
#my $cutoff  = shift or die "fold";
my $len_cutoff  = shift or die "len_cutoff";

my $fold_cutoff = shift or die "fold";
my $fC_diff_cutoff = shift or die "fC_diff";

#my $CG_cutoff  = shift or die "CG_cutoff";
#my $CHG_cutoff  = shift or die "CHG_cutoff";
#my $CHH_cutoff  = shift or die "CHH_cutoff";

die unless (-e $input);

open(IN, $input) or die;

my $head = <IN>;
chomp $head;
my @pts = split "\t", $head;
#die  unless ($pts[$index] =~ /_perC$/ and $pts[$index + 1] =~ /_perC$/ );
#die  unless ($pts[$index] =~ /_C_per$/ and $pts[$index + 1] =~ /_C_per$/ );

my $wmC_index = -1;

for my $i(0..$#pts){
	if($pts[$i] =~ /^wmC_/ ){
	#	$CG_index  = $i;
		$wmC_index = $i;
		last;
	}
}

die if($wmC_index == -1);




my $wmC_index2 = $wmC_index + 4;

my $fmC_index = $wmC_index + 8;
my $fmC_index2 = $fmC_index + 4;




#die  unless ($pts[$CG_index] =~ /_CG_per$/ and $pts[$CG_index2] =~ /_CG_per$/ );
#die  unless ($pts[$CHG_index] =~ /_CHG_per$/ and $pts[$CHG_index2] =~ /_CHG_per$/ );
#die  unless ($pts[$CHH_index] =~ /_CHH_per$/ and $pts[$CHH_index2] =~ /_CHH_per$/ );

die  unless ($pts[$wmC_index] =~ /^wmC_/ and $pts[$wmC_index2] =~ /^wmC_/ );
die  unless ($pts[$fmC_index] =~ /^fmC_/ and $pts[$fmC_index2] =~ /^fmC_/ );

#print join("\t", (@pts, "p_CG", "p_CHG", "p_CHH")), "\n";

print join("\t", @pts), "\n";

while(<IN>){
	chomp;
	my @a = split "\t";
	my $length   = $a[2] -$a[1] +1;
	#my ($CG_diff,  $CHG_diff ,$CHH_diff) = (-100) x3;
	
#	$CG_diff  = $a[$CG_index] - $a[$CG_index2]  if( $a[$CG_index] ne "NA" );
#	$CHG_diff = $a[$CHG_index] - $a[$CHG_index2] if($a[$CHG_index] ne "NA" );
#	$CHH_diff = $a[$CHH_index] - $a[$CHH_index2] if($a[$CHH_index] ne "NA");

	my $fold = -100;
	my $fC_diff = -100;
	
	$fold = $a[$wmC_index] / ($a[$wmC_index2] + 0.00001)  unless ( $a [$wmC_index] =~ /NA|NONE/);
	$fC_diff = $a[$fmC_index ] - $a[$fmC_index2];
	
#	if( $a[$index] / ($a[$index + 1] + 0.00000001) >= $fold)
#	if( $a[$index] - $a[$index + 1]  >= $cutoff)
	#if ( ( $length >= $len_cutoff ) and  (  $CG_diff >= $CG_cutoff or $CHG_diff >= $CHG_cutoff or $CHH_diff >= $CHH_cutoff) ){
	if( ( $length >= $len_cutoff ) and  ( $fold >=  $fold_cutoff) and   ( $fC_diff >= $fC_diff_cutoff) ) {
		print join("\t", @a), "\n" ;
	}
}

close(IN);

exit;