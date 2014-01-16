#!/usr/bin/perl -w

use strict;

#@F[25..26]
#wt_MethLevel_perC       mut_MethLevel_perC

#my $usage = "$0 \n <input> <index4/7> STDOUT\n\n";
#die $usage unless (@ARGV == 2);

#my $index = 25;

#my $index = 23;

my $CG_index = -1;


my $usage = "$0 \n <input> <length_cutoff> <CG_cutoff> <CHG_cutoff> <CHH_cutoff> STDOUT\n\n";
die $usage unless (@ARGV == 5);

my $input = shift or die "input";
#my $cutoff  = shift or die "fold";
my $len_cutoff  = shift or die "len_cutoff";
my $CG_cutoff  = shift or die "CG_cutoff";
my $CHG_cutoff  = shift or die "CHG_cutoff";
my $CHH_cutoff  = shift or die "CHH_cutoff";

die unless (-e $input);

open(IN, $input) or die;

my $head = <IN>;
chomp $head;
my @pts = split "\t", $head;
#die  unless ($pts[$index] =~ /_perC$/ and $pts[$index + 1] =~ /_perC$/ );
#die  unless ($pts[$index] =~ /_C_per$/ and $pts[$index + 1] =~ /_C_per$/ );

for my $i(0..$#pts){
	if($pts[$i] =~ /_CG_per$/ ){
		$CG_index  = $i;
		last;
	}
}

die if($CG_index == -1);




my $CG_index2 = $CG_index + 6;

my $CHG_index = $CG_index + 2;
my $CHG_index2 = $CHG_index + 6;

my $CHH_index = $CG_index + 4;
my $CHH_index2 = $CHH_index + 6;


die  unless ($pts[$CG_index] =~ /_CG_per$/ and $pts[$CG_index2] =~ /_CG_per$/ );
die  unless ($pts[$CHG_index] =~ /_CHG_per$/ and $pts[$CHG_index2] =~ /_CHG_per$/ );
die  unless ($pts[$CHH_index] =~ /_CHH_per$/ and $pts[$CHH_index2] =~ /_CHH_per$/ );



#print join("\t", (@pts, "p_CG", "p_CHG", "p_CHH")), "\n";

print join("\t", @pts), "\n";

while(<IN>){
	chomp;
	my @a = split "\t";
	my $length   = $a[2] -$a[1] +1;
	my ($CG_diff,  $CHG_diff ,$CHH_diff) = (0) x3;
	
	$CG_diff  = $a[$CG_index] - $a[$CG_index2]  if( $a[$CG_index] ne "NA" );
	$CHG_diff = $a[$CHG_index] - $a[$CHG_index2] if($a[$CHG_index] ne "NA" );
	$CHH_diff = $a[$CHH_index] - $a[$CHH_index2] if($a[$CHH_index] ne "NA");
	
#	if( $a[$index] / ($a[$index + 1] + 0.00000001) >= $fold)
#	if( $a[$index] - $a[$index + 1]  >= $cutoff)
	if ( ( $length >= $len_cutoff ) and  (  $CG_diff >= $CG_cutoff or $CHG_diff >= $CHG_cutoff or $CHH_diff >= $CHH_cutoff) ){
		print join("\t", @a), "\n" ;
	}
	
}

close(IN);

exit;