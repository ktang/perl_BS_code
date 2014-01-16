#!/usr/bin/perl -w
 
use strict;

#my @a = (3,5,6,1,0);
#my @a = (3,5,6,1, , , , );
# my @a = ( "chr1_1_200",  "chr1_301_400", "chr1_502_600", "chr1_800_1200", "chr2_1_200");
 my @a = ( "chr1_1_200",  "chr1_301_400", "chr1_502_600", "chr1_800_1200" );
# my @a = ( "chr1_1_200",  "chr1_301_400", "chr1_402_600", "chr1_700_1200" );

#my $mean = cal_mean(\@a);

my @b = ();

merge_list(\@a, \@b, 100);

print join("\n", @a), "\n\n";
print "new:\n";
print join("\n", @b), "\n\n";

exit;

sub merge_list{
	my ($ref_raw, $ref_new, $gap) = @_;
	my ($last_chr, $last_start, $last_end ) = ("chr0", 0, 0);
	my $last_index = scalar( @{$ref_raw} ) - 1;
	
	print STDERR "before merge: ", $last_index + 1, "\n\n";
	#print STDERR "last_index = $last_index\n\n";
	
	for my $i(0..$last_index){
		my ($chr, $start, $end) = split "_", $ref_raw->[$i];
		if($chr ne $last_chr || $start > $last_end + $gap + 1 ){
			if($last_chr ne "chr0") {
				push @{$ref_new}, join("_", ( $last_chr, $last_start, $last_end)) ;
			}
			($last_chr, $last_start, $last_end ) =  ($chr, $start, $end);
		}
		else{
			$last_end = $end;
		}
	
	}
	if($last_chr ne "chr0") {
		push @{$ref_new}, join("_", ( $last_chr, $last_start, $last_end)) ;
	}
	print STDERR "after merge: ", scalar(@{$ref_new}), "\n\n";
}