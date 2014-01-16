#!/usr/bin/perl -w
 
use strict;



#my @a = (3,5,6,1,0);
#my @a = (3,5,6,1, , , , );
my @a = (0);

my $mean = cal_mean(\@a);


print join("\t", @a), "\n\n";
print "mean = $mean\n";


exit;
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
		return ($sum / $num);
	}
}
