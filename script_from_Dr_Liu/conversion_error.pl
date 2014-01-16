#!/usr/bin/perl -w
# calculate conversion error rate
use strict;

my %meth; #chrC->total
foreach my $file(@ARGV){
	open(IN, $file) or die "Can't open $file:$!";
while(<IN>){
	if(/chrC/){
		my @temp = split /\t/;
		$meth{"chrC"}->{"total"} += $temp[3];
		$meth{"chrC"}->{"unconv"}  += $temp[3]*$temp[4];
	}elsif(/chrM/){
		my @temp = split /\t/;
		$meth{"chrM"}->{"total"} += $temp[3];
		$meth{"chrM"}->{"unconv"}  += $temp[3]*$temp[4];
	}elsif(/pUC19/){
        my @temp = split /\t/;
		$meth{"pUC19"}->{"total"} += $temp[3];
		$meth{"pUC19"}->{"unconv"}  += $temp[3]*$temp[4];
	}
}

print "Conversion error rate in $file", "\n";
print "chr\tC\tC+T\tErrorRate\n";
foreach my $chr("chrC", "chrM", "pUC19"){
    
	print join("\t", ($chr, $meth{$chr}->{"unconv"}, $meth{$chr}->{"total"},
	           $meth{$chr}->{"unconv"}/$meth{$chr}->{"total"})), "\n";
}
}
