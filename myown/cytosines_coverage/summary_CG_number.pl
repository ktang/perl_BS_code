#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;

my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 <input> <output>";
die $usage unless(@ARGV == 2);

my $input  = shift or die;
my $output = shift or die;

die unless (-e $input);
die if(-e $output);

#Total mC position in 12_pooled
#CG	1715335
#CHG	969188
#CHH	1728047
#Total	4412570


#C:   42859753
#CG:  5567716	12.99%
#CHG: 6093658	14.22%
#CHH: 31198379	72.79%

my $CG  =  5567716;	
my $CHG =  6093658;
my $CHH = 31198379;

open(IN, $input) or die;

if(!$debug){
	open (OUT, ">>$output") or die;
	print OUT join("\t", ("label", "mCG_num", "mCHG_num", "mCHH_num", "Total_mC", "mCG/CG", "mCHG/CHG", "mCHH/CHH")) , "\n";
}else{
	print join("\t", ("label", "mCG_num", "mCHG_num", "mCHH_num", "Total_mC", "mCG/CG", "mCHG/CHG", "mCHH/CHH")) , "\n";
}

my @types = ("CG", "CHG", "CHH");

while(<IN>){
	chomp;
	next unless(/Total\smC\sposition\sin\s(\S+)/);
	my $label = $1;
	my %nums;
	my %nums_per;

	for my $i (1..4){
		my $l = <IN>;
		chomp $l;
		#print STDERR $l, "\t";
		my @a = split "\t", $l;
		#print STDERR join(":", @a), "\n";
		$nums{$a[0]} = $a[1];
		$nums_per{$a[0]} = $a[1];
	}
	
	foreach my $type (@types){
		my $per = sprintf ("%.1f", 100 * $nums{$type} / $nums{"Total"}) . '%';
		$nums_per{$type} .= "($per)";
	}
	
#	print STDERR "\n";
	if(!$debug){
		#print OUT join("\t", ($label, $nums{"CG"}, $nums{"CHG"}, $nums{"CHH"}, $nums{"Total"})) , "\n";
		print OUT join("\t", ($label, $nums_per{"CG"}, $nums_per{"CHG"}, $nums_per{"CHH"}, $nums{"Total"}, 
		 		  sprintf("%.1f", 100 *$nums{"CG"} /$CG) . "%",
		 		  sprintf("%.1f",100 * $nums{"CHG"} /$CHG) . "%",
		 		  sprintf("%.1f", 100 *$nums{"CHH"} /$CHH) . "%"
		 		   )) , "\n";
	}else{
		print join("\t", ($label, $nums_per{"CG"}, $nums_per{"CHG"}, $nums_per{"CHH"}, $nums{"Total"}, 
		 		  sprintf("%.1f", 100 *$nums{"CG"} /$CG) . "%",
		 		  sprintf("%.1f",100 * $nums{"CHG"} /$CHG) . "%",
		 		  sprintf("%.1f", 100 *$nums{"CHH"} /$CHH) . "%"
		 		   )) , "\n";
	}
}


exit;

sub round {
    my($number) = shift;
    #return int($number + .5);
    return int($number + .5 * ($number <=> 0)); # take care of negative numbers too

}
