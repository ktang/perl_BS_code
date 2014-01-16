#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;

my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 \n <in1> STDOUT";
die $usage unless(@ARGV == 1);

my $file = shift or die;
die unless (-e $file);

open(IN, $file) or die;
my $i = 0;

my %nums;

while(<IN>){
	$i++;
	chomp;
	my @a = split "\t";
	my $type = $a[3];
	$nums{$type}++;
}	
close(IN);

die $file unless ($i == $nums{CG} + $nums{CHG} + $nums{CHH});

print  "$file\n";
print  $i, " = " , join(" + ", ($nums{CG} , $nums{CHG},  $nums{CHH} )), "\n\n";


exit;