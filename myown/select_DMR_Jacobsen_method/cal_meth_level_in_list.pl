#!/usr/bin/perl -w

#chr	start	end	DMC_num	wt_CG	wt_CG_per	wt_CHG	wt_CHG_per	wt_CHH	wt_CHH_per	mut_CG	mut_CG_per	mut_CHG	mut_CHG_per	mut_CHH	mut_CHH_per	avge_meth_level_wt	avge_meth_level_mut	Gene	Strand	GeneType	GeneAnnot	TE	TEFamily	Intergenic	Promoter_1kb	CG_num	CHG_num	CHH_num	high_weigth_meth_level	low_weigth_meth_level
#chr1	42319	42670	20	5/9=	55.5556	2/35=	5.7143	68/451=	15.0776	0/8=	0.0000	0/39=	0.0000	4/324=	1.2346	14.1034	0.9181	NONE	NONE	NONE	NONE	NONE	NONE	AT1G01070-AT1G01073	NONE40	15.1448	1.0974
#0        1       2      3    4              6                 8             10              12               14              
use strict;

my $usage = "$0 \n <input> <index4/7> STDOUT\n\n";

die $usage unless (@ARGV == 2);


my $input = shift or die "input";

die unless (-e $input);

my $index = shift or dei; #7;
open(IN, $input) or die;

my $head = <IN>;

my @deps;
my @mCs;

while(<IN>){
	chomp;
	my @a = split "\t";
	
	for my $i(0..5){
		my ($mC, $dep) = get_two($a[$index + 2*$i]);
		$deps[$i] += $dep;
		$mCs[$i]+=  $mC;
	}
	
}

close(IN);

print STDERR "\n";

for my $i(0..5){
	print STDERR "$mCs[$i] / $deps[$i] = ",  sprintf("%.4f", 100 * $mCs[$i] /  $deps[$i]) , "%", "\n";
}

print STDERR "\n\n";

exit;

sub get_two{
	my ($temp) = @_;
	
	die unless ($temp =~/=/);
	$temp =~ s/=//;
	return split "\/", $temp;
}