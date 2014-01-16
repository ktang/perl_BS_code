#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;

#my $script = "/Users/tang58/Kai_BS/myown/call_methylation_CG_CHG_CHH/call_methylation_CG_CHG_CHH_separated_v0.2.pl";
my $script = "/Users/tang58/Kai_BS/myown/call_methylation_CG_CHG_CHH/call_methylation_CG_CHG_CHH_separated_v0.2_chrC.pl";

die unless (-e $script);

my $debug = 0;

#my $p_value = 0.01;

my $max_dep = 178;# 170;

if($debug){
	print STDERR "debug = 1\n\n";
}

my $usage = "$0 \n <input_summary_file> <indir> <outdir> \n\n";
die $usage unless(@ARGV == 3);

my $input_summary_file = shift or die;
my $indir = shift or die;
my $outdir = shift or die;

die unless (-d $indir);
die unless (-d $outdir);

my %errors;
read_summary($input_summary_file, \%errors);


foreach my $pre (sort keys %errors){
		
		my $error_rate = $errors{$pre};
	
#<indir> <pre> <outdir> <meth_max> <error_rate>
		my $cmd = "time perl $script $indir $pre $outdir $max_dep $error_rate";
		print STDERR $cmd, "\n\n";
		if(!$debug){
			`$cmd`;
		}
	
}

exit;

sub round {
    my($number) = shift;
    #return int($number + .5);
    return int($number + .5 * ($number <=> 0)); # take care of negative numbers too

}

sub read_summary{
	my ($file, $ref)  = @_;
	die unless (-e $file);
	open(IN, $file) or die;
	my $head = <IN>;
	while(<IN>){
		chomp;
		my @a = split "\t";
		$ref->{$a[0]} = $a[-1];
	}
	
	close(IN); 
}
