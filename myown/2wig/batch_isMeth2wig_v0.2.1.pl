#!/usr/bin/perl -w

use strict;
use File::Spec;


my $debug = 0;

my $dep_cutoff = 4;

print STDERR "\n\n nodupl tag\n\n";

my $script = "/Users/tang58/Kai_BS/myown/2wig/isMeth2wig_v0.2.1.pl";
die unless (-e $script);

print STDERR "used script \n\n $script \n\n";


my $usage = "$0\n <indir> <outdir> \n\n";
die $usage unless(@ARGV == 2);

my $indir = shift or die;
my $outdir = shift or die;
die unless (-d $indir);
die unless (-d $outdir);

opendir(DIR, $indir);
my @files = grep /isMeth.+\.txt$/, readdir DIR;
closedir DIR;


my %pres;

foreach my $file(@files){
	if($file =~ /(\S+)_isMeth/){
		my $pre = $1;
		$pres{$pre} = $file;
	}else{
		die $file;
	}
}

#my $usage = "$0 <isMeth_file> <outdir> <pre> <depth_cutoff>";

foreach my $pre(sort keys %pres){
	my $file = $pres{$pre};
	my $input = File::Spec->catfile( $indir, $file);
	die unless (-e $input);
	my $cmd = "perl $script $input $outdir $pre $dep_cutoff";
	print STDERR $cmd, "\n\n";
	if(!$debug){
		`$cmd`;
	}
}

exit;
