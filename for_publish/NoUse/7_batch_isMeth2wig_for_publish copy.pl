#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;


print STDERR "used script \n\n /Users/tang58/Kai_BS/myown/2wig/isMeth2wig_v0.2.pl \n\n";

my $debug = 0;

print STDERR "\n\n nodupl tag\n\n";


my $script = "/Users/tang58/Kai_BS/for_publish/isMeth2wig_for_publish_v0.0.pl";
die unless (-e $script);

my $usage = "$0\n <indir> <outdir> <postfix(XXXPaper)>\n\n";
die $usage unless(@ARGV == 3);

my $indir = shift or die;
my $outdir = shift or die;

my $postfix = shift or die;

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



foreach my $pre(sort keys %pres){
	my $file = $pres{$pre};

	my $isMeth_file = File::Spec->catfile($indir, $file);
	die unless (-e $isMeth_file);
	my $cmd = "perl $script $isMeth_file $pre $outdir $postfix";
	print STDERR $cmd, "\n\n";
	if(!$debug){
		`$cmd`;
	}
}

exit;
