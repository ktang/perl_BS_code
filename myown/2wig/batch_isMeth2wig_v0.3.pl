#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;


print STDERR "used script \n\n /Users/tang58/Kai_BS/myown/2wig/isMeth2wig_v0.3.pl \n\n";

my $debug = 0;

#	print STDERR "\n\n nodupl tag\n\n";


my $script = "/Users/tang58/Kai_BS/myown/2wig/isMeth2wig_v0.3.pl";
die unless (-e $script);

my $usage = "$0 <indir> <outdir>";
die $usage unless(@ARGV == 2);

my $indir = shift or die;
my $outdir = shift or die;
die unless (-d $indir);
die unless (-d $outdir);

opendir(DIR, $indir);
my @files = grep /isMeth.+\.txt$/, readdir DIR;
closedir (DIR);


my %pres;

foreach my $file(@files){
	if($file =~ /(\S+)_isMeth/){
		my $pre = $1;
		$pres{$pre} = $file;
	}else{
		die $file;
	}
}


#my $dir = "/Volumes/My_Book/20120427_ShangHai_data/handled_data/cytosine_coverage_info";

foreach my $pre(sort keys %pres){
	my $file = $pres{$pre};
#	my $cmd = "time perl $script $indir $outdir $pre $file";
	#my $cmd = "time perl $script $indir $outdir $pre" . "_nodupl". " $file";
	my $cmd = "perl $script $indir $file $outdir $pre" . "_NonZero";
	print STDERR "\n", $cmd, "\n\n";
	if(!$debug){
		`$cmd`;
	}
}

exit;
