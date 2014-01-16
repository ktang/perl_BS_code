#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;

my $script = "/Users/tang58/Kai_BS/myown/fisher_p_value/count_CG_CHG_CHH_in_DMC_db.pl";

die unless (-e $script);

print STDERR "script used:\n";
print STDERR $script , "\n\n\n";

my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 <indir> <output>";
die $usage unless(@ARGV == 2);

my $indir = shift or die;
die unless (-d $indir);

my $output = shift or die;
die if(-e $output);

opendir(DIR, $indir) or die;

my @files = grep /_both\.txt$/, readdir DIR;
closedir DIR;

foreach my $file (@files){
	
		my $input = File::Spec->catfile($indir, $file);
		die unless (-e $input);
		my $cmd = "perl $script $input >> $output";
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
