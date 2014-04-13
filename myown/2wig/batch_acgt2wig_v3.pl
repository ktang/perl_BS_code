#!/usr/bin/perl -w

use strict;
use File::Spec;

print STDERR "\nuse linked files are best \n\n";

my $script = "/Users/tang58/Kai_BS/myown/2wig/acgt2wig_v3.pl";
die unless (-e $script);

print STDERR "script used:\n";
print STDERR $script, "\n\n\n";

my $debug = 0;


my $usage = "$0 \n  <indir> <outdir> \n\n";
die $usage unless(@ARGV == 2);



#my $usage = "$0 \n <isMeth_file> <outdir> <sample_label_in_wig> <wig_pre_name> \n\n";

my $indir = shift or die;
my $outdir = shift or die;
die unless (-d $indir);
die unless (-d $outdir);

opendir (DIR, $indir) or die;
my @files = grep /_forw\.txt$/, readdir DIR;
closedir DIR;

print STDERR "isMeth files:\n";
print STDERR join("\n", @files), "\n\n";

foreach my $file (@files){
		
#		my  $isMeth_file = File::Spec->catfile($indir, $file);
	#	die $isMeth_file unless (-e $isMeth_file);
		my $pre = "";
		if ($file =~ /(\S+)_forw\.txt$/) {
			$pre = $1;#code
		}else{
			die "$file\n";
		}
		
#		my $cmd = "time perl $script  $isMeth_file  $outdir $pre 1 4";
	#	my $cmd = "time perl $script $isMeth_file $outdir  $pre $pre";
       	my $cmd = "time perl $script  $indir $pre $outdir $pre 1 4";

		print STDERR $cmd, "\n\n";
		if(!$debug){
			`$cmd`;
		}
	
}

exit;
