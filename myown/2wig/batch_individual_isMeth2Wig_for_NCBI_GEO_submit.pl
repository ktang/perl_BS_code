#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;

print STDERR "\nuse linked files are best \n\n";

my $script = "/Users/tang58/Kai_BS/myown/2wig/individual_isMeth2Wig_for_NCBI_GEO_submit.pl";
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
my @files = grep /isMeth/, readdir DIR;
closedir DIR;

print STDERR "isMeth files:\n";
print STDERR join("\n", @files), "\n\n";

foreach my $file (@files){
		
		my  $isMeth_file = File::Spec->catfile($indir, $file);
		die $isMeth_file unless (-e $isMeth_file);
		my $pre = "";
		if ($file =~ /(\S+)_isMeth/) {
			$pre = $1;#code
		}else{
			die "$file\n";
		}
		
#		my $cmd = "time perl $script  $isMeth_file  $outdir $pre 1 4";
		my $cmd = "$script $isMeth_file $outdir  $pre $pre";
                #<isMeth_file> <outdir> <sample_label_in_wig> <wig_pre_name>
		print STDERR $cmd, "\n\n";
		if(!$debug){
			`$cmd`;
		}
	
}

exit;