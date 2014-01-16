#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

#die "\n\nincomplete!!!\n\n";

use strict;
use File::Spec;

#my $script = "/Users/tang58/Kai_BS/myown/cytosines_coverage/cal_cytosines_coverage_for_acgt_count_output.pl";
my $script = "/Users/tang58/Kai_BS/for_publish/11_beCalled_count_DMR_number_using_sliding_window.pl";

die unless (-e $script);

my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}

my %flags = ("hyper"=>1, "hypo"=>1);

my $usage = "$0 \n <indir> <outdir> <hyper/hypo>  <postfix(XXXPaper)> \n\n";
die $usage unless(@ARGV == 4);

my $indir = shift or die;
die unless (-e $indir);

my $outdir = shift or die;
die unless (-e $outdir);

my $flag = shift or die;
die "hyper hypo" unless (defined $flags{$flag});

my $postfix = shift or die "postfix";

#my $output = shift or die;

opendir(DIR, $indir) or die;
my @files = grep /$flag\S+\.txt$/, readdir DIR;
closedir DIR;

print STDERR join("\n", @files), "\n\n";

foreach my $file (@files){
	if($file =~ /(\S+)_vs_\S*\.txt$/){
		my $pre = $1;
		my $input = File::Spec->catfile($indir, $file);
		die unless (-e $input);
	#	my $cmd = "perl $script  $input  $pre no >> $output";
	#my $usage = "$0 \n <DMR_list> <outdir> <outpre> <postfix(XXXPaper)>\n\n";

		my $cmd = " perl $script $input $outdir $pre  $postfix";
		print STDERR $cmd, "\n\n";
		if(!$debug){
			`$cmd`;
		}
	}
}

exit;