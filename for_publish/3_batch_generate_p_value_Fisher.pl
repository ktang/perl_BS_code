#!/usr/bin/perl -w

use strict;
use File::Spec;

# my $script = "/Users/tang58/Kai_BS/myown/fisher_p_value/generate_p_value_Fisher_for_brat_bw_acgt_count_output_v0.2.pl";

my $script = "/Users/tang58/Kai_BS/for_publish/generate_p_value_Fisher_for_brat_bw_acgt_count_output_v1.1.pl";
die unless (-e $script);

my $debug = 0;

print STDERR "\n\n   WT_label must be complete not part!!!!!!!!!!!!\n\n   ";

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 \n <indir> <outdir> <WT_label> <postfix(edm2_paper)> \n\n";
die $usage unless(@ARGV == 4);

my $indir = shift or die;
die unless (-d $indir);

my $outdir = shift or die;
die unless (-d $outdir);

my $wt_label = shift or die;

my $postfix = shift or die "postfix";


opendir(DIR, $indir) or die;
#my @files = grep /_rev\.txt$/, readdir DIR;
my @files = grep {!/$wt_label/} (grep /isMeth\S+\.txt$/, readdir DIR);
closedir DIR;

opendir(DIR, $indir) or die;
my ($wt_name) = grep {/$wt_label/}  readdir DIR;
closedir DIR;

print STDERR "mut:\n";
print STDERR join("\n", @files), "\n\n\n";

print STDERR "wt: $wt_name \n\n\n";

my $wt_file = File::Spec->catfile($indir, $wt_name);
die unless (-e $wt_file);

#<WT_isMeth> <WT_pre> <mut_isMeth> <mut_pre>  <outdir>
foreach my $file (@files){
	if($file =~ /(\S+)_isMeth\S+\.txt$/){
	
	
		my $pre = $1;
		#next if ($pre =~ /col/);
		my $mut_file = File::Spec->catfile($indir, $file);
		die unless (-e $mut_file);

		#my $usage = "$0 <dir> <WT_pre> <mut_pre> <outdir>";
		my $cmd = "time perl $script $wt_file $wt_label $mut_file $pre $outdir $postfix";
		print STDERR $cmd, "\n\n";
		if(!$debug){
			`$cmd`;
		}
	}
}

exit;

