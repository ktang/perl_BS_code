#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;

my $debug = 0;

my $usage = "$0 \n<indir> <outdir> <p_value_cutoff>\n\n";
die $usage unless(@ARGV == 3);

my $indir = shift or die "indir";
my $outdir = shift or die "outdir";
my $p_cutoff = shift or die "p_value_cutoff";

print STDERR "p-value cutoff: $p_cutoff\n\n";
die unless (-d $indir and -d $outdir);

opendir(DIR, $indir) or die;

#my @files = grep/colA_vs\S+diff_bases\S+\.txt/, readdir DIR;
my @files = grep/\S+diff_bases\S+\.txt/, readdir DIR;

closedir(DIR);

print STDERR join("\n", @files), "\n\n";
#colA_vs_ape1l_2A_diff_bases_P005_depth_4_100_both.txt


#JKZ131_Col0_vs_rdd_pool_diff_bases_P0.05_AllMinDep4_ROS1Paper.txt
foreach my $file(@files){
	if($file =~ /(\S+_diff_bases_P)0\.\d+(\S+\.txt)$/ ){
		my $pre = $1;
		my $post = $2;
		my $input = File::Spec->catfile($indir, $file);
		die unless (-e $input);
		
		my $output = File::Spec->catfile($outdir, $pre . $p_cutoff . $post);
		print STDERR $output, "\n\n";
		next if (-e $output);
	}
}
if($debug){
#	print STDERR "deubg: colA_vs_Num77A_diff_bases_P005_depth_4_100_both.txt\n\n";
#	@files = ();
#	@files = ("colA_vs_Num77A_diff_bases_P005_depth_4_100_both.txt");
	exit;
}


foreach my $file(@files){
	
	if($file =~ /(\S+_diff_bases_P)0\.\d+(\S+\.txt)$/ ){
		my $pre = $1;
		my $post = $2;
		my $input = File::Spec->catfile($indir, $file);
		die unless (-e $input);
		
		#my $output = File::Spec->catfile($outdir, $pre . "_P" . $p_cutoff . "_depth_4_100_both.txt");
		my $output = File::Spec->catfile($outdir, $pre . $p_cutoff . $post);

		print STDERR $output, "\n";
		next if (-e $output);
		
		open(IN, $input) or die;
		open(OUT, ">$output") or die;
		while(<IN>){
			chomp;
			my @a = split "\t";
			next if ($a[0] eq "chrC" or $a[0] eq "chrM" or $a[0] eq "pUC19");
			if($a[-1] <= $p_cutoff){
				print OUT $_, "\n";
			}
		}
		close(IN);
		close(OUT);
	}
}
exit;