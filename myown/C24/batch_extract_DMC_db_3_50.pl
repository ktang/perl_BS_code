#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;

my $debug = 1;

my $usage = "$0 <indir> <outdir> <p_value_cutoff>";
die $usage unless(@ARGV == 3);

my $indir = shift or die "indir";
my $outdir = shift or die "outdir";
my $p_cutoff = shift or die "p_value_cutoff";

print STDERR "p-value cutoff: $p_cutoff\n\n";
die unless (-d $indir and -d $outdir);

opendir(DIR, $indir) or die;

my @files = grep/\S+diff_bases\S+\.txt/, readdir DIR;
closedir(DIR);

print STDERR join("\n", @files), "\n\n";
#colA_vs_ape1l_2A_diff_bases_P005_depth_3_50_both.txt

if($debug){
#	print STDERR "deubg: colA_vs_Num77A_diff_bases_P005_depth_3_50_both.txt\n\n";
#	@files = ();
#	@files = ("colA_vs_Num77A_diff_bases_P005_depth_3_50_both.txt");
}

#print STDERR join("\n", @files), "\n\n";


foreach my $file(@files){
	if($file =~ /(\S+_diff_bases)_P0\.*\d+(\S+)\.txt/){
		my $pre = $1;
		my $post= $2;
		my $input = File::Spec->catfile($indir, $file);
		die unless (-e $input);
		
		my $output = File::Spec->catfile($outdir, $pre . "_P" . $p_cutoff . $post . ".txt");
		print STDERR $output, "\n";
		next if (-e $output);
	}
}

if($debug){
#	die;
}

foreach my $file(@files){
	if($file =~ /(\S+_diff_bases)_P0\.*\d+(\S+)\.txt/){
		my $pre = $1;
		my $post= $2;
		my $input = File::Spec->catfile($indir, $file);
		die unless (-e $input);
		
		my $output = File::Spec->catfile($outdir, $pre . "_P" . $p_cutoff . $post . ".txt");
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