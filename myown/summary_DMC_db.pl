#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;

my $debug = 0;

my $usage = "$0 <indir>  <output>";
die $usage unless(@ARGV == 2);

my $indir = shift or die "indir";
#my $outdir = shift or die "outdir";
#my $p_cutoff = shift or die "p_value_cutoff";

my $output = shift or die "output";
die if (-e $output);

die unless (-d $indir );

opendir(DIR, $indir) or die;

my @files = grep/colA_vs\S+diff_bases\S+\.txt/, readdir DIR;
closedir(DIR);

print STDERR join("\n", @files), "\n\n";
#colA_vs_ape1l_2A_diff_bases_P005_depth_4_100_both.txt

if($debug){
	print STDERR "deubg: colA_vs_Num77A_diff_bases_P005_depth_4_100_both.txt\n\n";
#	@files = ();
	@files = ("colA_vs_Num77A_diff_bases_P005_depth_4_100_both.txt");
}

#print STDERR join("\n", @files), "\n\n";


foreach my $file(@files){
	if($file =~ /(colA_vs_\S+_diff_bases)_P0\S+_depth/){
		my $pre = $1;
		my $input = File::Spec->catfile($indir, $file);
		die unless (-e $input);
	}else{
		die $file, "\n";
	}
}

open (OUT, ">$output") or die;
print OUT join("\t", ("sample" , "total_DMC", "hyper_DMC", "hyper_CG",  "hyper_CHG", "hyper_CHH",  "hypo_DMC", "hypo_CG",  "hypo_CHG", "hypo_CHH" )), "\n";

													#wt								mut
#chr1    108     +       CHG     9       11      0.818182        1       8       0.125   0.00547749464158133
#  0	  1		 2		  3		 4		  5			6			 7		 8			9			10
foreach my $file(@files){
	if($file =~ /colA_vs_(\S+)_diff_bases_P0\S+_depth/){
		my $pre = $1;
		my $input = File::Spec->catfile($indir, $file);
		die unless (-e $input);
		
		my (%hypers, %hypos) ;
		my ($hyper, $hypo) = (0,0);
		
		open(IN, $input) or die;
		while(<IN>){
			chomp;
			my @a = split "\t";
			next if ($a[0] eq "chrC" or $a[0] eq "chrM" or $a[0] eq "pUC19");
			if($a[6] < $a[9]){
				$hypers{$a[3]}++;
				$hyper++;
			}
			else{
				$hypos{$a[3]}++;
				$hypo++;
			}
		}
		close(IN);
		print OUT join("\t", ($pre,$hyper + $hypo ,$hyper, $hypers{CG}, $hypers{CHG}, $hypers{CHH}, $hypo, $hypos{CG}, $hypos{CHG}, $hypos{CHH})), "\n";
	}
}


close(OUT);
exit;