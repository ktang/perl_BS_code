#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);


#input outdir pre p-value
use strict;
use File::Spec;

my $debug = 0;

my $p_value = 0.01;

my $usage = "$0 <indir>  <outdir>";
die $usage unless(@ARGV == 2);

my $script = "/Users/tang58/Kai_BS/myown/select_DMR_output_list.pl";
die unless (-e $script);

my $indir = shift or die "indir";
#my $outdir = shift or die "outdir";
#my $p_cutoff = shift or die "p_value_cutoff";

my $outdir = shift or die "outdir";

die unless (-d $indir and -d $outdir);

opendir(DIR, $indir) or die;

my @files = grep/colA_vs\S+diff_bases\S+\.txt/, readdir DIR;
closedir(DIR);

if($debug){
	print STDERR join("\n", @files), "\n\n";
}
#colA_vs_ape1l_2A_diff_bases_P005_depth_3_50_both.txt

#if($debug){
#	print STDERR "deubg: colA_vs_Num77A_diff_bases_P005_depth_3_50_both.txt\n\n";
#	@files = ();
#	@files = ("colA_vs_Num77A_diff_bases_P005_depth_3_50_both.txt");
#}

#print STDERR join("\n", @files), "\n\n";


foreach my $file(@files){
	if($file =~ /colA_vs_(\S+)_diff_bases_P0\S+_depth/){
		my $pre = $1;
		my $input = File::Spec->catfile($indir, $file);
		die unless (-e $input);
		
		my $cmd = "perl $script $input $outdir $pre $p_value";
		
		print STDERR $cmd, "\n\n";
		if(!$debug){
			`$cmd`;
		}
		
	}else{
		die $file, "\n";
	}
}


exit;