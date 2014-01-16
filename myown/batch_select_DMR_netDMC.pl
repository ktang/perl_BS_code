#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);


#input outdir pre p-value
use strict;
use File::Spec;

my $debug = 0;

#my $p_value = 0.01;

my $usage = "$0 <indir>  <outdir>";
die $usage unless(@ARGV == 2);

my $script = "/Users/tang58/Kai_BS/myown/select_DMR_output_list_for_individual_v2.0.pl";
die unless (-e $script);

my $indir = shift or die "indir";
#my $outdir = shift or die "outdir";
#my $p_cutoff = shift or die "p_value_cutoff";

my $outdir = shift or die "outdir";

die unless (-d $indir and -d $outdir);

#opendir(DIR, $indir) or die;
#my @files = grep/colA_vs\S+diff_bases\S+\.txt/, readdir DIR;
#closedir(DIR);

my @pres = ( "037306005190A", "Num12A", "Num27A", 
			 "Num77A", "akn1_2A", "ape_1A", "ape1l_2A",
			 "arp_1A", "ibm1-4A", "Jmij-19A",
			 "nrpd1-3A", "nrpel-11A", "P261A", "P313A",
			 "P314A", "P31A", "pkl_1A", "shh1A"	); # no colA


#colA_vs_ape1l_2A_diff_bases_P005_depth_3_50_both.txt

#if($debug){
#	print STDERR "deubg: colA_vs_Num77A_diff_bases_P005_depth_3_50_both.txt\n\n";
#	@files = ();
#	@files = ("colA_vs_Num77A_diff_bases_P005_depth_3_50_both.txt");
#}

#print STDERR join("\n", @files), "\n\n";


foreach my $pre(@pres){
	my $file = File::Spec->catfile($indir, $pre. "_mC.wig");
	die unless (-e $file);
}


foreach my $pre(@pres){
	my $file = File::Spec->catfile($indir, $pre. "_mC.wig");
	die unless (-e $file);
	my $cmd = "perl $script $file $outdir $pre";
	print STDERR $cmd, "\n\n";
	if(!$debug){
		`$cmd`;
	}
}


exit;