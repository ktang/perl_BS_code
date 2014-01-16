#!/usr/bin/perl -w
#extract_commod_DMR_for_JacobsenCellMethod.pl diff mut
# mutants have duplicate and the same WT

#this overlap script is specific for Jacobsen method
use strict;
use File::Spec;

my $script = "/Users/tang58/Kai_BS/Jacobsen_Cell_method/extract_common_DMR_diffMut_vsWT_for_JacobsenCellMethod.pl";
die unless (-e $script);

my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 \n <indir1>  <indir2>  <outdir>\n\n";

die $usage unless(@ARGV == 3);

my $indir1    = shift or die;
#my $vs_label1 = shift or die;
my $indir2    = shift or die;
#my $vs_label2 = shift or die;
my $outdir    = shift or die;

#die unless (-d $indir);
#die if( -e $output);

die unless (-d $indir1);
die unless (-d $indir2);
die unless (-d $outdir);

opendir(DIR1, $indir1 ) or die;
opendir(DIR2, $indir2 ) or die;

my @fs1 = grep /CG_hyper/ , readdir DIR1;
my @fs2 = grep /CG_hyper/ , readdir DIR2;

closedir (DIR1);
closedir (DIR2);

die unless (@fs1 == @fs2);

for my $i(0..$#fs1){
	my $f1 = $fs1[$i];
	my $f2 = $fs2[$i];
	
	print join("\n", ($f1, $f2)), "\n\n" if ($debug);
}

if ($debug) {
	#exit;#code
}

for my $i(0..$#fs1){
	my $f1 = $fs1[$i];
	my $f2 = $fs2[$i];
	
	my $label1 = get_label($f1);
	my $label2 = get_label($f2);
	
	
	my $cmd = "perl $script $indir1 $label1 $indir2 $label2 $outdir ";
	print STDERR $cmd, "\n\n";
	if ( !$debug) {
		`$cmd`;#code
	}
	
}

exit;

sub get_label{
	my  ($file) = @_;
	if ($file =~ /(\S+)_vs_/) {
		return $1;
	}else{
		die $file;
	}
	
}