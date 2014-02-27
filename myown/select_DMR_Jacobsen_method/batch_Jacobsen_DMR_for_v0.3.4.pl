#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);
use strict;
use File::Spec;
my $debug = 0;

#	print STDERR "\n\n nodupl tag\n\n";


my $script = "/Users/tang58/Kai_BS/myown/select_DMR_Jacobsen_method/select_DMR_Jacobsen_method_v0.3.4.pl";
die unless (-e $script);

print STDERR "used script \n\n  $script \n\n";


my $usage = "$0 \n <indir> <outdir> <control_WT_file> <WT_label>\n\n";
die $usage unless(@ARGV == 4);

my $indir = shift or die;
my $outdir = shift or die;
die unless (-d $indir);
die unless (-d $outdir);

#my $cutoff = shift or die;

my $wt_file = shift or die;
die unless (-e $wt_file);

my $wt_label = shift or die;

opendir(DIR, $indir);
my @files = grep /isMeth.*\.txt$/, readdir DIR;
closedir (DIR);


my %pres;

foreach my $file(@files){
#	next if($file =~ /col/i);
#	next if($file =~ /ibm/i);
#	next if($file =~ /ros1/i);
#	next if($file =~ /rdd/i);
#	next if($file =~ /WT/i);
	if($file =~ /(\S+)_isMeth/){
		my $pre = $1;
		$pres{$pre} = $file;
	}else{
		die $file;
	}
}


#my $dir = "/Volumes/My_Book/20120427_ShangHai_data/handled_data/cytosine_coverage_info";

foreach my $pre(sort keys %pres){
	my $input =  File::Spec->catfile($indir, $pres{$pre});
	die unless (-e $input);
	
	my $cmd = "time perl $script $wt_file $input $outdir $pre" . "_vs_" . $wt_label;

	
	print STDERR "\n", $cmd, "\n\n";
	if(!$debug){
		`$cmd`;
	}
}

exit;
