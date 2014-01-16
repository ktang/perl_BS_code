#!/usr/bin/perl -w

#v0.1.1
#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);
use strict;
use File::Spec;

my $debug = 0;

#	print STDERR "\n\n nodupl tag\n\n";
my $script = "/Users/tang58/Kai_BS/myown/select_DMR_Jacobsen_method/filter_Sep15_v0.1.1.pl";
die unless (-e $script);

print STDERR "used script \n\n  $script \n\n";

my ( $CG_num_cutoff, $CHG_num_cutoff, $CHH_num_cutoff ) = (5, 5, 20);
my ( $CG_diff_cutoff, $CHG_diff_cutoff, $CHH_diff_cutoff ) = (40, 20 ,10);


my $usage = "$0 \n <indir> <outdir> \n\n";
die $usage unless(@ARGV == 2);

my $indir = shift or die;
my $outdir = shift or die;
die unless (-d $indir);
die unless (-d $outdir);

#my $cutoff = shift or die;

#my $isMeth_file = shift or die;
#die unless (-e $isMeth_file);


opendir(DIR, $indir);
#my @files = grep /vs.+\.txt$/, readdir DIR;
my @files = grep /\.txt$/, readdir DIR;
closedir (DIR);


foreach my $file(@files){
	my $input =  File::Spec->catfile($indir, $file);
	die unless (-e $input);
	
	my $outfile ;
	
	if($file =~ /(\S+)\.txt/){
		$outfile = $1 . "_filtered_" . $CG_num_cutoff ."_".  $CG_diff_cutoff ."_". $CHG_num_cutoff."_".  $CHG_diff_cutoff."_". $CHH_num_cutoff ."_".$CHH_diff_cutoff . ".txt" ;
	}else{
		die $file;
	}
	
	my $output = File::Spec->catfile($outdir, $outfile);
	
	die if(-e $output);
	
	my $cmd = "perl $script $CG_num_cutoff   $CG_diff_cutoff  $CHG_num_cutoff $CHG_diff_cutoff $CHH_num_cutoff $CHH_diff_cutoff  $input > $output";

	
	print STDERR "\n", $cmd, "\n\n";
	if(!$debug){
		`$cmd`;
	}
}

exit;
