#!/usr/bin/perl -w

use strict;
use File::Spec;
my $debug = 0;
my $usage = "$0 <indir> <outdir>";

die $usage unless (@ARGV == 2);

my ($indir, $outdir) = @ARGV[0..1];

die unless (-d $indir and -d $outdir);

opendir (DIR, $indir) or die "cannot open $indir";


my @files = grep /cytosines_depth\.txt$/, readdir DIR;

print STDERR join("\n", @files) , "\n\n";

foreach my $file(@files){
	if($file =~ /(\S+)_cytosines/){
		
		my $label = $1;
		#my $label_removed = $label;
			
		my $input  = File::Spec->catfile($indir, $file);
		my $output = File::Spec->catfile($outdir, $label."_cytosines_depth_histogram.tiff");
		print STDERR "output:\n";
		print STDERR $output , "\n\n";
		
		die unless (-e $input);
		die if(-e $output);
		
		
		 #         1    2        3       4       5      6         7
		my $cmd = "R --slave --vanilla --args  $input  $output $label ". "< /Users/tang58/Kai_BS/myown/Draw_depth_distribution.R";
		print STDERR $cmd , "\n\n";

		if(!$debug){
			`$cmd`;
		}
		
	}else{
		die "wrong name\n";
	}
}

exit;