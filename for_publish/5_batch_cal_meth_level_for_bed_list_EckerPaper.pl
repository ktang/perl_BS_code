#!/usr/bin/perl -w

use strict;
use File::Spec;

my $script = "/Users/tang58/Kai_BS/for_publish/cal_meth_level_for_bed_list_EckerPaper_v0.0.pl";
die unless (-e $script);

my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 \n <indir_isMeth> <indir_bed_list> <outdir> <WT_label(col)> \n\n";
die $usage unless(@ARGV == 4);

my $indir_isMeth = shift or die;
die unless (-d $indir_isMeth);

my $indir_bed = shift or die;
die unless (-d $indir_bed);

my $outdir = shift or die;
die unless (-d $outdir);

my $wt_label = shift or die;

opendir(METDIR, $indir_isMeth) or die;
my ($wt_isMet_file) = grep {/$wt_label/}  readdir METDIR;
closedir METDIR;


opendir(METDIR, $indir_isMeth) or die;
#my @files = grep /_rev\.txt$/, readdir DIR;
my @isMet_files = grep {!/$wt_label/} (grep /isMeth\S+\.txt$/, readdir METDIR);
closedir METDIR;

opendir(LISTDIR, $indir_bed) or die;
my @hyper_lists = grep {/hyper/}  readdir LISTDIR;
closedir LISTDIR;

opendir(LISTDIR, $indir_bed) or die;
my @hypo_lists = grep {/hypo/}  readdir LISTDIR;
closedir LISTDIR;

die unless ( @isMet_files == @hyper_lists );
die unless ( @isMet_files == @hypo_lists );


# JKZ18_rdd_isMeth
# JKZ18_rdd_vs

for my $i (0..$#isMet_files){
	my $isMet_file = $isMet_files[$i];
	my $hyper_file = $hyper_lists[$i];
	my $hypo_file  = $hypo_lists[$i];
	
	my ($pre_isMet, $pre_hyper, $pre_hypo) = (1,2,3);
	
	if($isMet_file =~ /(\S+)_isMeth/){ $pre_isMet = $1;	}
	if($hyper_file =~ /(\S+)_vs/){	$pre_hyper= $1;}
	if($hypo_file  =~ /(\S+)_vs/){	$pre_hypo= $1;}
	
	die $pre_hyper if( $pre_isMet ne $pre_hyper );
	die $pre_hypo if( $pre_isMet ne $pre_hypo );
}


print STDERR "isMet:\n";
print STDERR join("\n", @isMet_files ), "\n\n";

print STDERR "hyper:\n";
print STDERR join("\n", @hyper_lists), "\n\n";

print STDERR "hypo:\n";
print STDERR join("\n", @hypo_lists), "\n\n";


print STDERR "wt: $wt_isMet_file \n\n\n";

my $wt_file = File::Spec->catfile( $indir_isMeth, $wt_isMet_file );
die unless (-e $wt_file);

#<isMeht_file_wt> <isMeht_file_mut>  <bed_like_file> <hyper|hypo> <outdir>

foreach my $i (0..$#isMet_files){
	
	my $isMet_file = File::Spec->catfile($indir_isMeth ,$isMet_files[$i]);
	my $hyper_file = File::Spec->catfile($indir_bed ,$hyper_lists[$i]);
	my $hypo_file  = File::Spec->catfile($indir_bed ,$hypo_lists[$i]);
	
	my $hyper_cmd = "perl $script $wt_file  $isMet_file $hyper_file hyper $outdir";
	
	print STDERR $hyper_cmd, "\n\n";
	if(!$debug){
		`$hyper_cmd`;
	}

	my $hypo_cmd  = "perl $script $wt_file  $isMet_file $hypo_file hypo $outdir";
	print STDERR $hypo_cmd, "\n\n";
	if(!$debug){
		`$hypo_cmd`;
	}	
}

exit;

