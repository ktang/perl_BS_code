#!/usr/bin/perl -w

# modified on Oct 21, 2013.
# the old version require A_vs_B pattern.

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);
# my ($volume,$directories,$file) =          File::Spec->splitpath( $path );

use strict;
use File::Spec;

my $debug = 0 ;

#my $script = "cal_meth_level_for_bed_list_EckerPaper_v0.1_single_isMeth.pl";

my $script = "/Users/tang58/Kai_BS/for_publish/batch_cal_meth_level_original_mC_NotEckerSuggest/cal_meth_level_for_bed_list_single_isMeth_original_mC_NotEckerSuggest_v0.0.pl";

die unless (-e $script);

print STDERR "\n\n script used is\n";
print STDERR $script, "\n\n";

my $usage = "$0 \n <isMeth_dir> <input_bed> <outdir> <output_name>\n\n";
die $usage unless(@ARGV == 4);

#my $orignal_file_name = shift or die;
my $dir_isMeth = shift or die;
my $orignal_file = shift or die;

#my $indir = shift or die;
my $outdir = shift or die;
my $output_name = shift or die;


#my $orignal_file = File::Spec->catfile($indir, $orignal_file_name);
die unless (-e $orignal_file);

#my ($volume,$bed_dir,$file_name) = File::Spec->splitpath( $orignal_file );
#my ( $mut_label, $wt_label);
#if( $file_name =~ /(\S+)_vs_(\S+)_hyp/ ){
#	$mut_label = $1;
#	$wt_label  = $2;
#}else{
#	die $file_name, "\n\nnot A_vs_B_hyp\n\n";
#}



#die unless (-d $indir);
die unless (-d $dir_isMeth);
die unless (-d $outdir);

my $output = File::Spec->catfile($outdir, $output_name);
die if(-e $output);

opendir(DIR, $dir_isMeth);
#my @files_raw = grep { $_ !~ $wt_label} (grep /\.txt$/, readdir DIR) ;
my @files_raw = grep /\.txt$/, readdir DIR ;
closedir DIR;

my @files;
my @labels;

my $k = -1;
foreach my $file (@files_raw){
	my $curr_label = get_label($file);
	#next if ( $curr_label eq $wt_label or $curr_label eq $mut_label);
	$k++;
	$files[$k] = $file;
	$labels[$k] = $curr_label;
}



if($debug ){
	print STDERR join("\n", @files), "\n\n\n";
}

#my $i = -1;

my $f0 = $orignal_file . ".temp0";

my $cp_cmd = "cp $orignal_file $f0";
print STDERR $cp_cmd, "\n\n";
if(!$debug){
	`$cp_cmd`;
}

#foreach my $file (@files){
for my $i (0..$#files){
	my $file = $files[$i];
	my $label = $labels[$i];
	
	my $bench_file =  $orignal_file . ".temp" . $i ;
	if(!$debug){
		die unless (-e $bench_file);
	}
	my $isMeth_file  =  File::Spec->catfile($dir_isMeth, $file );
	die unless (-e $isMeth_file);
	my $out_file =  $orignal_file . ".temp" . ($i+1) ;
	
	die if(-e $out_file);
	#<isMeht_file> <sample_label> <bed_like_file> <output> 	
	my $cmd = "$script $isMeth_file $label $bench_file  $out_file";
	
	print STDERR $cmd, "\n\n";
	if(!$debug){
		`$cmd`;
	}
}

my $final = $orignal_file . ".temp" . ($#files+1);

my $cp_out_cmd = "cp $final $output";
print STDERR $cp_out_cmd, "\n\n";
if(!$debug){
	`$cp_out_cmd`;
}

#my $rm_cmd = "rm -f $bed_dir/*txt.temp* ";
my $rm_cmd = "rm -f " . $orignal_file . ".temp* ";
print STDERR $rm_cmd, "\n\n";
if(!$debug){
	`$rm_cmd`;
}


exit;

#get_label($file);
sub get_label{
	my ($f) = @_;
	if ($f =~ /(\S+)_isMeth/){
		return $1;
	}else{
		die $f, "no label isMeth\n\n";
	}
}
