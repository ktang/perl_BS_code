#!/usr/bin/perl -w

#one_step_DMRs_v0.0.pl

#v0.0
# do the DMR finding in one script

# first, manullay create ln in a raw_isMeth_dir
use strict;
use File::Spec;

my $debug = 0;

#my $depth_cutoff = 4;
#my $p_value_cutoff = 0.05;

if($debug){
	print STDERR "debug = 1\n\n";
}

my $usage = "$0 \n <raw_link_isMeth_file> <outdir> <postfix(ROS1Paper)> <depth_cutoff>\n\n";
die $usage unless(@ARGV == 4);

my $raw_isMeth_dir = shift or die;
my $outdir = shift or die;
my $postfix  = shift or die;
my $depth_cutoff = shift or die;

#my $usage = "$0 \n <work_dir> <postfix(ROS1Paper)> <raw_link_isMeth_file> <WT_label_complete> \n\n";
#die $usage unless(@ARGV == 4);

#my $input = shift or die;
#my $output = shift or die;
#die unless (-e $input);
#die if( -e $output);

#my $work_dir = shift or die;
#my $postfix  = shift or die;
#my $raw_isMeth_dir = shift or die;
#my $WT_label = shift or die;

#die unless (-d $work_dir);
die unless (-d $outdir);
die unless (-d $raw_isMeth_dir);

opendir (DIR, $raw_isMeth_dir) or die  $!;
my @raw_isMeth_name = grep /isM.+\.txt$/, readdir DIR;
closedir DIR;

print STDERR join("\n", @raw_isMeth_name), "\n\n";



########
#	step_1 1_retain_isMeth_with_dep_cutoff_in_all_lib_v0.0.pl
#	<postfix> <outdir> <cutoff> <number_of_files> <isMeth> <label> [<isMeth> <label> ...] 
############



my $isMeth_num = $#raw_isMeth_name + 1;

# 1_retain_isMeth_with_dep_cutoff_in_all_lib_v0.0.pl
# <postfix> <outdir> <cutoff> <number_of_files> <isMeth> <label> [<isMeth> <label> ...] 

my $cmd_1_retain = "1_retain_isMeth_with_dep_cutoff_in_all_lib_v0.0.pl $postfix $outdir $depth_cutoff $isMeth_num ";
for my $file_name (@raw_isMeth_name){
	my $label = get_label( $file_name );
	my $real_file = File::Spec->catfile($raw_isMeth_dir, $file_name);
	die $file_name unless (-e $real_file);
	
	$cmd_1_retain =  $cmd_1_retain . " $real_file $label ";
	
}

print STDERR $cmd_1_retain, "\n\n";
if(!$debug){
	`$cmd_1_retain`;
}



exit;
sub get_label{
	my ($f) = @_;
	if ($f =~ /(\S+)_isMeth/){
		return $1;
	}else{
		die $f, "no label isMeth\n\n";
	}
}
