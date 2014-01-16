#!/usr/bin/perl -w

#one_step_DMRs_v0.0.pl

#v0.0
# do the DMR finding in one script

# first, manullay create ln in a raw_isMeth_dir

# v0.1
# skip step_1 1_retain_isMeth_with_dep_cutoff_in_all_lib_v0.0.pl

use strict;
use File::Spec;

my $debug = 0;

my $depth_cutoff = 4;
my $p_value_cutoff = 0.05;

if($debug){
	print STDERR "debug = 1\n\n";
}
#my $usage = "$0 \n <work_dir> <postfix(ROS1Paper)> <raw_link_isMeth_file> <WT_label_complete> \n\n";
#die $usage unless(@ARGV == 4);

my $usage = "$0 \n <work_dir> <postfix(ROS1Paper)>  <WT_label_complete>  <isMeth_ln>\n\n";
die $usage unless(@ARGV == 4);



#my $input = shift or die;
#my $output = shift or die;
#die unless (-e $input);
#die if( -e $output);

my $work_dir = shift or die;
my $postfix  = shift or die;
#my $raw_isMeth_dir = shift or die;
my $WT_label = shift or die;

my $isMeth_dep	 = shift or die;

die unless (-d $isMeth_dep);

die unless (-d $work_dir);
#die unless (-d $raw_isMeth_dir);

#opendir (DIR, $raw_isMeth_dir) or die  $!;
#my @raw_isMeth_name = grep /isM.+\.txt$/, readdir DIR;
#closedir DIR;

#print STDERR join("\n", @raw_isMeth_name), "\n\n";

#chdir : change working directory

#my $isMeth_dep		= File::Spec->catfile( $work_dir, "isMeth_dep" . $depth_cutoff );
my $DMC_dir   		= File::Spec->catfile( $work_dir, "0_DMC_db");
my $raw_list_dir	= File::Spec->catfile( $work_dir, "1.0_raw_list");
my $Ecker_dir 		= File::Spec->catfile( $work_dir, "1.1_EckerDetail");
my $annotation_dir 	= File::Spec->catfile( $work_dir, "1.2_annotation");
my $filter_dir 		= File::Spec->catfile( $work_dir, "2_filter_len100_wC2fold_fCdiff10");



########
#	step_1 1_retain_isMeth_with_dep_cutoff_in_all_lib_v0.0.pl
#	<postfix> <outdir> <cutoff> <number_of_files> <isMeth> <label> [<isMeth> <label> ...] 
############

unless (-d $isMeth_dep){
	my $cmd_mkdir = "mkdir $isMeth_dep";
	print STDERR $cmd_mkdir, "\n";
	if(!$debug){
		`$cmd_mkdir`;
	}
}

#my $isMeth_num = $#raw_isMeth_name + 1;

# 1_retain_isMeth_with_dep_cutoff_in_all_lib_v0.0.pl
# <postfix> <outdir> <cutoff> <number_of_files> <isMeth> <label> [<isMeth> <label> ...] 

#my $cmd_1_retain = "1_retain_isMeth_with_dep_cutoff_in_all_lib_v0.0.pl $postfix $isMeth_dep $depth_cutoff $isMeth_num ";
#for my $file_name (@raw_isMeth_name){
#	my $label = get_label( $file_name );
#	my $real_file = File::Spec->catfile($raw_isMeth_dir, $file_name);
#	die $file_name unless (-e $real_file);
	
#	$cmd_1_retain =  $cmd_1_retain . " $real_file $label ";
	
#}

#print STDERR $cmd_1_retain, "\n\n";
#if(!$debug){
#	`$cmd_1_retain`;
#}



########
#	step_2  3_batch_generate_p_value_Fisher.pl
#	 <indir> <outdir> <WT_label> <postfix(edm2_paper)> 
############

unless (-d $DMC_dir){
	my $cmd_mkdir = "mkdir $DMC_dir";
	print STDERR $cmd_mkdir, "\n";
	if(!$debug){
		`$cmd_mkdir`;
	}
}

my $cmd_p_value = " 3_batch_generate_p_value_Fisher.pl $isMeth_dep $DMC_dir $WT_label $postfix ";

print STDERR $cmd_p_value, "\n\n";
if(!$debug){
	`$cmd_p_value`;
}

#######
#	step_3  select DMR
########################
# 4_batch_select_DMR_IDM1_paper_method_v1.0_sliding_win.pl
# <indir> <outdir> <p_value> <dep> <postfix(edm2_paper)>

unless (-d $raw_list_dir){
	my $cmd_mkdir = "mkdir $raw_list_dir";
	print STDERR $cmd_mkdir, "\n";
	if(!$debug){
		`$cmd_mkdir`;
	}
}

my $cmd_DMR = "4_batch_select_DMR_IDM1_paper_method_v1.0_sliding_win.pl $DMC_dir $raw_list_dir $p_value_cutoff $depth_cutoff $postfix";
print STDERR $cmd_DMR, "\n\n";
if(!$debug){
	`$cmd_DMR`;
}

#######
#	step_4 Ecker_detail
########################
# 5_batch_cal_meth_level_for_bed_list_EckerPaper.pl
# <indir_isMeth> <indir_bed_list> <outdir> <WT_label(col)> 

unless (-d $Ecker_dir){
	my $cmd_mkdir = "mkdir $Ecker_dir";
	print STDERR $cmd_mkdir, "\n";
	if(!$debug){
		`$cmd_mkdir`;
	}
}

my $cmd_Ecker = "5_batch_cal_meth_level_for_bed_list_EckerPaper.pl $isMeth_dep $raw_list_dir $Ecker_dir $WT_label";

print STDERR $cmd_Ecker, "\n\n";
if(!$debug){
	`$cmd_Ecker`;
}

########
# 5 annotation
###############
# 6_batch_add_annotation_TAIR10_v0.1.pl
#<indir> <outdir> <pattern>

unless (-d $annotation_dir){
	my $cmd_mkdir = "mkdir $annotation_dir";
	print STDERR $cmd_mkdir, "\n";
	if(!$debug){
		`$cmd_mkdir`;
	}
}

my $cmd_anno = "6_batch_add_annotation_TAIR10_v0.1.pl $Ecker_dir $annotation_dir txt";

print STDERR $cmd_anno, "\n\n";
if(!$debug){
	`$cmd_anno`;
}

########
# 6 filter
###############
# 8_batch_filter_len100_wC2fold_fCdiff10.pl
# <indir> <outdir>

unless (-d $filter_dir){
	my $cmd_mkdir = "mkdir $filter_dir";
	print STDERR $cmd_mkdir, "\n";
	if(!$debug){
		`$cmd_mkdir`;
	}
}

my $cmd_filter = "8_batch_filter_len100_wC2fold_fCdiff10.pl $annotation_dir $filter_dir";

print STDERR $cmd_filter, "\n\n";
if(!$debug){
	`$cmd_filter`;
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
=haed

unless (-d ){
	my $cmd_mkdir = "mkdir ";
	print STDERR $cmd_mkdir, "\n";
	if(!$debug){
		`$cmd_mkdir`;
	}
}


print STDERR , "\n\n";
if(!$debug){
	``;
}

=cut
