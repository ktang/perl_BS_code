#!/usr/bin/perl -w
use strict;
use File::Spec;

my $debug = 0;

my $depth_cutoff = 4;
my $p_value_cutoff = 0.05;# 0.01;

my $script_dir = "/Users/tang58/Kai_BS/IDM1_Science_DMR_one_step_Nov19_2013";

if($debug){
	print STDERR "debug = 1\n\n";
}

my $usage = "$0 \n <1_work_dir> <2_postfix(ROS1Paper or self)> <3_isMeth_ln>
 <4_WT_isMeth_file_name (file in ln_isMeth dir)>  <5_WT_label_complete>  \n\n";
die $usage unless(@ARGV == 5);


my $work_dir 		= shift or die;
my $postfix  		= shift or die;
my $isMeth_dir		= shift or die;

my $WT_file_name 	= shift or die;
my $WT_label 		= shift or die;

die unless (-d $isMeth_dir);
die unless (-d $work_dir);

my $DMC_dir   		= File::Spec->catfile( $work_dir, "0_DMC_db");
my $raw_list_dir	= File::Spec->catfile( $work_dir, "1.0_raw_list");
my $list_with_meth_dir	= File::Spec->catfile( $work_dir, "1.1_MethDetail_notE");
my $anno_dir	 	= File::Spec->catfile( $work_dir, "1.2_anno");


########
#	step1  
#	 <indir> <outdir> <WT_label> <postfix(edm2_paper)> 
############

unless (-d $DMC_dir){
	my $cmd_mkdir = "mkdir $DMC_dir";
	print STDERR $cmd_mkdir, "\n";
	if(!$debug){
		`$cmd_mkdir`;
	}
}

my $script_step1_Fisher_p_value = File::Spec-> catfile($script_dir, "step1_batch_generate_p_value_Fisher.pl"); 
my $cmd_p_value = " $script_step1_Fisher_p_value $isMeth_dir $DMC_dir $WT_file_name $WT_label $postfix ";

print STDERR "step1\ncal p-value...\n",  $cmd_p_value, "\n\n";
if(!$debug){
	`$cmd_p_value`;
}

print STDERR "#" x40, "\n\n\n";

#######
#	step2 select DMR
########################
my $step2_script  = "/Users/tang58/Kai_BS/IDM1_Science_DMR_one_step_Nov19_2013/step2_batch_select_DMR_IDM1_paper_method.pl";
die unless (-e $step2_script);

unless (-d $raw_list_dir){
	my $cmd_mkdir = "mkdir $raw_list_dir";
	print STDERR $cmd_mkdir, "\n";
	if(!$debug){
		`$cmd_mkdir`;
	}
}

my $cmd_DMR = "$step2_script $DMC_dir $raw_list_dir $p_value_cutoff $depth_cutoff $postfix";
print STDERR "step2\nselect DMR...\n", $cmd_DMR, "\n\n";
if(!$debug){
	`$cmd_DMR`;
}
print STDERR "#" x40, "\n\n\n";

#######
#	step_3 meth_level
########################
# 5_batch_cal_meth_level_for_bed_list_EckerPaper.pl
# <indir_isMeth> <indir_bed_list> <outdir> <WT_label(col)> 

unless (-d $list_with_meth_dir){
	my $cmd_mkdir = "mkdir $list_with_meth_dir";
	print STDERR $cmd_mkdir, "\n";
	if(!$debug){
		`$cmd_mkdir`;
	}
}

my $step3_script = "/Users/tang58/Kai_BS/for_publish/batch_cal_meth_level_original_mC_NotEckerSuggest_v0.1_bed_may_overlap/batch_cal_meth_level_for_bed_list_original_mC_NotEckerSuggest_multiple_isMeth_v0.1.pl";
#<dir_isMeth> <indir> <outdir>

my $cmd_Ecker = "$step3_script $isMeth_dir $raw_list_dir $list_with_meth_dir";

print STDERR "step3\ncal meth level...\n", $cmd_Ecker, "\n\n";
if(!$debug){
	`$cmd_Ecker`;
}
print STDERR "#" x40, "\n\n\n";

########
# step4: annotation
###############
# 6_batch_add_annotation_TAIR10_v0.1.pl
#<indir> <outdir> <pattern>

unless (-d $anno_dir){
	my $cmd_mkdir = "mkdir $anno_dir";
	print STDERR $cmd_mkdir, "\n";
	if(!$debug){
		`$cmd_mkdir`;
	}
}

my $cmd_anno = "6_batch_add_annotation_TAIR10_v0.1.pl $list_with_meth_dir $anno_dir txt";

print STDERR "step4\anno...\n", $cmd_anno, "\n\n";
if(!$debug){
	`$cmd_anno`;
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
