#!/usr/bin/perl -w

#R --slave --vanilla --args Gene_all Wenfeng_Gene_all Gene Gene/ < /Users/tang58/Kai_BS/useful_R_script/10_gene_TE_bin_files_R_cmd_for_all_command_v0.1.r
#R --slave --vanilla --args TE_all Wenfeng_TE_all TE TE/ < /Users/tang58/Kai_BS/useful_R_script/10_gene_TE_bin_files_R_cmd_for_all_command_v0.1.rA

#ALL
#Gene
#TE
## R --slave --vanilla --args  file_pattern(Gene_all TE_all) png_pre  type_TE_Gene dir < /Users/tang58/Kai_BS/useful_R_script/10_gene_TE_bin_files_R_cmd_for_all_command_v0.1.r 

#Gene_length
## R --slave --vanilla --args   png_pre   dir_Gene < /Users/tang58/Kai_BS/useful_R_script/10_gene_TE_bin_files_R_cmd_Gene_command_input.r

#TE_length
## R --slave --vanilla --args  png_pre indir_TE < /Users/tang58/Kai_BS/useful_R_script/10_gene_TE_bin_files_R_cmd_TE_command_input.r


use strict;
use File::Spec;

my $Gene_length_script = "/Users/tang58/Kai_BS/useful_R_script/10_gene_TE_bin_files_R_cmd_Gene_command_input.r";
my $TE_length_script   = "/Users/tang58/Kai_BS/useful_R_script/10_gene_TE_bin_files_R_cmd_TE_command_input.r";

my $all_script = "/Users/tang58/Kai_BS/useful_R_script/10_gene_TE_bin_files_R_cmd_for_all_command_v0.1.r";

die unless (-e $Gene_length_script);
die unless (-e $TE_length_script);
die unless (-e $all_script);

my $TE_dir = "TE";
my $Gene_dir = "Gene";
unless ( -d $TE_dir and -d $Gene_dir){
	print STDERR "must have dir: TE and Gene\n\n";
	exit;
}

my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 \n <pre>\n\n";
die $usage unless(@ARGV == 1);
my $pre = shift or die;

######
# 1 Gene_all
my $pre_all_gene = $pre . "_Gene_all_meth_level";

## R --slave --vanilla --args  file_pattern(Gene_all TE_all) png_pre  type_TE_Gene dir < /Users/tang58/Kai_BS/useful_R_script/10_gene_TE_bin_files_R_cmd_for_all_command_v0.1.r 
my $cmd_gene_all = "R --slave --vanilla --args Gene_all $pre_all_gene Gene $Gene_dir < $all_script";
print STDERR $cmd_gene_all, "\n\n";
unless($debug){
	`$cmd_gene_all`;
}

######
# 2 TE_all
my $pre_all_TE = $pre . "_TE_all_meth_level";
my $cmd_TE_all = "R --slave --vanilla --args TE_all $pre_all_TE TE $TE_dir < $all_script";
print STDERR $cmd_TE_all, "\n\n";
unless($debug){
	`$cmd_TE_all`;
}

#######
# 3 Gene_length
my $pre_gene_length = $pre . "_Gene_length_meth_level";
my $cmd_Gene_length = "R --slave --vanilla --args $pre_gene_length $Gene_dir < $Gene_length_script";
print STDERR $cmd_Gene_length, "\n\n";
unless($debug){
	`$cmd_Gene_length`;
}
#######
# 4 TE_length
my $pre_TE_length = $pre . "_TE_length_meth_level";
my $cmd_TE_length = "R --slave --vanilla --args $pre_TE_length $TE_dir < $TE_length_script";
print STDERR $cmd_TE_length, "\n\n";
unless($debug){
	`$cmd_TE_length`;
}

my $mv_cmd = "mv */*png .";
print STDERR $mv_cmd, "\n";
unless($debug){
	`$mv_cmd`;
}

exit;
=head
print STDERR , "\n\n";
unless($debug){
	``;
}
=end