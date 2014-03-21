#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);
BEGIN { push @INC, '/Users/tang58/scripts_all/perl_code/Modules' }
use Kai_Module;
use strict;
use File::Spec;

my $debug = 0;

print STDERR "\n\n please note: maybe need to change trim and trim parameter -L \n\n";

my $L = 33;

if($debug){
	print STDERR "debug = 1\n\n";
}

my $brat_bw = "/Users/tang58/Software/BRAT/brat_bw-2.0.1/brat_bw";
die unless (-e $brat_bw);

#my $trim_tool = "/Users/tang58/Software/BRAT/brat_bw-2.0.1/trim";
my $trim_tool = "/Users/tang58/Software/BRAT/brat-1.2.3/trim";
die unless (-e $trim_tool);

my $ref = "/Users/tang58/DataBase/TAIR_Col0_genome/index/BRAT_BW/brat_bw_arabdopsis_index";
die unless (-d $ref);

my $fas_ref = "/Users/tang58/DataBase/BRAT_reference_files/Col0/Col0_references_names.txt";
die unless (-e $fas_ref);

my $remove_dupl = "/Users/tang58/Software/BRAT/brat_bw-2.0.1/remove-dupl";
die unless (-e $remove_dupl);

my $acgt_count = "/Users/tang58/Software/BRAT/brat_bw-2.0.1/acgt-count";
die unless (-e $acgt_count);

my $usage = "$0 \n<indir> <pre>\n\n";
die $usage unless(@ARGV == 2);

my $indir = shift or die "indir";
my $pre = shift or die "pre";
#my $outdir = shift or die "outdir";

my $outpre = $pre;

#opendir(DIR, $indir) or die "dir";
#my @files = grep /fq\.gz$/, readdir DIR;
#die unless (@files == 2);
#die "1: $files[0]" unless ($files[0] =~ /1\.fq\.gz$/);
#die "2: $files[1]" unless ($files[1] =~ /2\.fq\.gz$/);
#my $gz1    = File::Spec->catfile($indir, $files[0]);
#my $gz2    = File::Spec->catfile($indir, $files[1]);


my $gz1    = get_fqs("1.fq.gz", $indir);
my $gz2    = get_fqs("2.fq.gz", $indir);


my $fq_in1 = File::Spec->catfile($indir, $pre . "_1.fq");
my $fq_in2 = File::Spec->catfile($indir, $pre . "_2.fq");



if(!$debug){
	#die unless (-e $gz1);
	#die unless (-e $gz2);
	
	die if (-e $fq_in1);
	die if (-e $fq_in2);
}

######################
# 0 unzip
##########################

my $cmd_unzip1 = "time zcat $gz1 > $fq_in1";
print STDERR $cmd_unzip1, "\n\n";
if(!$debug){
	`$cmd_unzip1`;
}

if(!$debug){
	die unless (-e $gz1);
	die unless (-e $gz2);
	
	die unless (-e $fq_in1);
	die if (-e $fq_in2);
}

my $cmd_unzip2 = "time zcat $gz2 > $fq_in2";
print STDERR $cmd_unzip2, "\n\n";
if(!$debug){
	`$cmd_unzip2`;
}

if(!$debug){
	die unless (-e $gz1);
	die unless (-e $gz2);
	
	die unless (-e $fq_in1);
	die unless (-e $fq_in2);
}


#####################
# 1 brat_trim
#####################
my $q = 20;
my $m_trim = 0;

my $trim_cmd = "time $trim_tool -1 $fq_in1 -2 $fq_in2 -P $indir/$pre -q $q -L $L -m $m_trim";
print STDERR $trim_cmd,"\n\n";

my $rm_cmd = "rm -f $indir/*q";
print STDERR  $rm_cmd, "\n\n";

if(!$debug){
	`$trim_cmd`;
	
	`$rm_cmd`;
}

#####################
# 2 brat_bw
#####################
my $i = 0;
my $a = 1000;
my $m_brat_bw = 2;

#PE reads 1,2
my $in_P1 = File::Spec->catfile($indir, $pre."_reads1.txt");
my $in_P2 = File::Spec->catfile($indir, $pre."_reads2.txt");

#SE mates 1,2
my $in_S1 = File::Spec->catfile($indir, $pre."_mates1.txt");
my $in_S2 = File::Spec->catfile($indir, $pre."_mates2.txt");


my $outdir = $indir;

my $out_p    = File::Spec->catfile($outdir, $pre. "_paired_brat_bw_out.txt");
my $out_S1   = File::Spec->catfile($outdir, $pre. "_1_brat_bw_out.txt");
my $out_S2   = File::Spec->catfile($outdir, $pre. "_2_brat_bw_out.txt");
my $out_left = File::Spec->catfile($outdir, $pre. "_paired_brat_bw_out.txt.single_mates");

if(!$debug){
	die "wrong input" unless (-e $in_P1 and -e $in_P2 and -e $in_S1 and -e $in_S2);
	die "output exists" if(-e $out_p or -e $out_S1 or -e $out_S2);
}

my $log = File::Spec->catfile($outdir, $pre."_brat_bw_log.txt");


my $cmd_S1 = "time $brat_bw -s $in_S1 -P $ref -o $out_S1 -m $m_brat_bw  >> $log 2>&1";
print STDERR $cmd_S1,"\n\n";
open(OUT, ">>$log") or die;
print OUT $cmd_S1, "\n\n";
close OUT;

if(!$debug){
	`$cmd_S1`;
}


my $cmd_S2 = "time $brat_bw -s $in_S2 -P $ref -o $out_S2 -m $m_brat_bw -A >> $log 2>&1";
print STDERR $cmd_S2,"\n\n";
open(OUT, ">>$log") or die;
print OUT "\n\n", $cmd_S2, "\n\n";
close OUT;
if(!$debug){
	`$cmd_S2`;
}

my $cmd_P = "time $brat_bw -1 $in_P1 -2 $in_P2 -pe -o $out_p -P $ref -i $i -a $a -m $m_brat_bw  >> $log 2>&1";
print STDERR $cmd_P,"\n\n";
open(OUT, ">>$log") or die;
print OUT "\n\n",  $cmd_P, "\n\n";
close OUT;
if(!$debug){
	`$cmd_P`;
}


################################
# 3 remv-dupl
################################
#time /Users/tang58/Software/BRAT/brat_bw-2.0.1/remove-dupl -r ~/DataBase/BRAT_reference_files/Col0/Col0_references_names.txt -p   ./JKZ4_s1_ros1_4_pair_list.txt -1 ./JKZ4_s1_ros1_4_single_list.txt >> remove-dupl_log.txt  2>&1

my $single_list = File::Spec->catfile($indir, $outpre."_single_list.txt");
die if (-e $single_list);
my $pair_list   = File::Spec->catfile($indir, $outpre."_pair_list.txt");
die if (-e $pair_list );

if(!$debug){

	die if (-e $single_list);
	open(SOUT, ">$single_list") or die "Can't open $single_list: $!";

	die if (-e $pair_list );
	open(POUT, ">$pair_list") or die "Can't open $pair_list: $!";

	print SOUT $out_S1,"\n";
	print SOUT $out_S2,"\n";
	print SOUT $out_left, "\n";
	close SOUT;

	print POUT $out_p, "\n";
	close POUT;
}


my $cmd_remove = "time $remove_dupl -r $fas_ref -p $pair_list -1 $single_list >> $log  2>&1";
print STDERR $cmd_remove, "\n\n";
if(!$debug){
	`$cmd_remove`;
}

################################
# 4 acgt-count
################################
#.nodupl
my $out_p_nodupl    = File::Spec->catfile($outdir, $outpre."_paired_brat_bw_out.txt.nodupl");
my $out_S1_nodupl   = File::Spec->catfile($outdir, $outpre . "_1_brat_bw_out.txt.nodupl");
my $out_S2_nodupl   = File::Spec->catfile($outdir, $outpre. "_2_brat_bw_out.txt.nodupl");
my $out_left_nodupl = File::Spec->catfile($outdir, $outpre . "_paired_brat_bw_out.txt.single_mates.nodupl");

my $single_list_nodupl = File::Spec->catfile($indir, $outpre."_single_list_nodupl.txt");
die if (-e $single_list_nodupl);
my $pair_list_nodupl   = File::Spec->catfile($indir, $outpre."_pair_list_nodupl.txt");
die if (-e $pair_list_nodupl );
if(!$debug){

	die if (-e $single_list_nodupl);
	open(SOUT, ">$single_list_nodupl") or die "Can't open $single_list_nodupl: $!";

	die if (-e $pair_list_nodupl );
	open(POUT, ">$pair_list_nodupl") or die "Can't open $pair_list_nodupl: $!";

	print SOUT $out_S1_nodupl,"\n";
	print SOUT $out_S2_nodupl,"\n";
	print SOUT $out_left_nodupl, "\n";
	close SOUT;

	print POUT $out_p_nodupl, "\n";
	close POUT;
}

my $cmd_acgt_count = "time $acgt_count -r $fas_ref -P $indir/$outpre -p $pair_list_nodupl -s $single_list_nodupl -B";

print STDERR $cmd_acgt_count, "\n\n";
if(!$debug){
	`$cmd_acgt_count`;
}

print STDERR "\n\n OK \n\n";


exit;

