#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;

my $debug = 1;

print STDERR "\n\n please note: maybe need to change trim and trim parameter -L \n\n";


#my $script_trim = "/Users/tang58/scripts_all/perl_code/Purdue/raw_fastq_to_clean_fastq_v1.1.pl";
my $script_trim = "/Users/tang58/scripts_all/perl_code/Purdue/raw_fastq_to_clean_fastq_SE.pl";
die unless (-e $script_trim);

my $ref = "/Users/tang58/DataBase/TAIR_Col0_genome/index/BRAT_BW/brat_bw_arabdopsis_index";
die unless (-d $ref);

my $brat_bw = "/Users/tang58/Software/BRAT/brat_bw-2.0.1/brat_bw";
die unless (-e $brat_bw);


my $remove_dupl = "/Users/tang58/Software/BRAT/brat_bw-2.0.1/remove-dupl";
die unless (-e $remove_dupl);
my $fas_ref = "/Users/tang58/DataBase/BRAT_reference_files/Col0/Col0_references_names.txt";
die unless (-e $fas_ref);


my $acgt_count = "/Users/tang58/Software/BRAT/brat_bw-2.0.1/acgt-count";
die unless (-e $acgt_count);

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 <indir> <inpre> <outpre>";
die $usage unless(@ARGV == 3);

my $indir  = shift or die;
my $inpre  = shift or die;
my $outpre = shift or die;

my $outdir = $indir;
die unless(-d $indir);

##########################
# 1 trim adapter and low quality reads
##########################


my $adaptor = "AGATCGGAAGAGC";

my $phred = 20;
my $minLen = 24;
#my $dir = ".";

#my $raw_fq_2 = File::Spec->catfile($indir, $inpre . "_2_sequence.txt");
my $raw_fq = File::Spec->catfile($indir, $inpre . ".fastq");
die unless (-e $raw_fq);

my $trim_log = File::Spec->catfile($outdir,  $outpre."_trim_log.txt");

my $cmd = "time perl $script_trim --adaptor $adaptor --input $raw_fq --output_pre $outpre --dir $indir --minLen $minLen --phred $phred >> $trim_log";

print STDERR $cmd, "\n\n";

if(!$debug){
	`$cmd`;
}

###################################
# 2 fastq2brat input
###################################
#my $fq_P1 = $outpre."_reads1.fastq"; 
#my $fq_P2 = $outpre."_reads2.fastq"; 
#my $fq_S1 = $outpre."_mates1.fastq"; 
#my $fq_S2 = $outpre."_mates2.fastq"; 

my $fq_clean = $outpre . "_clean.fastq";

if(!$debug){
	#die unless (-e $fq_clean);
}

if(!$debug){
	fastq2brat($indir, $fq_clean);
#	fastq2brat($indir, $fq_P1);
#	fastq2brat($indir, $fq_P2);
#	fastq2brat($indir, $fq_S1);
#	fastq2brat($indir, $fq_S2);
}


##################################
# 3 brat
##################################

my $i = 0;
my $a = 1000;
my $m = 2;

#SE mates 1,2
#my $in_S1 = File::Spec->catfile($indir, $outpre."_mates1.txt");
#my $in_S2 = File::Spec->catfile($indir, $outpre."_mates2.txt");

my $in_SE = File::Spec->catfile($indir, $outpre."_clean.txt");

if(!$debug){
#	die "wrong input" unless (-e $in_P1 and -e $in_P2 and -e $in_S1 and -e $in_S2);
	die unless (-e $in_SE);

}

#my $out_p  = File::Spec->catfile($outdir, $outpre."_paired_brat_bw_out.txt");
#my $out_S1 = File::Spec->catfile($outdir, $outpre . "_1_brat_bw_out.txt");
#my $out_S2 = File::Spec->catfile($outdir, $outpre. "_2_brat_bw_out.txt");
#my $out_left = File::Spec->catfile($outdir, $outpre . "_paired_brat_bw_out.txt.single_mates");

my $out_SE = File::Spec->catfile($outdir, $outpre . "_SE_brat_bw_out.txt");


print STDERR "output:\n";
#print STDERR join("\n", ($out_S1, $out_S2, $out_p, $out_left)), "\n\n";
print STDERR $out_SE, "\n\n";
#die "output exists" unless (!(-e "$outdir/$out_p") and !(-e "$outdir/$out_S1") and !(-e "$outdir/$out_S2") );
#die "output exists" if(-e $out_p or -e $out_S1 or -e $out_S2);
die "output exists" if(-e $out_SE);


my $brat_bw_log = File::Spec->catfile($outdir, $outpre."_brat_bw_log.txt");


my $cmd_SE = "time $brat_bw -s $in_SE -P $ref -o $out_SE -m $m  >> $brat_bw_log 2>&1";
print STDERR $cmd_SE,"\n\n";


if(!$debug){
	open(OUT, ">>$brat_bw_log") or die;
	print OUT $cmd_SE, "\n\n";
	close OUT;
	`$cmd_SE`;
}



################################
# 4 remv-dupl
################################
#time /Users/tang58/Software/BRAT/brat_bw-2.0.1/remove-dupl -r ~/DataBase/BRAT_reference_files/Col0/Col0_references_names.txt -p   ./JKZ4_s1_ros1_4_pair_list.txt -1 ./JKZ4_s1_ros1_4_single_list.txt >> remove-dupl_log.txt  2>&1

my $single_list = File::Spec->catfile($indir, $outpre."_single_list.txt");
die if (-e $single_list);
#my $pair_list   = File::Spec->catfile($indir, $outpre."_pair_list.txt");
#die if (-e $pair_list );

if(!$debug){

	die if (-e $single_list);
	open(SOUT, ">$single_list") or die "Can't open $single_list: $!";

#	die if (-e $pair_list );
#	open(POUT, ">$pair_list") or die "Can't open $pair_list: $!";

	print SOUT $out_SE,"\n";
#	print SOUT $out_S2,"\n";
#	print SOUT $out_left, "\n";
	close SOUT;

#	print POUT $out_p, "\n";
#	close POUT;
}

#./remove-dupl -r references_names.txt -s single_results.txt

#my $cmd_remove = "time $remove_dupl -r $fas_ref -p $pair_list -1 $single_list >> $brat_bw_log  2>&1";
my $cmd_remove = "time $remove_dupl -r $fas_ref -s $single_list >> $brat_bw_log  2>&1";
print STDERR $cmd_remove, "\n\n";
if(!$debug){
	`$cmd_remove`;
}

################################
# 5 acgt-count
################################
#.nodupl
my $out_SE_nodupl = File::Spec->catfile($outdir, $outpre.  "_SE_brat_bw_out.txt.nodupl");

my $single_list_nodupl = File::Spec->catfile($indir, $outpre."_single_list_nodupl.txt");
die if (-e $single_list_nodupl);
if(!$debug){

	die if (-e $single_list_nodupl);
	open(SOUT, ">$single_list_nodupl") or die "Can't open $single_list_nodupl: $!";


	print SOUT $out_SE_nodupl, "\n";
	close SOUT;

}

#./acgt-count -r references_names.txt -P prefix -s files_with_singles.txt

#my $cmd_acgt_count = "time $acgt_count -r $fas_ref -P $indir/$outpre -p $pair_list_nodupl -s $single_list_nodupl -B";

my $cmd_acgt_count = "time $acgt_count -r $fas_ref -P $indir/$outpre -s $single_list_nodupl -B";


print STDERR $cmd_acgt_count, "\n\n";
if(!$debug){
	`$cmd_acgt_count`;
}



exit;

sub round {
    my($number) = shift;
    #return int($number + .5);
    return int($number + .5 * ($number <=> 0)); # take care of negative numbers too

}

sub fastq2brat{
	my ($dir, $fq_file) = @_;
	my $script_fastq2brat = "/Users/tang58/Kai_BS/myown/fastq2brat_bw_input.pl";
	die unless (-e $script_fastq2brat);

	my $input = File::Spec->catfile($dir, $fq_file);
	if($fq_file =~ /(\S+)\.fastq$/){
		my $output = File::Spec->catfile($dir, $1. ".txt");
		die if(-e $output);
		my $cmd = "perl $script_fastq2brat $input $output";
		print STDERR $cmd, "\n\n";
		if(!$debug){
			`$cmd`;
		}
		
	}else{
		die $fq_file;
	}
	
}
