#!/usr/bin/perl -w

use utf8;#可以吗？
use strict;
use File::Spec;

my $debug = 0;

if($debug){
	print STDERR "debug = 1\n";
}
#my $genome_size = "1.2e8";
#my $keep_dup = "all";

#my $usage = "\n$0 \n\n <treatment_ChIP> <control_Input> <outdir> <outpre> \n\n";
#die $usage unless(@ARGV == 4);

#my %scripts = (33=>"/Users/tang58/Kai_BS/myown/brat_bw_one_step/brat_bw_one_step_one_dir_nodupl_Col0_L33.pl" ,
#	       64=>"/Users/tang58/Kai_BS/myown/brat_bw_one_step/brat_bw_one_step_one_dir_nodupl_Col0_L64.pl");

my %phred = ( 33 => 1,
	      64 =>1
	    );

my $script = "/Users/tang58/Kai_BS/myown/brat_bw_one_step/brat_bw_one_step_one_dir_nodupl_Col0_L33_64_mutiple_fq.pl";
die unless (-e $script);

print STDERR "all sub_dir should in working dir\n";

#inlist file should like
#dir	pre
#Col0_stage2_BS_rep1_S33/	Col0_stage2_BS_rep1
#Col0_stage2_BS_rep2_S17/	Col0_stage2_BS_rep2

my $usage = "\n$0 \n\n <input_list>  <L33/64> \n\n";
die $usage unless(@ARGV == 2);

my $inlist = shift or die;
my $L_num = shift or die;
#die "33/64 for L, $usage" unless (defined $scripts{$L_num} );
die "33/64 for L, $usage" unless (defined $phred{$L_num} );
#my $script = $scripts{$L_num};
#die unless (-e $script);

open(IN, $inlist) or die;
while (<IN>) {
	chomp;
	my @a = split "\t";
	next if ($a[0] eq "dir" );
	
	my ($dir, $pre ) = @a;
	print STDERR $dir, "\n";
	die $dir unless (-d $dir);
}

close IN;
print STDERR "\n\n";

open(IN, $inlist) or die;
while (<IN>) {
	chomp;
	my @a = split "\t";
	next if ($a[0] eq "dir" );
	
	my ($dir, $pre ) = @a;
	
	#my $cmd = "time $script $dir $pre";
	my $cmd = "time $script $dir $pre $L_num";

	print STDERR $cmd, "\n\n";
	unless($debug){
#		open(OUT, ">>$log_file") or die;
		#print OUT $cmd, "\n";
		#close OUT;
		`$cmd`;
	}
}

close IN;

exit;
