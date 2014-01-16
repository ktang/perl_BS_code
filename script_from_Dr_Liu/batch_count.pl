#!/usr/bin/perl -w
use strict;

my $acgt = "/home/renyi/archives/brat/brat-1.1.19/acgt-count";
my $ref = "/mnt/disk4/renyi/bis_seq/col_ros2_11_09_10/references_names.txt";

opendir(CDIR, ".");
my @dirs = readdir CDIR;
foreach my $dir(@dirs){
    if((-d $dir) && ($dir =~ /lane/)){
		print STDERR "calculating coverage under $dir\n";
		chdir($dir);
		opendir(SDIR, ".");
		my @s_files = grep {/_brat\.txt/} readdir SDIR;
		#my @p_files = grep {/paired_brat\.txt/} readdir SDIR;
		my $pre = "NONE";
		if($s_files[0] =~ /(s_\d+)/){
			$pre = $1;
		}
		my ($single_list, $pair_list) = ($pre . "_single_res.txt", $pre . "_paired_res.txt");
		open(SOUT, ">$single_list") or die "Can't open $single_list: $!";
		open(POUT, ">$pair_list") or die "Can't open $pair_list: $!";
		
		foreach my $f(@s_files){
			if($f =~ /single/){
			    print SOUT $f, "\n";
			}elsif($f =~ /paired/){
				print POUT $f, "\n";
			}
		}
		close SOUT;
		close POUT;
		
		`$acgt -r $ref -P $pre -p $pair_list -s $single_list -B`;
	    chdir('..');
	    
	 
	}
}
