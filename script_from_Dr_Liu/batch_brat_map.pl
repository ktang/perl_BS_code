#!/usr/bin/perl -w
use strict;

my $debug = 0;
my $brat = "/home/renyi/archives/brat/brat-1.1.18/brat";
my $ref = "/mnt/disk4/renyi/bis_seq/col_ros2_11_09_10/references_names.txt";
opendir(CDIR, ".");
my @dirs = readdir CDIR;
my $log_file = "brat_map.log";
foreach my $dir(@dirs){
    if((-d $dir) && ($dir =~ /lane/)){
		#next if ($dir eq '.' || $dir eq '..');
        print STDERR "Map paired reads under $dir ...\n";
		opendir(SDIR, $dir) or die "Can't open $dir: $!";
		my @files1 = grep {/1_pairs\.txt$/} readdir SDIR;
		if($debug){
			print STDERR "files found: ", join(" ", @files1), "\n";
		}
        foreach my $file1(@files1){
			if($debug){
				print STDERR "doing file $file1\n";
			}
			my ($file2, $out);
	        if($file1 =~ /(s_\d+_)1_pairs\.txt/){
				$file2 = $1 . "2_pairs.txt";
				$out = $1 . "paired_brat.txt";
			}else{
				die "file ", $file1, " name does not match pattern";
			}
	        
	
            `$brat -r $ref -1 $dir/$file1 -2 $dir/$file2 -u -i 0 -a 1000 -m 2 -pe -bs -o $dir/$out >> $log_file 2>&1`;
		}
	}
}
