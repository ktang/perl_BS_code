#!/usr/bin/perl -w
use strict;

my $brat = "/home/renyi/archives/brat/brat-1.1.18/brat";
my $ref = "/mnt/disk4/renyi/bis_seq/col_ros2_11_09_10/references_names.txt";
opendir(CDIR, ".");
my @dirs = readdir CDIR;
my $log_file = "brat_single.log";
foreach my $dir(@dirs){
	
    if((-d $dir) && ($dir =~ /lane/)){
		#next if ($dir eq '.' || $dir eq '..');
        print STDERR "Map single reads under $dir ...\n";
		opendir(SDIR, $dir);
		my @files = grep {/_single\.txt/} readdir SDIR;
        foreach my $file(@files){
			my $out;
	        if($file =~ /(\S+_single)\.txt/){
				$out = $1 . "_brat.txt";
			}else{
				die "file ", $file, " name does not match pattern";
			}
	        if($file =~ /1_single\.txt/){

            `$brat -r $ref -s $dir/$file -m 2 -bs -o $dir/$out >> $log_file 2>&1`;
			}elsif($file =~ /2_single\.txt/){
				`$brat -r $ref -A -s $dir/$file -m 2 -bs -o $dir/$out >> $log_file 2>&1`;
}else{
	die "file $file does not match pattern";
}
            }
	}

}
