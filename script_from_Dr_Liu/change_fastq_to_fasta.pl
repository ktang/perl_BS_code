#!/usr/bin/perl -w
use strict;

opendir(CDIR, ".");
my @dirs = readdir CDIR;
my $log_file = "fastq_to_fasta.log";
foreach my $dir(@dirs){
	
	#if($dir =~ /Lane1_JKZ1_Col0/){
	#	next;
	#}
    if((-d $dir) && ($dir =~ /lane/)){
		#next if ($dir eq '.' || $dir eq '..');
        print STDERR "change fastq files under $dir ...\n";

	    #chdir($dir);
	    #print STDERR "checking quality under $dir ...\n";
	    opendir(SDIR, $dir);
		my @files = grep {/trimmed$/} readdir SDIR;
		foreach my $file(@files){
			my @pre = split /\./, $file;
			my $outfile = $pre[0] . "_trimmed.fa";
			`perl fastq2fasta.pl $dir/$file $dir/$outfile >> $log_file`;
		}
	}
}
