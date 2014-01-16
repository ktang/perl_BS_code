#!/usr/bin/perl -w

my $qual = 20;
if(@ARGV > 0){
	$qual = $ARGV[0];
}
opendir(CDIR, ".");
my @dirs = readdir CDIR;
foreach my $dir(@dirs){
	
	#if($dir =~ /Lane1_JKZ1_Col0/){
	#	next;
	#}
    if((-d $dir) && ($dir =~ /lane/)){
		next if ($dir eq '.' || $dir eq '..');
        print STDERR "checking quality under $dir ...\n";

	    chdir($dir);
	    #print STDERR "checking quality under $dir ...\n";
	
	    `perl ~/archives/SolexaQA/SolexaQA.pl *.txt -h $qual`;
	    chdir('..');
	}
}
