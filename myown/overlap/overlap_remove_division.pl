#!/usr/bin/perl -w

#v0.0
# length overlap

use strict;
use File::Spec;


my $debug = 0;


#############################
my $p = 1;
if($debug){
	print STDERR "debug = $debug\n\n";
}

#my $script = "/Users/tang58/Weiqiang_idm1_1_Nature_paper/add_annotation_Kai_v1.7_TAIR10.pl";
#my $script  = "./overlap_bed_first_in_second.pl";
#die unless (-e $script);

my $usage = "$0 \n <input>  \n\n";

die $usage unless (@ARGV == 1);

my $input = shift or die;
die unless (-e $input);
die unless ( $input =~ /overlap/);


my $output ;

if (  $input =~ /(\S+)\.txt$/ ){
	$output = $1 . "_simple.txt";
}

die if(-e $output);
open(IN, $input) or die;
open(OUT, ">>$output") or die;

my $head = <IN>;
print OUT $head;

while (<IN>){
	chomp;
	my @a = split "\t";
	for my $i(1..$#a){
		$a[$i] = remove($a[$i]);
	}
	print OUT join("\t", @a), "\n";
}

close(IN);
close(OUT);


exit;

sub remove {
	my ($foo) = @_;
	
	if ($foo =~ /\d+\s*\/\s*\d+\s*=\s*(\S*)/){
		return $1;
	}elsif($foo eq "--"){
		return $foo;
	}else{
		die $foo;
	}
}