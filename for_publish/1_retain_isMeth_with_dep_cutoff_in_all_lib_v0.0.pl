#!/usr/bin/perl -w

#raw name: multiple_isMeth2wig_with_depth_cutoff_v0.2.pl

#v0.1 
# output wig file and total number of cytosines meet the requirement.

#v0.2
# output isMet file itself but not wig file

#

use strict;
use File::Spec;

my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}
#my $usage = "$0 \n <outdir> <cutoff> <number_of_files> <isMeth> <label> [<isMeth> <label> ...] \n\n";
#die $usage unless(@ARGV >= 7);

my $usage = "$0 \n <postfix> <outdir> <cutoff> <number_of_files> <isMeth> <label> [<isMeth> <label> ...] \n\n";
die $usage unless(@ARGV >= 8);


my $postfix = shift or die;
my $outdir = shift or die;
my $dep_cutoff = shift or die;
my $num_of_files = shift or die;
my $last_index = $num_of_files - 1;

die unless (-d $outdir);

die unless (@ARGV == 2 * $num_of_files);

my @labels;
my @files;
my @outputs;

for my $i(0..$last_index){
	$files[$i]  = shift or die;
	$labels[$i] = shift or die;
	$outputs[$i] = File::Spec->catfile($outdir,  $labels[$i] . "_isMeth_depth" . $dep_cutoff . "_" . $postfix . ".txt" );
	die $files[$i] unless (-e $files[$i]);
	die if (-e $outputs[$i]);
}

print STDERR "depth_cutoff is $dep_cutoff \n\n";
#if($debug){
	print STDERR join("\n", @labels),  "\n\n";
	print STDERR join("\n", @files),   "\n\n";
	print STDERR join("\n", @outputs), "\n\n";
#}

if($debug){
	exit;
}
my @fhr;
my @fhw;
for my $i(0..$last_index){
	open($fhr[$i], "<" , $files[$i]) or die "$files[$i]";
	open($fhw[$i], ">" , $outputs[$i]) or die $outputs[$i], ": $!";

	my $head = readline($fhr[$i]);
	
	chomp $head;
	
	if($debug){
#		print STDERR $head, "\n\n";
	}
	
	my @h = split "\t", $head;
	die unless ($h[3] eq "type" and $h[5] eq "depth" and $h[7] eq "isMeth" );
	
	print {$fhw[$i]} $head, "\n";
}

my $l;
my %wig_pos;

my $total_num = 0;


#0		 1			2		3		4		5		6				7
#chr     pos     strand  type    num_C   depth   percentage      isMeth
#chr1    1       +       CHH     0       0       0       0
while($l = readline($fhr[0]) ){
	
	my $flag_less = 0;
	
	my @a;
	
	chomp $l;
	@{$a[0]} = split "\t", $l;
	
	my ($chr, $pos) = ( $a[0]->[0], $a[0]->[1] );
#	my $strand = ($a[0]->[2] eq "+") ? 1 : -1; 
	
	if($a[0]->[5] < $dep_cutoff){
		$flag_less = 1;
	}
	
	for my $i(1..$last_index){
		$l = readline($fhr[$i]);
		chomp $l;
		@{$a[$i]} = split "\t", $l;
		die $i,"\n", $l, "\n" unless ($a[$i]->[1] == $pos);
		
		if($a[$i]->[5] < $dep_cutoff){
			$flag_less = 1;
		}

	}
	
	if($flag_less == 0){
		$total_num++;
		for my $j(0..$last_index){
			#if($a[$j]->[7] == 1){
			#	$wig_pos{$j}->{$chr}->{$pos} = $strand * $a[$j]->[6];
			#}
			print {$fhw[$j]} join("\t", @{$a[$j]}) , "\n";
		}
	}
}

for my $i(0..$last_index){
	close($fhr[$i]);
	close($fhw[$i]);
}

print STDERR "\nnumber of cytosines meet the requirement all depth >= $dep_cutoff:\n";
print STDERR $total_num, "\n\n";

exit;
