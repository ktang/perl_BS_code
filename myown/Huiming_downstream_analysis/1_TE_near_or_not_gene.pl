#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;

my $debug = 1;

if($debug){
	print STDERR "debug = 1\n\n";
}

my $gap = 1000;

my $ref_file = "/Users/tang58/DataBase/TAIR10/TE/TAIR10_TE_closest_protein_coding_gene_nodupl.txt";
die unless (-e $ref_file);

#0		 1		  			2		 3			  4			 5				6		       7			  8				9		   10	           11   		  12
##chr    TE_start        TE_end  TE_length       TE_ID   TE_strand       TE_Family       TE_Super_Family gene    gene_start      gene_end        gene_strand     distance
#chr1    11897   11976   80      AT1TE00010      +       ATCOPIA24       LTR/Copia       AT1G01030       11649   13714   -       0
#chr1    16883   17009   127     AT1TE00020      -       ATREP4  RC/Helitron     AT1G01030       11649   13714   -       3169
#0		 1		  2		 3			4			 5			6		7				8				9		10	 11		  12
my $usage = "$0 <input> <output> STDOUT";
die $usage unless(@ARGV == 2);

my $input = shift or die;
my $output = shift or die;

die unless(-e $input);
die if (-e $output and $output ne "/dev/null");

my %TEs;
read_input($input, \%TEs);


my ($total_TE, $A_near, $A_far, $No_near, $No_far) = (0) x 5;

open(IN, $ref_file) or die;
open(OUT, ">$output") or die;

my $head = <IN>;

print OUT $head;

while(<IN>){
	$total_TE ++;
	chomp;
	my @a = split "\t";
	if(defined $TEs{$a[4]}){
		print OUT $_, "\n";
		if($a[-1] <= $gap){
			$A_near++;
		}else{
			$A_far++;
		}
	}else{
		if($a[-1] <= $gap){
			$No_near++;
		}else{
			$No_far++;
		}
	}
}

close(IN);
close(OUT);

print "total_TE: ", $total_TE, "\n";
print "associated_near_" . $gap. "bp:	", $A_near, "\n";
print "associated_not_near_" . $gap. "bp:	", $A_far, "\n";
print "NON_associated_near_" . $gap. "bp:	", $No_near, "\n";
print "NON_associated_not_near_" . $gap. "bp:	", $No_far, "\n\n";

exit;

sub round {
    my ($number) = shift;
    #return int($number + .5);
    return int($number + .5 * ($number <=> 0)); # take care of negative numbers too
}


sub read_input{
	my ($file, $ref) = @_;
	die unless (-e $file);
	open(IN, $file) or die;
	
	while(<IN>){
		chomp;
		$ref->{$_} = 1;
	}
	
	close(IN);	
}