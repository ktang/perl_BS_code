#!/usr/bin/perl -w

#v1.2 output percentage.

# count TE number in annotated file
#/Users/kaitang/Desktop/Nimblegen/Sep_20_new_analysis/result/send/label_ZDP_peaks_with_annotation.xls

#v1.1
# 8 kind of numbers:
# TE 	gene intergeneic
# 0/1   0/1 	0/1
#wrong
#chr     start   end     mut_num WT_num  Gene    Strand  GeneType        GeneAnnot       TE      TEFamily        Intergenic
#										 -7	(6)		-6		-5				-4			 -3			-2				-1

#Chr	start	end	number	Gene	Strand	GeneType	GeneAnnot	TE	TEFamily	Intergenic	Repeats	overlap_rdd	DT-id
#0        1       2     3    4        5       6           7          8    9            10         11       12        13
# inter     gene     TE
#	0		 0		  0		a0	NONE
#   0		 0		  1		a1	TE
#	0		 1		  0		a2	gene
#	0		 1		  1		a3	TE
#   1		 0		  0		a4	inter
#	1		 0		  1		a5	TE
#	1		 1		  0		a6	gene
#	1		 1		  1		a7	TE

use strict;




my $usage = "$0 <input> STDOUT";

die $usage unless(@ARGV == 1);

my $input = $ARGV[0];
die "wrong input" unless (-e $input);

open (IN, $input) or die "cannot open $input: $!";

my $head = <IN>;
chomp $head;
my @hs = split "\t", $head;

my $gene_label_index = -1;
for my $i(0..$#hs){
	if($hs[$i] eq "Gene" ){
		$gene_label_index  = $i;
		last;
	}
}
die if($gene_label_index == -1);

my $gene_type_index = $gene_label_index + 2;
my $TE_label_index  = $gene_label_index + 4;
my $inter_label_index = $gene_label_index + 6;


die unless ($hs[$gene_label_index] eq "Gene");
die unless ($hs[$gene_type_index] eq "GeneType");
die unless ($hs[$TE_label_index] eq "TE");
die unless ($hs[$inter_label_index] eq "Intergenic");
  
my $total = 0;

my ($te_num, $gene_num, $inter_num) = (0) x 3;

my $other_num = 0;
my $pseudogene_num = 0;
my $promoter_num = 0;
my $all_inter = 0;

# my @nums = (0, 0, 0, 0, 0, 0, 0, 0);
#my $gene_label_index = 4;
#my $gene_type_index = $gene_label_index + 1;
#my $TE_label_index  = $gene_label_index + 4;
#my $inter_label_index = $gene_label_index + 6;

my $te_pro_num = 0;
my $te_inter_num = 0;

while( <IN> ){
	chomp;
	$total ++;

	my @a = split /\t/;

	if( $a[$gene_type_index] =~ /protein/  &&   $a[$inter_label_index] ne "NONE" ){		$promoter_num++;	}
	if($a[$inter_label_index] ne "NONE" ) {$all_inter++;}
	
	if($a[$TE_label_index] ne "NONE"){
		$te_num ++;
		if($a[$gene_type_index] =~ /protein/){
			$te_pro_num ++;
		}
		
		if($a[$inter_label_index] ne "NONE" ){
			$te_inter_num ++;
		}
		
		
		
	}
	elsif( $a[$gene_label_index] ne "NONE"){
		if($a[$gene_type_index] =~ /transposable/){
			$te_num++;
			if($a[$gene_type_index] =~ /protein/){
				$te_pro_num ++;
			}
			if($a[$inter_label_index] ne "NONE" ){
				$te_inter_num ++;
			}
		}elsif(  $a[$gene_type_index] =~ /protein/){
			$gene_num++;
		}elsif( $a[$gene_type_index] =~ /pseudogene/ ){
			$pseudogene_num++;
		}else{
			$other_num ++;
			print STDERR $_, "\n";
		}
		
		
	}elsif( $a[$inter_label_index] ne "NONE"  ){
		$inter_num++;
	}else{
		die $_;
	}
	
}
close (IN);


  
print "total: $total\n";
#my $te_num = $nums[1] + $nums[3] + $nums[5] + $nums[7];
#my $gene_num = $nums[2] + $nums[6];
#my $inter_num = $nums[4];

if ($other_num + $te_num + $gene_num + $inter_num + $pseudogene_num != $total) {print STDERR  "wrong:\n\n sum != total"}
  
my $te_per = sprintf("%.1f", 100 * $te_num / $total);
my $gene_per = sprintf("%.1f", 100 * $gene_num / $total);
my $other_per = sprintf("%.1f", 100 * $other_num / $total);

my $pseudogene_per = sprintf("%.1f", 100 * $pseudogene_num / $total);

my $inter_per = 100 - $te_per - $gene_per - $other_per;



print join("\t", ("TE", $te_num, $te_per."%")), "\n";
print join("\t", ("gene", $gene_num, $gene_per."%")), "\n";
print join("\t", ("intergenic", $inter_num, $inter_per."%")), "\n";
print join("\t", ("other", $other_num, $other_per."%")), "\n";
print join("\t", ("pseudogene", $pseudogene_num, $pseudogene_per."%")), "\n\n";

print "TE+inter:\n", $te_inter_num, "\n";
print "protein+TE:\n", $te_pro_num, "\n";
print "all_intergenic_number:\n", $all_inter, "\n";
print "potential_promoter:\n", $promoter_num, "\n\n";


#  print STDERR "total = $total\nTE=$te\nNone=$none\ngene=$gene\ninter=$intergenic\n";


  exit;
