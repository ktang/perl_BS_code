#!/usr/bin/perl -w



use strict;

my $usage = "$0\n <input> <label> <head_yes/no> STDOUT \n\n";

die $usage unless(@ARGV == 3);

my $input = shift or die;
die "wrong input" unless (-e $input);

my $label = shift or die;

my $head_label = shift or die;
die $usage, "\nyes or no\n\n" unless ($head_label =~ /yes|no/i);




if($head_label =~ /yes/i){
	print join("\t", ("sample", "total", "TE", "gene","intergenic", "pseudogene", "other", "TE+protein","TE+inter", "all_intergenic_number"  )), "\n";
}


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



my $te_pro_num = 0;
my $te_inter_num = 0;

while( <IN> ){
	chomp;
	$total ++;

	my @a = split /\t/;

#	if( $a[$gene_type_index] =~ /protein/  &&   $a[$inter_label_index] ne "NONE" ){		$promoter_num++;	}
	if($a[$inter_label_index] ne "NONE" ) {$all_inter++;}
	
	if($a[$TE_label_index] ne "NONE"){
		$te_num ++;
		if($a[$gene_type_index] =~ /protein/){
			$te_pro_num ++;
		}
		
		if($a[$inter_label_index] ne "NONE" ){
			$te_inter_num ++;
		}
	}elsif( $a[$gene_label_index] ne "NONE"){
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
		#	print STDERR $_, "\n";
		}
	}elsif( $a[$inter_label_index] ne "NONE"  ){
		$inter_num++;
	}else{
		die $_;
	}
}
close (IN);



if ($other_num + $te_num + $gene_num + $inter_num + $pseudogene_num != $total) {print STDERR  "wrong:\n\n sum != total"}
  
my $te_per = sprintf("%.1f", 100 * $te_num / $total);
my $gene_per = sprintf("%.1f", 100 * $gene_num / $total);
my $other_per = sprintf("%.1f", 100 * $other_num / $total);

my $pseudogene_per = sprintf("%.1f", 100 * $pseudogene_num / $total);

my $inter_per = sprintf("%.1f", 100 - $te_per - $gene_per - $other_per -  $pseudogene_per) ;



my $TE 	 		=   $te_num  . " (" .  $te_per."%)";
my $gene 		=  $gene_num . " (" . $gene_per."%)"; 
my $intergenic  = $inter_num .  " (". $inter_per."%)" ;
my $pseudogene  = $pseudogene_num.  " (".  $pseudogene_per."%)";
my $other       = $other_num  .  " (".      $other_per. "%)";

print join("\t", ( $total,  $TE, $gene, $intergenic, $pseudogene, $other, $te_pro_num, $te_inter_num, $all_inter, $label)), "\n";



exit;
