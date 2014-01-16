#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;

my $debug = 1;

if($debug){
	print STDERR "debug = 1\n\n";
}

#my $te_file = "/Users/tang58/DataBase/TAIR10/TE/TAIR10_Transposable_Elements_sorted_bed.txt";
#die unless (-e $te_file);

my $dir = "/Users/tang58/DataBase/TAIR10/TE/";

my $usage = "$0 \n <1kb_1.5kb_2kb> <input> <output> \n\n";
die $usage unless(@ARGV == 3);

my $label = shift or die;
my $input = shift or die;
my $output = shift or die;

my $ref_file = File::Spec->catfile($dir, "TAIR10_GFF3_protein_coding_gene_loci_5chr_sorted_" . $label . "_promoter_with_TE.txt");
die unless (-e $ref_file);

die unless (-e $input);
die if (-e $output);

open(IN, $input) or die;
open(OUT, ">$output") or die;

#Transposon_Name orientation_is_5prime   Transposon_min_Start    Transposon_max_End      Transposon_Family       Transposon_Super_Family
#AT1TE00010      true    11897   11976   ATCOPIA24       LTR/Copia
#	0				1		2		3		4				5

my $bed_head = <IN>;
chomp $bed_head;
my %descs;# record the whole line

#read input file
my %mat_regions; # chr->start=\@
while(<IN>){
	next if(/^browser/);
	next if(/^track/);
	next unless(/\w+/);
	chomp;
	my ($chr, $start, $end) = split /\t/;
	$chr = lc $chr;
	if($chr eq "chloroplast"){
		$chr = "chrc";
	}
	if($chr eq "mitochondria"){
		$chr = "chrm";
	}
	$mat_regions{$chr}->{$start} = [$start, $end];
	$descs {$chr}->{$start} = $_;
}

close(IN);

open(TE, $ref_file) or die;

my %promoters;
my %TE;
#open(TF, $TE_frag) or die "Can't open $TE_frag:$!";

# 0			1				2			3					4				5				6			7
#chr     pro_start       pro_end protein_coding_gene     gene_start      gene_end        strand  TE_in_Promoter
#chr1    1631    3630    AT1G01010       3631    5899    +       NONE
my $l = <TE>;
while(<TE>){
	chomp;
	my @temp = split /\t/;
	my $chr = $temp[0];
	
	foreach my $st(keys %{$mat_regions{$chr}}){
		my ($start, $end) = @{$mat_regions{$chr}->{$st}};
		if(!($end < $temp[1] || $start > $temp[2])){
			if(defined $promoters{$chr}->{$start}){
				$promoters{$chr}->{$start} .= ":" . $temp[3];
				$TE{$chr}->{$start} .= ":" . $temp[7];
			}else{
				$promoters{$chr}->{$start} = $temp[3];
				$TE{$chr}->{$start} = $temp[7];
			}
		}
	}

}
close(TE);

print OUT join("\t", ($bed_head, "Promoter_". $label, "TE_in_Promoter")), "\n";
foreach my $chr(sort keys %mat_regions){
	foreach my $start(sort {$a <=> $b} keys %{$mat_regions{$chr}} ){
		my ($start, $end) = @{$mat_regions{$chr}->{$start}};
		my $TE = "NONE";
		my $pro = "NONE";
		if(defined $TE{$chr}->{$start}){
			$TE = $TE{$chr}->{$start};
		}
		if(defined $promoters{$chr}->{$start}){
			$pro = $promoters{$chr}->{$start};
		}
		print OUT join("\t", ($descs {$chr}->{$start},$pro, $TE)) , "\n";
	}
}

exit;