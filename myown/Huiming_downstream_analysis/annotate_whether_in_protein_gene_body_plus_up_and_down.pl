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

#my $dir = "/Users/tang58/DataBase/TAIR10/TE/";

# my $usage = "$0 \n <1kb_1.5kb_2kb> <input> <output> \n\n";
# die $usage unless(@ARGV == 3);

#my $usage = "$0 \n <up_bp> <down_bp> <input> <output> \n\n";
#die $usage unless(@ARGV == 4);

my $usage = "$0 \n <up_bp> <down_bp> <input>  \n\n";
die $usage unless(@ARGV == 3);

#my $label = shift or die;

my $up_bp = shift ;#or (die unless ($up_bp == 0) );
my $down_bp = shift;# or (die unless ($down_bp == 0) );


my $input = shift or die;
#my $output = shift or die;

#my $ref_file = File::Spec->catfile($dir, "TAIR10_GFF3_protein_coding_gene_loci_5chr_sorted_" . $label . "_promoter_with_TE.txt");
my $ref_file = "/Users/tang58/DataBase/TAIR10/GFF/TAIR10_GFF3_protein_coding_gene_loci_5chr_sorted.gff";

die unless (-e $ref_file);

die unless (-e $input);
#die if (-e $output);

open(IN, $input) or die;
#open(OUT, ">$output") or die;

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

open(REF, $ref_file) or die;

my %promoters;
my %TE;
#open(TF, $TE_frag) or die "Can't open $TE_frag:$!";

# 0			1				2			3					4				5				6			7
#chr     pro_start       pro_end protein_coding_gene     gene_start      gene_end        strand  TE_in_Promoter
#chr1    1631    3630    AT1G01010       3631    5899    +       NONE

#chr1    TAIR10  gene    3631    5899    .       +       .       AT1G01010
#chr1    TAIR10  gene    5928    8737    .       -       .       AT1G01020
#chr1    TAIR10  gene    11649   13714   .       -       .       AT1G01030
# 0			1		2		3		4	 5		 6		 7			8

my $l = <REF>;
while(<REF>){
	chomp;
	my @a = split /\t/;
	my $chr = $a[0];
	
	my ($ID, $gene_start, $gene_end, $strand) = ( $a[8], $a[3], $a[4], $a[6] );
	
	my ($start_ref, $end_ref);
	
	if( $strand eq "+"){
		$start_ref = $gene_start - $up_bp;
		$end_ref   = $gene_end   + $down_bp;
	}elsif($strand eq "-"){
		$start_ref = $gene_start - $down_bp;
		$end_ref   = $gene_end   + $up_bp;
	}else{
		die $_;
	}
	
	foreach my $st(keys %{$mat_regions{$chr}}){
		my ($start, $end) = @{$mat_regions{$chr}->{$st}};
		if(!($end < $start_ref || $start > $end_ref)){
			if(defined $promoters{$chr}->{$start}){
				$promoters{$chr}->{$start} .= ";" . $ID;
			}else{
				$promoters{$chr}->{$start} = $ID;
			}
		}
	}

}
close(REF);

my $num = 0;

#print  join("\t", ($bed_head, "Gene_up" . $up_bp . "_down" . $down_bp)), "\n";
foreach my $chr(sort keys %mat_regions){
	foreach my $start(sort {$a <=> $b} keys %{$mat_regions{$chr}} ){
		my ($start, $end) = @{$mat_regions{$chr}->{$start}};
		#my $TE = "NONE";
		my $pro = "NONE";
	#	if(defined $TE{$chr}->{$start}){
	#		$TE = $TE{$chr}->{$start};
	#	}
		if(defined $promoters{$chr}->{$start}){
			$pro = $promoters{$chr}->{$start};
			$num ++;
		}
	#	print OUT join("\t", ($descs {$chr}->{$start},$pro, $TE)) , "\n";
		#print  join("\t", ($descs {$chr}->{$start},$pro)) , "\n";
	}
}

print STDERR "Gene_up" . $up_bp . "_down" . $down_bp, "\n";

print STDERR $num, "\n\n";
exit;