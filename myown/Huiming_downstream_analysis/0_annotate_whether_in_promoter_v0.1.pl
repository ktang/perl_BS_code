#!/usr/bin/perl -w

#just count

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;
use Statistics::R;


my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}

#my $te_file = "/Users/tang58/DataBase/TAIR10/TE/TAIR10_Transposable_Elements_sorted_bed.txt";
#die unless (-e $te_file);


#my $usage = "$0 \n <1kb_1.5kb_2kb> <input> <output> \n\n";
my $usage = "$0 \n <1kb_1.5kb_2kb> <input_hyper_hypo_list> STDOUT \n\n";
die $usage unless(@ARGV == 2);

my $label = shift or die;
my $input = shift or die;
#my $output = shift or die;

my $dir = "/Users/tang58/DataBase/TAIR10/TE/";
die unless (-d $dir);
my $ref_file = File::Spec->catfile($dir, "TAIR10_GFF3_protein_coding_gene_loci_5chr_sorted_" . $label . "_promoter_with_TE.txt");

#my $dir = "/Users/tang58/DataBase/TAIR10/TE/set2";
#die unless (-d $dir);
#my $ref_file = File::Spec->catfile($dir, "TAIR10_GFF3_protein_coding_gene_loci_5chr_sorted_" . $label . "_promoter_with_TE_greater100.txt");



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
	
	$mat_regions{$chr}->{$start} = [$start, $end];
	#$descs {$chr}->{$start} = $_;
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

my ($total, $P_TE, $P_NoTE, $NoP_TE, $NoP_NoTE) = (0) x 5;

while(<TE>){
	$total++;
	chomp;
	my @temp = split /\t/;
	my $chr = $temp[0];
	
	my $flag = 0;
	
	foreach my $st(keys %{$mat_regions{$chr}}){
		my ($start, $end) = @{$mat_regions{$chr}->{$st}};
		if(!($end < $temp[1] || $start > $temp[2])){
			$flag = 1;
			last;
		}
	}
	
	if($flag){
		if($temp[-1] eq "NONE"){
			$P_NoTE++;
		}else{
			$P_TE++;
		}
	}else{
		if($temp[-1] eq "NONE"){
			$NoP_NoTE++;
		}else{
			$NoP_TE++;
		}
	}
}
close(TE);

print "total_gene: ", $total, "\n";
print "dependent_promoter_with_TE:\t", $P_TE , "\n";
print "dependent_promoter_without_TE:\t", $P_NoTE , "\n";
print "independent_promoter_with_TE:\t", $NoP_TE , "\n";
print "independent_promoter_without_TE:\t", $NoP_NoTE , "\n\n";

my $R = Statistics::R->new();
my $p_value = 1;
			
		    $R->set('a', $P_TE);
    		$R->set('b', $P_NoTE);
   			$R->set('c', $NoP_TE);
   			$R->set('d', $NoP_NoTE);
			$R->run(q`p = fisher.test(matrix( c(a,b,c,d ), ncol=2))$p.value`);
			
			$p_value = $R->get('p');


$R->stop();

print "p-value = $p_value\n\n";


exit;