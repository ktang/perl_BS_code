#!/usr/bin/perl -w
#outformat
#chr start end CG_num CG_total mCG CG_level CHG . . . CHH . . . 
use strict;
use File::Spec;

my $debug = 1;
if($debug){
	print STDERR "debug = 1\n\n";
}
my $bin_size     = 200;
my $sliding_size = 50;
my $const = 4000000;

my $covered_num = int ($bin_size / $sliding_size);


my $dep_cutoff = 2;

my %chr_len = ("chr1"=>30427671, "chr2"=>19698289, "chr3"=>23459830, "chr4"=>18585056, "chr5"=>26975502);

my %bin_last_index;
# bin start from 0
for my $chr (sort keys %chr_len){
	my $len = $chr_len{$chr};
	my $last = int( ($len - 1) / $sliding_size);
	$bin_last_index{$chr} = $last;
}
if($debug){
	foreach my $chr (sort keys %bin_last_index){
		print STDERR join("\t", ($chr, $bin_last_index{$chr})), "\n";
	}
}


my @chrs = ("chr1", "chr2", "chr3", "chr4", "chr5");

#my $usage = "$0 <input> <prefix_output> <outdir>\n/Volumes/Macintosh_HD_2/idm1_new_met_data/UCR_server_data/Lane1_JKZ1_Col0/s_1_isMeth.txt";
#label can be ros1,ros2,ros3,ros4,rdd,hsp,phd,zdp
#die $usage unless (@ARGV == 3);

my $usage = "$0 <WT_isMeth> <mut_isMeth> <outpre> <label>";
die $usage unless (@ARGV == 4);

my $wt_in  = shift or die;
my $mut_in = shift or die;
my $outpre = shift or die;
my $label  = shift or die;

die $wt_in  unless (-e $wt_in);
die $mut_in unless (-e $mut_in);

my $output = $outpre."_BinSize". $bin_size . "_sliding" . $sliding_size . ".txt";
print STDERR "output: $output\n\n";

die "output $output exist:$!" if (-e $output);

my (%CG_num, %CHG_num, %CHH_num);
my (%CG_depth, %CHG_depth, %CHH_depth);
my (%mCG, %mCHG, %mCHH);


if ($debug){
	$wt_in =  "/Volumes/My_Book/20120427_ShangHai_data/src/method5/debug/C24_WT_isMeth_chrC_error_separately_called_10000.txt";
	$mut_in = "/Volumes/My_Book/20120427_ShangHai_data/src/method5/debug/1-5_isMeth_chrC_error_separately_called_10000.txt";
 #	$output = "detail_debug_output.txt"; 
	$bin_size = 50;
	$sliding_size = 25;
	%chr_len = ("chr1"=> 58344);
	%bin_last_index = ('chr1' => int( ($chr_len{'chr1'} - 1) / $sliding_size));
}

open (OUT, ">$output") or die "cannot open $output: $!";
open (WT,  $wt_in )  or die "cannot open  $wt_in: $!";
open (MUT, $mut_in ) or die "cannot open $mut_in: $!";

my (%dep_num_wt, %dep_num_mut);
my (%mC_num_wt, %mC_num_mut);
my %ref_C_num;

my ($l_wt, $l_mut);
#chr     pos     strand  type    num_C   depth   percentage      isMeth
#chr1    1       +       CHH      0       0       0			       0
# 0      1       2         3      4      5       6                 7

my $h_wt  = <WT>;
my $h_mut = <MUT>;



#my $bin_size     = 200;
#my $sliding_size = 50;
#my $covered_num = int ($bin_size / $sliding_size);

my $line = 0;

while($l_wt = <WT>, $l_mut = <MUT>){
	$line++;
	if($line % $const == 0){
		print STDERR $line, "\tdone\n";
	}
	chomp $l_wt;
	chomp $l_mut;
	my @a_wt  = split "\t", $l_wt;
	my @a_mut = split "\t", $l_mut;
	
	die unless ($a_wt[1] == $a_mut[1]);
	my $chr = $a_wt[0];
	my $type = $a_wt[3];

	
	my $dep_wt = $a_wt[5];
	my $dep_mut = $a_mut[5];
	
	next if($chr =~ /chrC|chrM|pUC19/);
		
	if($dep_wt >= $dep_cutoff and $dep_mut >= $dep_cutoff){
		my $pos = $a_wt[1];
		my $mC_wt = $a_wt[4];
		my $mC_mut = $a_mut[4];
		
		my $base_index = int( ($pos - 1) / $sliding_size);
		
		for my $i (($base_index - $covered_num)..($base_index + $covered_num)){
			if($i >= 0){
				my $interval_start = $i * $sliding_size + 1;
				my $interval_end   = $interval_start + $bin_size - 1;
				if($pos >= $interval_start and $pos <= $interval_end){
					#$ref_den->{$chr}->[$i]++;
					$dep_num_wt{$type}->{$chr}->[$i]  += $dep_wt ;
					$dep_num_mut{$type}->{$chr}->[$i] += $dep_mut;
					
					$mC_num_wt{$type}->{$chr}->[$i]   += $mC_wt;
					$mC_num_mut{$type}->{$chr}->[$i]  += $mC_mut;
					$ref_C_num{$type}->{$chr}->[$i] ++;
				}
			}
		}
	}
}

#my (%dep_num_wt, %dep_num_mut);
#my (%mC_num_wt, %mC_num_mut);
#my %ref_C_num;


print OUT join("\t", ("chr", 'start', 'end', 
'CG_num', 'CG_wt', 'mCG_level_wt', 'CG_'.$label, 'mCG_level_'.$label, 
'CHG_num', 'CHG_wt', 'mCHG_level_wt', 'CHG_'.$label, 'mCHG_level_'.$label, 
'CHH_num', 'CHH_wt', 'mCHH_level_wt', 'CHH_'.$label, 'mCHH_level_'.$label, 
'C_num', 'C_wt', 'mC_level_wt', 'C_'.$label, 'mC_level_'.$label, 
)),"\n";

my @types = ("CG", "CHG", "CHH");

foreach my $chr(sort keys %chr_len){
	my $last = $bin_last_index{$chr};
	
	for (my $i = 0; $i < $bin_last_index{$chr}; $i++){
		
		my $start = $i * $sliding_size + 1;;
		my $end   = $start + $bin_size - 1;
		print OUT join("\t", ($chr, $start, $end)) , "\t";
		
		my ($total_C_num, $total_mC_wt, $total_dep_wt) = (0) x 3;
		my ( $total_mC_mut, $total_dep_mut ) = (0) x 2;
 		
		foreach my $type (@types){
			my ($ref_num, $level_wt, $level_mut) = (0) x 3;
			my ($mC_wt , $dep_wt ) = (0) x 2;
			my ($mC_mut, $dep_mut) = (0) x 2;
			if(defined $ref_C_num{$type}->{$chr}->[$i])  { $ref_num = $ref_C_num{$type}->{$chr}->[$i] }
			if(defined $dep_num_wt{$type}->{$chr}->[$i]) { $dep_wt = $dep_num_wt{$type}->{$chr}->[$i] }
			if(defined $dep_num_mut{$type}->{$chr}->[$i]){ $dep_mut = $dep_num_mut{$type}->{$chr}->[$i] }
			if(defined $mC_num_wt{$type}->{$chr}->[$i])  { $mC_wt = $mC_num_wt{$type}->{$chr}->[$i] }
			if(defined $mC_num_mut{$type}->{$chr}->[$i]) { $mC_mut = $mC_num_mut{$type}->{$chr}->[$i] }
			$total_C_num   += $ref_num;
			$total_mC_wt   += $mC_wt ;
			$total_dep_wt  += $dep_wt;
			$total_mC_mut  += $mC_mut;
			$total_dep_mut += $dep_mut;
			
			my $per_wt = 0;
			my $per_mut = 0;
			if($dep_wt !=0){$per_wt = sprintf("%.2f", 100 * $mC_wt / $dep_wt) }
			
			if($dep_mut != 0){ $per_mut = sprintf("%.2f", 100 * $mC_mut / $dep_mut)}
			
			print OUT join("\t", ($ref_num, "$mC_wt/$dep_wt=", $per_wt, "$mC_mut/$dep_mut=", $per_mut)), "\t";
		}
#	if ($chh_depth != 0) {$chh_level = sprintf("%.2f",100 * $mchh / $chh_depth) }
			my $total_per_wt = 0;
			my $total_per_mut = 0;
			if($total_dep_wt !=0){$total_per_wt = sprintf("%.2f", 100 * $total_mC_wt / $total_dep_wt) }
			
			if($total_dep_mut != 0){ $total_per_mut = sprintf("%.2f", 100 * $total_mC_mut / $total_dep_mut)}
			
			print OUT join("\t", ($total_C_num, "$total_mC_wt/$total_dep_wt=", $total_per_wt, "$total_mC_mut/$total_dep_mut=", $total_per_mut)), "\n";
 
	#print OUT join("\t", ($chr, $i * $bin_size + 1, ($i+1) * $bin_size , $cg_num, $cg_depth, $mcg, $cg_level,  $chg_num, $chg_depth, $mchg, $chg_level,  $chh_num, $chh_depth, $mchh, $chh_level)), "\n";
	}
	
	
}
close(OUT);
exit;


sub round {
    my($number) = shift;
    #return int($number + .5);
    return int($number + .5 * ($number <=> 0)); # take care of negative numbers too
}
