#!/usr/bin/perl -w
 
use strict;
use File::Spec;
use Text::NSP::Measures::2D::Fisher2::twotailed;

my $debug = 0;
if($debug){
	print STDERR "debug = 1\n\n";
}
my $bin_size     = 200;
my $sliding_size = 50;
my $const = 4000000;
my $allowed_gap  = 200;
my $covered_num = int ($bin_size / $sliding_size);

my $DMC_cutoff  = 10;
my $depth_cutoff = 2;
my $p_value_cutoff = 0.05;
my %chr_len = ("chr1"=>30427671, "chr2"=>19698289, "chr3"=>23459830, "chr4"=>18585056, "chr5"=>26975502);
my $R_script = "/Users/tang58/Kai_BS/myown/select_DMR_Jacobsen_method/p_adjust_BH.R";
die unless (-e $R_script);

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

#my $usage = "$0 <WT_isMeth> <mut_isMeth> <outdir> <outpre> <label>";
#die $usage unless (@ARGV == 5);

my $usage = "$0 <WT_isMeth> <mut_isMeth> <outdir> <outpre>";
die $usage unless (@ARGV == 4);

my $wt_in  = shift or die;
my $mut_in = shift or die;
my $outdir = shift or die;
my $outpre = shift or die;
#my $label  = shift or die;
die unless (-d $outdir);
die $wt_in  unless (-e $wt_in);
die $mut_in unless (-e $mut_in);

my $output = File::Spec->catfile($outdir, $outpre."_BinSize". $bin_size . "_sliding" . $sliding_size . "_BH_adj_p_value". $p_value_cutoff . "_depth". $depth_cutoff. ".txt");
print STDERR "output: $output\n\n";
my $temp_p_file =File::Spec->catfile($outdir, $outpre."_BinSize". $bin_size . "_sliding" . $sliding_size . "_BH_adj_p_value". $p_value_cutoff . "_depth". $depth_cutoff. "_p.txt"); 

my $temp_p_adj_file =File::Spec->catfile($outdir, $outpre."_BinSize". $bin_size . "_sliding" . $sliding_size . "_BH_adj_p_value". $p_value_cutoff . "_depth". $depth_cutoff. "_p_adj.txt"); 
die "output $output exist:$!" if (-e $output);
die if(-e $temp_p_file);
die if(-e $temp_p_adj_file);
my (%CG_num, %CHG_num, %CHH_num); # number of CXX that meet the requirement that both depth >= $depth_cutoff
my (%CG_depth, %CHG_depth, %CHH_depth);#depth sum
my (%mCG, %mCHG, %mCHH);               #mCXX sum


if ($debug){
	$wt_in =  "/Volumes/My_Book/20120427_ShangHai_data/src/method5/debug/C24_WT_isMeth_chrC_error_separately_called_10000.txt";
	$mut_in = "/Volumes/My_Book/20120427_ShangHai_data/src/method5/debug/1-5_isMeth_chrC_error_separately_called_10000.txt";
 #	$output = "detail_debug_output.txt"; 
	$bin_size = 50;
	$sliding_size = 25;
	%chr_len = ("chr1"=> 58344);
	%bin_last_index = ('chr1' => int( ($chr_len{'chr1'} - 1) / $sliding_size));
}

open (WT,  $wt_in )  or die "cannot open  $wt_in: $!";
open (MUT, $mut_in ) or die "cannot open $mut_in: $!";

die "output $output exist:$!" if (-e $output);
open (OUT, ">$output") or die "cannot open $output: $!";
die if(-e $temp_p_file);
open (P, ">$temp_p_file ") or die "cannot open $temp_p_file: $!";
die if(-e $temp_p_adj_file);
#open (P_ADJ, ">$temp_p_adj_file ") or die "cannot open $temp_p_adj_file :$!";
print OUT join("\t", ( "chr", "start", "end", "DMC_num","mut_CG", "mut_CG_per","mut_CHG", "mut_CHG_per","mut_CHH", "mut_CHH_per",
				"wt_CG", "wt_CG_per",	"wt_CHG", "wt_CHG_per","wt_CHH", "wt_CHH_per")), "\n";
print P join("\t", ("chr", "start", "end", "p_CG", "p_CHG", "p_CHH")), "\n";


my (%dep_num_wt, %dep_num_mut);
my (%mC_num_wt, %mC_num_mut);
my %ref_C_num;

my %DMCs;

my ($l_wt, $l_mut);
#chr     pos     strand  type    num_C   depth   percentage      isMeth
#chr1    1       +       CHH      0       0       0	           0
# 0      1       2         3      4       5       6                7

my $h_wt  = <WT>;
my $h_mut = <MUT>;

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
		
	if($dep_wt >= $depth_cutoff and $dep_mut >= $depth_cutoff){
		my $pos = $a_wt[1];
		my $mC_wt = $a_wt[4];
		my $mC_mut = $a_mut[4];
	
		my $per_wt  = $a_wt[6];
		my $per_mut = $a_mut[6];
		if($per_mut != 0){
			if($per_wt / $per_mut <= 0.5){
				$DMCs{$chr}->[$pos] = $type;
			}
		}
			
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
close(WT);
close(MUT);

my @types = ("CG", "CHG", "CHH");

foreach my $chr(sort keys %chr_len){
	my $last = $bin_last_index{$chr};
	
	for (my $i = 0; $i <= $bin_last_index{$chr}; $i++){
		
		my $start = $i * $sliding_size + 1;;
		my $end   = $start + $bin_size - 1;
		
		my ($total_C_num, $total_mC_wt, $total_dep_wt) = (0) x 3;
		my ( $total_mC_mut, $total_dep_mut ) = (0) x 2;
 		my @p_array = ();
		foreach my $type (@types){
			my $ref_num = 0; 
			my ($mC_wt , $dep_wt ) = (0) x 2;
			my ($mC_mut, $dep_mut) = (0) x 2;
			if(defined $ref_C_num{$type}->{$chr}->[$i])  { $ref_num = $ref_C_num{$type}->{$chr}->[$i] }
			if(defined $dep_num_wt{$type}->{$chr}->[$i]) { $dep_wt = $dep_num_wt{$type}->{$chr}->[$i] }
			if(defined $dep_num_mut{$type}->{$chr}->[$i]){ $dep_mut = $dep_num_mut{$type}->{$chr}->[$i] }
			if(defined $mC_num_wt{$type}->{$chr}->[$i])  { $mC_wt = $mC_num_wt{$type}->{$chr}->[$i] }
			if(defined $mC_num_mut{$type}->{$chr}->[$i]) { $mC_mut = $mC_num_mut{$type}->{$chr}->[$i] }
			
			my $n11 = $mC_wt;
			my $n1p = $mC_wt + $mC_mut;
			my $np1 = $dep_wt;
			my $npp = $dep_wt + $dep_mut;
			my $p_val = 1;
			$p_val = cal_Fisher_exact_test($n11, $n1p, $np1, $npp);
			if($p_val eq "error"){
				#print STDERR join("\t", ($chr, $start, $end, $n11, $n1p, $np1, $npp)) , "\n\n" if($debug);
				$p_val = 1;
			}
			push @p_array, $p_val;
			
		}
		print P join ("\t", ($chr, $start, $end, @p_array)), "\n";
	}	
}
close(P);

#####################
#	R cmd
#####################
#my $cmd = "R --slave --vanilla --args $indir $input $outdir $label < ask4_draw_50kb_bin_along_chr_modify_for_yCoor.R";
#	    1   2        3        4     5    6
my $R_cmd = "R --slave --vanilla --args $temp_p_file $temp_p_adj_file < $R_script ";

print STDERR "\n", $R_cmd, "\n\n";
`$R_cmd`;

die unless (-e $temp_p_adj_file);

my ($last_chr, $last_start, $last_end ) = ("chr0", 0, 0);

my @raw_list;

open(ADJ, $temp_p_adj_file) or die "cannot open $temp_p_adj_file";
my $h = <ADJ>;
while(<ADJ>){
	chomp;
	my @a = split "\t";
	my ($chr, $start, $end) = @a[0..2];
	if(  $a[3] <= $p_value_cutoff or  $a[4] <= $p_value_cutoff or $a[5] <= $p_value_cutoff ){
		if($chr ne $last_chr || $start > $last_end + $allowed_gap + 1 ){
			if($last_chr ne "chr0") {
				push @raw_list ,    join("_", ( $last_chr, $last_start, $last_end)) ;
			}
			($last_chr, $last_start, $last_end ) =  ($chr, $start, $end);
		}else{
			$last_end = $end
		}
	} 
}

close(ADJ);
if($debug){
	print STDERR join("\n", @raw_list), "\n";
	print STDERR "\n\n", scalar (@raw_list), "\n\n";
}
my %regions_in_list; 
my @border_list;
foreach my $i ( 0..$#raw_list ){
	my ($chr, $start, $end) = split "_", $raw_list[$i];
	my $DMC_num = 0;
	for my $j ($start..$end){
		if(defined $DMCs{$chr}->[$j]){
			$DMC_num++;
		}
	}
	if($DMC_num >= $DMC_cutoff ){
		my ($real_start, $real_end) = ($start, $end) ;
		for my $j($start..$end){
			if(defined $DMCs{$chr}->[$j]){
				$real_start = $j;
				last;
			}
		}
		for (my $j = $end; $j >= $real_start; $j--){
			if(defined $DMCs{$chr}->[$j]){
				$real_end = $j;
				last;
			}
		}
		push @border_list, join("_" , ($chr, $real_start, $real_end, $DMC_num));
	}
}

foreach my $i (0..$#border_list){
	my ($chr, $start, $end, $num) = split "_", $border_list[$i];
	for my $j($start..$end){
		$regions_in_list{$chr}->[$j] = $i;
	}
}

#chr     pos     strand  type    num_C   depth   percentage      isMeth
#chr1    1       +       CHH      0       0       0			       0
# 0      1       2         3      4      5       6                 7


open (WT,  $wt_in )  or die "cannot open  $wt_in: $!";
open (MUT, $mut_in ) or die "cannot open $mut_in: $!";

$h_wt  = <WT>;
$h_mut = <MUT>;

 $line = 0;
my (%mC_list, %dep_list);
while($l_wt = <WT>, $l_mut = <MUT>){
	$line++;
	if($line % $const == 0){
		print STDERR $line, "\tdone\n";
	}
	chomp $l_wt;
	chomp $l_mut;
	my @a_wt  = split "\t", $l_wt;
	my @a_mut = split "\t", $l_mut;
	my $pos = $a_wt[1];

#	die unless ($a_wt[1] == $a_mut[1]);
	my $chr = $a_wt[0];
	my $type = $a_wt[3];
	
	next unless (defined $regions_in_list{$chr}->[$pos]);
	my $order = $regions_in_list{$chr}->[$pos];	
	my $dep_wt = $a_wt[5];
	my $dep_mut = $a_mut[5];

	if($dep_wt >= $depth_cutoff and $dep_mut >= $depth_cutoff){
		my $mC_wt = $a_wt[4];
		my $mC_mut = $a_mut[4];
		$mC_list{"WT"}->{$type}->[$order] += $mC_wt;
		$mC_list{"mut"}->{$type}->[$order] += $mC_mut;

		$dep_list{"WT"}->{$type}->[$order] += $dep_wt;
		$dep_list{"mut"}->{$type}->[$order] += $dep_mut;

	}
}
close(WT);
close(MUT);


foreach my $i (0..$#border_list){
	my ($chr, $start, $end, $num) = split "_", $border_list[$i];
	my (%divide_mut, %divide_wt);
	my (%per_mut, %per_wt);

	foreach my $type (@types){
		my ($mC_wt , $dep_wt ) = (0) x 2;
		my ($mC_mut, $dep_mut) = (0) x 2;
		
		if(defined $dep_list{"WT"} ->{$type}->[$i] ) { $dep_wt  = $dep_list{"WT"} ->{$type}->[$i] }
		if(defined $dep_list{"mut"}->{$type}->[$i] ) { $dep_mut = $dep_list{"mut"}->{$type}->[$i] }
		if(defined $mC_list {"WT"} ->{$type}->[$i] ) { $mC_wt   = $mC_list {"WT"} ->{$type}->[$i] }
		if(defined $mC_list {"mut"}->{$type}->[$i] ) { $mC_mut  = $mC_list {"mut"}->{$type}->[$i] }
		$divide_mut{$type} = "$mC_mut/$dep_mut=";
		$divide_wt {$type} = "$mC_wt/$dep_wt=";
		if($dep_mut != 0) {
			$per_mut{$type} = sprintf("%.4f",  100 * $mC_mut/$dep_mut);
		}else{
			$per_mut{$type} = "NA";
		}	
		if($dep_wt != 0) {
			$per_wt{$type} = sprintf("%.4f",  100 * $mC_wt/$dep_wt);
		}else{
			$per_wt{$type} = "NA";
		}	
	}
	print OUT join("\t", ($chr, $start, $end, $num, 
							 $divide_mut{"CG"}, $per_mut{"CG"},  $divide_mut{"CHG"}, $per_mut{"CHG"}, $divide_mut{"CHH"}, $per_mut{"CHH"},
					 		 $divide_wt{"CG"}, $per_wt{"CG"},    $divide_wt{"CHG"}, $per_wt{"CHG"},   $divide_wt{"CHH"}, $per_wt{"CHH"}  )), "\n";
}


close(OUT);
exit;


sub round {
    my($number) = shift;
    #return int($number + .5);
    return int($number + .5 * ($number <=> 0)); # take care of negative numbers too
}

#		           WT	    mut
#		           word2   ~word2
#mC	     word1    n11      n12 | n1p
#No.T	~word1    n21      n22 | n2p
#      	          --------------
#		          np1      np2   npp

#$n11 = $wt_mC;
#$n1p = $wt_mC + $mut_mC;
#$np1 = $wt_depth;
#$npp = $wt_depth + $mut_depth;

sub cal_Fisher_exact_test{
	my($n11, $n1p ,$np1, $npp) = @_;
	my  $p_value = calculateStatistic( n11=>$n11,
					   n1p=>$n1p,
					   np1=>$np1,
					   npp=>$npp);
	my $errorCode;
	if( ($errorCode = getErrorCode())){
#		print STDERR $errorCode." - ".getErrorMessage(),"\n\n";
  		return "error";
	}
  	else{
		return $p_value;
	}
}
