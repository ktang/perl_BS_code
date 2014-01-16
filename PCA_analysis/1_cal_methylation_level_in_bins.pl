#!/usr/bin/perl -w

#/Users/tang58/Kai_BS/PCA_analysis/1_cal_methylation_level_in_bins.pl
#outformat
#coordinate	XX_mCXX_level
use strict;
use File::Spec;

my $debug = 0;
if($debug){
	print STDERR "debug: $debug\n";
}
my $bin_size = 5000;#300;
my $dep_cutoff = 4;

#my $debug_chr1_plus = 1;

#my %chr_len = ("chr1"=>30427671, "chr2"=>19698289, "chr3"=>23459830, "chr4"=>18585056, "chr5"=>26975502);

my %chr_len_simple = ("1"=>30427671, "2"=>19698289, "3"=>23459830, "4"=>18585056, "5"=>26975502);

my %bin_last_index;
for my $chr (sort keys %chr_len_simple){
	my $len = $chr_len_simple{$chr};
	my $last = int( ($len - 1) / $bin_size);
	$bin_last_index{$chr} = $last;
}

#my @chrs = ("chr1", "chr2", "chr3", "chr4", "chr5");
#my $usage = "$0 <input> <prefix_output> <outdir>\n/Volumes/Macintosh_HD_2/idm1_new_met_data/UCR_server_data/Lane1_JKZ1_Col0/s_1_isMeth.txt";
#label can be ros1,ros2,ros3,ros4,rdd,hsp,phd,zdp
#die $usage unless (@ARGV == 3);

my $usage = "$0 \n <input_isMeth> <sample_label> <outdir> \n\n";
die $usage unless (@ARGV == 3);

my $input_isMeth = shift or die;
my $sample_label = shift or die;
my $outdir	 = shift or die;
#my $outpre	 = shift or die;


#my $input = $ARGV[0];
#my $out_prefix = $ARGV[1];
#my $outdir = $ARGV[2];
#die "wrong input" unless (-e $input);

die "wrong input" unless (-e $input_isMeth);

die "wrong outdir" unless (-d $outdir);
#my $output = $out_prefix."_BinSize$bin_size". ".txt";
#if($debug){print STDERR "output: $output"}
#if (-e "$outdir/$output") {die "output $output exist:$!"}

# XXX_mCX_level_for_PCA_Bin$bin_size.txt

my $output_CG  = File::Spec->catfile($outdir, $sample_label . "_mCG_level_for_PCA_Bin"  . $bin_size . "_dep$dep_cutoff" . ".txt");
my $output_CHG = File::Spec->catfile($outdir, $sample_label . "_mCHG_level_for_PCA_Bin" . $bin_size . "_dep$dep_cutoff" . ".txt");
my $output_CHH = File::Spec->catfile($outdir, $sample_label . "_mCHH_level_for_PCA_Bin" . $bin_size . "_dep$dep_cutoff" . ".txt");

die if (-e $output_CG);
die if (-e $output_CHG);
die if (-e $output_CHH);


my %mC_h;
my %dep_h;

if ($debug){
#	$input = "/Users/tang58/try/debug/Weiqiang_idea/s1_isMeth_1000.txt";
#	$output = "detail_debug_output.txt"; 
#	$outdir = "/Users/tang58/try/debug/Weiqiang_idea/";
#	$bin_size = 50;
#	%chr_len = ("chr1"=>5576);
#	%bin_last_index = ('chr1' => 111);
	print STDERR join("\n", ( $output_CG, $output_CHG, $output_CHH )), "\n\n";
	print STDERR "OK\n\n";
	exit;
}

#open (OUT, ">$outdir/$output") or die "cannot open $output: $!";
open (IN, $input_isMeth ) or die "cannot open $input_isMeth: $!";

#chr     pos     strand  type    num_C   depth   percentage      isMeth
#chr1    1       +       CHH     0       0       0       0
#chr1    2       +       CHH     0       0       0       0
my $head = <IN>;
while(<IN>){
	my @a = split "\t";
	
#	if( $debug_chr1_plus ){
#		last if ($a[0] eq "chr2");
#	}
	
	next unless ( $a[5] >= $dep_cutoff );
	
	#my ($chr, $pos, $strand, $type, $mC_num, $dep_num, $per, $isMeth) = @_;
	my $chr = $a[0];
	$chr = simple_chr($chr);
	
	my $pos = $a[1];
	my $index = int ( ($pos -1 ) / $bin_size);
	
	my $type = $a[3];
	my $mC_num = $a[4];
	my $dep_num = $a[5];
	$mC_h{$type}->{$chr}->[$index] += $mC_num;
	$dep_h{$type}->{$chr}->[$index]+= $dep_num;
}
close IN;

output_file($output_CG,  \%{$mC_h{CG}},  \%{$dep_h{CG}},  $sample_label, "CG");
output_file($output_CHG, \%{$mC_h{CHG}}, \%{$dep_h{CHG}}, $sample_label, "CHG");
output_file($output_CHH, \%{$mC_h{CHH}}, \%{$dep_h{CHH}}, $sample_label, "CHH");


exit;

sub output_file{
	my ($output, $mC_ref, $dep_ref, $label, $type) = @_;
	die if(-e $output);
	open(OUT, ">>$output") or die;
	print OUT join("\t", ("coordinate", $label . "_m" . $type. "_level")), "\n";
	for my $chr (sort keys %{$mC_ref}){
	#	my $last_index = scalar(@{$mC_ref->{$chr}}) - 1;
		my $last_index  = $bin_last_index{$chr};
		for my $i (0..($last_index - 1)){
			my $s = $i * $bin_size + 1;
			my $e = ($i + 1 ) * $bin_size;
			my $cor = "$chr:$s-$e";
			my $per = "NA";
			my $mC_sub = 0;
			my $dep_sub = 0;
			if (defined $mC_ref->{$chr}->[$i]){
				$mC_sub = $mC_ref->{$chr}->[$i]
			}
			if(defined $dep_ref->{$chr}->[$i]){
				$dep_sub = $dep_ref->{$chr}->[$i]
			}
			if($dep_sub != 0){
				$per = sprintf ("%.2f", 100 * $mC_sub / $dep_sub)
			}
			
			print OUT join("\t", ($cor,$per)), "\n";
		}
		
		my $i = $last_index;
		my $s = $i * $bin_size + 1;
		my $e = $chr_len_simple{$chr};
		my $cor = "$chr:$s-$e";
		my $per = "NA";
		my $mC_sub = 0;
		my $dep_sub = 0;
		if (defined $mC_ref->{$chr}->[$i]){ $mC_sub = $mC_ref->{$chr}->[$i] }
		if(defined $dep_ref->{$chr}->[$i]){ $dep_sub = $dep_ref->{$chr}->[$i] }
		if($dep_sub != 0){
			$per = sprintf ("%.2f", 100 * $mC_sub / $dep_sub)
		}
		print OUT join("\t", ($cor,$per)), "\n";
	}
	
	close OUT;
}

sub simple_chr{
	my ($chr) = @_;
	if( $chr =~ /chr/i){
		$chr =~  s/chr//i;
	}
	if($chr eq "M" ){
		$chr = "Mt";
	}elsif( $chr eq "C"){
		$chr = "Pt";
	}
	return $chr;
}