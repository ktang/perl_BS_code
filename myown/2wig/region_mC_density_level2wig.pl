#!/usr/bin/perl -w
#input file is like 
# /Users/tang58/misc/Zhaobo_Lang/6_mbd7_91-2/20140312_email/Fig5C_6F_wig_for_mCG_density_level/TAIR10_WinSize1000bp_sliding500bp_in_colA_dep4_meth_level.txt
# indicate CG/CHG/CHH/C to output 2 wigs, mCX_density and mCX_meth_level

# 1 ID
# 2 chr
# 3 start
# 4 end
# 5 strand
# 6 num_CG
# 7 num_CHG
# 8 num_CHH
# 9 num_C

# 10 colA_covered_num_CG
# 11 colA_covered_num_CHG
# 12 colA_covered_num_CHH
# 13 colA_covered_num_C

# 14 colA_covered_per_CG
# 15 colA_covered_per_CHG
# 16 colA_covered_per_CHH
# 17 colA_covered_per_C

# 18 colA_nE_wmCG
# 19 colA_nE_per_wmCG
# 20 colA_nE_wmCHG
# 21 colA_nE_per_wmCHG
# 22 colA_nE_wmCHH
# 23 colA_nE_per_wmCHH
# 24 colA_nE_wmC
# 25 colA_nE_per_wmC

#mm useful
# 26 colA_nE_mmCG
# 27 colA_nE_per_mmCG
# 28 colA_nE_mmCHG
# 29 colA_nE_per_mmCHG
# 30 colA_nE_mmCHH
# 31 colA_nE_per_mmCHH
# 32 colA_nE_mmC
# 33 colA_nE_per_mmC

# 34 colA_nE_fmCG
# 35 colA_nE_per_fmCG
# 36 colA_nE_fmCHG
# 37 colA_nE_per_fmCHG
# 38 colA_nE_fmCHH
# 39 colA_nE_per_fmCHH
# 40 colA_nE_fmC
# 41 colA_nE_per_fmC
BEGIN { push @INC, '/Users/tang58/scripts_all/perl_code/Modules' }
use Kai_Module;

use strict;
use File::Spec;
my $debug = 1;
#my $const = 3000000;
if($debug){
	#$const = 500;
}

print STDERR "display highest point ";
my $highest = 15;
my $usage = "$0 \n<input> <CXX> <sample_label> <cover_cutoff (0-100)> <outdir> <pre (TAIR10_WinSize1000bp_sliding500bp)> <chr_col_number> <normalized_factor(1000)>\n\n";

die $usage unless(@ARGV == 8);

#my $indir = shift or die "indir";
my $input = shift or die "input";
my $CXX_label = shift or die "CXX_label";
my $sample_label = shift or die "sample_label";
my $cover_cutoff = shift or die "cover_cutoff";

my $outdir = shift or die "outdir";
my $pre = shift or die "shift";

my $chr_col_number = shift or die;

my $factor = shift or die;

die unless (-e $input);

print STDERR "input:\n", $input, "\n\n";

#my $mC_wig   = File::Spec->catfile($outdir, $pre . "_mC.wig");
#my $mCG_wig  = File::Spec->catfile($outdir, $pre . "_mCG.wig");
#my $mCHG_wig = File::Spec->catfile($outdir, $pre . "_mCHG.wig");
#my $mCHH_wig = File::Spec->catfile($outdir, $pre . "_mCHH.wig");

my $wig_den   = File::Spec->catfile($outdir, $pre . "_" . $sample_label . "_m" . $CXX_label . "_density" .  "_cover" . $cover_cutoff  . ".wig");
my $wig_level = File::Spec->catfile($outdir, $pre . "_" . $sample_label . "_m" . $CXX_label . "_wmlevel" .  "_cover" . $cover_cutoff .  ".wig");

#die "output exits \n" if(-e $mC_wig or -e $mCG_wig or -e $mCHG_wig or -e $mCHH_wig);
die "output exits \n" if(-e $wig_den);
die "output exits \n" if(-e $wig_level);
print STDERR "output:\n";
#print STDERR join("\n", ($mC_wig, $mCG_wig, $mCHG_wig, $mCHH_wig)), "\n\n";
print STDERR join("\n", ($wig_den, $wig_level )), "\n\n";

open(IN, $input) or die "cannot open $input";
my $h = <IN>;
my @h_a = split "\t", $h;
my $h_label = $sample_label . "_covered_num_" . $CXX_label;
my $first_index = -1;
for my $i(0..$#h_a){
	if ($h_a[$i] eq $h_label){
		$first_index = $i;
		last;
	}
}
if ($first_index == -1){
	die "head not find, $h_label ";
}
my ( $cover_num_i, $cover_per_i, $wm_per_i, $mm_per_i ) = ( $first_index ,  $first_index + 4 ,   $first_index + 9,   $first_index + 17);
die unless ( $h_a[$cover_per_i ]=~ /covered_per_$CXX_label/);
die unless ( $h_a[$wm_per_i ]=~ /per_wm$CXX_label/);
die unless ( $h_a[$mm_per_i ]=~ /per_mm$CXX_label/);

#close IN;
open(DEN, ">$wig_den") or die;
open(WM,  ">$wig_level") or die;
#open(CHG, ">$mCHG_wig") or die;
#open(CHH, ">$mCHH_wig") or die;
		
		
#print OUT join("\t", ("track", "type=wiggle_0","name=")), "\"$label_removed\"\t","description=\"$label_removed"," methlation\"\t", join("\t", ("viewLimits=-1.0:1.0", "color=255,102,0")), "\n";
#print CG  join("\t", ("track", "type=wiggle_0","name=")), "\"$label\"\t","description=\"$label"," CG methlation\"\t", join("\t", ("viewLimits=-1.0:1.0", "color=205,205,0")), "\n";
#print CHG join("\t", ("track", "type=wiggle_0","name=")), "\"$label\"\t","description=\"$label"," CHG methlation\"\t", join("\t", ("viewLimits=-1.0:1.0", "color=100,149,237")), "\n";
#print CHH join("\t", ("track", "type=wiggle_0","name=")), "\"$label\"\t","description=\"$label"," CHH methlation\"\t", join("\t", ("viewLimits=-1.0:1.0", "color=205,155,155")), "\n";
print DEN join("\t", ("track", "type=wiggle_0","name=")), "\"$sample_label\"\t","description=\"$sample_label m$CXX_label"," density\"\t", join("\t", ("viewLimits=0:$highest", "color=255,102,0")), "\n";
print WM join("\t", ("track", "type=wiggle_0","name=")), "\"$sample_label\"\t","description=\"$sample_label m$CXX_label" ," methlation\"\t", join("\t", ("viewLimits=0:100", "color=255,102,0")), "\n";
#less TAIR10_WinSize1000bp_sliding500bp_in_colA_dep4_meth_level.txt | cut -f 1-6,10,14,18,19,26,27 | cut -f11 | less| grep '/' | perl -lane '@a=split "/"; print $a[0]' | less
# less TAIR10_WinSize1000bp_sliding500bp_in_colA_dep4_meth_level.txt | cut -f 1-6,10,14,18,19,26,27 | cut -f11 | less| grep '/' | perl -lane '@a=split "/"; print $a[0]' | sort -g | awk ' BEGIN{OFS = "\t" } { a[i++]=$1;  s+=$1 }   END { mean = s/i; print i, mean, a[int(i*0.99)], a[int(i*0.95)], a[int(i*0.9)], a[int(i*0.75)], a[int(i/2)], a[int(i*0.25)], a[int(i*0.1)] , a[int(i*0.05)], a[int(i*0.01)];}'
# 234982	10.6413	72.573108	43.468225	31.669772	16.040954	2.3270201	0.308171	0.0833333	0	0

# close(OUT);
# close(CG);
# close(CHG);
# close(CHH);
close(DEN);
close(WM);
#my (%mC, %mCGs, %mCHGs, %mCHHs);

my %den_h;
my %wm_h;

if($debug){
	print "debug STOP\n\n";
	exit;
}

#$cover_num_i, $cover_per_i, $wm_per_i, $mm_per_i
while (<IN>){
#	$i++;
	chomp;
	my @a = split "\t";
	my $chr = $a[$chr_col_number - 1];
	$chr = Kai_Module::simple_chr($chr);
	my $pos = $a[$chr_col_number];
	
	if($a[$cover_per_i] eq "NA"){
		$den_h  {$chr}->{$pos} = 0;
		$wm_h{$chr}->{$pos} = 0;
		next;
	}
	
	if ($a[$cover_per_i] < $cover_cutoff){
		next;
	}
	$wm_h{$chr}->{$pos} = $a[$wm_per_i] ;
	my $den_val = eval sprintf("%.2f",  $a[$mm_per_i]*$a[$cover_num_i] / $factor); 
	$den_h {$chr}->{$pos} = $den_val;
	
}

#output_wig($mC_wig, \%mC);
#output_wig($mCG_wig, \%mCGs);
#output_wig($mCHG_wig, \%mCHGs);
#output_wig($mCHH_wig, \%mCHHs);
output_wig($wig_den , \%den_h);
output_wig($wig_level , \%wm_h);
exit;

sub output_wig{
	my ($file, $ref) = @_;
	open(FILE, ">>$file") or die;
	
	foreach my $chr (sort keys %{$ref}){
		print FILE "variableStep\tchrom=$chr\n";
		foreach my $pos(sort {$a<=>$b} keys %{$ref->{$chr}}){
			print FILE join("\t", ($pos, $ref->{$chr}->{$pos})), "\n";
		}
	}
	
	close(FILE);
}