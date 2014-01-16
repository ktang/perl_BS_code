#!/usr/bin/perl -w

#v0.5.1 Dec 5, 2013;
#
# in this version of script, the hash record is different from
# previous,
# as in previous version, there are four hashes to record the info
# now use only one hash to record
# val_type;

# val = eval sprintf("%.3f", $val)

#perl -e '$a= 0; print eval sprintf("%.3f", $a), "\n\n"'
#0
#perl -e '$a= 0.222522; print eval sprintf("%.3f", $a), "\n\n"'
#0.223


#v0.5
#input is original isMeth, so need depth cutoff
# and mC cutoff
# for edm2 paper dep>=4 and mC >= 1

#v0.4
# mC is not determined by last column.
# just output mC_number >= cutoff

# all C use the same artifical cutoff usually 2



#v0.1
# input hava a col of strand

#v0.2 
#input parameter is diff from v0.1

#v0.3
# mC is not determined by last column.
# just output mC_number > 0

use strict;
use File::Spec;
my $debug = 0;





#my $const = 3000000;
if($debug){
#	$const = 500;
	print STDERR "\n\n debug = 1 \n\n";
}

my $usage = "$0 \n <input_isMeth_name> <outdir> <out_pre> <mC_cutoff> <depth_cutoff>\n\n";

die $usage unless(@ARGV == 5);

#my $indir = shift or die "indir";
#my $infile = shift or die;
my $input = shift or die;
my $outdir = shift or die "outdir";
my $pre = shift or die "shift";

my $mC_cutoff = shift or die "mC_cutoff";
my $depth_cutoff = shift or die "depth_cutoff";

my $label = $pre;
my $label_removed = $pre;

#my $input = File::Spec->catfile($indir,$pre. "_isMeth.txt");
#my $input = File::Spec->catfile($indir,$infile);
die unless (-e $input);

#print STDERR "input:\n", $input, "\n\n";

my $post = "_mC" . $mC_cutoff . "_dep" . $depth_cutoff;

my $mC_wig   = File::Spec->catfile($outdir, $pre . $post.  "_mC.wig");
my $mCG_wig  = File::Spec->catfile($outdir, $pre . $post. "_mCG.wig");
my $mCHG_wig = File::Spec->catfile($outdir, $pre . $post. "_mCHG.wig");
my $mCHH_wig = File::Spec->catfile($outdir, $pre . $post. "_mCHH.wig");

die "output exits \n" if(-e $mC_wig or -e $mCG_wig or -e $mCHG_wig or -e $mCHH_wig);

if($debug){
	print STDERR "output:\n";
	print STDERR join("\n", ($mC_wig, $mCG_wig, $mCHG_wig, $mCHH_wig)), "\n\n";
	
	print STDERR "OK\n\n";
	
	exit;
}

open(IN, $input) or die "cannot open $input";
open(OUT, ">$mC_wig") or die;
open(CG,  ">$mCG_wig") or die;
open(CHG, ">$mCHG_wig") or die;
open(CHH, ">$mCHH_wig") or die;
		
		
print OUT join("\t", ("track", "type=wiggle_0","name=")), "\"$label_removed\"\t","description=\"$label_removed"," methlation\"\t", join("\t", ("viewLimits=-1.0:1.0", "color=255,102,0")), "\n";

print CG  join("\t", ("track", "type=wiggle_0","name=")), "\"$label\"\t","description=\"$label"," CG methlation\"\t", join("\t", ("viewLimits=-1.0:1.0", "color=205,205,0")), "\n";
print CHG join("\t", ("track", "type=wiggle_0","name=")), "\"$label\"\t","description=\"$label"," CHG methlation\"\t", join("\t", ("viewLimits=-1.0:1.0", "color=100,149,237")), "\n";
print CHH join("\t", ("track", "type=wiggle_0","name=")), "\"$label\"\t","description=\"$label"," CHH methlation\"\t", join("\t", ("viewLimits=-1.0:1.0", "color=205,155,155")), "\n";

#close(OUT);
#close(CG);
#close(CHG);
#close(CHH);
		
#my $last_chr = "chr0";

# 0      1        2      3       4        5      6                 7
#chr     pos     strand  type    num_C   depth   percentage      isMeth
#chr1    1       +       CHH     0       0       0       0
#chr1    2       +       CHH     0       0       0       0
# 0      1        2      3       4        5      6       7


#my (%mC, %mCGs, %mCHGs, %mCHHs);
my %m_h;
my $i = 0;
my $head = <IN>;#skip head;

my ($C_num, $CG_num, $CHG_num, $CHH_num ) = (0) x 4; 

while (<IN>){
	$i++;
	chomp;
	my @a = split "\t";
	
	#if($i % $const == 0){
	#	print STDERR $i, "...done\n";
	#}
	my ($mC, $dep) = @a[4..5];
	
	if($mC >= $mC_cutoff and  $dep >= $depth_cutoff ){
#	if($a[4] >= $cutoff){
	#if($a[4] > 0){
	#if ($a[7] == 1){
	#	my ($chr,$pos, $strand) = @{$records[$a[1]]};
		my $strand = ($a[2] eq "+") ? 1 : -1;
		
		$a[6] = eval sprintf("%.3f", $a[6]);
		
		my ($chr, $pos, $val, $type)  = ($a[0], $a[1], $a[6] * $strand, $a[3]);
		my $tmp = join("\t", ($type, $val));
		$m_h {$chr}->{$pos} = $tmp;

		#$mC{$chr}->{$pos} = $val;
		#$C_num ++;
		#if($type eq "CG"){
		#	$mCGs{$chr}->{$pos} = $val;
		#	$CG_num++;
		#}elsif($type eq "CHG"){
		#	$mCHGs{$chr}->{$pos} = $val;
		#	$CHG_num++;
		#}elsif($type eq "CHH"){
		#	$mCHHs{$chr}->{$pos} = $val;
		#	$CHH_num++;
		#}else{
		#	die "$_\n";
		#}
	#}
	}

}
	

#output_wig($mC_wig, \%mC);
#output_wig($mCG_wig, \%mCGs);
#output_wig($mCHG_wig, \%mCHGs);
#output_wig($mCHH_wig, \%mCHHs);

foreach my $chr (sort keys %m_h ){
#	print FILE "variableStep\tchrom=$chr\n";
	print OUT "variableStep\tchrom=$chr\n";
	print CG "variableStep\tchrom=$chr\n";
	print CHG "variableStep\tchrom=$chr\n";
	print CHH "variableStep\tchrom=$chr\n";
	
	foreach my $pos(sort {$a<=>$b} keys %{ $m_h{$chr} } ){
		#print FILE join("\t", ($pos, $ref->{$chr}->{$pos})), "\n";
		my $tmp = $m_h{$chr}->{$pos};
		my ($type, $val ) = split "\t",  $tmp;
		print OUT join("\t", ($pos, $val) ), "\n";
		$C_num ++;
		if($type eq "CG"){
			print CG join("\t", ($pos, $val) ), "\n";
			$CG_num++;
		}elsif($type eq "CHG"){
			print CHG join("\t", ($pos, $val) ), "\n";
			$CHG_num++;
		}elsif($type eq "CHH"){
			print CHH join("\t", ($pos, $val) ), "\n";
			$CHH_num++;
		}else{
			die "$chr, $pos, $tmp\n\n";
		}
	}
	
}
#print C join("\t", ($type, $val) ), "\n";

close(OUT);
close(CG);
close(CHG);
close(CHH);

print STDERR $pre,"\n";
print STDERR join("\t", ("C",   $C_num)), "\n";
print STDERR join("\t", ("CG",  $CG_num)), "\n";
print STDERR join("\t", ("CHG", $CHG_num)), "\n";
print STDERR join("\t", ("CHH", $CHH_num)), "\n\n\n";

exit;
