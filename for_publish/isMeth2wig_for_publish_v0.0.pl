#!/usr/bin/perl -w

#v0.1
# input hava a col of strand

#v0.2 
#input parameter is diff from v0.1

use strict;
use File::Spec;

my $debug = 0;


my $const = 3000000;
if($debug){
#	$const = 500;
}

#	my $cmd = "perl $script  $isMeth_file $pre $outdir $postfix";

# my $usage = "$0 <indir> <outdir> <pre> <file_name>";
 my $usage = "$0 \n <indir>  <pre>   <outdir> <file_name>\n\n";
die $usage unless(@ARGV == 4);

my $input = shift or die;
my $pre = shift or die "pre";
my $outdir = shift or die "outdir";
my $postfix = shift or die;

die unless (-e $input);
die unless (-d $outdir);

my $label = $pre;
#my $label_removed = $pre . "_for_" . $postfix;

#my $input = File::Spec->catfile($indir,$pre. "_isMeth.txt");
#my $input = File::Spec->catfile($indir,$infile);
die unless (-e $input);

print STDERR "input:\n", $input, "\n\n";

my $mC_wig   = File::Spec->catfile($outdir, $pre . "_for_" . $postfix. "_mC.wig");
my $mCG_wig  = File::Spec->catfile($outdir, $pre . "_for_" . $postfix. "_mCG.wig");
my $mCHG_wig = File::Spec->catfile($outdir, $pre . "_for_" . $postfix. "_mCHG.wig");
my $mCHH_wig = File::Spec->catfile($outdir, $pre . "_for_" . $postfix. "_mCHH.wig");

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
		
		
#print OUT join("\t", ("track", "type=wiggle_0","name=")), "\"$label_removed\"\t","description=\"$label_removed"," methlation\"\t", join("\t", ("viewLimits=-1.0:1.0", "color=255,102,0")), "\n";
print OUT join("\t", ("track", "type=wiggle_0","name=")), "\"$label\"\t","description=\"$label"," methlation\"\t", join("\t", ("viewLimits=-1.0:1.0", "color=255,102,0")), "\n";

print CG  join("\t", ("track", "type=wiggle_0","name=")), "\"$label\"\t","description=\"$label"," CG methlation\"\t", join("\t", ("viewLimits=-1.0:1.0", "color=205,205,0")), "\n";
print CHG join("\t", ("track", "type=wiggle_0","name=")), "\"$label\"\t","description=\"$label"," CHG methlation\"\t", join("\t", ("viewLimits=-1.0:1.0", "color=100,149,237")), "\n";
print CHH join("\t", ("track", "type=wiggle_0","name=")), "\"$label\"\t","description=\"$label"," CHH methlation\"\t", join("\t", ("viewLimits=-1.0:1.0", "color=205,155,155")), "\n";

close(OUT);
close(CG);
close(CHG);
close(CHH);
		
#my $last_chr = "chr0";

# 0      1        2      3       4        5      6                 7
#chr     pos     strand  type    num_C   depth   percentage      isMeth
#chr1    1       +       CHH     0       0       0       0
#chr1    2       +       CHH     0       0       0       0
# 0      1        2      3       4        5      6       7


my (%mC, %mCGs, %mCHGs, %mCHHs);
my $i = 0;
my $head = <IN>;#skip head;
while (<IN>){
	$i++;
	chomp;
	my @a = split "\t";
	
	if($i % $const == 0){
		print STDERR $i, "...done\n";
	}
	
	if ($a[7] == 1){
	#	my ($chr,$pos, $strand) = @{$records[$a[1]]};
		my $strand = ($a[2] eq "+") ? 1 : -1; 
		my ($chr, $pos, $val, $type)  = ($a[0], $a[1], $a[6] * $strand, $a[3]);
		$mC{$chr}->{$pos} = $val;
		
		if($type eq "CG"){
			$mCGs{$chr}->{$pos} = $val;
		}elsif($type eq "CHG"){
			$mCHGs{$chr}->{$pos} = $val;
		}elsif($type eq "CHH"){
			$mCHHs{$chr}->{$pos} = $val;
		}else{
			die "$_\n";
		}
	}
}

output_wig($mC_wig, \%mC);
output_wig($mCG_wig, \%mCGs);
output_wig($mCHG_wig, \%mCHGs);
output_wig($mCHH_wig, \%mCHHs);

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