#!/usr/bin/perl -w


use strict;
use File::Spec;

my $debug = 0;
#my $const = 3000000;
if($debug){
#	$const = 500;
	print STDERR "\n\n debug = 1 \n\n";
}

# my $usage = "$0 \n <indir> <input_isMeth_name> <outdir> <out_pre> <cutoff>\n\n";
# die $usage unless(@ARGV == 5);

my $usage = "$0 \n <input> <outdir> <outpre> \n\n";
die $usage unless(@ARGV == 3);

#my $indir = shift or die "indir";
my $input = shift or die "input";
my $outdir = shift or die "outdir";
my $pre = shift or die "shift";

#my $cutoff = shift or die "cutoff";


#my $input = File::Spec->catfile($indir,$pre. "_isMeth.txt");
#my $input = File::Spec->catfile($indir,$infile);
die unless (-e $input);

#print STDERR "input:\n", $input, "\n\n";

#my $mC_wig   = File::Spec->catfile($outdir, $pre . "_mC.wig");
#my $mCG_wig  = File::Spec->catfile($outdir, $pre . "_mCG.wig");
#my $mCHG_wig = File::Spec->catfile($outdir, $pre . "_mCHG.wig");
#my $mCHH_wig = File::Spec->catfile($outdir, $pre . "_mCHH.wig");

my $CG_hyper  = File::Spec->catfile($outdir, $pre . "_hyper_CG.txt");
my $CHG_hyper = File::Spec->catfile($outdir, $pre . "_hyper_CHG.txt");
my $CHH_hyper = File::Spec->catfile($outdir, $pre . "_hyper_CHH.txt");

my $CG_hypo   = File::Spec->catfile($outdir, $pre . "_hypo_CG.txt");
my $CHG_hypo  = File::Spec->catfile($outdir, $pre . "_hypo_CHG.txt");
my $CHH_hypo  = File::Spec->catfile($outdir, $pre . "_hypo_CHH.txt");

# die "output exits \n" if(-e $mC_wig or -e $mCG_wig or -e $mCHG_wig or -e $mCHH_wig);
die "output exits \n" if(-e $CG_hyper);
die "output exits \n" if(-e $CHG_hyper);
die "output exits \n" if(-e $CHH_hyper);
die "output exits \n" if(-e $CG_hypo);
die "output exits \n" if(-e $CHG_hypo);
die "output exits \n" if(-e $CHH_hypo);


if($debug){
	print STDERR "output:\n";
#	print STDERR join("\n", ($mC_wig, $mCG_wig, $mCHG_wig, $mCHH_wig)), "\n\n";
	print STDERR join("\n", (  $CG_hyper, $CHG_hyper, $CHH_hyper, $CG_hypo, $CHG_hypo, $CHH_hypo )), "\n\n";
	
	print STDERR "OK\n\n";
	
	exit;
}

my (%hyper, %hypo);

open(IN, $input) or die "cannot open $input";

# 0		  1		 2		  3		 4		 5			6			 7		  8			9					10
#chr1    310     +       CG      22      39      0.564103        3       18      0.166667        0.00876778042255472
#chr1    5449    +       CG      10      17      0.588235        1       20      0.05    0.000665849386650277

#my (%mC, %mCGs, %mCHGs, %mCHHs);


while (<IN>){
#	$i++;
	chomp;
	my @a = split "\t";
	my ($wt_per, $mut_per) = ($a[6], $a[9]);
	die $_, $wt_per, $mut_per if( $wt_per > 1);
	die $_, $wt_per, $mut_per if( $mut_per > 1);
	
	my $type = $a[3];
	my ($chr, $pos) = @a[0..1];
	
	if(  $mut_per > $wt_per){
		$hyper{$type}->{$chr}->{$pos} = 1;
	}elsif( $mut_per < $wt_per ){
		$hypo{$type}->{$chr}->{$pos} = 1;
	}else{
		die $_, "\n";
	}
}
output($CG_hyper, \%{$hyper{CG}});
output($CHG_hyper, \%{$hyper{CHG}});
output($CHH_hyper, \%{$hyper{CHH}});
output($CG_hypo, \%{$hypo{CG}}	);
output($CHG_hypo, \%{$hypo{CHG}});
output($CHH_hypo, \%{$hypo{CHH}});



exit;

sub output{
	my ($file, $ref) = @_;
	die "$file exists!!!" if(-e $file) ;
	open(OUT, ">$file") or die;
	foreach my $chr (sort keys %{$ref}){
		foreach my $pos(sort {$a<=>$b} keys %{$ref->{$chr}}){
			print OUT join("\t", ( $chr, $pos )), "\n";
		}
	}
	close(OUT);
}