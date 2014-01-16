#!/usr/bin/perl -w

#main purpose :
# display wig for 35S-SUC2 promoter
# use chrM

use strict;
use File::Spec;

my $debug = 0;
print STDERR "\n\n", "change desc for each paper", "\n\n\n";
#my $desc  = "Huang et al";
my $desc  = "35S-SUC2";

#my $usage = "$0 \n <isMeth_file> <outdir> <sample_label_in_wig> <wig_pre_name> \n\n";
#die $usage unless(@ARGV == 4);

my $usage = "$0 \n <acgt_dir> <acgt_pre> <outdir> <sample_label_in_wig> <wig_pre_name> \n\n";
die $usage unless(@ARGV == 5);

my $depth_cutoff = 4;
my $indir = shift or die;
my $inpre = shift or die;

#my $input = shift or die;
my $outdir = shift or die;
my $label = shift or die;
my $pre =   shift or die;

#die unless (-e $input and -d $outdir);
die unless ( -d $outdir);
die unless ( -d $indir);

#my $forw_file = File::Spec->catfile( $indir, $inpre . "_forw.txt" ) ;
#my $rev_file  = File::Spec->catfile( $indir, $inpre . "_rev.txt" ) ;

my $forw_file = File::Spec->catfile( $indir, $inpre . "_plus_meth.txt" ) ;
my $rev_file  = File::Spec->catfile( $indir, $inpre . "_minus_meth.txt" ) ;



die unless ( -e $forw_file  and -e $rev_file);

my $mCG_wig  = File::Spec->catfile($outdir, $pre . "_NCBI_depth" . $depth_cutoff ."_mCG.wig");
my $mCHG_wig = File::Spec->catfile($outdir, $pre . "_NCBI_depth" . $depth_cutoff ."_mCHG.wig");
my $mCHH_wig = File::Spec->catfile($outdir, $pre . "_NCBI_depth" . $depth_cutoff ."_mCHH.wig");

my $all_wig  = File::Spec->catfile($outdir, $pre . "_NCBI_depth" . $depth_cutoff ."_mC.wig");

die "output exits \n" if( -e $mCG_wig or -e $mCHG_wig or -e $mCHH_wig);
die "output exits \n" if( -e $all_wig );

open(CG,  ">$mCG_wig") or die;
open(CHG, ">$mCHG_wig") or die;
open(CHH, ">$mCHH_wig") or die;
open(ALL, ">$all_wig") or die;
# track type=wiggle_0 name="WT rep2 mCG" description="WT replicate 2 mCG Stroud et al" visibility=full color=0,100,0 graphType=bar
#CHG: color=0,0,255
#CHH: color=255,0,0

print ALL  join("\t", ("track", "type=wiggle_0","name=")), "\"$label mC\"\t", "description=\"$label"," C methylation $desc\"\t",  join("\t", ( "visibility=full", "color=255,102,0",   "graphType=bar", "viewLimits=-1.0:1.0" )), "\n";
print CG  join("\t", ("track", "type=wiggle_0","name=")), "\"$label mCG\"\t", "description=\"$label"," CG methylation $desc\"\t",  join("\t", ( "visibility=full", "color=205,205,0",   "graphType=bar", "viewLimits=-1.0:1.0" )), "\n";
print CHG join("\t", ("track", "type=wiggle_0","name=")), "\"$label mCHG\"\t","description=\"$label"," CHG methylation $desc\"\t", join("\t", ( "visibility=full", "color=100,149,237", "graphType=bar", "viewLimits=-1.0:1.0" )), "\n";
print CHH join("\t", ("track", "type=wiggle_0","name=")), "\"$label mCHH\"\t","description=\"$label"," CHH methylation $desc\"\t", join("\t", ( "visibility=full", "color=205,155,155", "graphType=bar", "viewLimits=-1.0:1.0" )), "\n";

close(CG);
close(CHG);
close(CHH);
close(ALL);

if($debug){
	print STDERR "output:\n";
	print STDERR join("\n", ($all_wig, $mCG_wig, $mCHG_wig, $mCHH_wig)), "\n\n";
	
	print STDERR "OK\n\n";
	exit;
}

my ( %mCs, %mCGs, %mCHGs, %mCHHs  );


read_acgt_file (  $forw_file, 1, \%mCs, \%mCGs, \%mCHGs, \%mCHHs );
read_acgt_file (  $rev_file ,-1, \%mCs, \%mCGs, \%mCHGs, \%mCHHs );
#open(IN, $input) or die "cannot open $input";

=head
my $head = <IN>;#skip head;
while (<IN>){
	chomp;
	my @a = split "\t";
	

	my $strand = ($a[2] eq "+") ? 1 : -1; 
	my ($chr, $pos, $val, $type)  = ($a[0], $a[1], $a[6] * $strand, $a[3]);
	my ( $mC_num, $dep  )  = @a[4..5];
		
	die $_ if ( $mC_num > $dep ) ;
		
	next if (  $dep < $depth_cutoff);
		
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
=cut

output_wig($mCG_wig, \%mCGs);
output_wig($mCHG_wig, \%mCHGs);
output_wig($mCHH_wig, \%mCHHs);
output_wig($all_wig, \%mCs);


exit;
#seq_35S_SUC2_UTR_MG     6       6       CG:2    0       -
#seq_35S_SUC2_UTR_MG     20      20      CHH:3   0       -
#seq_35S_SUC2_UTR_MG     21      21      CHH:3   0       -
#seq_35S_SUC2_UTR_MG     29      29      CG:3    0.333333        -

#read_acgt_file (  $forw_file, 1, \%mCs, \%mCGs, \%mCHGs, \%mCHHs );
sub read_acgt_file{
	my ( $file, $stand_num, $mCs_ref, $mCGs_ref, $mCHGs_ref, $mCHHs_ref ) = @_;
	
	die unless (-e $file);
	
	open(IN, $file) or die;
	while (<IN>) {
		chomp;
		my @a = split "\t";
		my $chr = "chrM";
		if( $a[0] =~ /chr./){
			$chr = $a[0];
		}
		my $pos = $a[1] + 1;
		my ($type, $dep) = split ":", $a[3];
		my $val = $a[4] * $stand_num;
		
		next if (  $dep < $depth_cutoff);
		
		$mCs_ref->{$chr}->{$pos} = $val;
		
		if($type eq "CG"){
			$mCGs_ref ->{$chr}->{$pos} = $val;
		}elsif($type eq "CHG"){
			$mCHGs_ref->{$chr}->{$pos} = $val;
		}elsif($type eq "CHH"){
			$mCHHs_ref->{$chr}->{$pos} = $val;
		}else{
			die "$_\n";
		}
	}
	close IN;
}
	
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