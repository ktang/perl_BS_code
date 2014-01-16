#!/usr/bin/perl -w

use strict;
use File::Spec;

my $debug = 0;
print STDERR "\n\n", "change desc for each paper", "\n\n\n";
my $desc  ="";# "Huang et al";

my $usage = "$0 \n <isMeth_file> <outdir> <sample_label_in_wig> <wig_pre_name> \n\n";
die $usage unless(@ARGV == 4);

my $depth_cutoff = 4;

my $input = shift or die;
my $outdir = shift or die;
my $label = shift or die;
my $pre =   shift or die;

die unless (-e $input and -d $outdir);

my $mCG_wig  = File::Spec->catfile($outdir, $pre . "_NCBI_depth" . $depth_cutoff ."_mCG.wig");
my $mCHG_wig = File::Spec->catfile($outdir, $pre . "_NCBI_depth" . $depth_cutoff ."_mCHG.wig");
my $mCHH_wig = File::Spec->catfile($outdir, $pre . "_NCBI_depth" . $depth_cutoff ."_mCHH.wig");

die "output exits \n" if( -e $mCG_wig or -e $mCHG_wig or -e $mCHH_wig);


open(CG,  ">$mCG_wig") or die;
open(CHG, ">$mCHG_wig") or die;
open(CHH, ">$mCHH_wig") or die;

# track type=wiggle_0 name="WT rep2 mCG" description="WT replicate 2 mCG Stroud et al" visibility=full color=0,100,0 graphType=bar
#CHG: color=0,0,255
#CHH: color=255,0,0

print CG  join("\t", ("track", "type=wiggle_0","name=")), "\"$label mCG\"\t", "description=\"$label"," CG methylation $desc\"\t",  join("\t", ( "visibility=full", "color=205,205,0",   "graphType=bar", "viewLimits=-1.0:1.0" )), "\n";
print CHG join("\t", ("track", "type=wiggle_0","name=")), "\"$label mCHG\"\t","description=\"$label"," CHG methylation $desc\"\t", join("\t", ( "visibility=full", "color=100,149,237", "graphType=bar", "viewLimits=-1.0:1.0" )), "\n";
print CHH join("\t", ("track", "type=wiggle_0","name=")), "\"$label mCHH\"\t","description=\"$label"," CHH methylation $desc\"\t", join("\t", ( "visibility=full", "color=205,155,155", "graphType=bar", "viewLimits=-1.0:1.0" )), "\n";

close(CG);
close(CHG);
close(CHH);

if($debug){
	print STDERR "output:\n";
	print STDERR join("\n", ( $mCG_wig, $mCHG_wig, $mCHH_wig)), "\n\n";
	
	print STDERR "OK\n\n";
	exit;
}

my ( %mCGs, %mCHGs, %mCHHs  ); 
open(IN, $input) or die "cannot open $input";


#new purude
#chr     pos     strand  type    num_C   depth   percentage      isMeth
#0	1	 2	 3	 4	 5	 6		 7


#UCR old
#0	1	2	 3	  4	 5	 6	  7	
#chr     pos     depth   num_C   percent strand  type    isMeth
#2 chr1    1       0       0       0       +       CHH     0
#3 chr1    2       0       0       0       +       CHH     0


my $head = <IN>;#skip head;
while (<IN>){
	chomp;
	my @a = split "\t";
	

#	my $strand = ($a[2] eq "+") ? 1 : -1; 
	my $strand = ($a[5] eq "+") ? 1 : -1; 
	#my ($chr, $pos, $val, $type)  = ($a[0], $a[1], $a[6] * $strand, $a[3]);
	#my ( $mC_num, $dep  )  = @a[4..5];
	
	my ($chr, $pos, $dep ,$mC_num, $val, $type) = ( @a[0..3], $a[4] * $strand, $a[6] );
		
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
