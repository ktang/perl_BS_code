#!/usr/bin/perl -w
# given acgt-count output, calculate the depth of coverage at each cytosine position
# and percentage of total cytosines that was covered by at least one or two reads.

#chr1    33      33      CHG:13  0.153846        -
#chr1    79      79      CHH:27  0.148148        -
# 0		  1		 2			3		4			 5

#v0.1 output depth file

use strict;
use File::Spec;

#my %chr_len = ("chr1"=>30427671, "chr2"=>19698289, "chr3"=>23459830, "chr4"=>18585056, "chr5"=>26975502, "chrM"=>366924, "chrC"=>154478, "pUC19"=>2686);


my $debug = 0;

my $dep_cutoff = 4;

my $constant = 2000000;
if ($debug) {$constant = 300 ;}

if($debug){
	print  STDERR "debug = $debug\n\n";
}

my $usage = "$0 \n <isMeth_file> <outdir> <pre> \n\n";
die $usage unless (@ARGV == 3);

#my $indir = shift or die "indir";
#my $infile = shift or die "infile";
my $isMeth_file = shift or die "isMeth";
my $outdir = shift or die "outdir";
my $pre = shift or die "pre";

#die "wrong indir" unless (-d $indir);
die "wrong outdir" unless (-d $outdir);

#my $in_forw = File::Spec->catfile($indir, $pre . "_forw.txt");
#my $in_rev  = File::Spec->catfile($indir, $pre . "_rev.txt");
#my $input = File::Spec->catfile($indir, $infile);

my $input = $isMeth_file;
die unless (-e $input);



if($debug){
	print  STDERR "input files:\n";
#	print  STDERR join("\n", ($in_forw, $in_rev)), "\n\n";
	print  join $input, "\n\n";
}
#die "wrong input" unless (-e $in_forw and -e $in_rev );

my $output = File::Spec->catfile($outdir,$pre. "_cytosines_coverage_info_DepCutoff". $dep_cutoff . ".txt");
die "$output exists" if (-e $output);

if($debug){
	print STDERR "output:\n$output\n";
}

#my $depth_file = File::Spec->catfile($outdir,$pre. "_depth_distribution_DepthCutoff.txt" );
#die "depth_file exists" if (-e $depth_file);

#my $mC_depth_file = File::Spec->catfile($outdir,$pre. "_depth_mC_DepthCutoff.txt");
#die if (-e $mC_depth_file);

if($debug){
#	print STDERR $depth_file, "\n\n";
}
if($debug){
	print STDERR "\nOK\n\n";
	exit;
}

#die "die debug\n\n" if ($debug);

die "$output exists\n\n" if(-e $output);

open (OUT, ">>$output") or die;
#open (DEP, ">>$depth_file") or die;
#print DEP join("\t", ("depth", "number")), "\n";

my %deps;

#my $max_dep = 0;

#my (%num_forw, %num_rev); #record cytosine number in reference %num{$chr}->{CG}
#my (%cover_1W, %cover_2W, %cover_1C, %cover_2C, %cover_depth_forw, %cover_depth_rev);
my (%cover_depth_forw, %cover_depth_rev);
my (%cytos_forw, %cytos_rev); #record read cytosine depth in each cytosine refrence

my ( @meth_level_CG, @meth_level_CHG, @meth_level_CHH  );



open(IN, $input) or die;
# 0		 1			2	  3			4		5		6				7
#chr     pos     strand  type    num_C   depth   percentage      isMeth
#chr1    334     +       CHH     0       20      0       0
#chr1    336     +       CHH     0       20      0       0

my $head = <IN>;

while(<IN>){
	chomp;
	my @a = split "\t";
	my ($type, $depth) = ( $a[3], $a[5]);
	
	next unless ( $depth >= $dep_cutoff );

	my $chr = $a[0];
	my $per = $a[6];
	my $C_num = $a[4];
	
	my $strand = $a[2];
	
	$deps{$depth} ++;
#	if($depth > $max_dep){
#		$max_dep = $depth;
#	}
	
	
	if($type eq "CG"){
		$meth_level_CG[$depth]->[$C_num]++;
	}
	elsif( $type eq "CHG" ){
		$meth_level_CHG[$depth]->[$C_num]++;
	}
	elsif( $type eq "CHH" ){
		$meth_level_CHH[$depth]->[$C_num]++;
	}else{
		die $_;
	}
	
	
	if( $strand eq "+" ){
		$cytos_forw{$chr}->{$type} += $C_num;
		$cover_depth_forw{$chr}->{$type} += $depth;

	}elsif( $strand eq "-" ){
		$cytos_rev{$chr}->{$type} += $C_num;
		$cover_depth_rev{$chr}->{$type} += $depth;

	}else{
		die $_;
	}
}

close (IN);


#for my $d (0..$max_dep){
#	my $num = 0;
#	if(defined $deps{$d}){
#		$num = $deps{$d};
#	}
#	print DEP join("\t", ($d, $num)), "\n";
#}


#output_mC_depth($mC_depth_file, \@meth_level_CG, \@meth_level_CHG, \@meth_level_CHH , $max_dep );


my @types = ("CG", "CHG", "CHH");


print OUT "Methylation_level:\n";
print OUT join("\t", ("chr", "CG+", "CG-", "CHG+", "CHG-", "CHH+", "CHH-", "forward", "reverse", "total")), "\n";
foreach my $chr(sort keys %cytos_forw){
	print OUT $chr;
	# server as denominator
	my $total_num = 0;
	my $total_num_f = 0;
	my $total_num_r = 0;
	
	#as numerator
	my $total = 0;
	my $total_f = 0;
	my $total_r = 0;	
	
	foreach my $type(@types){
		my ($f, $r) = (0, 0);
		
		if(defined $cytos_forw{$chr}->{$type}){$f = $cytos_forw{$chr}->{$type}}
		if(defined $cytos_rev{$chr}->{$type}){$r = $cytos_rev{$chr}->{$type}}
		
		$total += ($f + $r);
		$total_f += $f;
		$total_r += $r;
		
		unless(defined $cover_depth_forw{$chr}->{$type}) { $cover_depth_forw{$chr}->{$type} = 0}
		unless(defined $cover_depth_rev{$chr}->{$type} )  { $cover_depth_rev{$chr}->{$type} = 0}
		
		$total_num   += ($cover_depth_forw{$chr}->{$type} + $cover_depth_rev{$chr}->{$type} );
		$total_num_f += $cover_depth_forw{$chr}->{$type} ;
		$total_num_r += $cover_depth_rev{$chr}->{$type};
		
		
	#	print OUT "\t", $f, "/", $cover_depth_forw{$chr}->{$type}, "=", sprintf("%.3f",  $f / $cover_depth_forw{$chr}->{$type} * 100), "%";
	#	print OUT "\t", $r, "/", $cover_depth_rev{$chr}->{$type}, "=",  sprintf("%.3f", $r / $cover_depth_rev{$chr}->{$type} * 100), "%";
		my $per_forw = 0;
		my $per_rev	 = 0;
		if( $cover_depth_forw{$chr}->{$type} != 0 ){
			$per_forw = sprintf("%.3f",  $f / $cover_depth_forw{$chr}->{$type} * 100 );
		}
		
		if( $cover_depth_rev{$chr}->{$type} != 0 ){
			$per_rev = sprintf("%.3f", $r / $cover_depth_rev{$chr}->{$type} * 100);
		}
	
		print OUT "\t", $f, "/", $cover_depth_forw{$chr}->{$type}, "=", $per_forw , "%";
		print OUT "\t", $r, "/", $cover_depth_rev{$chr}->{$type}, "=",  $per_rev  , "%";
		
		
	}
	print OUT "\t";
	
	#print OUT $total_f, "/", $total_num_f, "=", sprintf("%.3f", $total_f  / $total_num_f * 100), "%\t";
	#print OUT $total_r, "/", $total_num_r, "=", sprintf("%.3f", $total_r  / $total_num_r * 100), "%\t";
	#print OUT $total,   "/", $total_num,   "=", sprintf("%.3f", $total / $total_num * 100)     , "%\n";
	
	my ( $per_f,  $per_r,  $per_total ) = (0, 0, 0);
	if( $total_num_f != 0) {$per_f     =  sprintf("%.3f", $total_f  / $total_num_f * 100);}
	if( $total_num_r != 0) {$per_r     =  sprintf("%.3f", $total_r  / $total_num_r * 100);}
	if( $total_num != 0) {$per_total =  sprintf("%.3f", $total / $total_num  * 100);}
	
	print OUT $total_f, "/", $total_num_f, "=", $per_f, "%\t";
	print OUT $total_r, "/", $total_num_r, "=", $per_r, "%\t";
	print OUT $total,   "/", $total_num,   "=", $per_total    , "%\n";

	#print OUT join("\t", ($chr,  , $cover_1W{$chr}/$chr_len{$chr} * 100, "%\n";
}



close(OUT);

exit;

sub output_mC_depth{
	my ($file, $CG_ref, $CHG_ref, $CHH_ref, $max) = @_;
	
	die if(-e $file);
	open(DEP, ">>$file") or die;
	print DEP join("\t", ("depth", "mC", "CG", "CHG", "CHH")), "\n";
	for my $i(0..$max){
		for my $j(0..$i){
			my ($num_cg, $num_chg, $num_chh) = (0, 0, 0);
			if(defined $CG_ref ->[$i]->[$j])  {$num_cg = $CG_ref->[$i]->[$j]}
			if(defined $CHG_ref->[$i]->[$j])  {$num_chg = $CHG_ref->[$i]->[$j]}
			if(defined $CHH_ref->[$i]->[$j])  {$num_chh = $CHH_ref->[$i]->[$j]}
			if( $num_cg + $num_chg + $num_chh  > 0){
				print DEP join("\t", ( $i, $j, $num_cg, $num_chg, $num_chh )), "\n";
			}
		}
	}
	
	close(DEP);
}



sub round {
    my ($number) = shift;
    #return int($number + .5);
    return int($number + 0.5 * ($number <=> 0)); # take care of negative numbers too
}