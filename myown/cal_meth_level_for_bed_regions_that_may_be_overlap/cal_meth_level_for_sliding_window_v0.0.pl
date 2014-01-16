#!/usr/bin/perl -w


# Nov_29 v0.0
# use Ecker suggested, mC number of  isMeth = 0 will be set to 0

#output head
#chr start end cg_num chg chh CXX(A/B=) CXX_per_label |||  mean_CXX


use strict;
use File::Spec;

my $debug = 0;

my $run_debug = 0;

if($run_debug){
	$debug = 0;
}

if($debug){
	print STDERR "debug:$debug\n";
	#print STDERR "input is special,read the script\n";
}

#my $window_size = 1000;
#my $allowed_gap = 1000;

my $bin_size     = 200;
my $sliding_size = 50;
my $const = 4000000;
#my $allowed_gap  = 100;


my %chr_len = ("chr1"=>30427671, "chr2"=>19698289, "chr3"=>23459830, "chr4"=>18585056, "chr5"=>26975502);

if($run_debug){
	#$debug = 0;
	print STDERR "run_debug = 1 \n\n";
	%chr_len = ("chr1"=> 62966 ); 
}

foreach my $chr (sort keys %chr_len){
	
	print join("\t", ($chr, $chr_len{$chr})) , "\n";
}

my $usage = "$0\n\n <input_isMeth> <outdir> <sample_lable> <dep_cutoff> \n\n";
die $usage unless (@ARGV == 4);

my $isMeth_file = shift or die "input";
my $outdir = shift or die "outdir";
my $pre = shift or die "pre";

my $dep_cutoff = shift or die "dep";


my $covered_num = int ($bin_size / $sliding_size);

die unless (-e $isMeth_file and -d $outdir);

my $output = File::Spec->catfile($outdir, $pre . "_WinSize" . $bin_size . "_sliding" . $sliding_size . "_MethLevel_dep" .  $dep_cutoff . "_EckerSuggest.txt" );

print STDERR "output:\n $output\n\n";


if ($debug){
	print STDERR "\nOK\n\n";
	exit;
}

die "output exists!!!\n\n" if(-e $output);
open(OUT, ">>$output") or die "$output exists!!!\n\n";

#perl -e '@a = ("CG","CHG"); $b = "col"; @c = map{$b."_". $_} @a; print join("\n", @c), "\n\n" ' 
#col_CG
#col_CHG

my @types = ("CG", "CHG", "CHH", "C");

print STDERR "\n\n $pre \n\n";

#	map {my $_ . "_num" } @types,
#	map {$_ . "_num_isMet_" . $pre } @types,
#	map {$_ . "_" . $pre } @types,
#	map {$_ . "_per_" . $pre } @types,
#	map { "mean_" . $_ . "_" . $pre } @types
		
		#	map {my $_ . "_num_isMet_"  } @types,
		#	map {my $_ . "_"  } @types,
		#	map {my $_ . "_per_" } @types,
		#	map { "mean_" . my $_ . "_"  } @types
		
#	"CG_", "CHG_", "CHH_", "C_",

print OUT join("\t",
	       ( "chr", "start", "end",
			"CG_num", "CHG_num", "CHH_num", "C_num",
			"mCG_num"  . "_" . $pre, "mCHG_num". "_" . $pre, "mCHH_num". "_" . $pre, "mC_num". "_" . $pre,
			"CG_". $pre, "CHG_". $pre, "CHH_". $pre, "C_". $pre,
			"CG_per" . "_" . $pre, "CHG_per" . "_" . $pre, "CHH_per" . "_" . $pre, "C_per" . "_" . $pre,
			"mean_CG" .  "_" . $pre, "mean_CHG" .  "_" . $pre, "mean_CHH" .  "_" . $pre, "mean_C" .  "_" . $pre, 
		)
	       ) , "\n";

#close OUT;


for my $i_chr (1..5){
	print STDERR  "reading chr$i_chr\n";
	my ( @types_temp,  @mCs_temp,  @depths_temp, @per_temp );
	
	open(DB, $isMeth_file) or die "db";
	my $curr_chr = "chr" . $i_chr;
	my $line = 0;
	my $const = 3000000;
	my $db_head = <DB>;
#0	chr
#1	pos
#2	strand
#3	type
#4	num_C
#5	depth
#6	percentage
#7	isMeth
	
	while(<DB>){
		$line++;
		if($line % $const == 0){
			print STDERR $line, "\t";
		}
		next unless(/$curr_chr/i);
		my @a = split "\t";
		
		next unless ( $a[5] >= $dep_cutoff );
		die unless (lc($a[0]) eq $curr_chr);

		my $pos = $a[1];
		my $type = $a[3];
		
		my $isMet = $a[7];
		
						
		$types_temp[$pos]   = $type;
		$depths_temp[$pos] = $a[5];
		
		if( $isMet == 1){
			$mCs_temp[$pos] = $a[4];
			$per_temp[$pos] = $a[6];
		}elsif( $isMet == 0 ){
			next;			
		}else{
			die "wrong isMeth\n\n", $_;
		}
	}
	close(DB);
	print STDERR "\n\n";
	
	
	my $chr_length = $chr_len{$curr_chr};
	
#	my $covered_tiles = int( $chr_length / $sliding_size) + 1;
	my $covered_tiles = int( ( $chr_length - $bin_size)/ $sliding_size) + 1;
	#foreach my $i (0..($covered_tiles - 2)){
	foreach my $i (0..($covered_tiles - 1)){
		my $start = $i * $sliding_size + 1;
		my $end   = $start + $bin_size - 1;
		
		my @out = cal_everything_in_interval( $start, $end,   \@types_temp,  \@mCs_temp,  \@depths_temp, \@per_temp );
		
		print OUT join("\t", ($curr_chr, $start, $end, @out)), "\n";
	}
	
	#my $last_index = $covered_tiles - 1;
	my $last_index = $covered_tiles ;
	
	my $last_start = $last_index * $sliding_size + 1;
	my $last_end   = $chr_length;
	
	my @last_out = cal_everything_in_interval( $last_start, $last_end,   \@types_temp,  \@mCs_temp,  \@depths_temp, \@per_temp );
	
	print OUT join("\t", ($curr_chr, $last_start, $last_end, @last_out)), "\n";
}
close OUT;
exit;

#my $line = cal_everything_in_interval( $start, $end,   \@types_temp,  \@mCs_temp,  \@depths_temp, \@per_temp );
sub cal_everything_in_interval{
	my ( $start_sub, $end_sub,   $types_ref,  $mCs_ref,  $depths_ref, $per_ref) = @_;
	
	my (%mC_h, %dep_h, %nums_h, %isMet_h, %mean_h);
	
	my @types_in_sub = ("CG", "CHG", "CHH", "C");

	foreach my $t (  @types_in_sub ){
		$nums_h{$t} = 0;
		$isMet_h{$t} = 0;
		$mC_h{$t} = 0 ;
		$dep_h{$t} = 0 ;
		
		@{$mean_h{$t}} = ();
		
	}
	
	for my $i ( $start_sub..$end_sub){
		
		if(defined $types_ref->[$i]){
			my $type_sub = $types_ref->[$i];
			my $depth_sub = $depths_ref->[$i];
			
			my $mC_sub = 0;
			my $per_sub = 0;
			
			if(defined $mCs_ref->[$i] ) {
				$mC_sub  =  $mCs_ref->[$i];
				$per_sub =  $per_ref->[$i];
				
				$isMet_h{$type_sub}++;
				$isMet_h{"C"}++;
			}
			
			$nums_h{$type_sub} ++;
			$nums_h{"C"} ++;
			
			$mC_h{$type_sub} += $mC_sub;
			$mC_h{"C"} += $mC_sub;
			
			$dep_h{$type_sub} += $depth_sub;
			$dep_h{"C"} += $depth_sub;
			
			push @{$mean_h{$type_sub}} , $per_sub;
			push @{$mean_h{"C"}} , $per_sub;
			
		}
		
	}
	
	my @results;
	
	foreach my $t (  @types_in_sub ){
		push @results, $nums_h{$t};
	}
	foreach my $t (  @types_in_sub ){
		push @results, $isMet_h{$t};
	}
	
	foreach my $t (  @types_in_sub ){
		push @results,  $mC_h{$t} . "/" . $dep_h{$t}  . "=";
	}
	
	foreach my $t (  @types_in_sub ){
		my $weight_level = "NA";
		if ( $dep_h{$t}  != 0){
			$weight_level = sprintf("%.4f",  100 * $mC_h{$t} / $dep_h{$t}   );
		}
		
		push @results, $weight_level;
		
	}
	
	foreach my $t (  @types_in_sub ){
		my $mean_level = cal_mean(\@{$mean_h{$t}});
		push @results, $mean_level;
	}
	
	return @results;
	
}

sub cal_mean{
	my ($ref) = @_;
	my $sum = 0;
	my $num = 0;
	
	foreach my $item (@{$ref}){
		$sum+=$item;
		$num++;
	}
	if($num == 0){
		return "NONE";
	}else{
		return ( sprintf ( "%.4f", 100 * $sum / $num));
	}
}
