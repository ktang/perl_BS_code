#!/usr/bin/perl -w

use strict;


my $usage ="$0 <indir> <outdir>";

die $usage unless (@ARGV == 2);

my ($indir, $outdir) = @ARGV[0..1];

die "wrong indir" unless (-d $indir);

die "wrong outdir" unless (-d $outdir);

opendir(DIR, $indir);

my @files = grep /\.fastq$/ , readdir DIR;

print join ("\n", @files), "\n\n";

open (STAT, ">>$indir/base_composition_stats.txt") or die "cannot open base_composition.txt:$!";

my ($total_t, $A_t, $T_t, $C_t, $G_t, $line_t) =  (0, 0, 0, 0, 0, 0);

foreach my $input (@files){
	my $output = "NONE";
	if($input =~ /(\S+)\.fastq$/){
		$output = $1."_brat_input.txt";
	}
	else{
		die "wrong name";
	}
	
	die "wrong input or output" unless (-e "$indir/$input" and !(-e "$outdir/$output"));
	open (IN, "$indir/$input")  or die "cannot open $input:$!";
	open (OUT, ">$outdir/$output") or die "cannot open $output:$!";
	
	print STDERR "handling $input...\n";
	
	my $i = 0;
	
	my ($total , $A, $T, $C, $G) = (0, 0, 0, 0, 0);
	
	while(my $l = <IN>){
		$i++;
		if ($i % 5000000 == 0){
			print STDERR "$i\tDONE\n";
		}
		if ($i % 4 == 2){
			print OUT $l;
			chomp $l;
			$total += (length($l));
			$A += ( $l =~ tr/A/A/);
			$T += ( $l =~ tr/T/T/);
			$C += ( $l =~ tr/C/C/);
			$G += ( $l =~ tr/G/G/);
		}
	}
	close (IN);
	close (OUT);
	print STDERR "handling $input...DONE\n";

	my $num = $i / 4;
	my $mean = sprintf("%.2f",$total / $num );
	print STAT "$input\n";
	my $A_p = sprintf("%.2f",100 * $A / $total );
	my $T_p = sprintf("%.2f",100 * $T / $total );
	my $C_p = sprintf("%.2f",100 * $C / $total );
	my $G_p = sprintf("%.2f",100 * $G / $total );
	print STAT "total_num:$num\tmean_length:$mean\n";
	print STAT "total_base:$total\tA:$A($A_p%)\tC:$C($C_p%)\tG:$G($G_p%)\tT:$T($T_p%)\n";
	$total_t += $total;
	$A_t += $A;
	$T_t += $T;
	$C_t += $C;
	$G_t += $G;
	$line_t += $num;	
}

my $mean = sprintf("%.2f",$total_t / $line_t );

my $A_p = sprintf("%.2f",100 * $A_t / $total_t );
my $T_p = sprintf("%.2f",100 * $T_t / $total_t );
my $C_p = sprintf("%.2f",100 * $C_t / $total_t );
my $G_p = sprintf("%.2f",100 * $G_t / $total_t );
print STAT "all above together:\n";
print STAT "total_reads:$line_t\tmean_length:$mean\n";
print STAT "total_base:$total_t\tA:$A_t($A_p%)\tC:$C_t($C_p%)\tG:$G_t($G_p%)\tT:$T_t($T_p%)\n";

close(STAT);