#!/usr/bin/perl -w
# generate input for meth db
use strict;
#use Math::CDF qw(:all);

#my $conver_error = 0.03;
#my $sig_level = 0.05;

my $debug = 0; 

print STDERR "debug:$debug\n";

my $usage = "$0 <pos_file> <indir> <outdir> <out_prefix> [<sample_id>]";

die $usage unless (@ARGV == 4 or @ARGV == 5);

my ($pos_file, $indir, $outdir, $pre) = @ARGV[0..3];
my $sample_id = 0;

if (@ARGV == 5) {$sample_id = $ARGV[4]}
#my $pos_file = $ARGV[0];
#my $d = $ARGV[]; #prefix

#my $pos_file = shift @ARGV;

die "wrong file" unless (-d $indir and -d $outdir );

if ($debug) {$pos_file = "/Users/tang58/try/perl_open/ath_col0_position_debug.txt"}
open(PF, $pos_file) or die "Can't open $pos_file: $!";
my %poss;
while(<PF>){
	next if(/pos_id/);
	chomp;
	my @temp = split /\t/;
	$poss{$temp[2]}->{$temp[3]} = [$temp[0], $temp[5]];
}

close PF;

my $output = $pre . "_meth_db_input.txt";

if (-e "$outdir/$output") {die "output exists"}

open(OUT, ">$outdir/$output") or die "Can't open $output: $!";

print STDERR "output: $outdir/$output\n";

print OUT join("\t", ("SampleID", "pos_id", "depth", "num_C", "percent", "type")), "\n";

opendir(DIR, $indir);
my @files = grep {/(forw|rev)\.txt/} readdir DIR;

if(@files > 2){
	die "more than two forw/rev files: ", join("\t", @files);
}

print STDERR join("\t", @files), "\n";
		   
#if ($debug) { die}

foreach my $file(@files){
	$file = $indir . "/" . $file;
	open(IN, $file) or die "Can't open $file: $!";
	while(<IN>){
		chomp;
		my ($chr, $pos1, $pos2, $t, $percent, $strand) = split /\t/;
		$pos1 += 1;
		my ($pos_id, $type) = (0, "NONE");
		if(defined $poss{$chr}->{$pos1}){
			$pos_id = $poss{$chr}->{$pos1}->[0];
			$type = $poss{$chr}->{$pos1}->[1];
		}else{
			die "chr $chr, pos $pos1, strand $strand, ", 
				  "not defined in position file";
		}
		#my $isMeth = 0;
        my $num_C = round($t * $percent);
                    
		print OUT join("\t", ($sample_id, $pos_id, $t, $num_C, $percent, $type)), "\n";
	}
}

close OUT;

print STDERR "DONE\n";

sub round {
    my($number) = shift;
    #return int($number + .5);
    return int($number + .5 * ($number <=> 0)); # take care of negative numbers too

}

