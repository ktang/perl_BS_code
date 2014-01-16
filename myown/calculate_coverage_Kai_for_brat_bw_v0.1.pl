#!/usr/bin/perl -w
# given brat output, calculate the depth of coverage at each genome position
# and percentage of genome that was covered by at least one or two reads.

#v0.2 paired-end reads represent same strand
#pair
#0       AGGTTTGATTGATAA      AACCCTCTTCT      chr3    -       2323351 2323179 0       0       2323351 2323179
#0				1					2			3		4			5		6  7	   8          9      10

#single-end
#7       ATTGGTATATGGATA      chrC    +       55143   0       55143
#0				1				2	  3			4	  5			6
use strict;
use File::Spec;
my %cover;
my %chr_len = ("chr1"=>30427671, "chr2"=>19698289, "chr3"=>23459830, "chr4"=>18585056, "chr5"=>26975502, "chrM"=>366924, "chrC"=>154478, "pUC19"=>2686);

my $debug = 0;
my $constant = 2000000;
if ($debug) {$constant = 300 ;}

if($debug){
	print  STDERR "debug = $debug\n\n";
}

my $usage = "$0 <dir>  <pre>";
die $usage unless (@ARGV == 2);

#my ($indir, $outdir, $pre, $flag) = @ARGV[0..3];

my $dir = shift or die "dir";
#my $outfile = shift or die "outfile";
my $pre = shift or die "pre";

die "wrong dir" unless (-d $dir);

my $in_P    = File::Spec->catfile($dir, $pre ."_paired_brat_bw_out.txt");

my $in_S1   = File::Spec->catfile($dir, $pre . "_1_brat_bw_out.txt"   );
my $in_S2   = File::Spec->catfile($dir, $pre . "_2_brat_bw_out.txt"    );
my $in_left = File::Spec->catfile($dir, $pre . "_paired_brat_bw_out.txt.single_mates");

print  STDERR "brat_file:\n";
print  STDERR join("\n", ($in_P, $in_S1, $in_S2, $in_left)), "\n\n";

die "wrong input" unless (-e $in_P and -e $in_S1 and -e $in_S2 and -e $in_left);

my $output = File::Spec->catfile($dir,$pre. "_coverage_info.txt");
die if (-e $output);
print STDERR "output:\t$output\n\n";

die "die debug\n\n" if ($debug);

open (OUT, ">>$output") or die;



handle_brat_bw_paired_end_output($in_P);

handle_brat_bw_single_end_output($in_S1);
handle_brat_bw_single_end_output($in_S2);
handle_brat_bw_single_end_output($in_left);



my (%cover_1W, %cover_2W, %cover_1C, %cover_2C, %cover_depth);
my $total_len = 0;
foreach my $chr(keys %chr_len){
	$total_len += $chr_len{$chr};
	foreach my $i(1..$chr_len{$chr}){
		if(defined $cover{$chr . "W"}->[$i]){ 
			$cover_depth{$chr . "W"} += $cover{$chr . "W"}->[$i];
			if($cover{$chr . "W"}->[$i] >= 1){
			    $cover_1W{$chr}++;
			    $cover_1W{"total"}++;
			}
			if($cover{$chr . "W"}->[$i] >= 2){
				$cover_2W{$chr}++;
				$cover_2W{"total"}++;
			}
		}
		if(defined $cover{$chr . "C"}->[$i]){
			$cover_depth{$chr . "C"} += $cover{$chr . "C"}->[$i]; 
			if($cover{$chr . "C"}->[$i] >= 1){
			    $cover_1C{$chr}++;
			    $cover_1C{"total"}++;
			}
			if($cover{$chr . "C"}->[$i] >= 2){
				$cover_2C{$chr}++;
				$cover_2C{"total"}++;
			}
		}
	}
}


print OUT "Coverage >= 1:\n";
print OUT "On Watson strand:\n";
foreach my $chr(sort keys %chr_len){
	print OUT $chr, "\t", $cover_1W{$chr}/$chr_len{$chr} * 100, "%\n";
}
print OUT "Whole genome", "\t", $cover_1W{"total"}/$total_len * 100, "%\n";
print OUT "On Crick strand:\n";
foreach my $chr(sort keys %chr_len){
	print OUT $chr, "\t", $cover_1C{$chr}/$chr_len{$chr} * 100, "%\n";
}
print OUT "Whole genome", "\t", $cover_1C{"total"}/$total_len * 100, "%\n";
print OUT "Coverage >= 2:\n";
print OUT "On Watson strand:\n";
foreach my $chr(sort keys %chr_len){
	print OUT $chr, "\t", $cover_2W{$chr}/$chr_len{$chr} * 100, "%\n";
}
print OUT "Whole genome", "\t", $cover_2W{"total"}/$total_len * 100, "%\n";
print OUT "On Crick strand:\n";
foreach my $chr(sort keys %chr_len){
	print OUT $chr, "\t", $cover_2C{$chr}/$chr_len{$chr} * 100, "%\n";
}
print OUT "Whole genome", "\t", $cover_2C{"total"}/$total_len * 100, "%\n";		
print OUT "Depth coverage:\n";
print OUT "On Watson strand:\n";
foreach my $chr(sort keys %chr_len){
	print OUT $chr, "\t", $cover_depth{$chr . "W"} / $chr_len{$chr}, "\n";
}
print OUT "On Crick strand:\n";
foreach my $chr(sort keys %chr_len){
	print OUT $chr, "\t", $cover_depth{$chr . "C"} /$chr_len{$chr}, "\n";
}

exit;

sub handle_brat_bw_paired_end_output{
	my ($file) = @_;
	print  STDERR "paired_end_output:\t$file\n";
	die unless (-e $file);
	open(IN, $file) or die;
	my $count = 0;
	while(<IN>){
		$count++;
		if ($count % $constant == 0) {print  STDERR $count,"\n"}
		chomp;
	#	my @temp = split "\t"; 
		my($id, $seq1, $seq2, $chr, $strand, $start1, $start2, $mis1, $mis2, $ori1, $ori2) = split "\t";
		
		my ($str, $first_start, $first_end, $second_start, $second_end) = ("", "", "", "", "");
		
		if($strand eq "+"){
			($str, $first_start, $first_end, $second_start, $second_end) = ("W", $start1 + 1, $start1 + length($seq1), $start2 + 1, $start2 + length($seq2) );
		}elsif($strand eq "-"){
			($str, $first_start, $first_end, $second_start, $second_end) = ("C", $start2 + 1, $start2 + length($seq2), $start1 + 1, $start1 + length($seq1));
		}else{
			die "strand direction is not + or - for read $id";
		}
		
		
		if($first_start > $second_start){print  STDERR "maybe wrong:\n", $_, "\n"};
		
		if($first_end >= ($second_start)){ #overlap
			if($first_end <= $second_end){ #part overlap
				for my $i( $first_start..$second_end){
					$cover{$chr . $str}->[$i]++;
				}
			}else{# one read in the other read
			#	print  STDERR "second in the first: $_\n";
				for my $i($first_start..$first_end){
					$cover{$chr . $str}->[$i]++;
				}
			}
		}
		else{
			for my $i($first_start..$first_end){
				$cover{$chr . $str}->[$i]++;
			}
				
			for my $i ($second_start..$second_end){
				$cover{$chr . $str}->[$i]++;
			}
		}
	}
	
	close(IN);
}

sub handle_brat_bw_single_end_output{
	my ($file) = @_;
	print  STDERR "single-end output:\t$file\n";
	die unless (-e $file);
	open(IN, $file) or die;
	while(<IN>){
		chomp;
		my ($id, $seq, $chr, $strand, $start, $mismatches, $ori) = split "\t";
		my $str = "";
		if($strand eq "+"){
			$str = "W";
		}elsif($strand eq "-"){
			$str = "C";
		}else{
			die "strand direction is not + or - for read $id";
		}
		
		for my $i( ($start+1)..($start+length($seq)) ){
			    $cover{$chr . $str}->[$i]++;
		}
	}
	close(IN)
}