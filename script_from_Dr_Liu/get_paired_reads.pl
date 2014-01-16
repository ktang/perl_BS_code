#!/usr/bin/perl -w
# pick up paired reads (reads from both direction exist and have length >= 20nt).
# output is one sequence per line, without sequence name
# singletons will be put into seperate files
use strict;
use Bio::SeqIO;

my $min_seq_len = 20;  # min seq len for a seq to keep
my $usage = "$0 <input dir> <output dir>";
die $usage unless(@ARGV >= 2);
my ($indir, $outdir) = @ARGV[0..1];
opendir(INDIR, $indir) or die "Can't open $indir: $!";
#my @files1 = readdir INDIR;
#if(@files1 < 1){
#	die "number of files is 0";
#}
my @files1 = grep {/_1_sequence_trimmed_noadapt\.fa/} readdir INDIR;
if(@files1 < 1){
   die "number of files is 0";
}
foreach my $file1(@files1){
	#if($file1 =~ /2\.fasta/){
	#	next;
	#}
	my ($n_short_1, $n_short_2) = (0, 0);
	my ($n_pair, $n_single1, $n_single2) = (0, 0, 0);
	my ($pre1, $pre2) = ("", "");
	if($file1 =~ /(s_\d+_)1(_sequence_trimmed_noadapt\.fa)/){
		$pre1 = $1;
		$pre2 = $2;
	}else{
		die "$file1 does not match name pattern";
	}
	my $file2 = $pre1 . "2" . $pre2;
	my $seqin1 = Bio::SeqIO->new(-file=>"$indir/$file1", -format=>'fasta');
	my $seqin2 = Bio::SeqIO->new(-file=>"$indir/$file2", -format=>'fasta');
	my %seqs2;
	while(my $seq2 = $seqin2->next_seq){
		if(defined $seqs2{$seq2->id}){
			die "seq name ", $seq2->id, " is not unique in file ", $file2;
		}else{
			if($seq2->length >= $min_seq_len){
		        $seqs2{$seq2->id} = $seq2->seq;
			}else{
				$n_short_2++;
			}
		}
	}
	$seqin2->close;
	my $pair1 = $pre1 . "1_pairs.txt";
	my $pair2 = $pre1 . "2_pairs.txt";
	my $single1 = $pre1 . "1_single.txt";
	my $single2 = $pre1 . "2_single.txt";
	open(P1, ">$outdir/$pair1") or die "Can't open $pair1: $!";
	open(P2, ">$outdir/$pair2") or die "Can't open $pair2: $!";
	open(S1, ">$outdir/$single1") or die "Can't open $single1: $!";
	open(S2, ">$outdir/$single2") or die "Can't open $single2: $!";
	while(my $seq1 = $seqin1->next_seq){
		if($seq1->length < $min_seq_len){
			$n_short_1++;
			next;
		}
		if(defined $seqs2{$seq1->id}){
			$n_pair++;
			print P2 $seqs2{$seq1->id}, "\n";
			print P1 $seq1->seq, "\n";
			delete $seqs2{$seq1->id};
		}else{
			print S1 $seq1->seq, "\n";
			$n_single1++;
		}
	}
	foreach my $k(keys %seqs2){
		print S2 $seqs2{$k}, "\n";
		$n_single2++;
	}
    close P1;
	close P2;
	close S2;
	close S1;
	print "Number of short reads in $file1: $n_short_1, $file2: $n_short_2\n";
	print "Number of paired reads: $n_pair, singletons in $file1: $n_single1, $file2: $n_single2\n";
}
