#!/usr/bin/perl -w

# remove 3' adaptor sequence. If the length of the result seq is less than 20nt,
# remove the short seqs.
# added 8/19/09: make sure there is no 'N' in the seqs and count number of seqs
# that did not pass the Illumina filter but has good adaptor
use strict;
use File::Spec;
use Bio::SeqIO;

my $debug = 0;
my $usage = "$0 <seqs_input> <seqs_adaptor_removed>";

my $adaptor_3 = "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG";
my $min_sRNA_len = 20;

my $adaptor_1st_1nt = substr($adaptor_3, 0, 1);
my $adaptor_1st_2nt = substr($adaptor_3, 0, 2);
my $adaptor_1st_3nt = substr($adaptor_3, 0, 3);
my $adaptor_1st_4nt = substr($adaptor_3, 0, 4);
my $adaptor_1st_5nt = substr($adaptor_3, 0, 5);
my $adaptor_1st_6nt = substr($adaptor_3, 0, 6);
my $adaptor_1st_7nt = substr($adaptor_3, 0, 7);
my $adaptor_1st_8nt = substr($adaptor_3, 0, 8);
my $adaptor_1st_9nt = substr($adaptor_3, 0, 9);
my $adaptor_1st_6nt_rev = reverse $adaptor_1st_6nt;
my $adaptor_1st_9nt_rev = reverse $adaptor_1st_9nt;
if($debug){
	print STDERR "adaptor_1st_9nt_rev: $adaptor_1st_9nt_rev\n";
}

die $usage unless (@ARGV >= 2);
my ($inDir, $outDir) = @ARGV[0..1];
opendir(INDIR, $inDir) or die "Cannot open dir $inDir: $!";
my @inFiles = grep {/_trimmed\.fa$/} readdir INDIR;

foreach my $inFile(@inFiles){
    my ($trim_1,$trim_2, $trim_3, $trim_4, $trim_5, $trim_6, $trim_7, 
	$trim_8, $trim_9_over, $no_trim) = 
	(0, 0, 0, 0, 0, 0, 0, 0, 0, 0); # number of seqs subject to different trim
    my ($no_short_seqs, $no_good_seqs) = (0, 0);
    my ($no_seqs_ge9) = (0);
	my ($no_seqs_ls9, $no_seqs_no) = (0,0);
    print  "Trim 3' adaptor for file $inFile\n";
	my @parts = split /\./, $inFile;
	my $ext = pop @parts;
	my $outFile = File::Spec->catfile($outDir, $parts[0] . "_noadapt"  . "." . $ext);
	$inFile = File::Spec->catfile($inDir, $inFile);
	my $seqin = Bio::SeqIO->new(-file=>$inFile, -format=>'fasta');
	#my $seqout = Bio::SeqIO->new(-file=>">$outFile", -format=>'fasta');
	my $seqout;
	open($seqout, ">$outFile") or die "Can't open $outFile:$!";
	while(my $seq = $seqin->next_seq){
		my $seqstr = $seq->seq;
		
		my $new_seqstr = $seqstr;
		my $flag_no	= 0;	my $flag_ge9	= 0;
		if($debug){
			#print STDERR "Doing seq ", $seq->id, " ", $seq->seq, "\n";
		}
		if($seqstr =~ /(\w*$adaptor_1st_9nt)\w*/){
			$trim_9_over++;
			$new_seqstr = $1;
            if($debug){
				print STDERR "match >= 9nt; new_seqstr: ", $new_seqstr, "\n";
			}

			$new_seqstr = reverse $new_seqstr;
			if($debug){
				print STDERR "rev new_seqstr: ", $new_seqstr, "\n";
			}
            
			$new_seqstr =~ s/$adaptor_1st_9nt_rev/lol/;
			if($debug){
				print STDERR "adaptor_rev replaced: ", $new_seqstr, "\n";
			}
			$new_seqstr = reverse $new_seqstr;
			$new_seqstr =~ s/lol//;			$flag_ge9	= 1;
		}elsif(substr($seqstr, length($seqstr) - 8) eq $adaptor_1st_8nt){
			$trim_8++;
			$new_seqstr = substr($seqstr, 0, length($seqstr) - 8);
		}elsif(substr($seqstr, length($seqstr) - 7) eq $adaptor_1st_7nt){
			$trim_7++;
			$new_seqstr = substr($seqstr, 0, length($seqstr) - 7);
		}elsif(substr($seqstr, length($seqstr) - 6) eq $adaptor_1st_6nt){
			$trim_6++;
			$new_seqstr = substr($seqstr, 0, length($seqstr) - 6);
		}elsif(substr($seqstr, length($seqstr) - 5) eq $adaptor_1st_5nt){
			$trim_5++;
			$new_seqstr = substr($seqstr, 0, length($seqstr) - 5);
		}elsif(substr($seqstr, length($seqstr) - 4) eq $adaptor_1st_4nt){
			$trim_4++;
			$new_seqstr = substr($seqstr, 0, length($seqstr) - 4);
		
		}elsif(substr($seqstr, length($seqstr) - 3) eq $adaptor_1st_3nt){
			$trim_3++;
			$new_seqstr = substr($seqstr, 0, length($seqstr) - 3);
		}elsif(substr($seqstr, length($seqstr) - 2) eq $adaptor_1st_2nt){
			$trim_2++;
			$new_seqstr = substr($seqstr, 0, length($seqstr) - 2);
		}elsif(substr($seqstr, length($seqstr) - 1) eq $adaptor_1st_1nt){
			$trim_1++;
			$new_seqstr = substr($seqstr, 0, length($seqstr) - 1);
		}else{
			$no_trim++;	$flag_no	= 1;
		}
		
		if(length($new_seqstr) < $min_sRNA_len){
			$no_short_seqs++;
		}elsif ($flag_no == 0) {
			$no_good_seqs++;
			#$seq->seq($new_seqstr);
			if ($flag_ge9	== 1) {
                    $no_seqs_ge9++;
					#$seqout->write_seq($seq);
					print $seqout ">",$seq->id, "\n";
					print $seqout $new_seqstr, "\n";
			}else{
					$no_seqs_ls9++;
				    #$seqout->write_seq($seq);
					print $seqout ">", $seq->id, "\n";
					print $seqout $new_seqstr, "\n";
			}
		}elsif ($flag_no == 1) {
			$no_good_seqs++;
			#$seq->seq($new_seqstr);
			$no_seqs_no++;
			#$seqout->write_seq($seq);
			print $seqout ">", $seq->id, "\n";
			print $seqout $new_seqstr, "\n";
		}
	}
    $seqin->close;
	close $seqout;


    print  "Number of seqs that are trimmed by different length:\n";
    print join("\t", (">=9", "8", "7", "6", "5", "4", "3", "2", "1", "0")), "\n";
    print join("\t", ($trim_9_over, $trim_8, $trim_7, $trim_6, $trim_5, $trim_4, $trim_3, $trim_2, $trim_1, $no_trim)), "\n";
    print "Number of short seqs (< $min_sRNA_len nt): $no_short_seqs", " Number of good seqs: $no_good_seqs\n";

	print ">= 9nt: $no_seqs_ge9 ", "<9nt: $no_seqs_ls9 ","no adaptor: $no_seqs_no\n";

}
