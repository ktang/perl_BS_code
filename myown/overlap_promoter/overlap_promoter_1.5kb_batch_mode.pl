#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;

my $debug = 1;

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 <indir> <outdir> <pattern>";
die $usage unless(@ARGV == 3);

my $indir  = shift or die;
my $outdir = shift or die;

my $pattern = shift or die;

die unless (-d $indir  );
die unless (-d $outdir );

opendir(DIR, $indir) or die;

my @files = grep /$pattern/ , readdir DIR;


my $promoter_file_1500bp = "/Users/tang58/DataBase/TAIR10/GFF/TAIR10_GFF3_protein_coding_gene_loci_5chr_sorted_1.5kb_promoter_list.txt";
die unless (-e $promoter_file_1500bp);

my @records;
my %promoter_pos;
read_promoter_file(\@records, \%promoter_pos, $promoter_file_1500bp);

#my $input = shift or die "input";
#my $output = shift or die "output";
#die unless (-e $input);
#die if(-e $output);


if ($debug){
        print STDERR join("\n", @files),"\n";
}

foreach my $file(@files){
        if( $file =~ /(\S+)\.txt$/){
                my $output = File::Spec->catfile($outdir,  $1 . "_1.5kb_promoter.txt");
				my $input = File::Spec->catfile ($indir, $file);
				die unless (-e $input);
				die if(-e $output);
				
				open(IN, $input) or die;
				open(OUT, ">>$output") or die;
				my $head = <IN>;
				chomp $head;
				print OUT join("\t", ($head, "promoter", "gene_annotation")), "\n";

				while(<IN>){
					chomp;
					my @a = split "\t";
					my ($chr, $start, $end) = @a[0..2];
					my %index = ();
					for my $i($start..$end){
						if(defined $promoter_pos{$chr}->[$i]){
							my @b = @{$promoter_pos{$chr}->[$i]};
							foreach my $j (@b){
								$index{$j} = 1;
							}
						}
					}
	
					my $pro  = "NONE";
					my $anno = "NONE";
	
					if( scalar(keys %index) == 0){
					#	print OUT join("\t", (@a, $pro, $anno)), "\n";
					}else{
						foreach my $i (sort {$a<=>$b} keys %index){
							my @b = split "\t", $records[$i];
							my $gene = $b[3];
							my $des  = $b[7];
			
							if($pro eq "NONE"){
								$pro = $gene;
								$anno = $des;
							}else{
								$pro  .= ";" . $gene;
								$anno .= ";" . $des;
							}
						}
					}
					print OUT join("\t", (@a, $pro, $anno)), "\n";
				}

				close(IN);


		}else{
                die $file;
        }
}

exit;

sub round {
    my($number) = shift;
    #return int($number + .5);
    return int($number + .5 * ($number <=> 0)); # take care of negative numbers too

}
# 0			1			2			3				4			5			6			7
#chr	pro_start	pro_end	protein_coding_gene	gene_start	gene_end	strand	short_desciption
sub read_promoter_file{
	my ($ref_annotation_array, $ref_pos, $file) = @_;
	die unless (-e $file);
	open(IN, $file) or die;
	my $head = <IN>;
	chomp $head;
	$ref_annotation_array->[0] = $head;
	my $i = 0;
	while(<IN>){
		$i++;
		chomp;
		my @a = split "\t";
		my ($chr, $start, $end) = @a[0..2];
		$ref_annotation_array->[$i] = $_;
		for my $j($start..$end){
			push @{$ref_pos->{$chr}->[$j]}, $i;
		}
	}
	
	close IN;
	print STDERR "read_promoter_file DONE \n\n";	
}