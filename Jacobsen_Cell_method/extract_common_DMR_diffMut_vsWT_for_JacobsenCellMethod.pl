#!/usr/bin/perl -w
#extract_commod_DMR_for_JacobsenCellMethod.pl diff mut
# mutants have duplicate and the same WT

#this overlap script is specific for Jacobsen method
use strict;
use File::Spec;

my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 \n <indir1> <mut_lable_rep1> <indir2> <mut_lable_rep2> <outdir>\n\n";
#1:vs_Kaust; (JKZ131_Col0)	2:vs_SH	(Col0_SH_pool)	3: outdir: common_both
die $usage unless(@ARGV == 5);

my $indir1    = shift or die;
my $vs_label1 = shift or die;
my $indir2    = shift or die;
my $vs_label2 = shift or die;
my $outdir    = shift or die;

#die unless (-d $indir);
#die if( -e $output);

die unless (-d $indir1);
die unless (-d $indir2);
die unless (-d $outdir);


my @types = ( "CG", "CHG", "CHH" );
my @hyps   = ("hyper", "hypo");

foreach my $type (@types){
	foreach my $hyp ( @hyps ){
		my @files1 = get_files ( $indir1, $type, $hyp, $vs_label1);
		my @files2 = get_files ( $indir2, $type, $hyp, $vs_label2);
		
		check_files ( \@files1, $vs_label1, \@files2, $vs_label2 );
		
		if($debug){
			print STDERR "\n", join("\t",  ($type, $hyp)), "\n";
			print STDERR join("\n", @files1), "\n\n";
			print STDERR join("\n", @files2), "\n\n";
		}
		
		#cal_overlap($output, $indir,$type, $hyp , \@files );
		for my $i (0..$#files1){
			my $f1 = File::Spec->catfile( $indir1, $files1[$i]);
			my $f2 = File::Spec->catfile( $indir2, $files2[$i]);
			die unless (-e $f1 and -e $f2);
			#my $output = generate_output_full_name ($outdir, $files1[$i], $vs_label1);
			my $output = generate_output_full_name_same_WT ($outdir, $files1[$i], $vs_label1, $vs_label2);
			die if(-e $output);
			if($debug){ print STDERR join("\n", ( $f1, $f2, $output) ), "\n\n" }
			
			if(!$debug){
				extract_DMR ($f1, $f2, $output);
			}
			
		}
		
	}
}

exit;

#extract_DMR ($f1, $f2, $output);
sub extract_DMR{
	my ($in1, $in2, $out) = @_;
	die unless (-e $in1);
	die unless (-e $in2);
	die if(-e $out);
	
	my %records;
	
	open(OUT, ">>$out") or die;
	open(IN1, $in1 ) or die;
	open(IN2, $in2 ) or die;
	
	my $h1 = <IN1>;
	chomp $h1;
	my $h2 = <IN2>;
	chomp $h2;
	my @a_h2 = split "\t", $h2;
	print OUT join("\t", ($h1, @a_h2[2..$#a_h2] )), "\n";
	
	while(<IN1>){
		chomp;
		my @a = split "\t";
		$records{$a[0]} = $_;
	}
	close IN1;
	
	while(<IN2>){
		chomp;
		my @a = split "\t";
		if(defined  $records{$a[0]}  ){
			print OUT join("\t", ( $records{$a[0]}, @a[2..$#a] )), "\n";
		}
	}
	
	close IN2;
	close OUT;
}
# my $output = generate_output_full_name ($outdir, $files1[$i], $vs_label1);
sub generate_output_full_name{
	my ( $dir, $file, $label ) = @_;
	if($file =~ /(\S+)_vs_$label\_(\S+)\.txt$/){
		return File::Spec->catfile( $dir, $1 . "_vs_both_WT_commonDMRs_" . $2 . ".txt" );
	}else{
		die $file;
	}
	
}

sub generate_output_full_name_same_WT{
	my ( $dir, $file, $label, $label2 ) = @_;
	if($file =~ /(\S+)\.txt$/){
		return File::Spec->catfile( $dir, $1 . "_mergedWith_" . $label2 . ".txt" );
	}else{
		die $file;
	}
}


#check_files ( \@files1, $vs_label1, \@files2, $vs_label2 );
sub check_files{
	my ($ref1, $label1, $ref2, $label2 ) = @_;
	my @array_temp1 = @{$ref1};
	my @array_temp2 = @{$ref2};
	my $num1 = $#array_temp1;
	my $num2 = $#array_temp2;
	
	die unless ( $num1 == $num2);
	
	for my $i(0..$num1){
		my $f1 = $array_temp1[$i];
		my $f2 = $array_temp2[$i];
		die "$f1 \n $label1"  unless ( $f1 =~ /$label1/);
		die $f2 unless ( $f2 =~ /$label2/);
		
		$f1 =~ s/$label1//;
		$f2 =~ s/$label2//;
		
		if( $f1 ne $f2){
			print STDERR "f1:$f1\n",  "f2:$f2\n\n";
			die;
		}
		
	}
	
	
}

sub get_files{
	my ( $indir_sub, $type_sub, $hyp_sub, $l_sub  ) = @_;
	
	die unless (-d $indir_sub);
	opendir(DIR, $indir_sub) or die;
	
	my $label = $type_sub . "_" . $hyp_sub;
	
	my @temp = grep /$l_sub/, grep /$label/, readdir DIR;
	
	closedir DIR;
	return @temp;
}