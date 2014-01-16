#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;

my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}

my $R_script = "/Users/tang58/Kai_BS/for_publish/10_gene_TE_bin_files_R_cmd_for_all_command.r";
die unless (-e $R_script );



#TAIR10_Transposable_Elements_sorted_bed.txt
#TE_super_family_name.txt
#   my $usage = "$0 \n<input_TE> <input_name> <outdir>\n\n";
#die $usage unless(@ARGV == 3);
#my $input = shift or die;
#die unless (-e $input);
#my $name_file = shift or die;
#die unless (-e $name_file);
#my $outdir = shift or die;
#die unless (-d $outdir);

my $name_file = "/Users/tang58/DataBase/TAIR10/TE/TE_super_family_name.txt";
die unless (-e $name_file);

my $usage = "$0 \n<do>\n\n";

die $usage unless (@ARGV ==1 and $ARGV[0] eq "do");

my %records ;
read_name($name_file, \%records);

foreach my $k (sort keys %records){
	print join("\t", ($k, $records{$k})), "\n";
}


foreach my $k (sort keys %records){
	#print join("\t", ($k, $records{$k})), "\n";
	my $label = $records{$k};
#	(\S+)_coordinate_bed.txt$
#	my $cmd0 = "head -1 $input > $output";
#	my $cmd = " less $input  | grep \'$k\$\' | perl -lane \' print join(\"\\t\", \@F\[0..2\], \$F\[5\], \@F\[3..4\], \@F\[6..7\]) \' >> $output";
	my $pattern = $label . "_in";
	my $pre = " ROS1Paper_TE_" . $label;
	my $cmd = " R --slave --vanilla --args $pattern  $pre < $R_script";
	print STDERR $cmd, "\n\n";
	if(!$debug){
	#	`$cmd0`;
		`$cmd`;
	}
}



exit;

#read_name($name_file, \%records);
sub read_name {
	my ($file, $ref) = @_;
	die unless (-e $file);
	open(IN, $file) or die;
	
	while(<IN>){
		chomp;
		my $ori = $_;
		next if (/TE_Super_Family/);
		if (/\//){
			my @b = split /\//, $ori;
			my $label = join("_", @b);
			$ref->{$ori} = $label;
		}elsif(/\?/){
			my $label = "LINE_like";
			$ref->{$ori} = $label;
			#print STDERR "\n\nhere\n\n";
		}else{
			$ref->{$ori} = $ori;
		}
	}
	close IN;
	
}
