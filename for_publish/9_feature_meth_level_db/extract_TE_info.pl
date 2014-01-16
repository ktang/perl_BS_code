#!/usr/bin/perl -w
#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);
# my ($volume,$directories,$file) =          File::Spec->splitpath( $path );

use strict;
use File::Spec;

my $TE_db_file = "/Users/tang58/Kai_BS/for_publish/9_feature_meth_level_db/TE/TE_all_coordinate_bed.txt";
die unless ( -e $TE_db_file );
#chr	start	end	strand	ID
#chr1	15827287	15838845	-	AT1TE52125

my $debug = 0;

print STDERR "\n\n\t\tlist has no head\n\n";

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 \n <input> <output>\n\n";
die $usage unless(@ARGV == 2);

my $input = shift or die;
my $output = shift or die;

die unless (-e $input);
die if(-e $output);

my %db_records;
my $head;
read_db($TE_db_file, \%db_records, \$head);

my %TEs;
read_list($input, \%TEs, \%db_records );
output_list($output, \%TEs, $head);

exit;

sub read_db{
	my ($file , $ref, $head_ref) = @_;
	die unless ( -e $file);
	
	open( IN, $file ) or die;
	my $h = <IN>;
	chomp $h;
	${$head_ref} = $h;
	while( <IN> ){
		chomp;
		my @a = split "\t";
		my $id = $a[-1];
		$ref->{$id} = $_;
	}
	close IN;
}

sub read_list{
	my ($file , $ref, $ref_db) = @_;
	die unless ( -e $file);
	
	open( IN, $file ) or die;
	while( <IN> ){
		chomp;
		my @a = split "\t";
		my $TE = $a[0];
		if( !defined $ref_db->{$TE} ){
			print STDERR $TE, "\n";
		}else{
			$ref->{$TE} = $ref_db->{$TE};
		}
	}
	close IN;
}

sub output_list{
	my ($file, $ref, $h) = @_;
	die if ( -e $file);
	
	open ( OUT, ">>$output" ) or die;
	print OUT $h, "\n";
	foreach my $k (sort keys %{$ref}){
		print OUT $ref->{$k}, "\n";
	}
	close OUT;
}