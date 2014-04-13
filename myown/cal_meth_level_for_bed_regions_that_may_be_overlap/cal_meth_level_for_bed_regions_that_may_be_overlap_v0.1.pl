#!/usr/bin/perl -w

#v0.1
#add more column to output more information
#output col
# C_num   -> [1,2, 3, ...]
# covered_C_num
# wmC wmC_per
# mmC mmC_per
# fmC fmC_per

#input_liet format
# chr start end strand

BEGIN { push @INC, '/Users/tang58/scripts_all/perl_code/Modules' }
use Kai_Module;

use strict;
use File::Spec;

my $debug = 0;
if($debug){
	print STDERR "debug = 1\n\n";
}

print STDERR "\n input bed must have head \n\n";

my $detail = 0;
if($debug == 0){
	$detail = 0;
}

my $E_or_nE = "nE";


#my  $usage = "$0 \n<dep_cutoff> <sample_label> <isMeth_file> <feature_bin_num> <flaking_bp> <flaking_bin_num> <input_bed_list> <outdir> <feature_pre>\n\n";
#die $usage unless(@ARGV == 9);
#time /Users/tang58/Kai_BS/myown/cal_meth_level_for_bed_regions_that_may_be_overlap/cal_meth_level_for_bed_regions_that_may_be_overlap_v0.1.pl 
#4 ../colA_isMeth_chrC_error_separately_called.txt  colA Gene_all_coordinate_SAF_head.txt 2 TAIR10_all_Gene .

my  $usage = "$0 \n<dep_cutoff> <isMeth_file>  <sample_label>  <input_bed_list>  <chr_col> <outpre> <outdir>\n\n";
die $usage unless(@ARGV == 7);

my $dep_cutoff = shift or die;
my $isMeth_file = shift or die "isMeth_file";
my $sample_label = shift or die;


my $input_list_file = shift or die "input_list_file";
my $chr_col = shift or die; # "chr_col records which column  is the chr "

my $outpre = shift or die;

my $outdir   = shift or die "outdir";
#my $feature_pre   = shift or die "outpre";

 
die unless (-e $input_list_file);
die unless (-d $outdir);

die unless (-e $isMeth_file and -e $input_list_file);

my $output = File::Spec->catfile($outdir, $outpre . "_in_" . $sample_label .  "_dep" . $dep_cutoff ."_meth_level.txt");

print STDERR "output:\n$output\n\n";
die if (-e $output);



if ($debug){
	#print STDERR $output, "\n\nOK\n\n";
	exit;
}

open(OUT, ">$output") or die;
 
#my (@target_list , @left_list);
#read_list($target_file, \@target_list);
#read_list($left_file,   \@left_list);
#read_list($input_list_file, \@input_list);
my @input_list;
my %pos_in_list_h;

Kai_Module::read_list_recored_pos_and_list($input_list_file, \%pos_in_list_h, \@input_list, $chr_col );

#h{1}->{CG/CHG/CHH/C}
my %seqed_dep_h;
my %seqed_mC_h;
my %seqed_per_h;
my %seqed_isMeth_h;

my %C_num_h ;
my %covered_C_num_h ;

=head
$seqed_dep_h{}->{} +=     ;
$seqed_mC_h{}->{} +=       ;
$seqed_per_h{}->{} +=      ;
$C_num_h{}->{} +=          ;
$covered_C_num_h{}->{} +=  ;
=cut

#chr     pos     strand  type    num_C   depth   percentage      isMeth
#chr1    109     +       CG      10      11      0.909091        1
#chr1    114     +       CHG     2       12      0.166667        1
print STDERR "Begin read DB\n";
open(DB, $isMeth_file) or die "db";
my $db_head = <DB>;

while (<DB>) {
	chomp;
	my @a = split "\t";
	my ($chr, $pos, $strand, $type, $mC, $dep, $per, $isMeth) = @a;
	$chr = Kai_Module::simple_chr($chr);
	next unless ( defined $pos_in_list_h{$chr}->[$pos]  );
	
	if ( $isMeth == 0 and  $E_or_nE eq "E") {
		$mC = 0;
		$per = 0;
	}
	
	my @indexes = split "," , $pos_in_list_h{$chr}->[$pos];
	for my $i(@indexes){
		
		$C_num_h{$i}->{$type} += 1;
		$C_num_h{$i}->{"C"}   += 1;
		
		if ( $dep >= $dep_cutoff) {
			$seqed_dep_h{$i}->{$type}     += $dep;
			$seqed_mC_h{$i}->{$type}      += $mC;
			$seqed_per_h{$i}->{$type}     += $per;
			$covered_C_num_h{$i}->{$type} += 1;
			
			$seqed_dep_h{$i}->{"C"}     += $dep;
			$seqed_mC_h{$i}->{"C"}      += $mC;
			$seqed_per_h{$i}->{"C"}     += $per;
			$covered_C_num_h{$i}->{"C"} += 1;
			
			if ( $isMeth == 1) {
				$seqed_isMeth_h{$i}->{"C"} += 1;
				$seqed_isMeth_h{$i}->{$type} += 1;
			}
			
			
		}
	}
}

my @types = ("CG", "CHG", "CHH", "C" );

my $h = $input_list[0];

my @C_num_head = map { "num_" . $_  } @types;
my @covered_C_num_head = map { $sample_label .  "_covered_num_" . $_  } @types;
my @covered_per_head  = map {$sample_label ."_covered_per_" . $_} @types;
my @wmC_head = map { $sample_label . "_" . $E_or_nE. "_wm" . $_ ,  $sample_label . "_" . $E_or_nE. "_per_wm" . $_ } @types;
my @mmC_head = map { $sample_label . "_" . $E_or_nE. "_mm" . $_ ,  $sample_label . "_" . $E_or_nE. "_per_mm" . $_ } @types;
my @fmC_head = map { $sample_label . "_" . $E_or_nE. "_fm" . $_ ,  $sample_label . "_" . $E_or_nE. "_per_fm" . $_ } @types;
#perl -e ' @types = ("CG", "CHG", "CHH", "C" ); @a = map { "A_". $_, "B_". $_} @types; print join(" ", @a), "\n\n" '
#A_CG B_CG A_CHG B_CHG A_CHH B_CHH A_C B_C


print OUT join("\t", ( $h, @C_num_head, @covered_C_num_head, @covered_per_head, @wmC_head, @mmC_head, @fmC_head  )) ,"\n";

for my $i(1..$#input_list){
	my $this_bed = $input_list[$i];

	my @C_num_print = (0) x 4;
	my @covered_C_num_print = (0) x 4;
	my @covered_per_print  = ( "NA" ) x 4;
	
	my @wmC_print = ( "NA" ) x 8;	
	my @mmC_print = ( "NA" ) x 8;
	my @fmC_print = ( "NA" ) x 8;
	
	Kai_Module::fill_C_num(\@C_num_print, \%C_num_h, $i);
	Kai_Module::fill_C_num(\@covered_C_num_print, \%covered_C_num_h, $i);
	
	Kai_Module::cal_covered_per( \@covered_per_print, \@covered_C_num_print, \@C_num_print  );
	
	Kai_Module::cal_meth_level ( \@wmC_print, \%seqed_mC_h, \%seqed_dep_h, $i);
	Kai_Module::cal_meth_level ( \@mmC_print, \%seqed_per_h, \%covered_C_num_h, $i);
	Kai_Module::cal_meth_level ( \@fmC_print, \%seqed_isMeth_h, \%covered_C_num_h, $i);
	
	print OUT join("\t", ( $this_bed,  @C_num_print, @covered_C_num_print, @covered_per_print, @wmC_print, @mmC_print, @fmC_print )), "\n";
	
}
=head
my %seqed_dep_h;
my %seqed_mC_h;
my %seqed_per_h;

my %C_num_h ;
my %covered_C_num_h ;

%seqed_isMeth_h

=cut
# C_num   -> [1,2, 3, ...]
# covered_C_num

# %covered

# wmC wmC_per
# mmC mmC_per
# fmC fmC_per

#output_list($output, \@input_list, \%total_depths, \%total_mCs, \%numbers );

close OUT;
exit;

