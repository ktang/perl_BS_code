#!/usr/bin/perl -w
# two bed files as input
# foreach region in first file, check whether there is a region in second file overlap with it
# output the region in last col for first file

# v1.2 output the list
use strict;
use File::Spec;

my $debug = 0;
if($debug){
	print STDERR "debug = 1\n\n";
}

my $usage = "$0 <annotated_file> <bed_file2> <symble_rdd_hyper> <output>";

#my $usage = "$0 <indir> <annotated_file> <bed_file2> <symble_rdd_hyper> <output>";
#my $usage = "$0 <indir> <bed_file1> <bed_file2> <gap> <outdir> <out_pre>";
#/Volumes/Macintosh_HD_2/idm1_new_met_data/script/systematically_parameter_for_region/list_output/v1.4/500_1000_5_10/annotated/idm1_1_hyper_v1.4_WinSize500_gap1000_initialCutoff5_reportCutoff10_annotated_TAIR10_repeat_id.txt
die $usage unless (@ARGV == 4);

print STDERR "\n\nfiles must have headers\n\n";


#my $dir = shift or die "dir";
my $main_file = shift or die "main";
my $bench_file = shift or die "bench";
my $label = shift or die "label";
my $output = shift or die "output";

#die unless (-d $dir);

my $main_input = $main_file;
#my $main_input = File::Spec->catfile($dir, $main_file);
#my $bench_input = File::Spec->catfile($dir, $bench_file);
my $bench_input = $bench_file;

die unless (-e $main_input and -e $bench_input);

die if (-e $output);

if($debug){
	print STDERR "\n\nOK\n\n\n";
	exit;
	
}

open(OUT, ">$output") or die "out";

open (IN2, $bench_input) or die "cannot open $bench_input: $!";

my %records; #record whole line;
my %coors;
my $num_f2 = 0;

#my ($head1, $head2);

my $h2 = <IN2>;

while (<IN2>){
#	next if (/start/);
	$num_f2++;
	chomp;
	my @a = split /\t/ ;
	my $chr = lc $a[0];
	my $start = $a[1] + 0;
	my $end = $a[2] + 0;
	
	#$records{$chr}->{$start} = $_;
	$coors{$chr}->{$start} = $end;
 		 
	for my $i ($start..$end){
		#$beds{$chr}->[$i] = [$chr, $start];
		$records{$chr}->[$i] = $start;
	} 
}

close(IN2);


open (IN1, $main_input) or die "cannot open $main_input: $!";

my $num_f1 = 0;
my $num_overlap = 0;

my $head = <IN1>;
chomp $head;
#my @a_head = split "\t" , $head;
#my $num_head = scalar (@a_head);
my $sym = "overlap_". $label;
#print join ("\t", (@a_head[0..($num_head-2)], $sym, $a_head[$num_head - 1])), "\n";

print OUT $head, "\t", $sym, "\n";

while (<IN1>){

	$num_f1 ++;
	chomp;
	my @a = split /\t/ ;
	my $chr = lc $a[0];
	my $start = $a[1];
	my $end = $a[2];
	my $flag = 0;
	#my ($former_chr, $former_start);
	my ($overlap_start, $overlap_end);
	for my $i ($start..$end){
		if(defined $records{$chr}->[$i]){
			$flag = 1;
			$overlap_start = $records{$chr}->[$i];
			$overlap_end  = $coors{$chr}->{$overlap_start};
		#	($former_chr, $former_start) = @{$beds{$chr}->[$i]};
			last;
		}
	}
	
	if ($flag == 1){
	#	print  join("\t",(@a[0..($num_head-2)], "OVERLAP", $a[$num_head - 1]) ) , "\n";
		print OUT join("\t",(@a, join("_", ($chr,$overlap_start, $overlap_end ))) ) , "\n";
	}
		
	else{
	#	print  join("\t", (@a[0..($num_head-2)], "NOT", $a[$num_head - 1])), "\n";
		print OUT join("\t", (@a, "NOT")), "\n";
	}
}

exit;