#!/usr/bin/perl -w

#for ShangHai use
# difference from v0.0
# the grep regular expression
#system("date");
use strict;
use File::Spec;


my $debug = 1;
my $p = 0;
if($debug){
	print STDERR "debug = $debug\n\n";
}

print STDERR "\n\nfiles must have head\n\n";

#my $script = "/Users/tang58/Weiqiang_idm1_1_Nature_paper/add_annotation_Kai_v1.7_TAIR10.pl";
#my $script  = "./overlap_bed_first_in_second.pl";
#die unless (-e $script);

my $usage = "$0 <indir> STDOUT";

die $usage , "\n\n" unless (@ARGV == 1);
my $dir = shift or die $usage;

die "wrong dir" unless (-d $dir );

opendir (DIR, $dir) or die "cannot open $dir: $!";

my @files = grep /\.txt$/, readdir DIR;

my @labels;
foreach my $i(0..$#files){
	my $file = $files[$i];
	if($file =~ /(\S+)_P0\S+_reduce\S+\.txt$/){
		#$labels[$i] = join("_", ($1, $2, $3));
		$labels[$i] = $1;
	}
}

foreach my $i(0..$#files){
	print STDERR join("\t", ($i + 1, $labels[$i], $files[$i])), "\n";
}

if ($debug){
#	print STDERR join("\n", @files),"\n";
#	print STDERR join("\n", @labels),"\n";
	
	die;
}



my @nums;

for my $i(0.. ($#files)){
	for my $j ( 0..$#files){
		$nums[$i][$j] = "--";
	}
}

for my $i(0.. ($#files)){
	for my $j ( 0..$#files){
		my $first = File::Spec->catfile($dir, $files[$i]);
		my $second = File::Spec->catfile($dir, $files[$j]);
		next if($first eq $second);
		if($p){
			print STDERR join("\t", ($i, $j)), "\n";
		}
		if(!$debug){
			print STDERR "handling ", $labels[$i], " and ", $labels[$j], "\n";
			$nums[$i][$j] = overlap($first, $second);
		}
		else{
			print STDERR join("\t", ($first, $second)), "\n";
		}
	}
}

if(!$debug){
	print join("\t", ("label", @labels)), "\n";
	for my $i(0.. ($#files)){
		print $labels[$i], "\t";
		for my $j ( 0..$#files){
			print $nums[$i][$j];
			if($j != $#files){print "\t"}
			else{
				print "\n";
			}
		}
	}
}

#print "\n\n";

exit;

sub overlap{
	my ($first, $second) = @_;
	open (IN1, $first)  or die "cannot open $first: $!";
	open (IN2, $second) or die "cannot open $second: $!";

	my %beds;
	my $num_f2 = 0;
	my $h = <IN2>;
	while (<IN2>){
#	next if (/Start/);
		$num_f2++;
		chomp;
		my @a = split /\t/ ;
		my $chr = lc $a[0];
		my $start = $a[1] + 0;
		my $end = $a[2] + 0;
	#	my $new_start = $start - $allowed_gap;
	
	#	my $new_start = $start;
	#	if($new_start <= 0 ){$new_start = 1}
	
		#for my $i ( $new_start..($end + $allowed_gap))
		for my $i ($start..$end){
			$beds{$chr}->[$i] = 1;
		}
	}

	my $num_f1 = 0;
	my $num_overlap = 0;
	$h = <IN1>;
	while (<IN1>){
		$num_f1 ++;
		chomp;
		my @a = split /\t/ ;
		my $chr = lc $a[0];
		my $start = $a[1] + 0;
		my $end = $a[2] + 0;
		my $flag = 0;
		for my $i ($start..$end){
			if(defined $beds{$chr}->[$i]){
				$flag = 1;
				last;
			}
		}
	
		if ($flag == 1){$num_overlap++}
	
	}
	my $per = sprintf ("%.1f",100 * $num_overlap / $num_f1);
	return "$num_overlap/$num_f1 = $per" . "%";
}