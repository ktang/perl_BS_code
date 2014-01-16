#!/usr/bin/perl -w

#this overlap script is  for Jacobsen method
#but cal overlap for all files
use strict;
use File::Spec;

my $debug = 0;

my $p = 1;
if($debug){
	print STDERR "debug = $debug\n\n";
}


my $usage = "$0 \n<indir> <output>\n\n";

die $usage unless (@ARGV == 2);
my $dir = shift or die $usage;
my $output = shift or die "output";
die if(-e $output);

die "wrong dir" unless (-d $dir );

opendir (DIR, $dir) or die "cannot open $dir: $!";

#my @files = grep /_Cnum4\.txt$/, readdir DIR;
my @files = grep /_Cnum4.*\.txt$/, readdir DIR;

my @labels;
foreach my $i(0..$#files){
	my $file = $files[$i];

	if($file =~ /(\S+)_JacobsenCell_Bin100_FDR0.01_sDep4_lDep4_Cnum4/){
		$labels[$i] = $1;
	}else{
		die $file;
	}
}

if ($debug){
	print STDERR join("\n", @files),"\n\n";
	print STDERR join("\n", @labels),"\n";
	
	exit;
}

my @nums;

for my $i(0.. ($#files)){
	for my $j ( 0..$#files){
		$nums[$i][$j] = "--";
	}
}

if(!$debug){
	open( OUT, ">>$output" ) or die;
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
			$nums[$i][$j] = overlap($first, $second);
		}
		else{
			print STDERR join("\t", ($first, $second)), "\n";
		}
	}
}

if(!$debug){

	print OUT join("\t", (1..scalar(@labels), "label" ) ), "\n";
	for my $i(0.. ($#files)){
#		print  OUT  $i+1 . "_" . $labels[$i], "\t";
		for my $j ( 0..$#files){
			print OUT $nums[$i][$j];
			if($j != $#files){print OUT "\t"}
			else{
				print OUT "\t",  $i+1 . "_" . $labels[$i],"\n";
			}
		}
#		print  OUT  $i+1 . "_" . $labels[$i], "\t";

	}
}

#print  OUT "\n";

exit;

sub overlap{
	my ($first, $second) = @_;
	open (IN1, $first)  or die "cannot open $first: $!";
	open (IN2, $second) or die "cannot open $second: $!";

	my %intervals;
	my $num_f2 = 0;
	my $h = <IN2>;
	while (<IN2>){
		$num_f2++;
		chomp;
		my @a = split /\t/ ;
		$intervals{$a[0]} = 1;
	}

	my $num_f1 = 0;
	my $num_overlap = 0;
	$h = <IN1>;
	while (<IN1>){
		$num_f1 ++;
		chomp;
		my @a = split /\t/ ;
		if ( defined $intervals{$a[0]} ){$num_overlap++}
	}
	my $per = sprintf ("%.1f",100 * $num_overlap / $num_f1);
	return "$num_overlap/$num_f1 = $per" . "%";
}
