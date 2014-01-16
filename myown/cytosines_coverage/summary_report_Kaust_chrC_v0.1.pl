#!/usr/bin/perl -w

# chrC as error rate
# v0.1
#

use strict;
use File::Spec;

my $debug = 0;

my $print = 0;

my $usage = "$0 \n <indir> <output>\n\n";
die $usage unless(@ARGV == 2);

my $indir  = shift or die "indir";
my $output = shift or die "output";
die unless (-d $indir);
die if(-e $output);

opendir(DIR, $indir) or die;
my @files = grep /_info\.txt$/, readdir (DIR);
closedir(DIR);

print STDERR join("\n", @files), "\n\n";

foreach my $file(@files){
	if($file =~ /(\S+)_cytosines/){
		my $file  = File::Spec->catfile($indir, $file);
		die unless (-e $file);
	}
	else{
		die;
	}
}

if(!$debug){
	die if (-e $output);
}
else{
	print STDERR "\nOK\n\n";
	exit;
}

open(OUT, ">>$output") or die;

print OUT join("\t", ("sample", "average_depth", "%_of_dep>=1",  "%_of_dep>=2", "CG_meth_level", "CHG_meth_level", "CHH_meth_level", "total_meth_level","error_rate_chrC")), "\n";
#print OUT join("\t", ("sample", "average_depth", "%_of_dep>=1",  "%_of_dep>=2", "CG_meth_level", "CHG_meth_level", "CHH_meth_level", "total_meth_level")), "\n";
my $l;


#my $gap_line_num = 6;

foreach my $infile(@files){
#	my $file = File::Spec->catfile($dir, $pre. "_cytosines_coverage_info.txt");
print STDERR $infile, "\n";
if($infile =~ /(\S+)_cytosines/){
	my $pre = $1;
	print STDERR $pre, "\n";
	my $file  = File::Spec->catfile($indir, $infile);
	die unless (-e $file);
	open(IN, $file) or die "cannot open $file:$!";
	
	my $dep1 = "NONE";
	
	while ( $l = <IN> ){
		next unless ( $l =~ /^chr1/);
		
		if($print ) {
			print STDERR "1:", $l , "\n\n";
		}		
		my ($first, $second) = (0, 0);
		
		chomp $l  ;
		my @a = split "\t", $l ;
		die unless ($a[0] eq "chr1");
		
		my ($t1, $t2) = get_two_value($a[-1]);
		$first += $t1;
		$second+= $t2;
		
		for my $i(2..5){
			$l = <IN>;
			chomp $l  ;
			my @a = split "\t", $l ;
			die unless ($a[0] eq "chr".$i);
			
			my ($t1, $t2) = get_two_value($a[-1]);
			$first += $t1;
			$second+= $t2;
		}
		if($print){
			print STDERR "$first / $second", "\n\n";
		}
		$dep1 = sprintf("%.3f", $first / $second * 100) . "%";
		last;
	}
	
	my $dep2 = "NONE";
	
	while ( $l = <IN> ){
		next unless ( $l =~ /^chr1/);
		my ($first, $second) = (0, 0);
		
		if($print ) {
			print STDERR "2:", $l , "\n\n";
		}
		
		chomp $l  ;
		my @a = split "\t", $l ;
		die unless ($a[0] eq "chr1");
		
		my ($t1, $t2) = get_two_value($a[-1]);
		$first += $t1;
		$second+= $t2;
		
		for my $i(2..5){
			$l = <IN>;
			chomp $l  ;
			my @a = split "\t", $l ;
			die unless ($a[0] eq "chr".$i);
			
			my ($t1, $t2) = get_two_value($a[-1]);
			$first += $t1;
			$second+= $t2;
		}
		if($print){
			print STDERR "$first / $second", "\n\n";
		}
		$dep2 = sprintf("%.3f", $first / $second * 100) . "%";
		last;
	}
	
	my $avg_dep = "NONE";
	while ( $l = <IN> ){
		next unless ( $l =~ /^chr1/);
		my ($first, $second) = (0, 0);
		
		if($print ) {
			print STDERR "3:", $l , "\n\n";
		}
		
		chomp $l  ;
		my @a = split "\t", $l ;
		die unless ($a[0] eq "chr1");
		
		my ($t1, $t2) = get_two_value($a[-1]);
		$first += $t1;
		$second+= $t2;
		
		for my $i(2..5){
			$l = <IN>;
			chomp $l  ;
			my @a = split "\t", $l ;
			die unless ($a[0] eq "chr".$i);
			
			my ($t1, $t2) = get_two_value($a[-1]);
			$first += $t1;
			$second+= $t2;
		}
		if($print){
			print STDERR "$first / $second", "\n\n";
		}
#		$avg_dep = sprintf("%.3f", $first / $second * 100) . "%";
		$avg_dep = sprintf("%.3f", $first / $second );
		last;
	}

#	for my $i(1..$gap_line_num){$l = <IN>}
	my ($first, $second) = (0, 0);
	#Methylation_level
	
	while($l = <IN>){
		next unless ($l =~ /Methylation_level/);
		last;
	}
	
	if($print ) {
		print STDERR $l , "\n\n";
		print STDERR join("\t", ($pre,$avg_dep, $dep1, $dep2) ), "\n\n";

	}
	
	$l = <IN> ; # skip the head
	
	my ($CG_first, $CG_second) = (0, 0);
	my ($CHG_first, $CHG_second) = (0, 0);
	my ($CHH_first, $CHH_second) = (0, 0);
	for my $i(1..5){
		$l = <IN>;
		chomp $l;
		my @a = split "\t",$l ;
		die unless ($a[0] eq "chr".$i);
		my ($t1, $t2) = get_two_value($a[-1]);
		$first += $t1;
		$second+= $t2;
		
		for my $j(1..2){
			my ($temp1, $temp2) = get_two_value($a[$j]);
			$CG_first += $temp1;
			$CG_second += $temp2;
		}
		for my $j(3..4){
			my ($temp1, $temp2) = get_two_value($a[$j]);
			$CHG_first += $temp1;
			$CHG_second += $temp2;
		}
		for my $j(5..6){
			my ($temp1, $temp2) = get_two_value($a[$j]);
			$CHH_first += $temp1;
			$CHH_second += $temp2;
		}
		
	}
	if($print ){
		print STDERR "$first / $second", "\n";
		print STDERR "$CG_first / $CG_second ", "\n";
		print STDERR "$CHG_first / $CHG_second ", "\n";
		print STDERR "$CHH_first / $CHH_second ", "\n";
	}
	my $cg = sprintf("%.3f", $CG_first / $CG_second * 100) . "%";
	my $chg = sprintf("%.3f", $CHG_first / $CHG_second * 100) . "%";
	my $chh = sprintf("%.3f", $CHH_first / $CHH_second * 100) . "%";
	my $meth_level = sprintf("%.3f", $first / $second * 100) . "%";
	$l = <IN>;
	die unless ($l =~ /chrC/);
	chomp $l;
	my @a = split "\t", $l;
	my ($t1, $t2) = get_two_value($a[-1]);
	my $error = $t1 / $t2;
	
	close(IN);
	
	print OUT join("\t", ($pre,$avg_dep, $dep1, $dep2, $cg, $chg, $chh, $meth_level, $error )), "\n";
}

else{
	die $infile;
}

}
close(OUT);
#print OUT join("\t", ("sample", "average_depth", "%_of_dep>=1",  "%_of_dep>=2", "CG_meth_level", "CHG_meth_level", "CHH_meth_level", "total_meth_level","error_rate")), "\n";

exit;

sub get_two_value{
	my ($item) = @_;
	my @a = split "=", $item;
	return (split /\//, $a[0]);
}