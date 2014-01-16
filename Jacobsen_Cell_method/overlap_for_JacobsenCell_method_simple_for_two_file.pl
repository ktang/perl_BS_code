#!/usr/bin/perl -w

#this overlap script is specific for Jacobsen method
use strict;
use File::Spec;

my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 \n <f1> <f2>\n\n";
die $usage unless(@ARGV == 2);

my $f1 = shift or die;
my $f2 = shift or die;

die unless ( -e $f1 );
die unless ( -e $f2 );

my %records;

open(IN1, $f1) or die;
my $h = <IN1>;

while(<IN1>){
	chomp;
	my @a = split "\t";
	$records{$a[0]} = 1;
}

close IN1;

open(IN2, $f2) or die;
$h = <IN2>;
my $n = 0;
while(<IN2>){
	chomp;
	my @a = split "\t";
	if( defined $records{$a[0]} ){
		$n++;
	}
}

close IN2;
print STDERR $n, "\n";

exit;

# cal_overlap($output, $indir,$type, $hyp , \@files );
sub cal_overlap{
	my ( $output_sub, $indir_sub, $type_sub, $hyp_sub, $ref ) = @_;
	my @files_sub = @{$ref};
	
	my @labels;
	foreach my $i(0..$#files_sub){
		my $file = $files_sub[$i];

		if($file =~ /(\S+)\.txt$/){
			$labels[$i] = $1;
		}else{
			die $file;
		}
	}

	
	my @nums;
	for my $i(0.. ($#files_sub)){
		for my $j ( 0..$#files_sub){
			$nums[$i][$j] = "--";
		}
	}

#	for my $i (0..( $#files_sub - 1 )){
#		for my $j ( ($i + 1).. $#files_sub){
	for my $i(0.. ($#files_sub)){
		for my $j ( 0..$#files_sub){
			my $file_i = File::Spec->catfile( $indir_sub, $files_sub[$i] );
			my $file_j = File::Spec->catfile( $indir_sub, $files_sub[$j] );
			next if ( $file_i eq  $file_j );
			if(! $debug){
				$nums[$i][$j] = overlap($file_i, $file_j);
				#$nums[$j][$i] = $nums[$i][$j];
			}
			
			
		}
	}
	
	if(!$debug){
		open(OUT, ">>$output_sub") or die;
		
		print OUT join("\t", ( $type_sub, $hyp_sub) ), "\n";
		print OUT join("\t", (1..scalar(@labels), "label" ) ), "\n";

		for my $i(0.. ($#files_sub)){
			for my $j ( 0..$#files_sub){
				print OUT $nums[$i][$j];
				if($j != $#files_sub){
					print OUT "\t";
				}
				else{
					print OUT "\t",  $i+1 . "_" . $labels[$i],"\n";
				}
			}
		}
		
		print OUT "\n\n";
		
		close OUT;
	}
}


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


# get_files ( $indir, $type, $hyp );
sub get_files{
	my ( $indir_sub, $type_sub, $hyp_sub  ) = @_;
	
	die unless (-d $indir_sub);
	opendir(DIR, $indir_sub) or die;
	
	my $label = $type_sub . "_" . $hyp_sub;
	
	my @temp = grep /$label/, readdir DIR;
	
	closedir DIR;
	return @temp;
}