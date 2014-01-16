#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;

my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 \n <outdir> <cutoff> <number_of_files> <isMeth> <label> [<isMeth> <label> ...] \n\n";
die $usage unless(@ARGV >= 7);

my $outdir = shift or die;
my $dep_cutoff = shift or die;
my $num_of_files = shift or die;
my $last_index = $num_of_files - 1;

die unless (-d $outdir);

die unless (@ARGV == 2 * $num_of_files);

my @labels;
my @files;
my @outputs;

for my $i(0..$last_index){
	$files[$i]  = shift or die;
	$labels[$i] = shift or die;
	$outputs[$i] = File::Spec->catfile($outdir,  $labels[$i] . "_depth" . $dep_cutoff . ".wig" );
	die unless (-e $files[$i]);
	die if (-e $outputs[$i]);
}

#if($debug){
	print STDERR join("\n", @labels),  "\n\n";
	print STDERR join("\n", @files),   "\n\n";
	print STDERR join("\n", @outputs), "\n\n";
#}

if($debug){
	exit;
}
my @fhr;
for my $i(0..$last_index){
	open($fhr[$i], "<" , $files[$i]) or die "$files[$i]";
#	my $fhs = $fhr[$i];
#	my $head = < $fhr[$i] >;
	my $head = readline($fhr[$i]);
	
	chomp $head;
	
	if($debug){
#		print STDERR $head, "\n\n";
	}
	
	my @h = split "\t", $head;
	die unless ($h[3] eq "type" and $h[5] eq "depth" and $h[7] eq "isMeth" );
}

my $l;
my %wig_pos;

#0		 1			2		3		4		5		6				7
#chr     pos     strand  type    num_C   depth   percentage      isMeth
#chr1    1       +       CHH     0       0       0       0
while($l = readline($fhr[0]) ){
	
	my $flag_less = 0;
	
	my @a;
	
	chomp $l;
	@{$a[0]} = split "\t", $l;
	
	my ($chr, $pos) = ( $a[0]->[0], $a[0]->[1] );
	my $strand = ($a[0]->[2] eq "+") ? 1 : -1; 
	
	if($a[0]->[5] < $dep_cutoff){
		$flag_less = 1;
	}
	
	for my $i(1..$last_index){
		$l = readline($fhr[$i]);
		chomp $l;
		@{$a[$i]} = split "\t", $l;
		die $i,"\n", $l, "\n" unless ($a[$i]->[1] == $pos);
		
		if($a[$i]->[5] < $dep_cutoff){
			$flag_less = 1;
		}

	}
	
	if($flag_less == 0){
		for my $j(0..$last_index){
			if($a[$j]->[7] == 1){
				$wig_pos{$j}->{$chr}->{$pos} = $strand * $a[$j]->[6];
			}
		}
	}
}

for my $i(0..$last_index){
	close($fhr[$i]);
}

for my $i(0..$last_index){
	output_wig($outputs[$i], \%{$wig_pos{$i}}, $labels[$i]);
}

exit;

sub output_wig{
	my ($file, $ref, $label) = @_;
	
	die if (-e $file);
		open(FILE, ">$file") or die;
	
	print FILE join("\t", ("track", "type=wiggle_0","name=")), "\"$label\"\t","description=\"$label"," methlation\"\t", join("\t", ("viewLimits=-1.0:1.0", "color=255,102,0")), "\n";
	
	foreach my $chr (sort keys %{$ref}){
		print FILE "variableStep\tchrom=$chr\n";
		foreach my $pos(sort {$a<=>$b} keys %{$ref->{$chr}}){
			print FILE join("\t", ($pos, $ref->{$chr}->{$pos})), "\n";
		}
	}
	
	close(FILE);
}