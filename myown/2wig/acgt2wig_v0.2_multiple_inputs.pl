#!/usr/bin/perl -w

# acgt2wig_v0.2_multiple_inputs.pl
#purpose: Huiming has GB03 BS-Seq, which the coverage is poor.
# just use the forw/rev file to get the wig file
# only output dep >= cutoff (4) for all files
# include 0/4 = 0
#chr5
#chrC
#chrM
#pUC19

use strict;
use File::Spec;

my $debug = 0;
my $p = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}

#GB03_BS_0320_rev.txt  GB03_BS_0320_forw.txt
# GB03_BS_0320
my $usage = "$0 \n <outdir> <cutoff> <indir> <number_of_labels>  <label1>  <label2> [ <label3>...] \n\n";
die $usage unless(@ARGV >= 6);

my $outdir = shift or die;
my $dep_cutoff = shift or die;
my $indir = shift or die;
my $num_of_files = shift or die;
my $last_index = $num_of_files - 1;

die unless (-d $outdir);

#die unless (@ARGV == 2 * $num_of_files);
die unless (@ARGV ==  $num_of_files);

my @labels;
my @files_forw;
my @files_rev;
my @wigs_C;
my @wigs_CG;
my @wigs_CHG;
my @wigs_CHH;

my @fhw_C;
my @fhw_CG;
my @fhw_CHG;
my @fhw_CHH;

for my $i(0..$last_index){
#	$files[$i]  = shift or die;
	$labels[$i] = shift or die;
	$files_forw[$i] = File::Spec->catfile($indir, $labels[$i] . "_forw.txt");
	$files_rev[$i]  = File::Spec->catfile($indir, $labels[$i] . "_rev.txt");
#	$outputs[$i] = File::Spec->catfile($outdir,  $labels[$i] . "_depth" . $dep_cutoff . ".wig" );
	$wigs_C[$i]   = File::Spec->catfile($outdir,  $labels[$i] . "_". $num_of_files . "f_dep" . $dep_cutoff . "_C.wig" );
	$wigs_CG[$i]  = File::Spec->catfile($outdir,  $labels[$i] . "_". $num_of_files . "f_dep" . $dep_cutoff . "_CG.wig" );
	$wigs_CHG[$i] = File::Spec->catfile($outdir,  $labels[$i] . "_". $num_of_files . "f_dep" . $dep_cutoff . "_CHG.wig" );
	$wigs_CHH[$i] = File::Spec->catfile($outdir,  $labels[$i] . "_". $num_of_files . "f_dep" . $dep_cutoff . "_CHH.wig" );
#	die unless (-e $files[$i]);
	die unless (-e $files_forw[$i]);
	die unless (-e $files_rev[$i]);
#	die if (-e $outputs[$i]);
	die if (-e $wigs_C[$i]);
	die if (-e $wigs_CG[$i]);
	die if (-e $wigs_CHG[$i]);
	die if (-e $wigs_CHH[$i]);
}

if($debug){
#	print STDERR join("\n", @labels),  "\n\n";
#	print STDERR join("\n", @files),   "\n\n";
#	print STDERR join("\n", @outputs), "\n\n";
	
	print STDERR join("\n", @labels), "\n\n";
	print STDERR join("\n", @files_forw), "\n\n";
	print STDERR join("\n", @files_rev), "\n\n";
	print STDERR join("\n", @wigs_C), "\n\n";
	print STDERR join("\n", @wigs_CG), "\n\n";
	print STDERR join("\n", @wigs_CHG), "\n\n";
	print STDERR join("\n", @wigs_CHH), "\n\n";
}

if($debug){
	exit;
}

for my $i(0..$last_index){
	open($fhw_C[$i], ">>" , $wigs_C[$i]) or die "$wigs_C[$i]";
	open($fhw_CG[$i], ">>" , $wigs_CG[$i]) or die "$wigs_CG[$i]";
	open($fhw_CHG[$i], ">>" , $wigs_CHG[$i]) or die "$wigs_CHG[$i]";
	open($fhw_CHH[$i], ">>" , $wigs_CHH[$i]) or die "$wigs_CHH[$i]";
	
	my $label = $labels[$i];
	my $label_removed = $label;

	print {$fhw_C[$i]}   join("\t", ("track", "type=wiggle_0","name=")), "\"$label_removed\"\t","description=\"$label_removed"," methlation\"\t", join("\t", ("viewLimits=-1.0:1.0", "color=255,102,0")), "\n";
	
	print {$fhw_CG[$i]}  join("\t", ("track", "type=wiggle_0","name=")), "\"$label\"\t","description=\"$label"," CG methlation\"\t", join("\t", ("viewLimits=-1.0:1.0", "color=205,205,0")), "\n";
	print {$fhw_CHG[$i]} join("\t", ("track", "type=wiggle_0","name=")), "\"$label\"\t","description=\"$label"," CHG methlation\"\t", join("\t", ("viewLimits=-1.0:1.0", "color=100,149,237")), "\n";
	print {$fhw_CHH[$i]} join("\t", ("track", "type=wiggle_0","name=")), "\"$label\"\t","description=\"$label"," CHH methlation\"\t", join("\t", ("viewLimits=-1.0:1.0", "color=205,155,155")), "\n";
	


}




my %wigs; #$wig{chr}->{pos} = {type, val1, val2...}
read_multiple_acgtFiles(\@files_forw, 1, \%wigs, $dep_cutoff);
read_multiple_acgtFiles(\@files_rev, -1, \%wigs, $dep_cutoff);


foreach my $chr (sort keys %wigs){
	for my $i(0..$last_index){
		print {$fhw_C[$i]} "variableStep\tchrom=$chr", "\n";
	
		print {$fhw_CG[$i]}  "variableStep\tchrom=$chr", "\n";
		print {$fhw_CHG[$i]} "variableStep\tchrom=$chr", "\n";
		print {$fhw_CHH[$i]} "variableStep\tchrom=$chr", "\n";
	
	}
	
	foreach my $pos (sort {$a <=> $b} keys %{$wigs{$chr}}){
		my ($type, @vals) = split ",", $wigs{$chr}->{$pos};
		for my $i(0..$last_index){
			print {$fhw_C[$i]} join("\t", ($pos, $vals[$i])), "\n";
			if($type eq "CG"){
				print {$fhw_CG[$i]}join("\t", ($pos, $vals[$i]) ), "\n";
			}elsif($type eq "CHG"){
				print {$fhw_CHG[$i]} join("\t", ($pos, $vals[$i]) ), "\n";
			}elsif($type eq "CHH"){
				print {$fhw_CHH[$i]} join("\t", ($pos, $vals[$i]) ), "\n";
			}else{
				die "$chr, $pos,\n\n";
			}
		}
	}

}


for my $i(0..$last_index){
	close ($fhw_C[$i]);
	close ($fhw_CG[$i]);
	close ($fhw_CHG[$i]);
	close ($fhw_CHH[$i]);
}
exit;

#0		 1		 2			3		4	 5
#chr1    116     116     CHG:16  0.375   -
#chr1    124     124     CHG:17  1       -
#chr1    125     125     CHH:17  0.470588        -
#eval sprintf("%.3f", $x)
# perl -e '@a= (1, -1,-0.3339, 0.375, 0.00000, 1/178); for $i (@a) {print eval sprintf("%.3f", $i), "\n"}'
# 1
# -1
# -0.334
# 0.375
# 0
# 0.006

# read_multiple_acgtFiles(\@files_rev, -1, \%wigs);
sub read_multiple_acgtFiles{
	my ($file_ref, $strand, $wig_ref, $dep_cutoff) = @_;
	my @fhr;
	my @files = @{$file_ref};
	my $last_index = $#files;
	
	
	for my $i (0..$last_index){
		open( $fhr[$i], "<" , $files[$i]) or die "$files[$i]";
		
	}
	
	my $l;
	while( $l = readline($fhr[0]) ){
		my $flag_less = 0;
		my @a = ();
		my @deps = ();
		my @vals = ();
		my $i = 0;
		

		
		chomp $l;
		@{$a[$i]} = split "\t", $l;
		my ($chr, $pos) = ( $a[$i]->[0], $a[$i]->[1] + 1 );
		last if ($chr =~ /chr[CM]/i);
		my ($type, $dep);
		
		
		if ( $a[$i]->[3] =~ /(C\S+):(\d+)/ ){
			($type, $dep) = ($1, $2);
		}
		else{
			die $l;
		}
		if ( $dep < $dep_cutoff){
			for $i(1..$last_index){
				$l = readline($fhr[$i]);
			}
			next;
		}
		else{
			push @vals, $type;
			my $tmp = eval sprintf("%.3f", $a[$i]->[4]);
			push @vals, $strand * $tmp;
			for $i(1..$last_index){
				$l = readline($fhr[$i]);
				chomp $l;
				@{$a[$i]} = split "\t", $l;
				die $i,"\n", $l, "\n" unless ($a[$i]->[1] + 1 == $pos);
				($type, $dep) = ("", "");
				my $tmp = eval sprintf("%.3f", $a[$i]->[4]);
				push @vals, $strand * $tmp;
				if ( $a[$i]->[3] =~ /(C\S+):(\d+)/ ){
					($type, $dep) = ($1, $2);
				}
				else{
					die $l;
				}
				if($dep < $dep_cutoff){
					$flag_less = 1;
				}
			}
			#my %wigs; #$wig{chr}->{pos}->{type}-> {val1, val2...}

			if ( $flag_less == 0){
				#$wig_ref->{$chr}->{$pos}->{$type} = join(",", @vals);
				my $sum = 0;
				for my $tmp (1..$#vals){
					$sum+= abs($vals[$tmp]);
				}
				if ($sum != 0){
					$wig_ref->{$chr}->{$pos} = join(",", @vals);
				}
			}
		}#else
	}#while
}

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

# print {$fhw_C[$i]} , "\n";
# 	
# print {$fhw_CG[$i]}  , "\n";
# print {$fhw_CHG[$i]} , "\n";
# print {$fhw_CHH[$i]} , "\n";
	