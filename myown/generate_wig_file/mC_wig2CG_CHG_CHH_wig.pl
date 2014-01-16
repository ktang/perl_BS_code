#!/usr/bin/perl -w

#v0.1
# input hava a col of strand

use strict;
use File::Spec;
my $debug = 1;
my $const = 3000000;
if($debug){
#	$const = 500;
}

my $ref_file = "/Volumes/My_Book/20120427_ShangHai_data/nodupl_July30/other19/nodupl_isMeth/colA_nodupl2_isMeth_chrC_error_separately_called.txt";
#chr     pos     strand  type    num_C   depth   percentage      isMeth
#chr1    1       +       CHH     0       0       0       0
#chr1    2       +       CHH     0       0       0       0

my $usage = "$0 <indir> <outdir>";

die $usage unless(@ARGV == 2);

my $indir = shift or die "indir";
my $outdir = shift or die "outdir";
#my $pre = shift or die "shift";

#my $infile = shift or die;

#my $label = $pre;
#my $label_removed = $pre;

#my $input = File::Spec->catfile($indir,$pre. "_isMeth.txt");

opendir(DIR, $indir) or die;
my @files = grep /\.wig$/, readdir DIR;
closedir(DIR);

my %pos_h;
unless($debug){
	read_ref($ref_file, \%pos_h);
}


foreach my $infile(@files){
	if($infile =~ /(S+)\.wig$/){
		
	#	my (%mCGs, %mCHGs, %mCHHs);
		
		my $pre = $1;
		my $input = File::Spec->catfile($indir,$infile);
		die unless (-e $input);
		my $label = "NA";
		if($pre =~ /(\S+)_nodupl/){
			$label = $1;
		}else{
			die $infile
		}

		print STDERR "input:\n", $input, "\n\n";

		#my $mC_wig   = File::Spec->catfile($outdir, $pre . "_mC.wig");
		my $mCG_wig  = File::Spec->catfile($outdir, $pre . "_mCG.wig");
		my $mCHG_wig = File::Spec->catfile($outdir, $pre . "_mCHG.wig");
		my $mCHH_wig = File::Spec->catfile($outdir, $pre . "_mCHH.wig");
	#	die "output exits \n" if(-e $mC_wig or -e $mCG_wig or -e $mCHG_wig or -e $mCHH_wig);
		die "output exits \n" if( -e $mCG_wig or -e $mCHG_wig or -e $mCHH_wig);
		
		print STDERR join("\t", ( $label, $pre)), "\n";
		
		unless ($debug){
			my %mC;
			open(CG,  ">$mCG_wig") or die;
			open(CHG, ">$mCHG_wig") or die;
			open(CHH, ">$mCHH_wig") or die;
			print CG  join("\t", ("track", "type=wiggle_0","name=")), "\"$label\"\t","description=\"$label"," CG methlation\"\t", join("\t", ("viewLimits=-1.0:1.0", "color=205,205,0")), "\n";
			print CHG join("\t", ("track", "type=wiggle_0","name=")), "\"$label\"\t","description=\"$label"," CHG methlation\"\t", join("\t", ("viewLimits=-1.0:1.0", "color=100,149,237")), "\n";
			print CHH join("\t", ("track", "type=wiggle_0","name=")), "\"$label\"\t","description=\"$label"," CHH methlation\"\t", join("\t", ("viewLimits=-1.0:1.0", "color=205,155,155")), "\n";

			close(CG);
			close(CHG);
			close(CHH);
			
			open (IN, $input) or die "cannot open $input:$!";
			my $chr = "";
			while (<IN>){
				next if (/^track/);
				if (/variableStep\s+chrom=(\w+)/){
					$chr = lc $1;
					next;
				}
				chomp;
				my ($pos, $val) = split /\t/;
				die $_ unless (defined $pos_h{$cht}->[$pos]);
				my $type = $pos_h{$cht}->[$pos];
				$mC{$type}->{$chr}->[$pos] = $val;
			}
			close(IN);

			
			output_wig($mCG_wig, \%);
			output_wig($mCHG_wig, \%mCHGs);
			output_wig($mCHH_wig, \%mCHHs);

		}
		
		
	}
	else{
		die $infile;
	}


}


exit;






#chr     pos     strand  type    num_C   depth   percentage      isMeth
#chr1    1       +       CHH     0       0       0       0
#chr1    2       +       CHH     0       0       0       0

#read_ref($ref_file, \%pos);
sub read_ref{
	my($file, $ref) = @_;
	die unless (-e $file);
	open(IN, $file) or die;
	my $head = <IN>;
	
	while(<IN>){
		chomp;
		my @a = split "\t";
		$ref->{$a[0]}->[$a[1]] = $a[3];
	}
	
	close(IN);
}

sub output_wig{
	my ($file, $ref) = @_;
	open(FILE, ">>$file") or die;
	
	foreach my $chr (sort keys %{$ref}){
		print FILE "variableStep\tchrom=$chr\n";
		foreach my $pos(sort {$a<=>$b} keys %{$ref->{$chr}}){
			print FILE join("\t", ($pos, $ref->{$chr}->{$pos})), "\n";
		}
	}
	
	close(FILE);
}