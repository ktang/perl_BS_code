#!/usr/bin/perl -w

# acgt2wig_v3.pl
#v3: input is the actg *_forw.txt and *_rev.txt
# output the wig files that has depth >= XX and mC>=1
# output three digit in the wig
#


#isMeth2wig v0.5.1 Dec 5, 2013;
#
# in this version of script, the hash record is different from
# previous,
# as in previous version, there are four hashes to record the info
# now use only one hash to record
# val_type;

# val = eval sprintf("%.3f", $val)

#perl -e '$a= 0; print eval sprintf("%.3f", $a), "\n\n"'
#0
#perl -e '$a= 0.222522; print eval sprintf("%.3f", $a), "\n\n"'
#0.223
BEGIN { push @INC, '/Users/tang58/scripts_all/perl_code/Modules' }
use Kai_Module;
use strict;
use File::Spec;
my $debug = 0;

#my $const = 3000000;
if($debug){
#	$const = 500;
	print STDERR "\n\n debug = 1 \n\n";
}

my $usage = "$0 \n <acgt_dir> <acgt_pre> <outdir> <out_pre> <mC_cutoff> <depth_cutoff>\n\n";

die $usage unless(@ARGV == 6);

my $indir = shift or die "indir";
#my $infile = shift or die;
#my $input = shift or die;
my $inpre = shift or die;
my $outdir = shift or die "outdir";
my $pre = shift or die "shift";

my $mC_cutoff = shift or die "mC_cutoff";
my $depth_cutoff = shift or die "depth_cutoff";

my $label = $pre;
my $label_removed = $pre;

#my $input = File::Spec->catfile($indir,$pre. "_isMeth.txt");
#my $input = File::Spec->catfile($indir,$infile);
#die unless (-e $input);

#print STDERR "input:\n", $input, "\n\n";

my $forw_file = File::Spec->catfile( $indir, $inpre . "_forw.txt" ) ;
my $rev_file  = File::Spec->catfile( $indir, $inpre . "_rev.txt" ) ;
die unless (-e $forw_file );
die unless (-e $rev_file );

my $post = "_mC" . $mC_cutoff . "_dep" . $depth_cutoff;

my $mC_wig   = File::Spec->catfile($outdir, $pre . $post.  "_mC.wig");
my $mCG_wig  = File::Spec->catfile($outdir, $pre . $post. "_mCG.wig");
my $mCHG_wig = File::Spec->catfile($outdir, $pre . $post. "_mCHG.wig");
my $mCHH_wig = File::Spec->catfile($outdir, $pre . $post. "_mCHH.wig");

die "output exits \n" if(-e $mC_wig or -e $mCG_wig or -e $mCHG_wig or -e $mCHH_wig);

if($debug){
	print STDERR "output:\n";
	print STDERR join("\n", ($mC_wig, $mCG_wig, $mCHG_wig, $mCHH_wig)), "\n\n";
	
	print STDERR "OK\n\n";
	
	exit;
}

my %m_h;
read_acgt_file (  $forw_file, 1, \%m_h);
read_acgt_file (  $rev_file ,-1, \%m_h);


open(OUT, ">$mC_wig") or die;
open(CG,  ">$mCG_wig") or die;
open(CHG, ">$mCHG_wig") or die;
open(CHH, ">$mCHH_wig") or die;
		
		
print OUT join("\t", ("track", "type=wiggle_0","name=")), "\"$label_removed\"\t","description=\"$label_removed"," methlation\"\t", join("\t", ("viewLimits=-1.0:1.0", "color=255,102,0")), "\n";

print CG  join("\t", ("track", "type=wiggle_0","name=")), "\"$label\"\t","description=\"$label"," CG methlation\"\t", join("\t", ("viewLimits=-1.0:1.0", "color=205,205,0")), "\n";
print CHG join("\t", ("track", "type=wiggle_0","name=")), "\"$label\"\t","description=\"$label"," CHG methlation\"\t", join("\t", ("viewLimits=-1.0:1.0", "color=100,149,237")), "\n";
print CHH join("\t", ("track", "type=wiggle_0","name=")), "\"$label\"\t","description=\"$label"," CHH methlation\"\t", join("\t", ("viewLimits=-1.0:1.0", "color=205,155,155")), "\n";



foreach my $chr (sort keys %m_h ){
#	print FILE "variableStep\tchrom=$chr\n";
	print OUT "variableStep\tchrom=$chr\n";
	print CG "variableStep\tchrom=$chr\n";
	print CHG "variableStep\tchrom=$chr\n";
	print CHH "variableStep\tchrom=$chr\n";
	
	foreach my $pos(sort {$a<=>$b} keys %{ $m_h{$chr} } ){
		#print FILE join("\t", ($pos, $ref->{$chr}->{$pos})), "\n";
		my $tmp = $m_h{$chr}->{$pos};
		my ($type, $val ) = split "\t",  $tmp;
		print OUT join("\t", ($pos, $val) ), "\n";

		if($type eq "CG"){
			print CG join("\t", ($pos, $val) ), "\n";
		}elsif($type eq "CHG"){
			print CHG join("\t", ($pos, $val) ), "\n";
			#$CHG_num++;
		}elsif($type eq "CHH"){
			print CHH join("\t", ($pos, $val) ), "\n";
		#	$CHH_num++;
		}else{
			die "$chr, $pos, $tmp\n\n";
		}
	}
	
}
#print C join("\t", ($type, $val) ), "\n";

close(OUT);
close(CG);
close(CHG);
close(CHH);
exit;


sub read_acgt_file{
	my ( $file, $stand_num, $mCs_ref) = @_;
	
	die unless (-e $file);
	
	open(IN, $file) or die;
	while (<IN>) {
		chomp;
		my @a = split "\t";
		my $chr = $a[0];
		last if ($chr eq "chrC");
		$chr = Kai_Module::simple_chr($chr);
		my $pos = $a[1] + 1;
		my ($type, $dep) = split ":", $a[3];
		#my $val = $a[4] * $stand_num;
		next if (  $dep < $depth_cutoff);
		
		my $val = eval sprintf ("%.3f", $a[4] );
		if ($val != 0){
			$val *= $stand_num;
			my $tmp = join("\t", ($type, $val));
			$mCs_ref->{$chr}->{$pos} = $tmp;
		}
	}
	close IN;
}
