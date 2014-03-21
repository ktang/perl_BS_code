#!/usr/bin/perl -w

# 2014 Feb 27
# this script is write for GSE39045
# the input is from STDIN.

#modified on Mar2
# as the input is TAIR8, so no strand


#v0.5.1 Dec 5, 2013;
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
use strict;
use File::Spec;
my $debug = 0;
my $next_flag = 0;

#my $const = 3000000;
if($debug){
#	$const = 500;
	print STDERR "\n\n debug = 1 \n\n";
}

#my $usage = "$0 \n <outdir> <out_pre> <mC_cutoff> <depth_cutoff>\n\n";
#die $usage unless(@ARGV == 4);

my $usage = "$0 \n <indir> <outdir> <mC_cutoff> <depth_cutoff>\n\n";
die $usage unless(@ARGV == 4);

my $indir = shift or die "indir";
#my $infile = shift or die;
#my $input = shift or die;
my $outdir = shift or die "outdir";
#my $pre = shift or die "shift";

my $mC_cutoff = shift ;
die "mC_cutoff" if ($mC_cutoff < 0);
my $depth_cutoff = shift or die "depth_cutoff";

#my $label = $pre;
#my $label_removed = $pre;

#my $input = File::Spec->catfile($indir,$pre. "_isMeth.txt");
#my $input = File::Spec->catfile($indir,$infile);
#die unless (-e $input);

#print STDERR "input:\n", $input, "\n\n";

opendir(DIR, $indir ) or die;
my @cg_files = grep /CG\.gff\.gz$/ , readdir DIR;
closedir(DIR);
	
print STDERR join("\n", @cg_files), "\n\n";
	
my $isMeth_file = "/Volumes/Macintosh_HD_2/GEO_data_download_from_NCBI/GSE39045/colA_isMeth_chrC_error_separately_called.txt";
die unless(-e $isMeth_file);
my %strands_h;

if ($debug) {
	$isMeth_file = "/Volumes/Macintosh_HD_2/GEO_data_download_from_NCBI/GSE39045/debug/isMeth.txt";
}
die unless(-e $isMeth_file);

#record_strand($isMeth_file, \%strands_h);

foreach my $cg_file(@cg_files){
	
	my $pre = "NA";
	
	if ($cg_file =~ /(\S+)_BS-chr.all.single/) {
		$pre = $1;
	}else{
		die $cg_file;
	}
	
	my $files_str = get_files($indir, $pre);
	
	print STDERR join("\n", ($pre, $files_str)),"\n\n";
	
	open(IN, "zcat $files_str |") or die;
	
	my $label = $pre;
	my $label_removed = $pre;

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
		if($next_flag){
		    next;
        }
	}
	
	#open(IN, $input) or die "cannot open $input";
	open(OUT, ">$mC_wig") or die;
	open(CG,  ">$mCG_wig") or die;
	open(CHG, ">$mCHG_wig") or die;
	open(CHH, ">$mCHH_wig") or die;
			
			
	print OUT join("\t", ("track", "type=wiggle_0","name=")), "\"$label_removed\"\t","description=\"$label_removed"," methlation\"\t", join("\t", ("viewLimits=-1.0:1.0", "color=255,102,0")), "\n";
	
	print CG  join("\t", ("track", "type=wiggle_0","name=")), "\"$label\"\t","description=\"$label"," CG methlation\"\t", join("\t", ("viewLimits=-1.0:1.0", "color=205,205,0")), "\n";
	print CHG join("\t", ("track", "type=wiggle_0","name=")), "\"$label\"\t","description=\"$label"," CHG methlation\"\t", join("\t", ("viewLimits=-1.0:1.0", "color=100,149,237")), "\n";
	print CHH join("\t", ("track", "type=wiggle_0","name=")), "\"$label\"\t","description=\"$label"," CHH methlation\"\t", join("\t", ("viewLimits=-1.0:1.0", "color=205,155,155")), "\n";
	
	#close(OUT);
	#close(CG);
	#close(CHG);
	#close(CHH);
			
	#my $last_chr = "chr0";
	
	# 0      1        2      3       4        5      6       7		8
	#chr1    .       CHH     93      93      1.0000  .       .       c=1;t=0;n=1
	#chr1    .       CHH     94      94      0.0000  .       .       c=0;t=1;n=1
	#chr1    .       CHH     100     100     0.5000  .       .       c=1;t=1;n=1
	
	#my (%mC, %mCGs, %mCHGs, %mCHHs);
	my %m_h;
	my $i = 0;
	#my $head = <IN>;#skip head;
	
	my ($C_num, $CG_num, $CHG_num, $CHH_num ) = (0) x 4; 
	
	#while (<IN>){
	while (<IN>) {
		
		$i++;
		chomp;
		my @a = split "\t";
		my ($mC, $dep) = (0,0);
		
		#if ($a[8] =~ /c=(\d+);t=(\d+);n/) {
		if ($a[8] =~ /c=(\d+);t=(\d+)/) {
			$mC = $1;
			$dep = $mC+ $2;
		}
		else{
			die $_;
		}
		
		#($mC, $dep) = @a[4..5];
		
		if($mC >= $mC_cutoff and  $dep >= $depth_cutoff ){
			
			my ($chr, $pos) = ($a[0], $a[3]);
			#my $strand = ($a[2] eq "+") ? 1 : -1;
			#my $strand = 1;
			
			#if (defined $strands_h{$chr}->[$pos]) {
			#	$strand =$strands_h{$chr}->[$pos];
			#}
			#$a[6] = eval sprintf("%.3f", $a[6]);
			
			#my $val = eval sprintf("%.3f", $a[5] * $strand);
			my $val = eval sprintf("%.3f", $a[5] );
			my $type = $a[2];
#			my ($chr, $pos, $val, $type)  = ($a[0], $a[1], $a[6] * $strand, $a[3]);
			my $tmp = join("\t", ($type, $val));
			$m_h {$chr}->{$pos} = $tmp;
	
		}
	
	}
		
	close IN;
	
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
			$C_num ++;
			if($type eq "CG"){
				print CG join("\t", ($pos, $val) ), "\n";
				$CG_num++;
			}elsif($type eq "CHG"){
				print CHG join("\t", ($pos, $val) ), "\n";
				$CHG_num++;
			}elsif($type eq "CHH"){
				print CHH join("\t", ($pos, $val) ), "\n";
				$CHH_num++;
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
	
	print  $pre,"\n";
	print  join("\t", ("C",   $C_num)), "\n";
	print  join("\t", ("CG",  $CG_num)), "\n";
	print  join("\t", ("CHG", $CHG_num)), "\n";
	print  join("\t", ("CHH", $CHH_num)), "\n\n\n";

}
exit;

#record_strand($isMeth_file, \%strands_h);
# 0      1        2      3       4        5      6      	 7
#chr     pos     strand  type    num_C   depth   percentage      isMeth
#chr1    1       +       CHH     0       0       0       0
#chr1    2       +       CHH     0       0       0       0

sub record_strand{
	my ($file, $ref ) = @_;
	die unless (-e $file);
	
	print STDERR $file, "\n\n";
	print STDERR "recording strands...\n";
	
	open(IN, $file)  or die;
	my $h = <IN>;
	while (<IN>) {
		my @a = split "\t";
		my $tmp = 1;
		if ($a[2] eq "-") {
			$tmp = -1;
		}
		
		$ref->{$a[0]}->[$a[1]] = $tmp;
	}
	
	close IN;
	print STDERR "Done\n";
}

#	my $files_str = get_files($indir, $pre);
sub get_files{
	my ($dir, $pre_sub) = @_;
	die unless(-d $dir);
	opendir(DIR, $dir) or die;
	my @names = grep /$pre_sub/, readdir DIR;
	closedir (DIR);
	
	my @files = map { File::Spec->catfile($dir, $_)} @names;
	foreach my $f(@files){
		die $f unless (-e $f);
	}
	
	return (join(" ", @files));
	
}
