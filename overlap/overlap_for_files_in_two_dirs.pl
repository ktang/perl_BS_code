#!/usr/bin/perl -w

#overlap_for_files_in_two_dirs.pl
# one dir called db_dir
#one is bench_dir, for each file in bench dir, cal overlap number of it with db_dir;
# each file has one output
use strict;
use File::Spec;

my $debug = 0;

if($debug){
	print STDERR "debug = $debug\n\n";
}

my %types = ( "bed"  =>1,
	      "coor" => 1	
);

my $usage = "$0 \n<indir_bench mutants_cared> <indir_db> <outdir> <coor type:bed or coor>\n\n";

die $usage unless (@ARGV == 4);
my $indir_bench = shift or die $usage;
my $db_dir =  shift or die $usage;
#my $output = shift or die "output";
my $outdir = shift or die "outdir";
#die if(-e $output);

my $coor_type = shift or die;
unless (defined $types{$coor_type}){
	die "coor type:bed or coor";
}

#die "wrong dir" unless (-d $dir );
die "wrong dir" unless (-d $indir_bench );
die "wrong dir" unless (-d $db_dir );
die "wrong dir" unless (-d $outdir );

opendir (DIR, $indir_bench) or die "cannot open $indir_bench: $!";
#my @files = grep /_Cnum4\.txt$/, readdir DIR;
my @files = grep /\.txt$/, readdir DIR;
closedir DIR;


opendir(DB, $db_dir) or die;
my @dbs = grep /\.txt$/, readdir DB;
closedir DB;

if ($debug){
	print STDERR join("\n", @files),"\n\n";
	#print STDERR join("\n", @dbs),"\n";
	
	#exit;
}



for my $i(0.. ($#files)){
	my $file = $files[$i];
	#get_output ($indir_bench, $)
	print STDERR "\n", $i, "...\n";
	
	if ($file =~ /(\S+)\.txt$/) {
		my $pre = $1;
		if ($debug) {
			print STDERR $pre, "\n\n";#code
		}
		
		my $input = File::Spec->catfile($indir_bench, $file);
		die unless (-e $input);
		
		my $self_num = get_line_number($input);
		
		my $output = File::Spec->catfile($outdir, "overlap_" . $pre . ".txt" );
		die if (-e $output);
		open(OUT, ">>$output") or die;
		
		print OUT join("\t", ("list_num", "overlap", "per_of_list", "per_of_self", "file_name")), "\n";
		
		print OUT join("\t", ($self_num, "--", "--", "--", $file)) , "\n";
		
		if (! $debug) {
			my $j = 0;
			foreach my $db(@dbs){
				$j++;
				print STDERR $j, " ";
				my $db_file = File::Spec->catfile($db_dir, $db);
				die unless (-e $db_file);
				my @tmp = get_overlap($input, $db_file, $coor_type, $self_num);
				print OUT join("\t", (@tmp, $db)), "\n"
			}
		}
		
		
		
		close OUT;
		
	}else{
		die $file;
	}
	
}

exit;

#my @tmp = get_overlap($input, $db_file, $coor_type);
sub get_overlap{
	my ($input_sub, $db_file_sub, $coor_type_sub,$self_num_sub) = @_;
	
	my %r;
	
	die unless (-e $input_sub);
	die unless (-e $db_file_sub);
	
	open(DB, $db_file_sub) or die;
	my $h = <DB>;
	my $n_db = 0;
	while (<DB>) {
		$n_db ++;
		chomp;
		my @a = split "\t";
		my ($chr, $start, $end);
		if ($coor_type_sub eq "bed" ) {
			($chr, $start, $end) = @a[0..2];
		}elsif($coor_type_sub eq "coor" ){
			if ($a[0] =~ /(\S+):(\d+)-(\d+)/) {
				($chr, $start, $end) = ($1, $2, $3);
			}else{
				die $a[0];
			}
			
		}else{
			die "coor";
		}
		
		$chr = simple_chr($chr);
		
		for my $i ($start..$end){
			$r{$chr}->[$i] = 1;
		}
	}
	close DB;
	
	open(IN, $input_sub) or die;
	$h = <IN>;
	my $n_list = 0;
	my $overlap = 0;
	while (<IN>) {
		$n_list ++;
		chomp;
		my @a = split "\t";
		my ($chr, $start, $end);
		if ($coor_type_sub eq "bed" ) {
			($chr, $start, $end) = @a[0..2];
		}elsif($coor_type_sub eq "coor" ){
			if ($a[0] =~ /(\S+):(\d+)-(\d+)/) {
				($chr, $start, $end) = ($1, $2, $3);
			}else{
				die $a[0];
			}
			
		}else{
			die "coor";
		}
		
		$chr = simple_chr($chr);
		my $flag = 0;
		for my $i ($start..$end){
			if ( defined $r{$chr}->[$i]) {
				$flag = 1;
				last;
			}
		}
		if ($flag) {
			$overlap++;
		}
		
	}
	close IN;

	die "$n_list != $self_num_sub " if ($n_list != $self_num_sub);
	
	my $p_self = sprintf ("%.2f", 100 * $overlap / $self_num_sub);
	my $p_db = sprintf ("%.2f", 100 * $overlap / $n_db);
	
	return ($n_db, $overlap, $p_db, $p_self);
	
# ("list_num", "overlap", "per_of_list", "per_of_self", "file_name")), "\n";

}

sub simple_chr{
	my ($chr) = @_;
	if( $chr =~ /chr/i){
		$chr =~  s/chr//i;
	}
	if($chr eq "M" ){
		$chr = "Mt";
	}elsif( $chr eq "C"){
		$chr = "Pt";
	}
	return $chr;
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

#		my $self_num = get_line_number($input);
sub get_line_number{
	my ($file) =  @_;
	die unless (-e $file);
	open(IN, $file) or die;
	my $i = -1;
	while (<IN>) {
		$i++;
	}
	
	close IN;
	return ($i);
}