#!/usr/bin/perl -w

# input two files with head and annotation,
# output overlap of each category(TE, Gene, IG)


use strict;


my $debug = 0;
if($debug){
	print STDERR "\n\n debug =1 \n\n";
}

# my $usage = "$0\n <input> <label> <head_yes/no> STDOUT \n\n";
# die $usage unless(@ARGV == 3);

my $usage = "$0\n <input1> <label1> <input2> <label2> STDOUT";
die $usage unless ( @ARGV == 4 );


my $input1 = shift or die;
die "wrong input" unless (-e $input1);
my $label1 = shift or die;


my $input2 = shift or die;
die "wrong input" unless (-e $input2);
my $label2 = shift or die;

my ( @TE_list1, @gene_list1, @IG_list1, @all_list1 );
my ( @TE_list2, @gene_list2, @IG_list2, @all_list2 );



record_file( $input1, \@TE_list1, \@gene_list1, \@IG_list1, \@all_list1);
record_file( $input2, \@TE_list2, \@gene_list2, \@IG_list2, \@all_list2);

if($debug){
	print STDERR scalar(@TE_list1) , "\n";
	print STDERR scalar(@gene_list1) , "\n";
	print STDERR scalar(@IG_list1) , "\n";
	print STDERR scalar(@all_list1) , "\n";
	print STDERR scalar(@TE_list2) , "\n";
	print STDERR scalar(@gene_list2) , "\n";
	print STDERR scalar(@IG_list2) , "\n";
	print STDERR scalar(@all_list2) , "\n\n\n";
}



my $all_overlap  = overlap(\@all_list1, \@all_list2);
my $te_overlap   = overlap(\@TE_list1, \@TE_list2);
my $gene_overlap = overlap(\@gene_list1, \@gene_list2);
my $IG_overlap   = overlap(\@IG_list1, \@IG_list2);

print join("\t", ("feature", $label1, $label2, "overlap_num", "percentage")), "\n";

print join("\t", ("all", $all_overlap  )), "\n";
print join("\t", ("TE",  $te_overlap )), "\n";
print join("\t", ("gene",$gene_overlap )), "\n";
print join("\t", ("IG",  $IG_overlap)), "\n";

exit;

sub overlap{
	my ($first, $second) = @_;
	
	my %beds;
	my $num_f2 = 0;
	
	my $last_index2 = scalar(@{$second}) - 1;
	if($debug){
		print STDERR "last_index2 = $last_index2 \n";
	}
	for my $j(0..$last_index2){
		$num_f2++;
		my $this = $second->[$j];
		chomp $this;
#		my @a = split /\t/, $this;
		my @a = @{$this};
		my $chr = lc $a[0];
		my $start = $a[1];
		my $end = $a[2];
		for my $i ($start..$end){
			$beds{$chr}->[$i] = 1;
		}
	}

	my $num_f1 = 0;
	my $num_overlap = 0;
#	$h = <IN1>;
#	while (<IN1>){
	my $last_index1 = scalar(@{$first}) - 1;
	for my $k(0.. $last_index1){
		$num_f1 ++;
		my $this = $first->[$k];
		chomp $this;
	#	my @a = split /\t/ , $this;
		my @a = @{$this};
		my $chr = lc $a[0];
		my $start = $a[1];
		my $end = $a[2];
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
	return join("\t", ($num_f1, $num_f2, $num_overlap, "$num_overlap/$num_f1=$per" . "%"));
}


# record_file ( $input1, \@TE_list1, \@gene_list1, \@IG_list1);
sub record_file{
	my ( $input, $TE_ref, $gene_ref, $IG_ref, $all_ref) = @_;
	open (IN, $input) or die "cannot open $input: $!";
	my $head = <IN>;
	chomp $head;
	my @hs = split "\t", $head;
	my $gene_label_index = -1;
	for my $i(0..$#hs){
		if($hs[$i] eq "Gene" ){
			$gene_label_index  = $i;
			last;
		}
	}
	die if($gene_label_index == -1);
	my $gene_type_index = $gene_label_index + 2;
	my $TE_label_index  = $gene_label_index + 4;
	my $inter_label_index = $gene_label_index + 6;
	die unless ($hs[$gene_label_index] eq "Gene");
	die unless ($hs[$gene_type_index] eq "GeneType");
	die unless ($hs[$TE_label_index] eq "TE");
	die unless ($hs[$inter_label_index] eq "Intergenic");

	while( <IN> ){
		chomp;
		my @a = split /\t/;
		push @{$all_ref}, [@a[0..2]];
#	if( $a[$gene_type_index] =~ /protein/  &&   $a[$inter_label_index] ne "NONE" ){		$promoter_num++;	}
#		if($a[$inter_label_index] ne "NONE" ) {$all_inter++;}
		if($a[$TE_label_index] ne "NONE"){
			#$te_num ++;
			push @{$TE_ref}, [@a[0..2]];
		}elsif( $a[$gene_label_index] ne "NONE"){
			if($a[$gene_type_index] =~ /transposable/){
				#$te_num++;
				push @{$TE_ref}, [@a[0..2]];

			}elsif(  $a[$gene_type_index] =~ /protein/){
				#$gene_num++;
				push @{$gene_ref}, [@a[0..2]];
			}
			#elsif( $a[$gene_type_index] =~ /pseudogene/ ){
			#	$pseudogene_num++;
			#}else{
			#	$other_num ++;
			#	print STDERR $_, "\n";
			#}
		}elsif( $a[$inter_label_index] ne "NONE"  ){
			#$inter_num++;
			push @{$IG_ref}, [@a[0..2]];

		}else{
			die $_;
		}
	}
	close (IN);
}