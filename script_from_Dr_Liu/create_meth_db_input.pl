#!/usr/bin/perl -w
# generate input for meth db
use strict;
#use Math::CDF qw(:all);

#my $conver_error = 0.03;
#my $sig_level = 0.05;


my $pos_file = shift @ARGV;
open(PF, $pos_file) or die "Can't open $pos_file: $!";
my %poss;
while(<PF>){
	next if(/pos_id/);
	chomp;
	my @temp = split /\t/;
	$poss{$temp[2]}->{$temp[3]} = [$temp[0], $temp[5]];
}

close PF;

foreach my $dir(@ARGV){
    opendir(CDIR, $dir);
    my @dirs = sort readdir CDIR;
    foreach my $d(@dirs){
				my $this_dir = "$dir/$d";
        if((-d $this_dir) && ($this_dir =~ /Lane/i)){
		my $sample_id = 0;
		if($d =~ /JKZ_?(\d+)/){
			$sample_id = $1;
		}else{
			 die "dir $d does not match pattern";
		}
           my $output = $d . "_meth_db_input.txt";
		open(OUT, ">$output") or die "Can't open $output: $!";
        print OUT join("\t", ("SampleID", "pos_id", "depth", "num_C", "percent", "type")), "\n";

		    print STDERR "Generating input under $this_dir\n";
		    opendir(SDIR, $this_dir);
		    my @files = grep {/(forw|rev)\.txt/} readdir SDIR;
			if(@files > 2){
				die "more than two forw/rev files: ", join("\t", @files);
			}
		    foreach my $file(@files){
				$file = $this_dir . "/" . $file;
				open(IN, $file) or die "Can't open $file: $!";
				while(<IN>){
					chomp;
					my ($chr, $pos1, $pos2, $t, $percent, $strand) = split /\t/;
					$pos1 += 1;
					my ($pos_id, $type) = (0, "NONE");
					if(defined $poss{$chr}->{$pos1}){
						$pos_id = $poss{$chr}->{$pos1}->[0];
						$type = $poss{$chr}->{$pos1}->[1];
					}else{
						die "chr $chr, pos $pos1, strand $strand, ", 
						  "not defined in position file";
					}
					#my $isMeth = 0;
                    my $num_C = round($t * $percent);
                    #if($t > 0){
					
                    #my $p = pbinom($num_C, $t, (1-$conver_error));
					#if($p >= $sig_level){
				    #		$isMeth = 1;
					#}
					#}
					print OUT join("\t", ($sample_id, $pos_id, $t, $num_C, $percent, $type)), "\n";
               }
			}
			
	}
close OUT;
}
}
sub round {
    my($number) = shift;
    #return int($number + .5);
    return int($number + .5 * ($number <=> 0)); # take care of negative numbers too

}

