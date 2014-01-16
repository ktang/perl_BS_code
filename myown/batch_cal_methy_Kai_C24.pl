#!/usr/bin/perl -w
#
use strict;
my $debug = 0;
print STDERR "debug: $debug\n";
my $max_id = 42858083; # C24

my $cal_mehty_script = "/Users/tang58/Kai_BS/myown/calculate_methy_status_notype_Kai.pl";

my $log = "cal_meth.log";
my %rates = ("JKZ19_C24_Luc_149_s2_meth_db_input.txt" => 0.05225119,
			 "JKZ20_C24_Luc_150_s3_meth_db_input.txt" => 0.03963806,
			 "JKZ21_ros1_1_s5_meth_db_input.txt" => 0.04458593,
			 "JKZ22_ros1_1_152_s6_meth_db_input.txt" => 0.04652252,
			 "JKZ23_ros3_1_153_s7_meth_db_input.txt" => 0.137921,
			 "JKZ24_ros3_1_154_s8_meth_db_input.txt" =>  0.03968124);
foreach my $file(sort keys %rates){
	die unless (-e $file);
	print STDERR "working on file $file ...\n";
	my ($pre, $ext) = split /\./, $file;
	my $outfile = $pre . "_onetype" . "." . $ext;
	my $cutoff_file = $pre . "_cutoff" . "." . $ext;
	die "$outfile exists" unless (!(-e $outfile ));
	my $r = $rates{$file};
#	`perl calculate_methy_status_notype.pl $file $r > $outfile 2>>$log`;
	my $cmd = "date; time perl $cal_mehty_script $file $r $cutoff_file $max_id > $outfile 2>>$log ";
	print STDERR "$cmd\n\n";
	if (!$debug) {
		system("$cmd");
	} 
}

