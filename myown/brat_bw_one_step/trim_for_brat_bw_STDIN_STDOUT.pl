#!/usr/bin/perl -w

use strict;

my $debug = 0;
if ($debug) {
	print STDERR "debug:$debug\n\n";
}



my $length = 101;
my $beginning_trimmed_bp = 13;
my $end_trimmed_bp = 18;

print STDERR "\n\nlength = $length, beginning_trimmed_bp = $beginning_trimmed_bp and end_trimmed_bp =$end_trimmed_bp may need to be Changed!!! \n\n\n";


my $left_len = $length - $beginning_trimmed_bp - $end_trimmed_bp;

my $i = 0;
while ( <> ) {
	$i++;
	if ($i%4 == 2) {
		chomp;
		my $seq = substr( $_,  $beginning_trimmed_bp,  $left_len);
		print join( "\t", ($seq,  $beginning_trimmed_bp,  $end_trimmed_bp)),"\n";
	}
	
}


exit;