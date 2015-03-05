#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 gmap.gff3.summary  maxDistance\n\n";

my $gmap_summary_file = $ARGV[0] or die $usage;
my $max_distance = $ARGV[1];
if (! defined $max_distance) {
    die $usage;
}

open (my $fh, $gmap_summary_file) or die $!;
while (<$fh>) {
    if (/^\#/) { # header
        print $_;
        next;
    }

    chomp;
    my $line = $_;
    my @x = split(/\t/);
    my $info = $x[3];
    
    if ($info ne '.') {
        my ($geneA, $distA, $breakpointA, $geneB, $distB, $breakpointB, $chim) = split(/;/, $info);
        if ($distA <= $max_distance && $distB <= $max_distance
            && $geneA ne $geneB
            
            ) {
            print "$line\t$chim\n";
        }
    }
    
}

exit(0);
