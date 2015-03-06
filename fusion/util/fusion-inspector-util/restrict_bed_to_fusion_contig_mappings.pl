#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 fusion_reads_mappings.bed\n\n";

my $bed_file = $ARGV[0] or die $usage;

open (my $fh, $bed_file) or die "Error, cannot open file $bed_file";
while (<$fh>) {
    my $line = $_;
    chomp;
    my @x = split(/\t/);
    my $contig = $x[0];
    my $read = $x[3];

    my ($read_contig, $read_name) = split(/\|/, $read);
    
    if ($read_contig eq $contig) {
        print $line;
    }
}
close $fh;


exit(0);


        
