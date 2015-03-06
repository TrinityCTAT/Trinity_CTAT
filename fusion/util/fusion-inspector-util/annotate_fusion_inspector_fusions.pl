#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../../");
use FusionAnnotator;

my $usage = "usage: $0 fusion_inspector_summary.all_samples.dat [max_dist=100000]\n\n";

my $summary_file = $ARGV[0] or die $usage;
my $MAX_DISTANCE = $ARGV[1] || 100000;


main: {

    open (my $fh, $summary_file) or die $!;
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my $geneA = $x[1];
        my $geneB = $x[3];

        my $coord_info_A = $x[2];
        my ($chrA, $coordA) = split(/:/, $coord_info_A);
        
        my $coord_info_B = $x[4];
        my ($chrB, $coordB) = split(/:/, $coord_info_B);

        my @annotations = &FusionAnnotator::get_annotations($geneA, $geneB);
        
        if ($chrA eq $chrB && abs($coordA-$coordB) <= $MAX_DISTANCE) {
            push (@annotations, "NEIGHBORS");
        }
        
        if (@annotations) {
            push (@x, join(",", @annotations));
        }

        print join("\t", @x) . "\n";
    }
    close $fh;

    exit(0);
}

