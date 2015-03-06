#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");

use __GLOBALS__;
use FusionAnnotator;

my $usage = "usage: $0 gmap.map.gff3.chims_described.D0\n\n";

my $chim_file = $ARGV[0] or die $usage;

main: {

    open (my $fh, $chim_file) or die "Error, cannot open file $chim_file";
    while (<$fh>) {
        if (/^\#/) { next; }
        unless (/\w/) { next; }
        chomp;
        my @x = split(/\t/);
        
        my $trinity_acc = $x[0];
        my $fusion_info = $x[3];
        my ($geneA, $deltaA, $breakpointA, $geneB, $deltaB, $breakpointB, $fusion_name) = split(/;/, $fusion_info);

        my @annotations = &FusionAnnotator::get_annotations($geneA, $geneB);

        my ($chrA, $coordA) = split(/:/, $breakpointA);
        my ($chrB, $coordB) = split(/:/, $breakpointB);

        print join("\t", $fusion_name, $trinity_acc, $geneA, $chrA, $coordA, $geneB, $chrB, $coordB, 
                   join(",", @annotations) ) . "\n";
    }
    
    close $fh;

    exit(0);
    
}
