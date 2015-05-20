#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use FusionAnnotator;


my $usage = "\n\n\tusage: $0 fusion_preds.txt\n\n";

my $fusion_preds = $ARGV[0] or die $usage;

main: {

    open (my $fh, $fusion_preds) or die "Error, cannot open file $fusion_preds";
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        if (/^\#/) {
            push (@x, "FusionAnnotations");
        }
        else {
            my $fusion = $x[0];
            my ($geneA, $geneB) = split(/--/, $fusion);
            if ($geneA eq $geneB) {
                next; # no selfies
            }
            if (my @annots = &FusionAnnotator::get_annotations($geneA, $geneB)) {
                push (@x, join(",", @annots));
            }
            else {
                push (@x, '.');
            }
        }
        print join("\t", @x) . "\n";
    }
    close $fh;

    exit(0);
}


                
